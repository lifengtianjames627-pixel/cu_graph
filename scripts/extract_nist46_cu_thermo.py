"""Extract overall Cu2+ complex thermodynamics from NIST SRD 46.

Source: Burgess, D. R. (2004), NIST SRD 46, https://doi.org/10.18434/M32154
"""
from __future__ import annotations

import hashlib
import json
import re
import sqlite3
import tempfile
import zipfile
from contextlib import closing
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data" / "SRD_46_SQL.zip"
DB = ROOT / "data" / "NIST_SRD_46.db"
OUT_DIR = ROOT / "data"

EXPECTED_SHA256 = "141269bb8c6d9e8271a5b4bff35f7c9fa913938bb50a62c096847de67e382d18"
BATCH_SIZE = 5000

# Ligands relevant to the cu_calc UV–vis study (name substrings, case-insensitive).
LIGAND_PATTERNS = [
    ("EDTA", ["edta", "ethylenediaminetetraacet"]),
    ("Cl-", ["chloride", " chloro", "cl-"]),
    ("NO3-", ["nitrate", " nitrato"]),
    ("SCN-", ["thiocyanate", " isothiocyanate", "scn"]),
    ("OH-", ["hydroxide", " hydroxo", "oh-"]),
    ("SO4 2-", ["sulfate", " sulphate", "sulfato"]),
    ("H2O (aqua)", ["water", "aqua"]),
]

CU_METAL = "Cu 2+ "


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def sqlite_type(mysql_type: str) -> str:
    value = mysql_type.lower()
    if value.startswith(("int", "tinyint", "smallint", "mediumint", "bigint")):
        return "INTEGER"
    if value.startswith(("float", "double", "decimal", "numeric", "real")):
        return "REAL"
    if value.startswith(("blob", "binary", "varbinary")):
        return "BLOB"
    return "TEXT"


def parse_mysql_schema(sql_text: str) -> dict:
    table_match = re.search(r"CREATE TABLE\s+`([^`]+)`\s*\(", sql_text)
    if not table_match:
        raise ValueError("CREATE TABLE statement not found")
    table_name = table_match.group(1)
    body_start = table_match.end()
    body_end = sql_text.find(") ENGINE=", body_start)
    if body_end < 0:
        raise ValueError(f"End of CREATE TABLE not found for {table_name}")
    body = sql_text[body_start:body_end]

    columns = []
    primary_key = None
    indexes = []
    for line in body.splitlines():
        column_match = re.match(r"\s*`([^`]+)`\s+([^\s,]+)", line)
        if column_match:
            columns.append(
                {
                    "name": column_match.group(1),
                    "mysql_type": column_match.group(2),
                    "sqlite_type": sqlite_type(column_match.group(2)),
                }
            )
            continue
        primary_match = re.match(r"\s*PRIMARY KEY\s*\(`([^`]+)`\)", line)
        if primary_match:
            primary_key = primary_match.group(1)
            continue
        index_match = re.match(r"\s*KEY\s+`([^`]+)`\s*\(([^)]+)\)", line)
        if index_match:
            index_columns = re.findall(r"`([^`]+)`", index_match.group(2))
            if index_columns:
                indexes.append({"name": index_match.group(1), "columns": index_columns})

    if not columns or primary_key is None:
        raise ValueError(f"Invalid schema for {table_name}")
    return {
        "table": table_name,
        "columns": columns,
        "primary_key": primary_key,
        "indexes": indexes,
    }


def mysql_unescape(value: str) -> str | None:
    if value == r"\N":
        return None
    mapping = {
        "0": "\0",
        "b": "\b",
        "n": "\n",
        "r": "\r",
        "t": "\t",
        "Z": "\x1a",
        "\\": "\\",
        "'": "'",
        '"': '"',
    }
    result = []
    index = 0
    while index < len(value):
        char = value[index]
        if char == "\\" and index + 1 < len(value):
            escaped = value[index + 1]
            result.append(mapping.get(escaped, escaped))
            index += 2
        else:
            result.append(char)
            index += 1
    return "".join(result)


def convert_value(value: str | None, column_type: str):
    if value is None:
        return None
    if column_type == "INTEGER" and value != "":
        return int(value)
    if column_type == "REAL" and value != "":
        return float(value)
    return value


def iter_rows(raw_data: bytes, schema: dict):
    text = raw_data.decode("cp1252")
    expected_columns = len(schema["columns"])
    types = [column["sqlite_type"] for column in schema["columns"]]
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != expected_columns:
            raise ValueError(
                f"{schema['table']}.txt line {line_number}: expected "
                f"{expected_columns} fields, found {len(fields)}"
            )
        yield tuple(
            convert_value(mysql_unescape(value), column_type)
            for value, column_type in zip(fields, types, strict=False)
        )


def create_table(connection: sqlite3.Connection, schema: dict) -> None:
    definitions = []
    for column in schema["columns"]:
        definition = f"{quote_identifier(column['name'])} {column['sqlite_type']}"
        if column["name"] == schema["primary_key"]:
            definition += " PRIMARY KEY"
        definitions.append(definition)
    connection.execute(
        f"CREATE TABLE {quote_identifier(schema['table'])} ({', '.join(definitions)})"
    )


def create_indexes(connection: sqlite3.Connection, schema: dict) -> None:
    for index in schema["indexes"]:
        index_name = f"idx_{schema['table']}_{index['name']}"
        columns = ", ".join(quote_identifier(value) for value in index["columns"])
        connection.execute(
            f"CREATE INDEX {quote_identifier(index_name)} "
            f"ON {quote_identifier(schema['table'])} ({columns})"
        )


def insert_table(connection: sqlite3.Connection, schema: dict, raw_data: bytes) -> int:
    column_names = [column["name"] for column in schema["columns"]]
    placeholders = ", ".join("?" for _ in column_names)
    quoted_columns = ", ".join(quote_identifier(value) for value in column_names)
    statement = (
        f"INSERT INTO {quote_identifier(schema['table'])} "
        f"({quoted_columns}) VALUES ({placeholders})"
    )
    count = 0
    batch = []
    for row in iter_rows(raw_data, schema):
        batch.append(row)
        if len(batch) >= BATCH_SIZE:
            connection.executemany(statement, batch)
            count += len(batch)
            batch.clear()
    if batch:
        connection.executemany(statement, batch)
        count += len(batch)
    return count


def build_database(source: Path, output: Path) -> None:
    archive_hash = sha256_file(source)
    if archive_hash != EXPECTED_SHA256:
        raise ValueError(f"Archive checksum mismatch: {archive_hash}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=output.stem + ".",
        suffix=".tmp.db",
        dir=output.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)

    try:
        with zipfile.ZipFile(source) as archive, closing(
            sqlite3.connect(temporary_path)
        ) as connection:
            sql_members = sorted(name for name in archive.namelist() if name.endswith(".sql"))
            connection.execute("PRAGMA journal_mode = DELETE")
            connection.execute("BEGIN")
            schemas = []
            for sql_member in sql_members:
                sql_data = archive.read(sql_member)
                schema = parse_mysql_schema(sql_data.decode("ascii", errors="replace"))
                txt_member = f"{schema['table']}.txt"
                raw_data = archive.read(txt_member)
                create_table(connection, schema)
                insert_table(connection, schema, raw_data)
                schemas.append(schema)
            for schema in schemas:
                create_indexes(connection, schema)
            connection.commit()
            connection.execute("VACUUM")

        temporary_path.replace(output)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def is_overall_beta(name: str) -> bool:
    """True for cumulative/overall beta_n, not stepwise K_n."""
    n = name.strip()
    if re.search(r"\bK\s*\d", n, re.I):
        return False
    if re.search(r"\bK\s*=", n, re.I):
        return False
    if re.search(r"stepwise|successive|sequential", n, re.I):
        return False
    return bool(re.search(r"\bbeta\b|\bβ\b|\*beta|\*β", n, re.I))


def ligand_group(name: str, formula: str) -> str | None:
    text = f"{name} {formula}".lower()
    for label, patterns in LIGAND_PATTERNS:
        if any(p in text for p in patterns):
            return label
    return None


def fetch_cu_records(connection: sqlite3.Connection) -> pd.DataFrame:
    query = """
    SELECT
        v.verkn_ligand_metalID AS record_id,
        m.name_metal AS metal,
        l.name_ligand AS ligand_name,
        l.formula AS ligand_formula,
        lc.name_ligandclass AS ligand_class,
        b.name_beta_definition AS equilibrium,
        v.temperature AS temp_C,
        v.ionicstrength AS ionic_strength,
        s.name_solvent AS solvent,
        ct.name_constanttyp AS value_type,
        v.constant AS value,
        fn.name_footnote AS footnote
    FROM verkn_ligand_metal AS v
    JOIN metal AS m ON m.metalID = v.metalNr
    JOIN liganden AS l ON l.ligandenID = v.ligandenNr
    LEFT JOIN ligand_class AS lc ON lc.ligand_classID = l.ligand_classNr
    LEFT JOIN beta_definition AS b ON b.beta_definitionID = v.beta_definitionNr
    JOIN constanttyp AS ct ON ct.constanttypID = v.constanttypNr
    LEFT JOIN solvent AS s ON s.solventID = v.solventNr
    LEFT JOIN footnote AS fn ON fn.footnoteID = v.footnoteNr
    WHERE m.name_metal = ?
      AND ct.name_constanttyp IN ('K', 'H', 'S')
    """
    return pd.read_sql_query(query, connection, params=(CU_METAL,))


def pivot_thermo(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["ligand_group"] = df.apply(
        lambda r: ligand_group(str(r["ligand_name"]), str(r["ligand_formula"])), axis=1
    )
    df = df[df["ligand_group"].notna()]
    df = df[df["equilibrium"].apply(lambda x: is_overall_beta(str(x)))]

    df["value_type"] = df["value_type"].map({"K": "logK", "H": "delta_H_kJ_mol", "S": "delta_S_J_mol_K"})
    df["match_key"] = (
        df["ligand_group"].astype(str)
        + "|"
        + df["equilibrium"].astype(str)
        + "|"
        + df["temp_C"].astype(str)
        + "|"
        + df["ionic_strength"].astype(str)
        + "|"
        + df["solvent"].fillna("").astype(str)
    )

    wide = df.pivot_table(
        index=[
            "ligand_group",
            "ligand_name",
            "ligand_formula",
            "equilibrium",
            "temp_C",
            "ionic_strength",
            "solvent",
            "footnote",
        ],
        columns="value_type",
        values="value",
        aggfunc="first",
    ).reset_index()

    for col in ("logK", "delta_H_kJ_mol", "delta_S_J_mol_K"):
        if col not in wide.columns:
            wide[col] = pd.NA

    wide = wide.sort_values(["ligand_group", "temp_C", "ionic_strength", "equilibrium"])
    return wide


def summarize_best(wide: pd.DataFrame) -> pd.DataFrame:
    """Pick one representative row per ligand group (25 C, I=0 or lowest I)."""

    def score(row: pd.Series) -> tuple:
        temp = row["temp_C"]
        ionic = row["ionic_strength"]
        temp_penalty = abs(float(temp) - 25.0) if pd.notna(temp) else 999.0
        ionic_penalty = float(ionic) if pd.notna(ionic) else 999.0
        completeness = sum(pd.notna(row[c]) for c in ("logK", "delta_H_kJ_mol", "delta_S_J_mol_K"))
        return (-completeness, temp_penalty, ionic_penalty)

    rows = []
    for group, chunk in wide.groupby("ligand_group"):
        chunk = chunk.copy()
        chunk["_score"] = chunk.apply(score, axis=1)
        best = chunk.sort_values("_score").iloc[0]
        rows.append(best.drop("_score"))
    summary = pd.DataFrame(rows)
    summary = summary[
        [
            "ligand_group",
            "ligand_name",
            "equilibrium",
            "temp_C",
            "ionic_strength",
            "solvent",
            "logK",
            "delta_H_kJ_mol",
            "delta_S_J_mol_K",
            "footnote",
        ]
    ]
    return summary.sort_values("ligand_group")


def main() -> None:
    if not ARCHIVE.is_file():
        raise FileNotFoundError(f"Missing archive: {ARCHIVE}")

    if not DB.is_file():
        print(f"Building SQLite database from {ARCHIVE.name} ...")
        build_database(ARCHIVE, DB)
        print(f"Built {DB}")

    with closing(sqlite3.connect(DB)) as connection:
        raw = fetch_cu_records(connection)
        beta_samples = pd.read_sql_query(
            """
            SELECT DISTINCT b.name_beta_definition
            FROM verkn_ligand_metal AS v
            JOIN metal AS m ON m.metalID = v.metalNr
            JOIN beta_definition AS b ON b.beta_definitionID = v.beta_definitionNr
            WHERE m.name_metal = ?
            ORDER BY 1
            LIMIT 80
            """,
            connection,
            params=(CU_METAL,),
        )

    OUT_DIR.mkdir(exist_ok=True)
    beta_samples.to_csv(OUT_DIR / "nist46_cu_beta_definitions_sample.csv", index=False)

    wide = pivot_thermo(raw)
    summary = summarize_best(wide)

    wide_path = OUT_DIR / "nist46_cu_overall_thermo_all.csv"
    summary_path = OUT_DIR / "nist46_cu_overall_thermo_summary.csv"
    wide.to_csv(wide_path, index=False)
    summary.to_csv(summary_path, index=False)

    print(f"Raw Cu2+ records (K/H/S): {len(raw)}")
    print(f"Overall beta rows (filtered ligands): {len(wide)}")
    print(f"Summary rows: {len(summary)}")
    print(f"wrote {wide_path}")
    print(f"wrote {summary_path}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
