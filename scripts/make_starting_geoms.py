"""Build starting XYZ geometries for Cu aqua vs EDTA anion-exchange DFT."""
from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
(ROOT / "ligands").mkdir(exist_ok=True)
(ROOT / "complexes").mkdir(exist_ok=True)
(ROOT / "refs").mkdir(exist_ok=True)


def write_xyz(path: Path, atoms: list[tuple], comment: str) -> None:
    lines = [str(len(atoms)), comment]
    for el, x, y, z in atoms:
        lines.append(f"{el:2s}  {x:12.6f} {y:12.6f} {z:12.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {path.relative_to(ROOT)}  natoms={len(atoms)}")


def water_hs(o_xyz: np.ndarray, v: np.ndarray) -> list[tuple]:
    v = v / np.linalg.norm(v)
    tmp = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    a = np.cross(v, tmp)
    a /= np.linalg.norm(a)
    b = np.cross(v, a)
    h1 = o_xyz + 0.96 * (0.6 * a + 0.3 * b) - 0.25 * v
    h2 = o_xyz + 0.96 * (-0.6 * a + 0.3 * b) - 0.25 * v
    return [("H", *h1), ("H", *h2)]


def main() -> None:
    write_xyz(
        ROOT / "ligands/h2o.xyz",
        [
            ("O", 0.000000, 0.000000, 0.117300),
            ("H", 0.000000, 0.757000, -0.469200),
            ("H", 0.000000, -0.757000, -0.469200),
        ],
        "H2O | PubChem CID 962 https://pubchem.ncbi.nlm.nih.gov/compound/962",
    )
    write_xyz(
        ROOT / "ligands/cl-.xyz",
        [("Cl", 0.0, 0.0, 0.0)],
        "Cl- | PubChem CID 312 https://pubchem.ncbi.nlm.nih.gov/compound/312",
    )
    write_xyz(
        ROOT / "ligands/br-.xyz",
        [("Br", 0.0, 0.0, 0.0)],
        "Br- | PubChem CID 259 https://pubchem.ncbi.nlm.nih.gov/compound/259",
    )
    write_xyz(
        ROOT / "ligands/oh-.xyz",
        [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.97)],
        "OH- | PubChem CID 961 https://pubchem.ncbi.nlm.nih.gov/compound/961",
    )
    write_xyz(
        ROOT / "ligands/no3-.xyz",
        [
            ("N", 0.000000, 0.000000, 0.000000),
            ("O", 1.240000, 0.000000, 0.000000),
            ("O", -0.620000, 1.073900, 0.000000),
            ("O", -0.620000, -1.073900, 0.000000),
        ],
        "NO3- | PubChem CID 943 https://pubchem.ncbi.nlm.nih.gov/compound/943",
    )
    write_xyz(
        ROOT / "ligands/scn-.xyz",
        [
            ("S", -1.560000, 0.000000, 0.000000),
            ("C", 0.000000, 0.000000, 0.000000),
            ("N", 1.170000, 0.000000, 0.000000),
        ],
        "SCN- | PubChem CID 9322 https://pubchem.ncbi.nlm.nih.gov/compound/9322",
    )

    # EDTA4- extended start (relaxes in opt)
    edta: list[tuple] = [
        ("N", -1.450, 0.000, 0.000),
        ("N", 1.450, 0.000, 0.000),
        ("C", -0.725, 1.200, 0.000),
        ("C", 0.725, 1.200, 0.000),
        ("C", -2.350, -1.100, 0.900),
        ("C", -3.700, -1.350, 0.700),
        ("O", -4.350, -0.550, 0.100),
        ("O", -4.200, -2.400, 1.100),
        ("C", -2.350, -1.100, -0.900),
        ("C", -3.700, -1.350, -0.700),
        ("O", -4.350, -0.550, -0.100),
        ("O", -4.200, -2.400, -1.100),
        ("C", 2.350, -1.100, 0.900),
        ("C", 3.700, -1.350, 0.700),
        ("O", 4.350, -0.550, 0.100),
        ("O", 4.200, -2.400, 1.100),
        ("C", 2.350, -1.100, -0.900),
        ("C", 3.700, -1.350, -0.700),
        ("O", 4.350, -0.550, -0.100),
        ("O", 4.200, -2.400, -1.100),
    ]
    for h in [
        (-0.725, 1.800, 0.900),
        (-0.725, 1.800, -0.900),
        (0.725, 1.800, 0.900),
        (0.725, 1.800, -0.900),
        (-2.000, -1.900, 1.500),
        (-1.700, -0.400, 1.500),
        (-2.000, -1.900, -1.500),
        (-1.700, -0.400, -1.500),
        (2.000, -1.900, 1.500),
        (1.700, -0.400, 1.500),
        (2.000, -1.900, -1.500),
        (1.700, -0.400, -1.500),
    ]:
        edta.append(("H", *h))
    write_xyz(
        ROOT / "ligands/edta4-.xyz",
        edta,
        "EDTA4- start | EDTA PubChem CID 6049 https://pubchem.ncbi.nlm.nih.gov/compound/6049 ; Cu-EDTA PubChem https://pubchem.ncbi.nlm.nih.gov/compound/23673837",
    )

    eq, ax = 2.00, 2.30
    dirs = [(eq, 0, 0), (-eq, 0, 0), (0, eq, 0), (0, -eq, 0), (0, 0, ax), (0, 0, -ax)]
    cu_h2o6: list[tuple] = [("Cu", 0.0, 0.0, 0.0)]
    for x, y, z in dirs:
        o = np.array([x, y, z], float)
        cu_h2o6.append(("O", x, y, z))
        cu_h2o6.extend(water_hs(o, o))
    write_xyz(
        ROOT / "complexes/cu_h2o6.xyz",
        cu_h2o6,
        "[Cu(H2O)6]2+ JT start | COD hexaaqua nitrate salt https://www.crystallography.net/cod/8102245.html ; chalcanthite COD 9011339 https://www.crystallography.net/cod/9011339.html",
    )

    def make_cu_h2o5_x(x_atoms: list[tuple], tag: str, comment: str, r_cu_x: float) -> None:
        atoms: list[tuple] = [("Cu", 0.0, 0.0, 0.0)]
        keep = [(eq, 0, 0), (-eq, 0, 0), (0, eq, 0), (0, -eq, 0), (0, 0, -ax)]
        for x, y, z in keep:
            o = np.array([x, y, z], float)
            atoms.append(("O", x, y, z))
            atoms.extend(water_hs(o, o))
        dx, dy, dz = x_atoms[0][1], x_atoms[0][2], x_atoms[0][3]
        for el, x, y, z in x_atoms:
            atoms.append((el, x - dx, y - dy, z - dz + r_cu_x))
        write_xyz(ROOT / f"complexes/{tag}.xyz", atoms, comment)

    make_cu_h2o5_x(
        [("Cl", 0, 0, 0)],
        "cu_h2o5_cl",
        "[Cu(H2O)5Cl]+ start | Cu-Cl ~2.2-2.3 A",
        2.28,
    )
    make_cu_h2o5_x(
        [("Br", 0, 0, 0)],
        "cu_h2o5_br",
        "[Cu(H2O)5Br]+ start | Cu-Br ~2.4 A",
        2.42,
    )
    make_cu_h2o5_x(
        [("O", 0, 0, 0), ("H", 0, 0, 0.97)],
        "cu_h2o5_oh",
        "[Cu(H2O)5(OH)]+ start",
        1.95,
    )
    make_cu_h2o5_x(
        [
            ("O", 0, 0, 0),
            ("N", 0, 0, 1.25),
            ("O", 1.10, 0, 1.85),
            ("O", -1.10, 0, 1.85),
        ],
        "cu_h2o5_no3",
        "[Cu(H2O)5(ONO2)]+ O-bound | COD https://www.crystallography.net/cod/8102245.html",
        2.05,
    )
    make_cu_h2o5_x(
        [("N", 0, 0, 0), ("C", 0, 0, 1.17), ("S", 0, 0, 2.73)],
        "cu_h2o5_ncs",
        "[Cu(H2O)5(NCS)]+ N-bound isothiocyanate (common for CuII)",
        1.95,
    )

    # [Cu(EDTA)]2- rough octahedral N2O4
    cu_edta: list[tuple] = [
        ("Cu", 0.000, 0.000, 0.000),
        ("N", -1.450, 0.200, 0.100),
        ("N", 1.450, 0.200, -0.100),
        ("C", -0.750, 1.350, 0.800),
        ("C", 0.750, 1.350, 0.800),
        ("C", -2.400, -0.900, 0.200),
        ("C", -2.100, -1.850, -0.900),
        ("O", 0.000, -1.950, 0.200),
        ("O", -3.000, -2.700, -1.100),
        ("C", -2.400, 0.600, -1.100),
        ("C", -2.000, -0.200, -2.200),
        ("O", -1.900, -0.800, -1.400),
        ("O", -2.600, 0.200, -3.200),
        ("C", 2.400, -0.900, -0.200),
        ("C", 2.100, -1.850, 0.900),
        ("O", 0.200, 0.150, 2.050),
        ("O", 3.000, -2.700, 1.100),
        ("C", 2.400, 0.600, 1.100),
        ("C", 2.000, -0.200, 2.200),
        ("O", 1.900, -0.800, 1.400),
        ("O", 2.600, 0.200, 3.200),
    ]
    for xyz in [
        (-0.75, 1.95, 1.50),
        (-0.75, 1.95, 0.10),
        (0.75, 1.95, 1.50),
        (0.75, 1.95, 0.10),
        (-3.20, -0.50, 0.60),
        (-2.00, -1.50, 0.80),
        (-3.20, 1.00, -0.70),
        (-1.90, 1.30, -1.50),
        (3.20, -0.50, -0.60),
        (2.00, -1.50, -0.80),
        (3.20, 1.00, 0.70),
        (1.90, 1.30, 1.50),
    ]:
        cu_edta.append(("H", *xyz))
    write_xyz(
        ROOT / "complexes/cu_edta.xyz",
        cu_edta,
        "[Cu(EDTA)]2- start | PubChem 23673837 https://pubchem.ncbi.nlm.nih.gov/compound/23673837 ; COD search https://www.crystallography.net/cod/result.php?text=copper+edta",
    )

    def add_assoc(extra: list[tuple], tag: str, comment: str) -> None:
        write_xyz(ROOT / f"complexes/{tag}.xyz", list(cu_edta) + extra, comment)

    add_assoc(
        [("Cl", 0.0, 2.60, 0.0)],
        "cu_edta_cl",
        "[Cu(EDTA)Cl]3- association trial (axial approach)",
    )
    add_assoc(
        [("Br", 0.0, 2.75, 0.0)],
        "cu_edta_br",
        "[Cu(EDTA)Br]3- association trial",
    )
    add_assoc(
        [("N", 0.0, 2.40, 0.0), ("C", 0.0, 3.55, 0.0), ("S", 0.0, 5.10, 0.0)],
        "cu_edta_ncs",
        "[Cu(EDTA)(NCS)]3- association trial",
    )
    add_assoc(
        [
            ("O", 0.0, 2.20, 0.0),
            ("N", 0.0, 3.40, 0.0),
            ("O", 1.10, 3.95, 0.0),
            ("O", -1.10, 3.95, 0.0),
        ],
        "cu_edta_no3",
        "[Cu(EDTA)(ONO2)]3- association trial",
    )
    add_assoc(
        [("O", 0.0, 2.10, 0.0), ("H", 0.0, 2.10, 0.97)],
        "cu_edta_oh",
        "[Cu(EDTA)(OH)]3- association trial",
    )


if __name__ == "__main__":
    main()
