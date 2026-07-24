# Cu(II) aqua vs EDTA — ligand exchange thermodynamics (ORCA 6.1.1)

## Research questions

1. **How do different ligands / anions affect the optical properties of Cu(II) complexes?**  
   (experiment: UV–vis ± FTIR; computation: optional TD-DFT after thermo)
2. **What are the thermodynamic conditions for ligand exchange?**  
   (computation: Opt + Freq → ZPE → RRHO **Gibbs free energy** ΔG)

Kinetics (barriers) are **out of scope** for this first DFT pass; we answer feasibility from **ΔG**, not rates.

---

## Hypothesis (poster one-liner)

Water ligands on `[Cu(H₂O)₆]²⁺` are **labile**, so anions (NO₃⁻, Cl⁻, Br⁻, SCN⁻/NCS⁻, OH⁻) can enter the **inner sphere**.  
EDTA **chelates** Cu(II) as `[Cu(EDTA)]²⁻`, so the same anions are **much less likely** to displace donors or form a stable inner-sphere adduct — informing when second-sphere / outer-anion chemistry is even possible.

---

## Thermodynamic models

### A. Aqua inner-sphere exchange

```text
[Cu(H₂O)₆]²⁺ + X⁻  →  [Cu(H₂O)₅X]⁺ + H₂O

ΔG_ex = G(Cu(H₂O)₅X) + G(H₂O) − G(Cu(H₂O)₆) − G(X⁻)
```

**More negative ΔG_ex** ⇒ anion entry thermodynamically favored.

### B. EDTA association (protection test)

```text
[Cu(EDTA)]²⁻ + X⁻  →  [Cu(EDTA)X]³⁻

ΔG_assoc = G(Cu(EDTA)X) − G(Cu(EDTA)) − G(X⁻)
```

Expect **ΔG_assoc ≫ ΔG_ex** (less favorable / positive). Also inspect Cu–X distance after opt (bonded vs outer approach).

X⁻ series: **NO₃⁻ < Cl⁻ ≈ Br⁻ < NCS⁻ < OH⁻** (coordinating strength; not identical to spectrochemical Δ).

---

## Method / basis

| Setting | Value | Why |
|---------|-------|-----|
| Engine | **ORCA 6.1.1** | |
| Functional | **PBE0-D3BJ** | hybrid + dispersion; OK for open-shell Cu(II) |
| Wavefunction | **UKS**, MULT=**2** | Cu(II) d⁹ doublet |
| Basis | light **def2-SVP**; **def2-TZVP** on Cu, Br, S, Cl | Br/S need TZVP; cheap screening on H/C/N/O |
| Solvent | **CPCM(Water)** | **required** for charged ΔG in water (gas-phase ion ΔG is meaningless) |
| gCP | `GCP(DFT/SVP)` | BSSE-ish correction with small basis |
| Thermo | `TightOpt Freq` @ **298.15 K** | ZPE + RRHO → **Final Gibbs free energy** |

Landmark XYZ sources: [`refs/SOURCES.md`](refs/SOURCES.md).

Regenerate starting geometries:

```bash
python scripts/make_starting_geoms.py
```

---

## Species & charges (`* xyz q m`)

### Ligands (`ligands/`)

| File | Species | q m |
|------|---------|-----|
| `h2o.xyz` | H₂O | 0 1 |
| `cl-.xyz` / `br-.xyz` / `oh-.xyz` / `no3-.xyz` / `scn-.xyz` | anions | −1 1 |
| `edta4-.xyz` | EDTA⁴⁻ | −4 1 |

Lab EDTA is usually Na₂H₂EDTA; DFT uses **EDTA⁴⁻** as the fully deprotonated limit — state this on the poster.

### Complexes (`complexes/`)

| File | Species | q m |
|------|---------|-----|
| `cu_h2o6.xyz` | `[Cu(H₂O)₆]²⁺` | **2 2** |
| `cu_h2o5_cl.xyz` … `_br` `_oh` `_no3` `_ncs` | `[Cu(H₂O)₅X]⁺` | **1 2** |
| `cu_edta.xyz` | `[Cu(EDTA)]²⁻` | **−2 2** |
| `cu_edta_cl.xyz` … | `[Cu(EDTA)X]³⁻` association trials | **−3 2** |

NCS is modeled **N-bound** (isothiocyanate), the usual Cu(II) linkage.

---

## How to run DFT

```bash
# WSL / Linux with ORCA on PATH
bash run_workflow.sh prep      # write inputs/*.inp
bash run_workflow.sh run       # Opt+Freq all species (skip finished)
bash run_workflow.sh analyze   # print G and ΔG_ex / ΔG_assoc table
# FORCE=1 bash run_workflow.sh run
```

`analyze` pulls ORCA’s **Final Gibbs free energy** (includes ZPE via the thermo ensemble) and prints kcal/mol reaction free energies.

Optional later: TD-DFT on optimized geometries for the optical RQ.

---

## Experiment (aligned with your materials)

**Stock (as planned):** 1.0 M CuSO₄·5H₂O, NaNO₃, KBr, NaCl, EDTA, SCN⁻, OH⁻.

**Critical UV–vis note:** do **not** put undiluted 1 M Cu mixes straight into a 1 cm cuvette — saturate the detector. Dilute to ~**0.01–0.05 M Cu** (or lower) for spectra; keep **anion excess** (e.g. 5–20 equiv) when testing substitution.

### Suggested workflow

1. Record UV–vis of diluted **CuSO₄(aq)** (blank / standard).  
2. Add each anion separately → watch λ_max / color (Cl⁻/Br⁻/SCN⁻/OH⁻ stronger perturbation than NO₃⁻; OH⁻ risks Cu(OH)₂ ppt).  
3. Prepare **Cu + EDTA** (pH ~4–6 or as instructed) until spectrum matches chelate (complete complexation).  
4. Add the **same anions** to Cu–EDTA → expect **little** inner-sphere change vs aqua series.  
5. FTIR: Cu–O / Cu–N / anion bands; EDTA carboxylate shifts on chelation.  
6. Poster: abstract, each experiment, ΔG table from DFT, shortcomings, next steps.

**Safety:** SCN⁻ — never mix with strong acid (HCN risk). OH⁻ corrosive. CuSO₄ toxic/environmental. EDTA irritant.

---

## Poster checklist

- Abstract  
- Hypothesis + reaction schemes (A/B)  
- Methods (UV–vis/FTIR + ORCA PBE0-D3BJ/CPCM)  
- Results: spectra + **ΔG_ex vs ΔG_assoc** table  
- Analysis / shortcomings (no kinetics TS; EDTA protonation state; CPCM limits)  
- Possible uses: when second-sphere anion engineering is viable  
- Where next: en chelate, TD-DFT optics, explicit water / AIMD if needed  

---

## Shortcomings (be honest)

- ΔG ≠ rate (no TS / kinetics yet).  
- Association model `[Cu(EDTA)X]³⁻` is a **bound trial**; opt may push X to outer sphere — that result supports the hypothesis.  
- High formal charges + CPCM: trends more trustworthy than absolute kcal/mol.  
- Starting XYZ are sketches; validate vs COD/PubChem links in `refs/SOURCES.md` after opt.
