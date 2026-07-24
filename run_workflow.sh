#!/usr/bin/env bash
# Cu(II) ligand-exchange thermodynamics -- ORCA 6.1.1
#
# Hypothesis:
#   [Cu(H2O)6]2+  : labile waters -> anions X- ENTER inner sphere (DeltaG_ex < 0)
#   [Cu(EDTA)]2-  : chelate lock   -> X- association UNFAVORABLE / outer only
#
# Reactions (aqueous CPCM):
#   Aqua exchange:
#     [Cu(H2O)6]2+ + X-  ->  [Cu(H2O)5X]+ + H2O
#     dG_ex = G(CuH2O5X) + G(H2O) - G(CuH2O6) - G(X)
#
#   EDTA association (does X attach?):
#     [Cu(EDTA)]2- + X-  ->  [Cu(EDTA)X]3-
#     dG_assoc = G(CuEDTA_X) - G(CuEDTA) - G(X)
#
# Workflow: Opt + Freq -> ZPE + RRHO thermo -> Gibbs G(298 K)
# Optical RQ (optional later): TD-DFT on optimized Cu species.
#
# Usage:
#   bash run_workflow.sh prep
#   bash run_workflow.sh run
#   bash run_workflow.sh analyze
#   FORCE=1 bash run_workflow.sh run

set -euo pipefail

# ---------------------------- USER CONFIG ------------------------------------
FUNCTIONAL="PBE0"
DISP="D3BJ"
BASIS_LIGHT="def2-SVP"
BASIS_HEAVY="def2-TZVP"   # Cu, Br, S (and Cl optional heavy)
GCP_LEVEL="DFT/SVP"
USE_GCP=1
USE_CPCM=1                # REQUIRED for charged aqueous DeltaG (set 0 only for debug)
SOLVENT="Water"
NPROCS=8
ORCA=orca
TEMP_K=298.15

# Species tags <-> xyz + charge + multiplicity
# Free ligands (closed shell unless noted)
LIGANDS=(h2o cl- br- oh- no3- scn- edta4-)
declare -A LIG_XYZ LIG_CHG LIG_MULT
LIG_XYZ[h2o]="ligands/h2o.xyz";     LIG_CHG[h2o]=0;   LIG_MULT[h2o]=1
LIG_XYZ[cl-]="ligands/cl-.xyz";     LIG_CHG[cl-]=-1;  LIG_MULT[cl-]=1
LIG_XYZ[br-]="ligands/br-.xyz";     LIG_CHG[br-]=-1;  LIG_MULT[br-]=1
LIG_XYZ[oh-]="ligands/oh-.xyz";     LIG_CHG[oh-]=-1;  LIG_MULT[oh-]=1
LIG_XYZ[no3-]="ligands/no3-.xyz";   LIG_CHG[no3-]=-1; LIG_MULT[no3-]=1
LIG_XYZ[scn-]="ligands/scn-.xyz";   LIG_CHG[scn-]=-1; LIG_MULT[scn-]=1
LIG_XYZ[edta4-]="ligands/edta4-.xyz"; LIG_CHG[edta4-]=-4; LIG_MULT[edta4-]=1

# Cu complexes (all Cu(II) doublet)
COMPLEXES=(cu_h2o6 cu_h2o5_cl cu_h2o5_br cu_h2o5_oh cu_h2o5_no3 cu_h2o5_ncs \
           cu_edta cu_edta_cl cu_edta_br cu_edta_oh cu_edta_no3 cu_edta_ncs)
declare -A CPLX_XYZ CPLX_CHG CPLX_MULT
for s in "${COMPLEXES[@]}"; do
  CPLX_XYZ[$s]="complexes/${s}.xyz"
  CPLX_MULT[$s]=2
done
CPLX_CHG[cu_h2o6]=2
CPLX_CHG[cu_h2o5_cl]=1
CPLX_CHG[cu_h2o5_br]=1
CPLX_CHG[cu_h2o5_oh]=1
CPLX_CHG[cu_h2o5_no3]=1
CPLX_CHG[cu_h2o5_ncs]=1
CPLX_CHG[cu_edta]=-2
CPLX_CHG[cu_edta_cl]=-3
CPLX_CHG[cu_edta_br]=-3
CPLX_CHG[cu_edta_oh]=-3
CPLX_CHG[cu_edta_no3]=-3
CPLX_CHG[cu_edta_ncs]=-3

# Map anion tag -> aqua product / edta product / free ligand key
ANIONS=(cl- br- oh- no3- ncs)
declare -A ANION_LIG AQUA_PROD EDTA_PROD
ANION_LIG[cl-]="cl-";   AQUA_PROD[cl-]="cu_h2o5_cl";   EDTA_PROD[cl-]="cu_edta_cl"
ANION_LIG[br-]="br-";   AQUA_PROD[br-]="cu_h2o5_br";   EDTA_PROD[br-]="cu_edta_br"
ANION_LIG[oh-]="oh-";   AQUA_PROD[oh-]="cu_h2o5_oh";   EDTA_PROD[oh-]="cu_edta_oh"
ANION_LIG[no3-]="no3-"; AQUA_PROD[no3-]="cu_h2o5_no3"; EDTA_PROD[no3-]="cu_edta_no3"
ANION_LIG[ncs]="scn-";  AQUA_PROD[ncs]="cu_h2o5_ncs";  EDTA_PROD[ncs]="cu_edta_ncs"

# ---------------------------- DERIVED ----------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p logs inputs

GCP_KW=""
[ "$USE_GCP" = "1" ] && GCP_KW="GCP(${GCP_LEVEL})"
CPCM_KW=""
[ "$USE_CPCM" = "1" ] && CPCM_KW="CPCM(${SOLVENT})"

_now()    { date '+%Y-%m-%d %H:%M:%S'; }
_fmt_dur() { local s=$1 h m; h=$((s/3600)); s=$((s%3600)); m=$((s/60)); s=$((s%60));
             if [ "$h" -gt 0 ]; then printf '%dh %dm %ds' "$h" "$m" "$s"
             elif [ "$m" -gt 0 ]; then printf '%dm %ds' "$m" "$s"
             else printf '%ds' "$s"; fi; }

_job_ok()  { grep -qiE 'ORCA TERMINATED NORMALLY|ORCA TERMINATED|FINAL SINGLE POINT ENERGY|The optimization converged' "$1"; }
is_done()  { [ -s "logs/$1.log" ] && _job_ok "logs/$1.log"; }

_append_xyz_block() {
  local out="$1" xyz="$2" chg="$3" mult="$4"
  {
    printf '* xyz %s %s\n' "$chg" "$mult"
    tail -n +3 "$xyz"
    printf '*\n'
  } >> "$out"
}

# Heavy-element NewGTO lines depending on file contents
_has_elem() { grep -qE "^[[:space:]]*$1([[:space:]]|$)" "$2"; }

_heavy_basis_block() {
  local xyz="$1"
  echo '  Basis "'"${BASIS_LIGHT}"'"'
  if _has_elem Cu "$xyz"; then
    echo '  NewGTO Cu "'"${BASIS_HEAVY}"'" end'
  fi
  if _has_elem Br "$xyz"; then
    echo '  NewGTO Br "'"${BASIS_HEAVY}"'" end'
  fi
  if _has_elem S "$xyz"; then
    echo '  NewGTO S "'"${BASIS_HEAVY}"'" end'
  fi
  if _has_elem Cl "$xyz"; then
    echo '  NewGTO Cl "'"${BASIS_HEAVY}"'" end'
  fi
}

# Opt + analytic Freq -> ZPE + Gibbs (RRHO) in ORCA thermo block
gen_optfreq() {
  local out="$1" xyz="$2" chg="$3" mult="$4" uks_kw=""
  [ "$mult" != "1" ] && uks_kw="UKS"
  local basis_tmp; basis_tmp="$(_heavy_basis_block "$xyz")"
  cat > "inputs/${out}.inp" <<EOF
! ${FUNCTIONAL} ${DISP} ${GCP_KW} ${CPCM_KW} ${uks_kw} TightOpt Freq TightSCF PModel
%basis
${basis_tmp}
end
%scf
  Convergence Tight
  end
%geom
  Calc_Hess true
  end
%freq
  Temp ${TEMP_K}
  end
%pal nprocs ${NPROCS} end
%maxcore 4000
EOF
  _append_xyz_block "inputs/${out}.inp" "$xyz" "$chg" "$mult"
}

gen_tddft_opt() {
  # Optional optical single point on OPTIMIZED geometry (after run)
  local out="$1" xyz="$2" chg="$3" mult="$4"
  local basis_tmp; basis_tmp="$(_heavy_basis_block "$xyz")"
  cat > "inputs/${out}.inp" <<EOF
! ${FUNCTIONAL} ${DISP} ${GCP_KW} ${CPCM_KW} UKS TightSCF def2/J
%basis
${basis_tmp}
end
%tddft
  nroots 40
  maxdim 3
  dosoc false
  end
%scf
  Convergence Tight
  end
%pal nprocs ${NPROCS} end
%maxcore 4000
EOF
  _append_xyz_block "inputs/${out}.inp" "$xyz" "$chg" "$mult"
}

# =============================================================================
JOB_IDX=0; JOB_OK=0; JOB_FAIL=0; JOB_SKIP=0; T_START=0

run_job() {
  local inp="$1" name; name="$(basename "$inp" .inp)"
  if [ -z "${FORCE:-}" ] && is_done "$name"; then
    JOB_SKIP=$((JOB_SKIP+1))
    printf '[%s] ... SKIP: %s\n' "$(_now)" "$name"
    return 0
  fi
  JOB_IDX=$((JOB_IDX+1))
  printf '\n[%s] >>> START %s\n' "$(_now)" "$name"
  local t0 t1 dt; t0=$(date +%s)
  if ${ORCA} "${inp}" > "logs/${name}.log" 2>&1; then
    t1=$(date +%s); dt=$((t1-t0))
    if _job_ok "logs/${name}.log"; then
      JOB_OK=$((JOB_OK+1))
      printf '[%s] <<< OK %s | %s\n' "$(_now)" "$name" "$(_fmt_dur $dt)"
    else
      JOB_FAIL=$((JOB_FAIL+1))
      printf '[%s] !!! FAIL %s | %s\n' "$(_now)" "$name" "$(_fmt_dur $dt)"
    fi
  else
    t1=$(date +%s); dt=$((t1-t0)); JOB_FAIL=$((JOB_FAIL+1))
    printf '[%s] XXX CRASH %s | %s\n' "$(_now)" "$name" "$(_fmt_dur $dt)"
  fi
}

prep_all() {
  echo "== generating Opt+Freq (ZPE/Gibbs) inputs =="
  echo "   CPCM=${USE_CPCM} (${SOLVENT})  GCP=${USE_GCP}  T=${TEMP_K} K"
  for s in "${LIGANDS[@]}"; do
    local xyz="${LIG_XYZ[$s]}"
    [ -s "$xyz" ] || { echo "  missing $xyz"; continue; }
    gen_optfreq "01_lig_${s}_optfreq" "$xyz" "${LIG_CHG[$s]}" "${LIG_MULT[$s]}"
  done
  for s in "${COMPLEXES[@]}"; do
    local xyz="${CPLX_XYZ[$s]}"
    [ -s "$xyz" ] || { echo "  missing $xyz"; continue; }
    gen_optfreq "02_cplx_${s}_optfreq" "$xyz" "${CPLX_CHG[$s]}" "${CPLX_MULT[$s]}"
  done
  echo "inputs in inputs/*.inp"
}

run_all() {
  prep_all
  echo "== executing =="
  T_START=$(date +%s)
  for s in "${LIGANDS[@]}"; do
    [ -s "inputs/01_lig_${s}_optfreq.inp" ] && run_job "inputs/01_lig_${s}_optfreq.inp"
  done
  for s in "${COMPLEXES[@]}"; do
    [ -s "inputs/02_cplx_${s}_optfreq.inp" ] && run_job "inputs/02_cplx_${s}_optfreq.inp"
  done
  local t_end total; t_end=$(date +%s); total=$((t_end-T_START))
  printf '\n[%s] DONE | %s | ok %d fail %d skip %d\n' \
    "$(_now)" "$(_fmt_dur $total)" "$JOB_OK" "$JOB_FAIL" "$JOB_SKIP"
  printf 'next: bash run_workflow.sh analyze\n'
}

# --- parsers ---
get_scf_energy() {
  grep -oiE 'FINAL SINGLE POINT ENERGY[[:space:]]+-?[0-9.]+' "$1" | tail -n1 | grep -oE '\-?[0-9.]+'
}
get_gibbs() {
  # ORCA: "Final Gibbs free energy         -xxx.xxxxxx Eh"
  local g
  g=$(grep -i 'Final Gibbs free energy' "$1" | tail -n1 | grep -oE '\-?[0-9]+\.[0-9]+' | tail -n1 || true)
  if [ -z "$g" ]; then
    g=$(grep -iE 'Gibbs free energy[[:space:]]+-?[0-9.]+' "$1" | tail -n1 | grep -oE '\-?[0-9]+\.[0-9]+' | tail -n1 || true)
  fi
  echo "$g"
}
get_zpe() {
  grep -i 'Zero point energy' "$1" | tail -n1 | grep -oE '\-?[0-9]+\.[0-9]+' | head -n1 || true
}
has_imag() { grep -qiE 'imaginary mode|imaginary frequency' "$1"; }
spin_s2()  { grep -oE '<S\*\*2>:[[:space:]]+[0-9.]+' "$1" | tail -n1 | grep -oE '[0-9.]+' || echo "NA"; }

analyze() {
  echo "================ Cu ligand-exchange thermo (ORCA) ================"
  echo "CPCM=${USE_CPCM}  GCP=${USE_GCP}  T=${TEMP_K} K"
  echo
  printf "%-22s %16s %16s %10s %8s\n" "species" "E_SCF/Eh" "G/Eh" "<S2>" "imag?"
  printf -- "--------------------------------------------------------------------\n"

  declare -A G E
  for s in "${LIGANDS[@]}"; do
    local f="logs/01_lig_${s}_optfreq.log"
    [ -s "$f" ] || continue
    E[$s]=$(get_scf_energy "$f"); G[$s]=$(get_gibbs "$f")
    printf "%-22s %16s %16s %10s %8s\n" "lig_$s" "${E[$s]:-}" "${G[$s]:-}" "$(spin_s2 "$f")" "$(has_imag "$f" && echo YES || echo no)"
  done
  for s in "${COMPLEXES[@]}"; do
    local f="logs/02_cplx_${s}_optfreq.log"
    [ -s "$f" ] || continue
    E[$s]=$(get_scf_energy "$f"); G[$s]=$(get_gibbs "$f")
    printf "%-22s %16s %16s %10s %8s\n" "cplx_$s" "${E[$s]:-}" "${G[$s]:-}" "$(spin_s2 "$f")" "$(has_imag "$f" && echo YES || echo no)"
  done

  echo
  echo "---- Aqua inner-sphere exchange ----"
  echo "  [Cu(H2O)6]2+ + X- -> [Cu(H2O)5X]+ + H2O"
  echo "  dG_ex (kcal/mol) = [G(prod)+G(H2O)-G(cu_h2o6)-G(X)] * 627.5095"
  echo "  more NEGATIVE => anion entry thermodynamically favored"
  echo
  printf "%-10s %12s %s\n" "X" "dG_ex" "note"
  printf -- "----------------------------------------\n"
  for a in "${ANIONS[@]}"; do
    local lig="${ANION_LIG[$a]}" prod="${AQUA_PROD[$a]}"
    if [ -z "${G[cu_h2o6]:-}" ] || [ -z "${G[$lig]:-}" ] || [ -z "${G[$prod]:-}" ] || [ -z "${G[h2o]:-}" ]; then
      printf "%-10s %12s %s\n" "$a" "--" "missing G"
      continue
    fi
    local dg
    dg=$(awk -v gp="${G[$prod]}" -v gw="${G[h2o]}" -v gc="${G[cu_h2o6]}" -v gx="${G[$lig]}" \
      'BEGIN{printf "%.2f", (gp+gw-gc-gx)*627.5095}')
    local note=""
    awk -v d="$dg" 'BEGIN{exit !(d<0)}' && note="favored (inner entry)"
    awk -v d="$dg" 'BEGIN{exit !(d>=0)}' && note="uphill"
    printf "%-10s %12s %s\n" "$a" "$dg" "$note"
  done

  echo
  echo "---- EDTA association (chelate protection) ----"
  echo "  [Cu(EDTA)]2- + X- -> [Cu(EDTA)X]3-"
  echo "  dG_assoc (kcal/mol) = [G(CuEDTA_X)-G(CuEDTA)-G(X)] * 627.5095"
  echo "  POSITIVE / much less negative than aqua dG_ex => anions stay out"
  echo
  printf "%-10s %12s %s\n" "X" "dG_assoc" "note"
  printf -- "----------------------------------------\n"
  for a in "${ANIONS[@]}"; do
    local lig="${ANION_LIG[$a]}" prod="${EDTA_PROD[$a]}"
    if [ -z "${G[cu_edta]:-}" ] || [ -z "${G[$lig]:-}" ] || [ -z "${G[$prod]:-}" ]; then
      printf "%-10s %12s %s\n" "$a" "--" "missing G"
      continue
    fi
    local dg
    dg=$(awk -v gp="${G[$prod]}" -v gc="${G[cu_edta]}" -v gx="${G[$lig]}" \
      'BEGIN{printf "%.2f", (gp-gc-gx)*627.5095}')
    printf "%-10s %12s\n" "$a" "$dg"
  done

  echo
  echo "Also check optimized Cu-X distances: short covalent vs long outer-sphere."
  echo "ZPE lines: grep -i 'Zero point energy' logs/*.log"
  echo "Refs: refs/SOURCES.md"
}

case "${1:-run}" in
  prep)    prep_all ;;
  run|resume) run_all ;;
  analyze) analyze ;;
  *) echo "usage: bash run_workflow.sh [prep|run|analyze]  FORCE=1 to re-run"; exit 1 ;;
esac
