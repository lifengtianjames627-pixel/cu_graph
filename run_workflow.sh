#!/usr/bin/env bash
# Cu(II) d9 open-shell workflow -- SECOND-SPHERE BINDING / ACETONE SELECTIVITY
# Engine: ORCA 6.1.1 (Linux/WSL)  -- keywords follow ORCA 6.1.1 manual
#
# Goal: find a Cu(II) host whose second sphere SELECTIVELY binds acetone over other VOCs.
# Three hosts tuned by anion H-bond strength (BF4 < Cl < TFA):
#   host_cl  : Cu(H2O)2Cl2          (Cl inner-sphere)
#   host_bf4 : [Cu(H2O)6](BF4)2     (BF4- weak outer / non-coordinating)
#   host_tfa : [Cu(H2O)6](TFA)2     (TFA- H-bonding outer sphere)
#
# Model: guest G ADDS into the second sphere (anion stays). Binding energy:
#   dE_bind(G)      = E(host.G) - E(host) - E(G)                 [kcal/mol]
#   dE_bind_BBCP(G) = E(host.G) - E(host@ghost) - E(guest@ghost) [Boys-Bernardi CP]
# With USE_GCP=1, FINAL SINGLE POINT ENERGY already includes gCP (ORCA 6.1.1).
# Selectivity = compare dE_bind across guests; want acetone most negative.
#
# Layout of xyz (atom ORDER matters for ghost-CP fragments):
#   host_<h>.xyz        : Cu, inner waters, 2 anions (NO outer water)
#   guests/<g>.xyz      : isolated guest (shared across hosts)
#   complex_<h>_<g>.xyz : host_<h> atoms FIRST, then guest G atoms
#
# Usage:
#   bash run_workflow.sh prep
#   bash run_workflow.sh run
#   bash run_workflow.sh analyze
#   FORCE=1 bash run_workflow.sh run
# Set HOSTS=(tfa); extend to (tfa bf4 cl) when ready.

set -euo pipefail

# ---------------------------- USER CONFIG ------------------------------------
CHARGE_HOST=0          # Cu2+ + 2 anions = 0 for tfa/bf4/cl salts
CHARGE_GUEST=0         # neutral VOC guests
MULT=2                 # Cu(II) d9 doublet. Keep 2 for all Cu-containing systems.
FUNCTIONAL="PBE0"
DISP="D3BJ"
BASIS_HEAVY="def2-TZVP"   # Cu (and optionally heavy anions if you bump later)
BASIS_LIGHT="def2-SVP"    # H C O F Cl B N S -- screening level
GCP_LEVEL="DFT/SVP"       # ORCA 6.1.1 gCP param set matching def2-SVP (manual Table 2.50)
USE_GCP=1                 # 1 = add !GCP(${GCP_LEVEL}) to DFT energy jobs
USE_GHOST_CP=1            # 1 = also run Boys-Bernardi ghost-atom CP (ORCA 6.1.1 official)
NPROCS=8
ORCA=orca
GUESTS=(acetone water methanol ethyl_acetate toluene)
HOSTS=(tfa)               # start with TFA; later: (tfa bf4 cl)
TARGET=acetone

# ---------------------------- DERIVED ----------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p logs

GCP_KW=""
[ "$USE_GCP" = "1" ] && GCP_KW="GCP(${GCP_LEVEL})"

_now()    { date '+%Y-%m-%d %H:%M:%S'; }
_fmt_dur() { local s=$1 h m; h=$((s/3600)); s=$((s%3600)); m=$((s/60)); s=$((s%60));
             if [ "$h" -gt 0 ]; then printf '%dh %dm %ds' "$h" "$m" "$s"
             elif [ "$m" -gt 0 ]; then printf '%dm %ds' "$m" "$s"
             else printf '%ds' "$s"; fi; }

_job_ok()  { rg -qi 'ORCA TERMINATED NORMALLY|ORCA TERMINATED|FINAL SINGLE POINT ENERGY|The optimization converged' "$1"; }
is_done()  { [ -s "logs/$1.log" ] && _job_ok "logs/$1.log"; }

# write * xyz ... * block from an xyz file; optional ghost_mode:
#   none | ghost_tail (last nguest atoms become ghosts) | ghost_head (first nhost atoms become ghosts)
_append_xyz_block() {
  local out="$1" xyz="$2" chg="$3" mult="$4" mode="${5:-none}" nhost="${6:-0}" nguest="${7:-0}"
  {
    printf '* xyz %s %s\n' "$chg" "$mult"
    case "$mode" in
      ghost_tail)
        # host real, guest ghost  -> E(host in complex basis)
        awk -v nh="$nhost" 'NR>2 {
          if (NR-2 <= nh) print $0;
          else { el=$1; $1=""; printf "%s :%s\n", el, $0 }
        }' "$xyz"
        ;;
      ghost_head)
        # host ghost, guest real  -> E(guest in complex basis)
        awk -v nh="$nhost" 'NR>2 {
          if (NR-2 <= nh) { el=$1; $1=""; printf "%s :%s\n", el, $0 }
          else print $0
        }' "$xyz"
        ;;
      *)
        tail -n +3 "$xyz"
        ;;
    esac
    printf '*\n'
  } >> "$out"
}

# =============================================================================
#  Input generators  (ORCA 6.1.1)
# =============================================================================
gen_optfreq() {
  local out="$1" xyz="$2" chg="$3" mult="$4"
  cat > "${out}.inp" <<EOF
! ${FUNCTIONAL} ${DISP} ${GCP_KW} UKS TightOpt Freq TightSCF PModel
%basis
  Cu ${BASIS_HEAVY}
  H C O F Cl B N S ${BASIS_LIGHT}
  end
%scf
  Convergence Tight
  end
%geom
  Calc_Hess true
  end
%pal nprocs ${NPROCS} end
%maxcore 4000
EOF
  _append_xyz_block "${out}.inp" "$xyz" "$chg" "$mult"
}

gen_tddft() {
  local out="$1" xyz="$2" chg="$3" mult="$4"
  # Single-point TD-DFT on given geometry. Mixed basis via %basis (not on ! line).
  # DoSOC false: ORCA 6.1.1 SOC from UKS reference is not supported (manual).
  cat > "${out}.inp" <<EOF
! ${FUNCTIONAL} ${DISP} ${GCP_KW} UKS TightSCF def2/J
%basis
  Cu ${BASIS_HEAVY}
  H C O F Cl B N S ${BASIS_LIGHT}
  end
%tddft
  nroots 60
  maxdim 3
  dosoc false
  end
%scf
  Convergence Tight
  end
%pal nprocs ${NPROCS} end
%maxcore 4000
EOF
  _append_xyz_block "${out}.inp" "$xyz" "$chg" "$mult"
}

# Boys-Bernardi CP (ORCA 6.1.1 official): two ghost-atom single points at complex geometry.
# dE_BBCP = E(complex) - E(A with B ghosts) - E(B with A ghosts)
gen_ghost_cp() {
  local tag="$1" xyz="$2" nhost="$3" nguest="$4" chg_cplx="$5" mult_cplx="$6"
  # A = host fragment in complex basis (guest atoms ghosted)
  cat > "04_bbcp_${tag}_A.inp" <<EOF
! ${FUNCTIONAL} ${DISP} ${GCP_KW} UKS TightSCF PModel
%basis
  Cu ${BASIS_HEAVY}
  H C O F Cl B N S ${BASIS_LIGHT}
  end
%scf
  Convergence Tight
  end
%pal nprocs ${NPROCS} end
%maxcore 4000
EOF
  _append_xyz_block "04_bbcp_${tag}_A.inp" "$xyz" "$CHARGE_HOST" "$MULT" ghost_tail "$nhost" "$nguest"

  # B = guest fragment in complex basis (host atoms ghosted); closed-shell guest
  cat > "04_bbcp_${tag}_B.inp" <<EOF
! ${FUNCTIONAL} ${DISP} ${GCP_KW} TightSCF PModel
%basis
  Cu ${BASIS_HEAVY}
  H C O F Cl B N S ${BASIS_LIGHT}
  end
%scf
  Convergence Tight
  end
%pal nprocs ${NPROCS} end
%maxcore 4000
EOF
  _append_xyz_block "04_bbcp_${tag}_B.inp" "$xyz" "$CHARGE_GUEST" "1" ghost_head "$nhost" "$nguest"
}

# =============================================================================
#  Runner
# =============================================================================
JOB_IDX=0; JOB_OK=0; JOB_FAIL=0; JOB_SKIP=0; T_START=0

run_job() {
  local inp="$1" name; name="$(basename "$inp" .inp)"
  if [ -z "${FORCE:-}" ] && is_done "$name"; then
    JOB_SKIP=$((JOB_SKIP+1))
    printf '[%s] ... mission SKIP (already done): %s\n' "$(_now)" "$name"
    return 0
  fi
  JOB_IDX=$((JOB_IDX+1))
  printf '\n==================================================================\n'
  printf '[%s] >>> mission %d: START  %s\n' "$(_now)" "$JOB_IDX" "$name"
  printf '      input: %s   log: logs/%s.log\n' "$inp" "$name"
  printf '==================================================================\n'
  local t0 t1 dt; t0=$(date +%s)
  if ${ORCA} "${inp}" > "logs/${name}.log" 2>&1; then
    t1=$(date +%s); dt=$((t1-t0))
    if _job_ok "logs/${name}.log"; then
      JOB_OK=$((JOB_OK+1))
      printf '[%s] <<< mission %d: COMPLETED  %s | time: %s\n' "$(_now)" "$JOB_IDX" "$name" "$(_fmt_dur $dt)"
    else
      JOB_FAIL=$((JOB_FAIL+1))
      printf '[%s] !!! mission %d: FINISHED BUT LIKELY FAILED  %s | time: %s\n' "$(_now)" "$JOB_IDX" "$name" "$(_fmt_dur $dt)"
      printf '      check tail logs/%s.log for SCF/opt errors\n' "$name"
    fi
  else
    t1=$(date +%s); dt=$((t1-t0)); JOB_FAIL=$((JOB_FAIL+1))
    printf '[%s] XXX mission %d: CRASHED  %s | time: %s\n' "$(_now)" "$JOB_IDX" "$name" "$(_fmt_dur $dt)"
    printf '      ORCA exited non-zero; see logs/%s.log\n' "$name"
  fi
  printf '      totals -> run %d | ok %d | failed %d | skipped %d\n' "$JOB_IDX" "$JOB_OK" "$JOB_FAIL" "$JOB_SKIP"
}

prep_all() {
  echo "== generating ORCA 6.1.1 inputs =="
  echo "   USE_GCP=${USE_GCP} (${GCP_KW:-off})   USE_GHOST_CP=${USE_GHOST_CP}"
  for h in "${HOSTS[@]}"; do
    local hxyz="host_${h}.xyz"
    if [ -s "$hxyz" ]; then
      gen_optfreq "01_host_${h}_optfreq" "$hxyz" "$CHARGE_HOST" "$MULT"
      gen_tddft    "05_host_${h}_tddft"    "$hxyz" "$CHARGE_HOST" "$MULT"
    else
      echo "  [skip host] $hxyz missing"
    fi
    for g in "${GUESTS[@]}"; do
      local cxyz="complex_${h}_${g}.xyz"
      if [ -s "$cxyz" ]; then
        local ng nc nhost
        ng=$(head -n1 "guests/${g}.xyz"); nc=$(head -n1 "$cxyz"); nhost=$((nc - ng))
        gen_optfreq "03_complex_${h}_${g}_optfreq" "$cxyz" "$((CHARGE_HOST+CHARGE_GUEST))" "$MULT"
        gen_tddft    "05_complex_${h}_${g}_tddft"   "$cxyz" "$((CHARGE_HOST+CHARGE_GUEST))" "$MULT"
        if [ "$USE_GHOST_CP" = "1" ]; then
          gen_ghost_cp "${h}_${g}" "$cxyz" "$nhost" "$ng" "$((CHARGE_HOST+CHARGE_GUEST))" "$MULT"
        fi
      else
        echo "  [skip complex] $cxyz not provided (host atoms first, then ${g})"
      fi
    done
  done
  for g in "${GUESTS[@]}"; do
    local gxyz="guests/${g}.xyz"
    [ -s "$gxyz" ] || { echo "  [skip] $gxyz missing"; continue; }
    gen_optfreq "02_guest_${g}_optfreq" "$gxyz" "$CHARGE_GUEST" "1"
  done
}

run_all() {
  prep_all
  echo "== executing =="
  T_START=$(date +%s)
  printf '[%s] workflow START | hosts: %s\n' "$(_now)" "${HOSTS[*]}"
  for g in "${GUESTS[@]}"; do
    [ -s "02_guest_${g}_optfreq.inp" ] && run_job "02_guest_${g}_optfreq.inp"
  done
  for h in "${HOSTS[@]}"; do
    [ -s "01_host_${h}_optfreq.inp" ] && run_job "01_host_${h}_optfreq.inp"
    for g in "${GUESTS[@]}"; do
      [ -s "03_complex_${h}_${g}_optfreq.inp" ] && run_job "03_complex_${h}_${g}_optfreq.inp"
      if [ "$USE_GHOST_CP" = "1" ]; then
        [ -s "04_bbcp_${h}_${g}_A.inp" ] && run_job "04_bbcp_${h}_${g}_A.inp"
        [ -s "04_bbcp_${h}_${g}_B.inp" ] && run_job "04_bbcp_${h}_${g}_B.inp"
      fi
      [ -s "05_complex_${h}_${g}_tddft.inp" ] && run_job "05_complex_${h}_${g}_tddft.inp"
    done
    [ -s "05_host_${h}_tddft.inp" ] && run_job "05_host_${h}_tddft.inp"
  done
  local t_end total; t_end=$(date +%s); total=$((t_end-T_START))
  printf '\n==================================================================\n'
  printf '[%s] WORKFLOW FINISHED | total time: %s\n' "$(_now)" "$(_fmt_dur $total)"
  printf '      run %d | ok %d | failed %d | skipped %d\n' "$JOB_IDX" "$JOB_OK" "$JOB_FAIL" "$JOB_SKIP"
  [ "$JOB_FAIL" -gt 0 ] && printf '      !! %d failed -- fix then re-run (finished jobs auto-skipped)\n' "$JOB_FAIL"
  printf '      next: bash run_workflow.sh analyze\n'
  printf '==================================================================\n'
}

# =============================================================================
#  Analysis
# =============================================================================
get_scf_energy() { rg -oi 'FINAL SINGLE POINT ENERGY\s+-?[0-9.]+' "$1" | tail -n1 | rg -o '\-?[0-9.]+'; }
has_imag()        { rg -qi 'imaginary frequency' "$1"; }
spin_s2()         { rg -o '<S\*\*2>:\s+[0-9.]+' "$1" | tail -n1 | rg -o '[0-9.]+' || echo "NA"; }

analyze() {
  echo "================ SUMMARY (ORCA 6.1.1) ================"
  echo "USE_GCP=${USE_GCP}  USE_GHOST_CP=${USE_GHOST_CP}  TARGET=${TARGET}"
  printf "%-26s %-14s %-10s %-8s\n" "system" "E_hartree" "<S^2>" "imag?"
  printf -- "----------------------------------------------------\n"
  declare -A Eg
  for g in "${GUESTS[@]}"; do
    local f="logs/02_guest_${g}_optfreq.log"
    [ -s "$f" ] || continue
    Eg[$g]=$(get_scf_energy "$f")
    printf "%-26s %-14s %-10s %-8s\n" "guest_$g" "${Eg[$g]}" "$(spin_s2 "$f")" "$(has_imag "$f" && echo YES || echo no)"
  done

  for h in "${HOSTS[@]}"; do
    local fh="logs/01_host_${h}_optfreq.log"
    [ -s "$fh" ] || { echo "host_${h} log missing -- run first"; continue; }
    local E_host; E_host=$(get_scf_energy "$fh")
    printf "%-26s %-14s %-10s %-8s\n" "host_${h}" "$E_host" "$(spin_s2 "$fh")" "$(has_imag "$fh" && echo YES || echo no)"

    # cache target binding for DD column
    local dE_target=""
    local fct="logs/03_complex_${h}_${TARGET}_optfreq.log"
    if [ -s "$fct" ] && [ -n "${Eg[$TARGET]:-}" ]; then
      local Ect; Ect=$(get_scf_energy "$fct")
      dE_target=$(awk "BEGIN{printf \"%.6f\", ($Ect - $E_host - ${Eg[$TARGET]})*627.5095}")
    fi

    echo; echo "---- Binding energy for host = ${h} ----"
    echo "dE_bind     = E(complex) - E(host) - E(G)                         [kcal/mol]"
    echo "dE_bind_BBCP= E(complex) - E(host@guest-ghosts) - E(guest@host-ghosts)"
    echo "(more negative = binds better; DD_vs_${TARGET}>0 => ${TARGET} stronger)"
    echo
    printf "%-16s %12s %12s %12s %s\n" "guest" "dE_bind" "dE_BBCP" "DD_vs_${TARGET}" "note"
    printf -- "------------------------------------------------------------------------\n"
    for g in "${GUESTS[@]}"; do
      local fc="logs/03_complex_${h}_${g}_optfreq.log"
      [ -s "$fc" ] || { printf "%-16s %12s %12s %12s %s\n" "$g" "--" "--" "--" "no complex log"; continue; }
      [ -z "${Eg[$g]:-}" ] && continue
      local Ec dE dE_cp dd note=""
      Ec=$(get_scf_energy "$fc")
      dE=$(awk "BEGIN{printf \"%.3f\", ($Ec - $E_host - ${Eg[$g]})*627.5095}")
      local fA="logs/04_bbcp_${h}_${g}_A.log" fB="logs/04_bbcp_${h}_${g}_B.log"
      if [ -s "$fA" ] && [ -s "$fB" ]; then
        local EA EB
        EA=$(get_scf_energy "$fA"); EB=$(get_scf_energy "$fB")
        if [ -n "$EA" ] && [ -n "$EB" ]; then
          dE_cp=$(awk "BEGIN{printf \"%.3f\", ($Ec - $EA - $EB)*627.5095}")
        else dE_cp="--"; fi
      else dE_cp="--"; fi
      if [ -n "$dE_target" ]; then
        dd=$(awk "BEGIN{printf \"%+.3f\", $dE - $dE_target}")
      else dd="--"; fi
      [ "$g" = "$TARGET" ] && note="<-- target"
      printf "%-16s %12s %12s %12s %s\n" "$g" "$dE" "$dE_cp" "$dd" "$note"
    done

    echo; echo "---- d-d / LMCT roots (host_${h} vs complexes, lowest 5) ----"
    for f in "logs/05_host_${h}_tddft.log" "logs/05_complex_${h}_"*"_tddft.log"; do
      [ -s "$f" ] || continue
      echo "### $f"
      rg -A2 'STATE\s+[0-9]+:' "$f" | head -n 20
    done
  done

  echo; echo "---- JT geometry check ----"
  echo "Measure 2 longest Cu-O (axial) vs 4 shorter (equatorial) in optimized .xyz files."
}

case "${1:-run}" in
  prep)    prep_all ;;
  run|resume) run_all ;;
  analyze) analyze ;;
  *) echo "usage: bash run_workflow.sh [prep|run|analyze]   FORCE=1 to re-run finished"; exit 1 ;;
esac
