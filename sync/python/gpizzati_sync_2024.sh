#!/usr/bin/env bash
# Full sync TAMU skim_v4 vs gpizzati hmm_skim v7, Run3_2024: every data dataset, and the main MC
# with all the analysis weights (compare_gpizzati_v7.py --weight analysis). One HTCondor job per
# dataset; then `summary` collects the per-dataset READMEs in one table.
#
#   bash sync/python/gpizzati_sync_2024.sh submit     # submit the missing datasets
#   bash sync/python/gpizzati_sync_2024.sh local DS   # run one dataset here
#   bash sync/python/gpizzati_sync_2024.sh summary    # SUMMARY.md once the jobs are done
set -u
REPO=/afs/cern.ch/work/v/vdamante/Hmumu_analysis
OUT=${SYNC_OUT:-/eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7_allweights}
PLOTS=${SYNC_PLOTS:-/eos/user/v/vdamante/www/H_mumu/skim_v4/studies/sync/gpizzati_v7_allweights}
DATA="Muon0_Run2024C Muon0_Run2024D Muon0_Run2024E Muon0_Run2024F Muon0_Run2024G Muon0_Run2024H Muon0_Run2024I_v1 Muon0_Run2024I_v2
      Muon1_Run2024C Muon1_Run2024D Muon1_Run2024E Muon1_Run2024F Muon1_Run2024G Muon1_Run2024H Muon1_Run2024I_v1 Muon1_Run2024I_v2"
MC="GluGluHto2Mu VBFHto2Mu_M125_powheg VBFHto2Mu_M125_amcatnlo DYto2Mu_MLL_105to160_amcatnloFXFX
    DYto2Mu_MLL_105to160_amcatnloFXFX_Fil_VBF EWK_2Mu2J_MLL_105to160_herwig DYto2Mu_M_50_amcatnloFXFX TTto2L2Nu"
# MC skims of several hundred GB: a subset of files (same lumisections on both sides are compared)
declare -A MAXFILES=([DYto2Mu_M_50_amcatnloFXFX]=60 [TTto2L2Nu]=60 [DYto2Mu_MLL_105to160_amcatnloFXFX]=80)

run_one() {
  local ds=$1 extra=""
  [[ $ds == Muon* ]] || extra="--weight analysis --threads 8"
  [ -n "${MAXFILES[$ds]:-}" ] && extra="$extra --max-files ${MAXFILES[$ds]}"
  cd $REPO && set -- && source env.sh >/dev/null 2>&1
  python3 sync/python/compare_gpizzati_v7.py --dataset $ds $extra --plots --plots-dir $PLOTS/$ds --output-dir $OUT
}

case "${1:-}" in
  local) run_one "$2" ;;
  submit)
    L=$REPO/htcondor/sync_gpizzati_2024; mkdir -p $L
    : > $L/datasets.txt
    running=$(condor_q -af Args 2>/dev/null)
    for ds in $DATA $MC; do [ -f $OUT/$ds/README.md ] || grep -qx "$ds" <<< "$running" || echo $ds >> $L/datasets.txt; done
    cat > $L/job.sh <<EOF
#!/bin/bash
bash $REPO/sync/python/gpizzati_sync_2024.sh local "\$1"
EOF
    chmod +x $L/job.sh
    cat > $L/submit.sub <<EOF
executable     = $L/job.sh
arguments      = \$(ds)
output         = $L/\$(ds).out
error          = $L/\$(ds).err
log            = $L/condor.log
request_cpus   = 8
request_memory = 40000
+JobFlavour    = "tomorrow"
getenv         = False
queue ds from $L/datasets.txt
EOF
    echo "$(wc -l < $L/datasets.txt) datasets"; condor_submit $L/submit.sub ;;
  summary)
    python3 - "$OUT" <<'EOF'
import sys, csv
from pathlib import Path
out = Path(sys.argv[1])
rows = ["# Sync TAMU skim_v4 vs gpizzati v7, Run3_2024 (data: all eras; MC: all analysis weights)", "",
        "Yields on the events present in both skims, TAMU / gpizzati. Details: `<dataset>/README.md`.", "",
        "| dataset | baseline Z | baseline H | baseline SR | VBF Z | VBF H | VBF SR | weight factors (baseline Z, ratio of means) |",
        "|---|---|---|---|---|---|---|---|"]
for d in sorted(p for p in out.iterdir() if (p / "yields.csv").exists()):
    y = {(r["category"], r["region"]): r for r in csv.DictReader(open(d / "yields.csv"))}
    cell = lambda c, r: (f"{float(y[c, r]['TAMU_shared']) / float(y[c, r]['gpizzati_shared']):.4f}"
                         if (c, r) in y and float(y[c, r]["gpizzati_shared"]) else "-")
    f = ""
    if (d / "weight_factors.csv").exists():
        f = ", ".join(f"{r['factor'][2:] if r['factor'].startswith('w_') else r['factor']} {float(r['ratio_of_means']):.3f}"
                      for r in csv.DictReader(open(d / "weight_factors.csv"))
                      if r["category"] == "baseline" and r["region"] == "Z_sideband")
    rows.append(f"| {d.name} | " + " | ".join(cell(c, r) for c in ("baseline", "VBF")
                                              for r in ("Z_sideband", "H_sideband", "Signal_Fit")) + f" | {f} |")
# all 2024 data together: sums of the per-dataset yields
data = [d for d in sorted(out.iterdir()) if d.name.startswith("Muon") and (d / "yields.csv").exists()]
total = {}
for d in data:
    for r in csv.DictReader(open(d / "yields.csv")):
        t = total.setdefault((r["category"], r["region"]), [0.0, 0.0])
        t[0] += float(r["TAMU_shared"]); t[1] += float(r["gpizzati_shared"])
cell = lambda c, r: f"{total[c, r][0] / total[c, r][1]:.5f}" if (c, r) in total and total[c, r][1] else "-"
rows[6:6] = [f"| **all 2024 data ({len(data)}/16 datasets)** | " + " | ".join(
    cell(c, r) for c in ("baseline", "VBF") for r in ("Z_sideband", "H_sideband", "Signal_Fit")) + " | |"]
rows += ["", "All 2024 data, events: TAMU / gpizzati / ratio", "", "| category | region | TAMU | gpizzati | ratio |", "|---|---|---|---|---|"]
rows += [f"| {c} | {r} | {t[0]:,.0f} | {t[1]:,.0f} | {t[0] / t[1]:.5f} |" for (c, r), t in total.items() if t[1]]
(out / "SUMMARY.md").write_text("\n".join(rows) + "\n")
print("\n".join(rows))
EOF
    # the markdown and CSV summaries also go next to the plots, on the web area
    cp $OUT/SUMMARY.md $PLOTS/; [ -f $OUT/DIAGNOSTICS.md ] && cp $OUT/DIAGNOSTICS.md $PLOTS/
    for d in $OUT/*/; do
      [ -f $d/README.md ] && mkdir -p $PLOTS/$(basename $d) && cp $d/README.md $d/*.csv $PLOTS/$(basename $d)/
    done
    ;;
  *) sed -n 2,9p "$0" ;;
esac
