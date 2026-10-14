# Commands to run every comparison

All commands are run from the repository root after `source env.sh` (this sets `$HMUMU_WWW`). Details of what each step computes are in `sync/GROUPS_COMPARISON.md`; the state of the work and the findings are in `sync/GROUPS_SYNC_STATUS_20261013.md`.

```bash
source env.sh
OUT=$HMUMU_WWW/skim_v4/studies/sync/groups_comparison_13Oct     # group-vs-group output (online area)
CACHE=/tmp/$USER/run3_sync_extracts                              # extracted columns, reused between runs
```

## 1. Group against group (all pairs, all samples, all regions and channels)

```bash
# what will be run (pairs x common samples); nothing is computed
python3 sync/python/compare_run3_sync_groups.py --list

# everything: every pair, the all-groups overlays and the indices
python3 sync/python/compare_run3_sync_groups.py --output $OUT --cache-dir $CACHE
```

The first read of a large file takes several minutes (EOS I/O); later runs use the cache. Groups: TAMU, Marmara, PKU, NNU, Pisa, Purdue. A pair only uses the samples both groups provide, and Pisa has no ggH category (those cells are empty).

### One pair, one sample, one cell

```bash
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 PKU --output $OUT
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 PKU --sample Data_Run2024H --output $OUT
python3 sync/python/compare_run3_sync_groups.py --group1 Marmara --group2 PKU --sample ggH \
    --region Signal_Fit --channel VBF --output $OUT
```

### Faster: pairs in parallel

Different pairs write to different folders, so they can run as separate processes. Warm the cache first, one process per (group, sample), so that two processes never write the same cache file:

```bash
python3 - <<'EOF'
import sys; sys.path.insert(0, "sync/python")
from pathlib import Path
import compare_run3_sync_groups as c
c.extract("PKU", "VBFH", Path("/tmp/$USER/run3_sync_extracts"))     # repeat for every group and sample, one process each
EOF

for pair in "TAMU PKU" "TAMU Marmara" "TAMU NNU" "TAMU Pisa" "TAMU Purdue" "Purdue PKU" "Purdue Pisa" "Purdue NNU" "Purdue Marmara" "Marmara PKU" "Marmara NNU" "Marmara Pisa" "PKU NNU" "PKU Pisa" "NNU Pisa"; do
  set -- $pair
  python3 sync/python/compare_run3_sync_groups.py --group1 $1 --group2 $2 --output $OUT --cache-dir $CACHE > log_$1_$2.txt 2>&1 &
done; wait

# finally the all-groups overlays and the indices (pairs that already have a summary.json are skipped)
python3 sync/python/compare_run3_sync_groups.py --resume --output $OUT --cache-dir $CACHE
```

### Rewrite only the markdown summaries and the index

```bash
python3 sync/python/compare_run3_sync_groups.py --summaries-only --output $OUT
```

### Output

`$OUT/<GROUP1>-<GROUP2>/SUMMARY.md` (quantitative tables and qualitative reading of the pair), `<sample>/<region>/<channel>/SUMMARY.md` with plots, and `$OUT/ALL_GROUPS/<sample>/...` with all groups on one plot. Use `--resume` to skip pair/sample folders that already have a `summary.json`; do not use it after changing a `treat_*` function or a variable expression, because the old plots would be kept.

## 2. New TAMU tuples (variant `ggfvbf_pisa_20261013`)

Jets stay at 35/25 GeV (TAMU decision 2026-10-14); Pisa's exchange tuples use 30/30 GeV, so jet differences with Pisa are expected.

```bash
OUTT=/eos/user/v/vdamante/H_mumu/Sync/ggfvbf_pisa_20261013

# submit (data)
python3 sync/python/special_sync_submit.py --variant ggfvbf_pisa_20261013 --tag ggfvbf_pisa_20261013 \
    --group-name TAMU --output-root $OUTT --flavour tomorrow --max-jobs 500
# submit (MC): NO --mc-lumi-split (the split keeps only 44% of the events, Pisa keeps a different random 50%);
# output in a separate root so the split version is not overwritten
OUTN=/eos/user/v/vdamante/H_mumu/Sync/ggfvbf_pisa_20261013_nosplit
python3 sync/python/special_sync_submit.py --variant ggfvbf_pisa_20261013 --tag ggfvbf_pisa_20261013 \
    --group-name TAMU --output-root $OUTN --flavour tomorrow \
    --dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg --dataset TTH_Hto2Mu --dataset WplusH_Hto2Mu --dataset ZH_Hto2Mu

# when the jobs are done: merge into TAMU-<dataset>/sync_tuple.root
python3 sync/python/special_sync_merge.py --tag ggfvbf_pisa_20261013 --group-name TAMU --output-root $OUTT
```

**Always check the merged entries** against the sum over the jobs of `inputs.txt`: `hadd` exits 0 but silently drops whole files when
the branch set or a leaf type differs between jobs (`sample_id` Int_t vs Long_t, `OrthogonalEraTag` present only in some jobs). This
happened for Muon0+Muon1, ttH, WplusH and VH; they were rebuilt with an RDataFrame Snapshot per (branch, type) signature, common
columns, `sample_id` cast to Long_t, then `hadd`. The combined `TAMU-VH` merge must be redone the same way. The queue is capped at
`HMUMU_MAX_QUEUED_JOBS=3000` jobs.

## 3. TAMU against any group on the unpublished tuples (not online)

`HMUMU_TAMU_STAGING` makes the TAMU treater read `<dir>/TAMU-<dataset>/sync_tuple.root` (datasets Run2024H, GluGluHto2Mu, VBFHto2Mu_M125_powheg, TTH_Hto2Mu, ZH_Hto2Mu, WplusH_Hto2Mu) instead of `run3_sync/TAMU/`. Use a separate cache, and write outside `$OUT` (the Pisa redo goes to `sync/results_newTamu_Pisa`):

```bash
export HMUMU_TAMU_STAGING=/eos/user/v/vdamante/H_mumu/Sync/ggfvbf_pisa_20261013     # or the _nosplit root for the MC
C=/tmp/$USER/run3_sync_extracts_newTamu
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 Pisa --output sync/results_newTamu_Pisa --cache-dir $C --no-overlay
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 PKU  --output /tmp/$USER/groups_newTamu --cache-dir $C --no-overlay
```

(The older dedicated diagnostic against Filippo's tuples, `compare_filippo_diagnostic.py`, compares with his diagnostic ntuples, which use 35/25 GeV jets, not with the exchange tuples.)

## 4. Publishing the new TAMU files to the exchange area

The new files replace the old ones in `run3_sync/TAMU/` (the old ones are not kept there). Do it only after the validation against Pisa and after checking the entries of every file; keep a private copy of the old files outside the shared area.

```bash
BASE=/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/TAMU
mkdir -p /eos/user/$USER/H_mumu/Sync/old_exchange_ggfvbf_20261015 && cp -n $BASE/TAMU-*.root /eos/user/$USER/H_mumu/Sync/old_exchange_ggfvbf_20261015/
# then overwrite with the merged files under the exchange names, e.g.
cp $OUTT/TAMU-Muon0_Run2024H/sync_tuple.root ...   # exchange names: TAMU-Run2024H.root (Muon0+Muon1), TAMU-GluGluHto2Mu.root, ...
```

After the replacement, rerun section 1: the TAMU cache is invalidated automatically when the file changes, and the online plots are **overwritten** (the pair folders are rewritten, not kept next to the old ones).

## 5. Adding a group

See "How groups are plugged in" in `sync/GROUPS_COMPARISON.md`: write `treat_<group>_ntuple()`, register it in `TREATERS`, add the `<group>=` expression to the variables in `VARIABLES`, check with `--list`.
