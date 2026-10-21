# sync — H→μμ Run 3 synchronisation tools

What this folder contains: (1) the scripts that produce the TAMU exchange tuples (`TAMU-<dataset>.root`) that the other groups compare against, and (2) the standalone tool that compares the published tuples of any two groups (`standalone/`, no framework needed).
Older sync blocks (gpizzati v7, Pisa of September, analysis variant, reproducibility test) are not part of this branch.
Commands run from the repository root after `source env.sh` (proxy: `voms-proxy-init` in a real terminal).

## 1. Producing the tuples for the other groups

The exchange tuple has the columns of `TAMU_SYNC_TUPLE_VARIABLES.md`; the request sent to the groups is `ggfvbf_sync_20261015.yaml`. One Condor job per NanoAOD file; the variant selects the selection yaml
(`config/Run3_2024/selections_sync_*.yaml`) and the wrapper.

| variant | selection | script / wrapper |
|---|---|---|
| `ggfvbf_20261015`, `ggfvbf_pisa_20261013`, `ggfvbf_iso_20261014` | ggF/VBF exchange tuple (50 < m < 200, raw pT > 20) | `special_sync_skim.py` / `run_special_sync_skim.sh` |
| `ggfvbf_update_20261016` | **current**: update after the Hmm of 15 Oct (golden JSON + MET, trigger match, OS, ID/iso, 2 good muons, jets |eta| < 4.7, data blinded 115–135 on any mass) | `special_sync_skim_update.py` / `run_special_sync_skim_update.sh` |

```bash
source env.sh
# 1. submit (rerun the same command: only missing/failed jobs are resubmitted); MC extension: --extension
python3 sync/python/special_sync_submit.py --variant ggfvbf_update_20261016 --tag update_after_Hmm_15Oct_v3 \
    --dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg --dataset TTH_Hto2Mu --dataset ZH_Hto2Mu --dataset WplusH_Hto2Mu \
    --dataset Muon0_Run2024H --dataset Muon1_Run2024H
python3 sync/python/special_sync_submit.py --variant ggfvbf_update_20261016 --tag update_after_Hmm_15Oct_v3 --dataset VBFHto2Mu_M125_powheg --extension
# 2. when all jobs have a job_N.root and job_N_report.json: scan, hadd, check entries == sum of jobs, copy to run3_sync/TAMU/
#    (SRC defaults to Sync/update_after_Hmm_15Oct_v3; override with SRC=...). Data: both streams in one file; ext1 in a separate file.
bash sync/python/merge_sync_update_local.sh <workdir> TAMU-GluGluHto2Mu TAMU-GluGluHto2Mu
bash sync/python/merge_sync_update_local.sh <workdir> TAMU-Run2024H TAMU-Muon0_Run2024H TAMU-Muon1_Run2024H
#    merge_sync_update.sh is the same as a Condor job.
```

Before publishing: check that the merged file has no data events in 115–135 on `m_mumu_raw`, `m_mumu_nano`, `m_mumu_nano_corr`, `m_mumu_FSR_corr`, and golden JSON/MET flags = 1 (hadd silently drops trees with a different schema, so the entries check is not optional;
cast `sample_id` consistently between Muon0 and Muon1). Back up the old published files first.

Older, non-Condor way to make a loose tuple from an existing skim: `prepare_tuple_forSync.py <skim dirs> -o out.root [--mc]` (loose 50 < m < 200, every region/category flag stored; raw jet pT added by `add_raw_jet_pt.py`, reading `Jet_rawFactor` from NanoAOD).
`special_sync_merge.py --tag <tag>` merges the jobs of the first-round variants (refuses if a report is missing). `collect_corrections.py` writes the inventory of the corrections wired into the skim (imported by `prepare_tuple_forSync.py`).

## 2. Standalone tool: compare the published tuples of two groups

`compare_run3_sync_groups.py` needs only numpy, awkward, uproot and matplotlib (`pip install -r sync/requirements_sync.txt`) and read access to `/eos/cms/store/group/phys_higgs/cmshmm/run3_sync`. No framework, no ROOT.
Each group has a `treat_<group>_ntuple` function that projects its schema onto a common description (mass, ggH/VBF channel, weights, flags, common variables); the comparison only sees that description.
Groups: TAMU, Marmara, PKU, NNU, Pisa, Purdue, MIT, Caltech, Iowa. Samples: `Data_Run2024H`, `ggH`, `VBFH`, `ttH`, `ZH`, `WplusH`. Regions: Z_sideband, H_sideband, Signal_Fit (data blinded). Channels: ggH, VBF, inclusive, VBF_onlyLeading, ggH_noVBF_onlyLeading.

```bash
python3 H_mumu_pub/sync/standalone/compare_run3_sync_groups.py --list                                     # groups, samples, pairs
python3 H_mumu_pub/sync/standalone/compare_run3_sync_groups.py --group1 TAMU --group2 Pisa --sample ggH --region Z_sideband --channel inclusive --output out/
python3 H_mumu_pub/sync/standalone/compare_run3_sync_groups.py --output out/                              # every pair (add --resume to keep finished ones)
python3 H_mumu_pub/sync/standalone/compare_run3_sync_groups.py --reference TAMU --output out/             # one group against all: overlay plots with ratios
python3 H_mumu_pub/sync/standalone/compare_run3_sync_groups.py --output out/ --summaries-only             # rebuild the summaries from the summary.json files
```

Output: `out/<A>-<B>/<sample>/<region>/<channel>/` with plots and `SUMMARY.md`, plus `summary.json` per sample and a root `SUMMARY.md`. `--cache-dir` (default `/tmp/$USER/run3_sync_extracts`) holds the extracted columns; `--no-overlay` skips the all-groups plots.
Adding a group: write `treat_<Group>_ntuple` and register it in `TREATERS`.


Current results: `results/sync/sync_0_v1/considerations/` (`SUMMARY_SHORT.md`, `PAIR_BY_PAIR.md`, `TODO.md`); regenerate with `bash results/sync/sync_0_v1/scripts/run_sync_0_v1.sh` (TAMU, Marmara, Pisa) and `run_caltech.sh`.
