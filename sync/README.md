# sync — H→μμ Run 3 synchronisation tools

Everything needed to compare our skim (TAMU, `skim_v4`) with other groups, and to produce the tuples the other groups compare against.
All commands run from the repository root after `source env.sh` (proxy: `voms-proxy-init` in a real terminal). `./hmumu sync <tool> ...` is a shortcut for `python3 sync/python/<tool>.py ...`.

## Who does what

| Block | Purpose | Scripts (`sync/python/`) |
|---|---|---|
| **1. gpizzati (hmm_skim v7)** | event-by-event sync of our skim with Giorgio's, Run3_2024 data and MC, all analysis weights | `compare_gpizzati_v7.py`, `gpizzati_sync_2024.sh`, `tamu_analysis_weights.py` |
| **2. Pisa, end of September** | VBF Z-sideband/H-sideband sync with Filippo's diagnostic ntuples (cutflow, per-event, MC weights) | `compare_cutflow_pisa.py`, `compare_filippo_diagnostic.py`, `plot_filippo_comparison.py`, `compare_filippo_weights.py`, `sync_mc_weights.py`, `tuple_weights_friend.py` |
| **3. Tuples for the other groups** | NanoAOD → exchange tuple (`TAMU-<dataset>.root`) published in `run3_sync/TAMU/` | `special_sync_skim.py`, `special_sync_skim_analysis.py`, `special_sync_skim_update.py`, `special_sync_submit.py`, `special_sync_merge.py`, `run_special_sync_skim*.sh`, `merge_sync_update.sh`, `merge_sync_update_local.sh`, `prepare_tuple_forSync.py`, `add_raw_jet_pt.py`, `collect_corrections.py` |
| **4. Standalone group-vs-group tool** | compares the published tuples of any two groups, no framework needed | `compare_run3_sync_groups.py`, `../requirements_sync.txt` (the standalone tool in `H_mumu_pub/sync/standalone/` is the one to run) |

Other files: `ggfvbf_sync_20261015.yaml` (description of the exchange tuple, sent to the groups), `TAMU_SYNC_TUPLE_VARIABLES.md` (every column of the TAMU tuple and how it is filled),
Results and status of all rounds: `results/sync/` (`SYNC_LOG.md`, `sync_0_v0` = before 15 Oct, `sync_0_v1` = after 15 Oct).
`sync_september23_Pisa.pdf` (slides of the Pisa sync).

## 1. Sync with gpizzati hmm_skim v7 (Run3_2024)

`compare_gpizzati_v7.py` applies the selection of `config/Run3_2024/selections.yaml` (OS, trigger, muon pT 26/20, ML, TT below 30 GeV, b-veto, `VBF_def`; categories baseline/ggF/VBF,
regions Z_sideband/H_sideband/Signal_Fit) to our skim and to his v7 skim, each with its own objects, and compares them event by event on the lumisections present in both.
MC events are also matched on the raw eta of the two muons (extensions reuse run:lumi:event).

```bash
# one dataset, tables only
python3 sync/python/compare_gpizzati_v7.py --dataset Muon0_Run2024C --output-dir $OUT
# MC with plots; several datasets in one call (data in Signal_Fit needs --unblind)
python3 sync/python/compare_gpizzati_v7.py --dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg --weight genWeight --plots --output-dir $OUT
# explicit inputs (one --dataset only)
python3 sync/python/compare_gpizzati_v7.py --dataset DYto2Mu_MLL_105to160_amcatnloFXFX --ours-input <folder> --theirs-input '<glob>' --output-dir $OUT --plots
```

Per dataset: `README.md`, `yields.csv`, `cutflow.csv`, `agreement.csv`, `events_one_side_only.csv`, `coverage.json`, and with `--plots` the distributions and ratios.

Full 2024 run with all the analysis weights (TAMU `weight__Central` from `tamu_analysis_weights.py`, DY reweights included, against his genWeight × PU × muon ID/iso × trigger SF × rescale):

```bash
bash sync/python/gpizzati_sync_2024.sh submit    # HTCondor, one job per missing dataset
bash sync/python/gpizzati_sync_2024.sh local DS  # one dataset here
bash sync/python/gpizzati_sync_2024.sh summary   # SUMMARY.md with all the datasets
```

Outputs: `/eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7_allweights/<dataset>/`, plots in `$HMUMU_WWW/skim_v4/studies/sync/gpizzati_v7_allweights/`.
The outputs of the 2024 run were deleted on 2026-10-20 (not needed any more): rerun the commands above to reproduce them.

## 2. Pisa (Filippo), end of September — variant `filippo_20260925`

The special sync skim reproduces Pisa's VBF selection (`config/Run3_2024/selections_sync_filippo_20260925.yaml`).

```bash
python3 sync/python/special_sync_submit.py --variant filippo_20260925          # HTCondor, one job per NanoAOD file
python3 sync/python/special_sync_merge.py  --tag filippo_20260925              # tuples + cutflow CSVs
python3 sync/python/compare_filippo_diagnostic.py --tag filippo_20260925 --output-dir $OUT   # per-event, where each side loses the other's events
python3 sync/python/plot_filippo_comparison.py    --tag filippo_20260925 --output-dir $WWW   # 17 distributions TAMU/Pisa with ratio
python3 sync/python/compare_cutflow_pisa.py                                     # our cutflow in their order + bridge of the definition differences (skim only)
```

MC: add `--dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg` to submit/merge, then `--filippo-input '<his ntuple>'` to the two compare scripts. For VBFH Signal_Fit use
`--region Signal_Fit --filippo-step HiggsMass`, `--weight genWeight --filippo-weight Generator_weight` in the plot script.
The outputs of this round (`confronto_filippo_*`) were deleted on 2026-10-20: rerun the commands above to reproduce them.
Weighted comparison (genWeight, PU, muon ID/iso, trigger): `sync_mc_weights.py` builds our weights for the tuple events (friend tuple), `tuple_weights_friend.py` does the same from an exchange tuple that already carries the weights,
`compare_filippo_weights.py` compares them factor by factor with his.

## 3. Producing the tuples for the other groups

The exchange tuple has the columns of `TAMU_SYNC_TUPLE_VARIABLES.md`; the request sent to the groups is `ggfvbf_sync_20261015.yaml`. One Condor job per NanoAOD file; the variant selects the selection yaml
(`config/Run3_2024/selections_sync_*.yaml`) and the wrapper.

| variant | selection | script / wrapper |
|---|---|---|
| `special_20260923`, `filippo_20260925` | first Pisa sync, data (and MC) | `special_sync_skim.py` / `run_special_sync_skim.sh` |
| `ggfvbf_20261015`, `ggfvbf_pisa_20261013`, `ggfvbf_iso_20261014` | ggF/VBF exchange tuple (50 < m < 200, raw pT > 20) | `special_sync_skim.py` / `run_special_sync_skim.sh` |
| `ggfvbf_analysis_20261014` | main-analysis selections as flags, Pisa definitions as flags only (private) | `special_sync_skim_analysis.py` / `run_special_sync_skim_analysis.sh` |
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

## 4. Standalone tool: compare the published tuples of two groups

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
