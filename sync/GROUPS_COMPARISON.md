# Group-by-group sync comparison

Tool: `sync/python/compare_run3_sync_groups.py`.
It compares the exchange ntuples of two groups published in
`/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/<GROUP>/`, sample by sample, region by region and channel by channel.

## Quick start

```bash
source env.sh
OUT=$HMUMU_WWW/skim_v4/studies/sync/groups_comparison_13Oct

# what would be run (pairs x common samples), nothing is computed
python3 sync/python/compare_run3_sync_groups.py --list

# every pair of groups, every common sample, all regions and channels
python3 sync/python/compare_run3_sync_groups.py --output $OUT

# one pair, one sample
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 PKU \
    --sample Data_Run2024H --output $OUT

# one cell only
python3 sync/python/compare_run3_sync_groups.py --group1 Marmara --group2 PKU \
    --sample ggH --region Signal_Fit --channel VBF --output $OUT
```

## Options

| option | meaning |
|---|---|
| `--output DIR` | output root (required unless `--list`) |
| `--group1`, `--group2` | the two groups, from `TAMU, Marmara, PKU, NNU, Pisa, Purdue`; must be given together and differ. Without them, all pairs are run |
| `--sample S` | restrict to a sample; repeatable. Only samples present in both groups are compared |
| `--region R` | `Z_sideband`, `H_sideband`, `Signal_Fit`; repeatable |
| `--channel C` | `ggH`, `VBF`; repeatable |
| `--resume` | skip pair/sample folders that already have `summary.json` |
| `--cache-dir DIR` | where the extracted columns are cached (default `/tmp/vdamante/run3_sync_extracts`) |
| `--list` | print `GROUP1-GROUP2/sample` and exit |

Samples: `Data_Run2024H`, `ggH`, `VBFH`, `ttH`, `ZH`, `WplusH` (whatever both groups provide).

## Output layout

```
<OUT>/<GROUP1>-<GROUP2>/<sample>/SUMMARY.md, summary.json
<OUT>/<GROUP1>-<GROUP2>/<sample>/<region>/SUMMARY.md
<OUT>/<GROUP1>-<GROUP2>/<sample>/<region>/<channel>/SUMMARY.md + one plot per common variable
```

e.g. `TAMU-PKU/Data_Run2024H/Z_sideband/VBF/`. The root of `<OUT>` gets a top-level summary.
Each leaf `SUMMARY.md` contains the event accounting (events absent from one tuple vs present but classified in another
region/channel, with the migration table), the event-by-event residuals on common events (relative tolerance 1e-3),
the efficiency of every flag, the yields, and an English explanation of the differences.

Regions (on the corrected dimuon mass): Z_sideband 70-110 GeV; H_sideband 110-115 and 135-150 GeV; Signal_Fit 115-135 GeV.
**Data in Signal_Fit is blinded**: only counts, no distributions.
MC is compared with two weight modes: `genWeights` and `all_common_weights`.

## How groups are plugged in

Every group has a `treat_<group>_ntuple()` function that returns a `GroupSpec`, the projection of its native schema onto
a common description: sample-to-file map, tree name, lumi branch, corrected mass, the ggH/VBF channel expressions,
the weight expressions, the flags to show as efficiencies, and notes on the interpretation. The comparison itself only
sees this common description, so any two groups can be paired.

To add a group (e.g. RICE):

1. In `compare_run3_sync_groups.py`, write `treat_<NEW>_ntuple()` returning a `GroupSpec`.
2. Register it in `TREATERS`.
3. In `VARIABLES`, add the `<NEW>=` expression to each variable `V(...)` it provides (expressions can use `abs`, `sqrt`, `pi`, `dphi`).
4. Run `--list` to check the new pairs, then `--group1 <NEW> --group2 <OTHER>`.

## Reproducibility: how it works, step by step

### Requirements

Python 3 with `numpy`, `awkward`, `uproot`, `matplotlib` (nothing from the TAMU framework is imported).
The script reads the files directly from `/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/` (`BASE` at the top of the script);
change `BASE` or the file paths in the `treat_*` functions to run on a local copy.
The only hard-coded paths besides `BASE` are the column cache (`--cache-dir`) and `--output`.

### 1. Projection of each native tuple (`treat_<group>_ntuple`)

Each group publishes a different schema. A `GroupSpec` states, in plain strings, how to read it:

| field | meaning |
|---|---|
| `files` | sample name -> file under `run3_sync/` (samples: `Data_Run2024H`, `ggH`, `VBFH`, `ttH`, `ZH`, `WplusH`) |
| `tree`, `lumi_branch` | TTree name and the luminosity-block branch (`run` and `event` are assumed) |
| `mass` | the corrected (FSR/calibrated) dimuon mass used to define the regions |
| `channels` | `ggH` and `VBF` as boolean expressions on branches, defined on the whole tuple (no implicit filtering) |
| `weights` | `genWeights` (generator weight only) and `all_common_weights` (the group's full per-event weight) |
| `flags` | boolean branches reported as efficiencies, never used for the selection |
| variables | taken from the table `VARIABLES`: one expression per group, `None` if the group lacks it |

Expressions are evaluated with `eval` on awkward arrays. Allowed syntax: branch names, `+ - * / & | ~ == < >`, `abs`, `sqrt`, `pi`,
`dphi(a,b)`, `branch[:,i]` (i-th element of a jagged branch, NaN if missing), and `branch[mask][:,i]` (i-th element after a boolean mask,
used for jets with `jet_passAnalysis`). **Expressions are literals written in the script; never feed it untrusted strings.**
A variable is silently skipped for a pair if one of the two groups does not provide it, or if its branches are missing in the file.

### 2. Event identity

Events are matched by a 64-bit key built from (run, lumi, event). For data: `(run-370000)<<48 | lumi<<34 | event`.
For MC: bit 63 set, `run<<59 | lumi<<39 | event`. Keys are unique within a tuple (the first occurrence is kept if not).
Nothing else is used to match events, so a group must publish run, lumi and event numbers.

### 3. Extraction and cache

For each (group, sample) the script reads only the needed branches in 200 MB chunks, evaluates mass, channels, variables, flags and weights,
and keeps the events with **70 < m < 150 GeV and in the ggH or VBF channel**. It also keeps the set of *all* event keys of the file
(`full_keys`), which is how it tells "this event is not in the other tuple" from "this event is there but is classified differently".
The result is cached in `<cache-dir>/<group>_<sample>.npz`; the cache is invalidated automatically when the file size/mtime or the
treatment (expressions, flags, weights) changes. The first read of a large file takes about 5 minutes (EOS I/O); later runs are fast.
Delete the cache directory to force a re-read.

### 4. Regions, channels, blinding

Regions are cuts on the corrected mass: `Z_sideband` 70-110, `H_sideband` 110-115 and 135-150, `Signal_Fit` 115-135 GeV.
The channel is whatever each group defines as ggH / VBF (see the treatment notes below); it is **not** re-derived from jets,
because the point is to compare the groups' own definitions. A sample's process is irrelevant: every sample is looked at in both channels.
`Data_Run2024H` in `Signal_Fit` is blinded: only the event accounting is written, no distributions.

### 5. What is computed in every (pair, sample, region, channel) cell

1. **Event accounting.** Let A and B be the two groups. Events selected by both; events selected only by A, split into
   *absent from the B tuple* (input coverage or preselection difference) and *present in the B tuple but not selected here*
   (definition difference); the same for B.
2. **Migration table.** For the events selected only by A and present in B, the first (channel / region) cell of B where they land
   (or "outside 70-150 GeV or ggH/VBF"). This says, for example, that A's VBF events are B's ggH events.
3. **Event-by-event residuals** on the events selected by both: for each common variable the fraction of events whose values differ by
   more than 0.1% (relative to max(|a|, |b|, 1)), the median difference and the 99th percentile of |a-b|.
   Jet variables of events with fewer jets than needed are treated as undefined and ignored.
4. **Distributions.** For every common variable and weight mode (`unweighted` for data; `genWeights` and `all_common_weights` for MC):
   unit-normalised histograms of both groups (density per unit of the variable), a ratio panel A/B
   (bins empty in either histogram have no ratio and are left out; the y range zooms on the differences, between +-2% and 0.5-1.5),
   means, standard deviations and the **total-variation distance** between the unit-normalised shapes (0 = identical, 1 = disjoint).
   Plots are saved as PNG and PDF. Because shapes are normalised, overall yield differences are only in the quoted sum of weights.
5. **Flags.** For each group, how many of the selected events pass each of its boolean flags.

### 6. Markdown and the English explanation

`SUMMARY.md` of a leaf is generated from those numbers: a rule-based explanation (how many events each group selects and share in common,
where the group-only events are lost or move to, the three variables with the largest shape distance and the type of difference they suggest
- muons, jets or MET -, the residual table), the event accounting, the tables, the plots, and the treatment notes of both groups.
The explanation is only a hint of where to look. It does not know the groups' definitions beyond the notes, and it never proves a cause:
that needs the object-level comparison (muon indices, first failing cut) of the one-sided events.
Each `<sample>/summary.json` has every number behind the markdown, and the root of the output has an index.

### 7. Reproducing the TAMU results from scratch

```bash
git clone git@github.com:valeriadamante/Hmumu_analysis.git && cd Hmumu_analysis
python3 sync/python/compare_run3_sync_groups.py --list
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 PKU --output my_output
```

The result depends only on the files in `run3_sync/` and on the `treat_*` functions; there is no random component.

### 8. Common mistakes

- Comparing the number of rows of two tuples: tuples are produced with different preselections, only the per-cell counts are meaningful.
- Reading `Marmara/data_*.root` or `Marmara/GluGluH-*.root`: they are NanoAOD, not analysis tuples.
- Comparing yields: weights are normalised differently by each group; compare shapes and the event accounting.
- Changing a `treat_*` function without deleting the cache is safe (the signature is checked), but changing the **value** of a branch without changing the file is not detected.

## Interpretation of each group's tuple (to be confirmed with the groups)

- **TAMU**: channels are the stored `ggF`/`VBF` flags (Pisa-style). Dimuon = the single good pair when it exists (`pisa_candidate`).
- **Marmara**: channel = `sel_full && cat_isGGH/VBF[0]`; the meaning of index 1 is unknown. `data_Run2024H_Muon*.root` and the `GluGluH-/VBFH-*.root` files are NanoAOD, not tuples, so Marmara data cannot be compared and MC uses `tree_*_loose_LAST.root`.
- **PKU**: one exclusive category per event (VBF, GGH, TTH, VH), only events passing the PKU selection are stored.
- **NNU**: only `passAnalysisSelection` is a selection (VBF); ggH is our reading (same common flags, category != 3). Only VBFH is available.

## Status (2026-10-14)

- Groups: TAMU, Marmara, PKU, NNU, Pisa, Purdue (RICE and `PKU-VH` are not compared). Full commands: `sync/RUN_ALL_COMPARISONS.md`.
- `HMUMU_TAMU_STAGING=<dir>` makes the TAMU treater read `<dir>/TAMU-<dataset>/sync_tuple.root` (datasets `Run2024H`, `GluGluHto2Mu`,
  `VBFHto2Mu_M125_powheg`, `TTH_Hto2Mu`, `ZH_Hto2Mu`, `WplusH_Hto2Mu`) instead of `run3_sync/TAMU/`, so the new TAMU tuples can be compared before they are published.
- TAMU keeps jets at pT > 35/25 GeV for VBF; the Pisa exchange tuples use 30/30 GeV, so jet differences TAMU-Pisa are expected.
- Old `TAMU-Pisa_25Sep2026` folders in the output come from an earlier study.
