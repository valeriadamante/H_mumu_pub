# skim_v4_sync

One exact NanoAOD input file per run, using the current skim corrections and
nominal analysis selection. The default is Muon0 Run2024C MINIv6NANOv15-v1,
first file in `config/Run3_2024/samples_withfiles.yaml`. Agree on that exact
input file with the other analysis before comparing counts.

```bash
source env.sh
# Print the input, output and command without running:
python3 sync/python/run_sync_skim.py
# Process the complete file and write to the production EOS area:
python3 sync/python/run_sync_skim.py --run
# Use an agreed file or an MC sample instead:
python3 sync/python/run_sync_skim.py --run --era Run3_2024 \
  --dataset DATASET_KEY --input-file /path/to/exact/nanoaod.root
```

Default output root:
`<skim>_sync` (`skim` in `config/paths.yaml`).
Subdirectories are `ERA/DATASET/INPUT_FILE_STEM/CATEGORY`.
Use `--output-dir` to change the root. Existing output directories are refused
so an earlier comparison cannot be silently overwritten. `--n-events N`
creates a distinct `test_N` subdirectory; omit it for an exact full-file comparison.

The default category is `baseline`; `--category VBF` and `--category ggF` are
also available. Definitions are read from the era's `selections.yaml`.
OS remains a flag in the general skim; the sync `baseline` requires OS through
its normal analysis definition. Trigger matching, muon 26/20 GeV cuts, muon ID,
isolation and b-tag veto therefore follow the current analysis selection.
Data output contains only `H_sideband`: 110 < m_mumu < 115 or
135 < m_mumu < 150 GeV. MC also includes `Signal_Fit`: 115 < m_mumu < 135 GeV.
The strict boundaries are inherited from the YAML. Systematic variations are
disabled for sync. MC retains the normal skim's sample/era splitting and weights;
synchronization counts and distributions are unweighted.

Each run writes:

- `events.root`: the selected sync events and skim branches, plus raw jet pT.
- `H_sideband_events.jsonl`: sorted per-event muon and jet details.
- `H_sideband_event_ids.txt`: sorted `run:luminosityBlock:event` identifiers.
- Corresponding `Signal_Fit_*` files for MC only.
- `summary.json`: exact input, region definitions, counts and duplicate IDs.
- `cutflow.json`: the skim and sync selection counts.
- `distributions.root`: unweighted mass, muon pT/eta, jet multiplicity and
  selected-jet pT/eta histograms, separated by region. Jet histograms include
  every selected jet, so their entries count jets rather than events.
- `manifest.json`, `configuration/` and `production.log`: input dataset,
  command, source/configuration hashes, configuration copies and correction log.

Besides the objects, the export also carries the dimuon variables and the
high-level ggH/VBF training inputs (`pt_mumu`, `eta_mumu`, `phi_mumu`, `y_mumu`,
`dR_mumu`, `cosTheta_CS`, `phi_CS`, `R_pt`, `minDeltaEtaSigned`, `minDeltaPhi`,
`Zeppenfeld_Var`, `pt_centrality`, `m_jj`, `m_jj_ls`, `delta_eta_jj`,
`delta_eta_jj_ls`, `pt_vbfj1j2`, `era_code` and the soft-activity counts), so a
disagreement can be traced to the DNN inputs and not only to the objects.

Muon details include raw/noCorr and corrected + FSR pT, charge, ID/isolation,
trigger matching and gen matching. Jet details include all jets and the indices
of selected jets, correction input pT, true raw pT, corrected pT, eta/phi,
horn/veto-map/selection flags. `Jet_pt_nocorr` is the NanoAOD pT before our
reapplication of JEC; true raw pT is explicitly stored as
`Sync_Jet_pt_raw = Jet_pt_nocorr * (1 - Jet_rawFactor)`.

Write the event list in the column order agreed with the other analysis:

```bash
python3 sync/python/export_event_list.py ours/H_sideband_events.jsonl \
  --output ours/H_sideband_event_list.txt
```

`--precision N` sets the number of decimals, `--missing STR` the placeholder for
absent values, `--extra` adds the supplementary columns, `--mc-weights` the MC
weights, and `--no-header` drops the header line.

Compare the event lists before comparing distributions:

```bash
python3 sync/python/compare_sync_events.py \
  ours/H_sideband_event_ids.txt theirs/H_sideband_event_ids.txt \
  --output comparison.json
```

Exit status is 0 only for identical event multisets. Differences and duplicate
counts are reported; no event deduplication is performed. Event IDs remain
integers throughout JSON export, including values above 2**53.

## Compact sync tuple and VBF cutflow

`sync/python/prepare_tuple_forSync.py` re-reads an existing skim and writes the
compact sync tuple plus its cutflow. The selections come from
`config/Run3_2024/selections_sync.yaml`; `--cutflow-key` chooses which
cumulative step list in that file drives the cutflow, `--selection` the filter
applied to the output tuple.

```bash
source env.sh
SKIM=$HMUMU_SKIM/Run3_2024
python3 sync/python/prepare_tuple_forSync.py $SKIM/Muon0_Run2024H $SKIM/Muon1_Run2024H \
  -o $HMUMU_SYNC/sync_tuple_Run2024H_VBF_Z_CR_H_SB.root \
  --selection "VBF && Z_CR_H_SB" --cutflow-key cutflow_vbf_Z_CR_H_SB --threads 12
```

Full instructions (Italian): [`LOOSE_TUPLE.md`](LOOSE_TUPLE.md).

Loose tuple for the sync with Pisa (default since 2026-10-09): `--selection` is
`mass_50_200` (50 < m_mumu < 200 GeV, the skim window), and every region,
category and muon requirement of `selections_sync.yaml` with `store: true` is a
boolean branch, so different definitions can be compared without re-skimming.
Weights and observables come from `common.prepare_rdf` (`--no-framework` skips it).
The tuple keeps run/lumi/event, `HLT_IsoMu24` and the trigger matching, muon pT
raw (`_raw_noCorr`), ScaRe (`_raw_scale`, `_raw_corr`) and ScaRe+FSR (`mu*_pt`),
ID/isolation, `m_mumu`/`pt_mumu` and their `_noCorr` versions, the `Jet_*`
(all jets) and `SelectedJet_*` collections with b-tag scores and flags, jet
ID/veto-map/horn flags, the leading and VBF dijets, and with `--mc` every nominal
weight (`weight__Central`, its factors, `inv_N_orig`).

Jet pT comes in three versions: `Jet_pt` (our JEC, + JER in MC), `Jet_pt_nocorr`
(NanoAOD pT, with the NanoAOD JEC) and `Jet_pt_raw` (no JEC, no JER), also as
`SelectedJet_pt_raw`, `leadingjet_pt_raw`, `subleadingjet_pt_raw`,
`thirdleadingjet_pt_raw`, `vbfjet1_pt_raw`, `vbfjet2_pt_raw`.
skim_v4 does not store `Jet_rawFactor`, so `add_raw_jet_pt.py` (run automatically
at the end, `--no-raw-jet-pt` skips it) reads it from the NanoAOD files of the
datasets (`samples_withfiles.yaml`), matches the events on run:lumi:event and on the
jet collection (NanoAOD `Jet_pt` must equal `Jet_pt_nocorr` jet by jet, which also
separates the MC extensions), stops if any tuple event has no unique match, and
sets `Jet_pt_raw = Jet_pt_nocorr * (1 - Jet_rawFactor)`. The NanoAOD reads are cached
in `<output stem>_nano_cache/` (rerunning resumes; delete it when the tuple is final).
On VBFH (190 NanoAOD files) the reading takes about 5 minutes with 24 workers.

`<output stem>_metadata.json` (also the TNamed `metadata` in the tuple) lists the
datasets, luminosity and cross sections, corrections and payloads, selections and
which pT is used for each cut and scale factor.

```bash
source env.sh
SKIM=$HMUMU_SKIM/Run3_2024
OUT=$HMUMU_SYNC/loose_tuples_20261009
python3 sync/python/prepare_tuple_forSync.py $SKIM/Muon0_Run2024H $SKIM/Muon1_Run2024H \
  -o $OUT/sync_tuple_loose_Run2024H.root --threads 12
for ds in VBFHto2Mu_M125_powheg GluGluHto2Mu; do
  python3 sync/python/prepare_tuple_forSync.py $SKIM/$ds --mc \
    -o $OUT/sync_tuple_loose_$ds.root --threads 12
done
```

The steps before `muon_pre_sel_check` are the skim `report_*.json` counts, so
the cutflow starts from the NanoAOD events before any cut. `cutflow_vbf_Z_CR_H_SB`
splits `HasVBF` into `vbf_jets_ge2`, `vbf_pair_pt`, `vbf_pair_mjj` and
`vbf_pair_deta` in the order `FindVBFJets` applies them; `vbf_pair_deta` must
reproduce `VBF_def` exactly. Steps that the skim already enforces stay in the
list on purpose and show zero rejected events.

Since 2026-10-09 the sync muon cuts are leading > 26 GeV and subleading > 20 GeV
on the raw pT (`mu1_rawpt_g26`, `mu2_rawpt_g20` on `mu*_pt_raw_noCorr`), as Pisa;
before they were 26/26 GeV on the corrected pT. `cutflow_vbf_Z_CR_H_SB_mu2_20` and
the `VBF_mu2_20` category repeat the same chain with 26/20 GeV on the corrected
(ScaRe + FSR) pT. `VBF_mu2_20` is not stored, so the tuple keeps the nominal
column set:

```bash
python3 sync/python/prepare_tuple_forSync.py $SKIM/Muon0_Run2024H $SKIM/Muon1_Run2024H \
  -o $HMUMU_SYNC/sync_tuple_Run2024H_VBF_Z_CR_H_SB_mu2pt20.root \
  --selection "VBF_mu2_20 && Z_CR_H_SB" \
  --cutflow-key cutflow_vbf_Z_CR_H_SB_mu2_20 --threads 12
```

`sync/python/compare_cutflow_pisa.py` (`hmumu sync compare_cutflow_pisa`) reconciles our VBF cutflow with the Pisa one.
It runs on our skim only and takes their numbers from the table hard-coded in
the file; it writes our cutflow in their step order, the same cutflow with
their jet and muon definitions, and a bridge that switches the definition
differences on one at a time.

## Sync VBF with Pisa (Filippo), Run2024H

The `filippo_20260925` variant of the special sync skim reproduces Pisa's VBF
selection (definitions in `config/Run3_2024/selections_sync_filippo_20260925.yaml`).

```bash
source env.sh
python3 sync/python/special_sync_submit.py --variant filippo_20260925   # HTCondor, one job per NanoAOD file
python3 sync/python/special_sync_merge.py --tag filippo_20260925        # tuples + cutflow CSVs
python3 sync/python/compare_filippo_diagnostic.py --tag filippo_20260925 \
    --output-dir /eos/user/v/vdamante/H_mumu/Sync/confronto_filippo_20260925
python3 sync/python/plot_filippo_comparison.py --tag filippo_20260925 \
    --output-dir /eos/user/v/vdamante/www/H_mumu/skim_v4/studies/sync/filippo_20260925
```

MC (ready, to run once Filippo's MC diagnostic ntuples exist): the same chain
takes any dataset of `config/Run3_2024/samples_withfiles.yaml`. MC tuples keep
`genWeight`; `--mc-lumi-split` applies the production Run3_2024 event split
(off by default: it hashes the job's file list, so it does not reproduce the
production skim's events, and the other analysis may not split at all).

```bash
python3 sync/python/special_sync_submit.py --variant filippo_20260925 --dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg
python3 sync/python/special_sync_merge.py --tag filippo_20260925 --dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg
python3 sync/python/compare_filippo_diagnostic.py --tag filippo_20260925 --dataset GluGluHto2Mu \
    --filippo-input '/eos/home-f/fcattafe/Hmm_sync_valeria/<his MC file, {ds} = dataset>' \
    --output-dir /eos/user/v/vdamante/H_mumu/Sync/confronto_filippo_mc
python3 sync/python/plot_filippo_comparison.py --tag filippo_20260925 --dataset GluGluHto2Mu \
    --filippo-input '<same>' --weight genWeight --filippo-weight <his weight branch> \
    --output-dir /eos/user/v/vdamante/www/H_mumu/skim_v4/studies/sync/filippo_mc
```

VBF signal region, VBFH 2024 (Filippo's `VBF_SignalRegion` ntuple, 2026-10-05): his
`HiggsMass` is 115-135 GeV, so ours is `VBF && Signal_Fit`; his last step `MC_selection`
is a random ~50% event split that we do not reproduce, so his side stops at `HiggsMass`.
His `*_cumulativeWeight` start from `Generator_weight` and are not divided by the sum of
weights (yield = sum / sum(genEventSumw) * xsec * lumi).

```bash
FIL='/eos/home-f/fcattafe/Hmm_sync_valeria/hmm_2024Summer24_VBF-Hto2Mu_M-125_powheg_*_VBF_SignalRegion_diagnostic_ntuples.root'
python3 sync/python/compare_filippo_diagnostic.py --dataset VBFHto2Mu_M125_powheg --region Signal_Fit \
    --filippo-step HiggsMass --filippo-input "$FIL" --output-dir $HMUMU_SYNC/confronto_filippo_mc_vbfh
python3 sync/python/plot_filippo_comparison.py --dataset VBFHto2Mu_M125_powheg --region Signal_Fit \
    --filippo-step HiggsMass --filippo-input "$FIL" --weight genWeight --filippo-weight Generator_weight \
    --output-dir $HMUMU_WWW/skim_v4/studies/sync/filippo_mc_vbfh
```

Weighted comparison with all the weights in common (genWeight, PU, muon ID/iso and
trigger SFs): `sync_mc_weights.py` (our analysis weights for the tuple events) and
`compare_filippo_weights.py`; results of 2026-10-07 in `sync/SYNC_PISA_MC_VBFH_2024.md`.

`compare_filippo_diagnostic.py` matches events on run:lumi:event and lists,
for the events selected by one side only, where the other side loses them.
`plot_filippo_comparison.py` draws 17 distributions TAMU (ours) vs Pisa with
the ratio; the ratio uncertainty counts only the non-shared events.

## Sync with gpizzati hmm_skim v7 (Run3_2024, data and MC)

`compare_gpizzati_v7.py` applies the selection of `config/Run3_2024/selections.yaml`
(OS, trigger, muon pT 26/20, ML, TT below 30 GeV, b-veto, VBF_def; categories
baseline/ggF/VBF, regions Z_sideband/H_sideband/Signal_Fit) to our skim_v4 and to
gpizzati's v7 skim, each with its own objects, and compares them event by event
on the lumisections present in both. MC events are matched also on the raw eta
of the two muons, because extensions reuse run:lumi:event.

```bash
source env.sh
# data, default inputs, CSV + markdown only
python3 sync/python/compare_gpizzati_v7.py --dataset Muon0_Run2024C \
    --output-dir /eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7
# MC with genWeight and plots; several datasets in one call
python3 sync/python/compare_gpizzati_v7.py --dataset GluGluHto2Mu --dataset VBFHto2Mu_M125_powheg \
    --weight genWeight --plots --output-dir /eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7
# explicit inputs: a file, a folder or a glob for each side (one --dataset only)
python3 sync/python/compare_gpizzati_v7.py --dataset DYto2Mu_MLL_105to160_amcatnloFXFX \
    --ours-input /path/to/our/skim/folder --theirs-input '/path/to/his/DY*_skimmed_*.root' \
    --output-dir /tmp/sync_dy --plots --plots-dir /eos/user/v/vdamante/www/H_mumu/...
```

Per dataset it writes `README.md` (input events, coverage, yields, yields on the
events present in both skims, where one-side events are lost, cutflows,
variable agreement), `yields.csv`, `cutflow.csv`, `agreement.csv`,
`events_one_side_only.csv`, `coverage.json`, and with `--plots` the
distributions and ratios in `plots/<category>_<region>/`. Data in Signal_Fit
needs `--unblind`. Our 2024 MC skims keep only part of the lumisections
(`Orthogonal_lumi_split_for_Run3_2024`) and include the extensions, so for MC
the event-level sync is the table restricted to the events present in both skims.

### Full 2024 sync with all the analysis weights

`sync/python/gpizzati_sync_2024.sh` runs `compare_gpizzati_v7.py` on every 2024 data dataset and on
the main MC with `--weight analysis`: TAMU `weight__Central` from the framework
(`sync/python/tamu_analysis_weights.py`, DY reweights included; `--no-dy-weights` to drop them)
against gpizzati genWeight x PU x muon ID/iso x trigger SF x rescale (his `results.json`).
Each README then has a table of the single weight factors on the common events.

```bash
bash sync/python/gpizzati_sync_2024.sh submit    # HTCondor, one job per missing dataset
bash sync/python/gpizzati_sync_2024.sh summary   # SUMMARY.md with every dataset
```

Outputs: `/eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7_allweights/<dataset>/`, plots in
`$HMUMU_WWW/skim_v4/studies/sync/gpizzati_v7_allweights/`. DY M-50, TT 2L2Nu and DY 105-160 use a
subset of files (60/60/80).
