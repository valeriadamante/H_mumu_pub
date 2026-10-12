# First ggF/VBF synchronization — 15 October

The `ggfvbf_20261015` variant of `sync/python/special_sync_skim.py` reads the
requested 2024 NanoAOD datasets directly.  It does not start from the analysis
skim, because that skim has already applied trigger matching, multiplicity and
lepton-veto selections.

## Ntuple preselection

The only event filters are:

```text
at least two NanoAOD muons with Muon_pt > 20 GeV (data, ggF, VBF)
at least two NanoAOD muons with Muon_pt > 10 GeV (ttH, ZH, W+H)
50 < mass(mu1, mu2) < 200 GeV
```

`mu1` and `mu2` are the two highest-`Muon_pt` muons passing the channel-specific threshold. The mass
is calculated from their uncorrected NanoAOD `(pt, eta, phi, mass)`.  ScaRe and
FSR are applied only afterwards to the corrected comparison branches.

The golden JSON, trigger, trigger matching, charge, muon ID/isolation,
electron veto, b veto, MET filters, jet quality, VBF and ggF definitions are
not event filters.  Their decisions are stored as flags (the golden JSON flag
is `pass_golden_json`).

## Submit and merge

```bash
source env.sh

OUT=/eos/cms/store/group/phys_higgs/cmshmm/run3_sync
TAG=ggfvbf_20261015
GROUP=vdamante

python3 sync/python/special_sync_submit.py \
  --variant "$TAG" --tag "$TAG" --group-name "$GROUP" \
  --output-root "$OUT" --flavour tomorrow

# Safe to repeat: completed jobs are skipped.
python3 sync/python/special_sync_submit.py \
  --variant "$TAG" --tag "$TAG" --group-name "$GROUP" \
  --output-root "$OUT" --flavour tomorrow

# Run after every job has produced both ROOT and report JSON files.
python3 sync/python/special_sync_merge.py \
  --tag "$TAG" --group-name "$GROUP" --output-root "$OUT"
```

The default dataset list is exactly the request: Muon0/Run2024H,
Muon1/Run2024H, VBF H, gluon-fusion H, ttH, ZH and W+H. For internal
sample entries containing an additional `ext1` campaign, only the explicitly
requested primary dataset is used.  Missing cached file lists (currently ttH
and VH) are obtained with a read-only `dasgoclient` query.

Each tuple is written as `{GROUPNAME}-{DATASETNAME}/sync_tuple.root`. The merge
also produces the combined data and VH directories:

- `vdamante-Run2024H/sync_tuple.root` (Muon0 + Muon1);
- `vdamante-VH/sync_tuple.root` (ZH + W+H).

## Stored information

- Event identifiers, string and integer sample identifiers, trigger and
  matching decisions, golden-JSON decision, mass regions and ggF/VBF flags.
- Both muons: raw NanoAOD and ScaRe/FSR kinematics, charge, ID, isolation,
  impact parameters and matching information; raw and corrected dimuon
  kinematics.
- All jets and selected jets: corrected, pre-JERC NanoAOD and raw pT,
  eta/phi/mass, b tagging, ID/veto/horn/cleaning flags, flavours for MC,
  selected indices and leading/VBF-pair observables.
- `MET_pt`, `MET_phi` and `MET_collection`, using `PuppiMET`.
- For ZH/W+H, the complete NanoAOD `FatJet_*` collection, `nFatJet` and
  `FatJet_pt_beforeJERC`.
- MC weights kept separately: `weight_generator`, `weight_pileup`, every
  central per-muon ID/isolation/reconstruction factor,
  `weight_muon_idiso`, and `weight_trigger`.  No DY or dimuon-pT analysis
  reweighting is evaluated.

Each per-file ROOT output has a matching report JSON.  The merger refuses to
run if any requested input is missing and writes cutflow and category/region
CSV summaries beside the final tuples.
