# Preparing training tuples from analysis skims

Use `common.prepare_rdf.prepare_rdf` to build the same event observables,
selection columns, and weights used by histogram production. The function
returns ROOT RDataFrames. You choose whether to make histograms, write a
training tuple with `Snapshot`, or export arrays with `AsNumpy`.

This guide starts from **the skim ntuples produced by this analysis**, with an
`Events` tree and their JSON bookkeeping reports. Raw NanoAOD files must first
go through `analysis/skim.py`: `prepare_rdf` expects the skim schema, including
selected objects and correction weights.

## 1. Understand the three steps

```text
Skim ROOT files + bookkeeping reports + era configuration
                          |
                     prepare_rdf
                          |
        RDF with observables, weights, and selection flags
                          |
           choose a region/category and output columns
                          |
              Snapshot / AsNumpy / histograms
```

An RDF is a computation graph. `Define` adds a column and `Filter` adds a
selection; both return a new node. Keep the returned node, for example
`rdf = rdf.Filter("VBF")`. Event processing usually happens when an action
needs its result, such as `Snapshot`, `AsNumpy`, or `Count().GetValue()`.
Optional DNN inference can also materialize inputs during preparation.

`prepare_rdf` returns a dictionary:

- `nodes["inclusive"]` is the prepared graph before a final region/category
  filter. **“Inclusive” does not mean that the VBF or signal-region selection
  has already been applied.** The original skim cuts still apply.
- With `split_jet_multiplicity=True`, additional entries contain component
  selections. Leave this disabled for a first training-tuple production.
- No usable ROOT inputs produce an empty dictionary, `{}`.

Call `prepare_rdf` once on a raw skim graph. Reusing its returned nodes for
several outputs is fine; passing an already prepared node back into the
function would repeat definitions and weight corrections.

## 2. Set up the environment and validate a dataset

Run from the repository root:

```bash
source env.sh

python3 analysis/validate_dataset.py \
  --era Run3_2025 \
  --dataset-name VBFHto2Mu_M125_amcatnlo \
  --root-input /path/to/skim/Run3_2025/VBFHto2Mu_M125_amcatnlo \
  --json-input /path/to/reports/Run3_2025/VBFHto2Mu_M125_amcatnlo \
  --output-manifest /path/to/manifests/Run3_2025/VBFHto2Mu_M125_amcatnlo.json
```

Replace the paths with your production locations. If a passed validation
manifest already exists for this production, reuse it. It contains the era,
exact dataset name, data/MC flag, accepted ROOT files, and accepted reports.
Use a separate manifest and output file for each dataset and era.

The era configurations have different purposes:

| Configuration | Used for |
|---|---|
| `config/<era>/maincfg.yaml` | Analysis settings, including `bTagAlgo` |
| `config/<era>/selections.yaml` | Mass regions, object selections, and categories |
| `config/<era>/systematics.yaml` | Central weight expression and variation definitions |
| `config/<era>/samples.yaml` | Dataset definitions; use the exact dataset name |

The example below uses `Signal_Fit` and `VBF`, both defined in the Run3_2025
selection configuration. These are example training selections, not a choice
made automatically by `prepare_rdf`. Choose the region appropriate to your
training study. For DY/EWK, the generated mass window must also match the
analysis region; the example checks the project's dataset routing rule.

## 3. Write a complete central training tuple

Save the following example as `make_train_tuple.py` in the repository root.
It writes one nominal tuple for one MC dataset. Run it separately for signal
and background samples, assigning the label appropriate to your classifier.
The example deliberately requires MC; using collision data for training is a
separate choice about labels and sample selection.

```python
import argparse
import json
from pathlib import Path

from common.prepare_rdf import prepare_rdf, prepare_region_dataframes
from common.utilities import (
    dataset_region_allowed,
    get_config,
    get_segmentation_dict,
    initialize_root_runtime,
    read_manifest,
)

parser = argparse.ArgumentParser()
parser.add_argument("--manifest", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--label", type=int, choices=(0, 1), required=True)
parser.add_argument("--region", default="Signal_Fit")
parser.add_argument("--category", default="VBF")
args = parser.parse_args()

initialize_root_runtime()
manifest = read_manifest(args.manifest, expected_stage="validation")
if manifest.get("status") != "passed":
    raise RuntimeError("Use a passed validation manifest")
if manifest["is_data"]:
    raise ValueError("This example builds labelled MC training tuples")
if manifest.get("tree", "Events") != "Events":
    raise ValueError("prepare_rdf expects the Events tree")

root_files = manifest["valid_root_files"]
report_files = manifest["valid_json_files"]
if not root_files or not report_files:
    raise RuntimeError("The manifest has no usable MC inputs")
era = manifest["era"]
dataset = manifest["dataset"]
if not dataset_region_allowed(dataset, args.region):
    raise ValueError("The dataset mass window does not match this region")

cfg_dir = Path("config") / era
main_cfg = get_config(cfg_dir / "maincfg.yaml")
sel_cfg = get_config(cfg_dir / "selections.yaml")
syst_cfg = get_config(cfg_dir / "systematics.yaml")
central = {"Central": {"weight": "weight__Central"}}

# Always use reports for the entire validated dataset, even for an input batch.
seg_dict = get_segmentation_dict(report_files)
if not seg_dict:
    raise RuntimeError("Missing generator normalization in the MC reports")

nodes = prepare_rdf(
    input_files=root_files,
    dataset_name=dataset,
    era=era,
    is_data=False,
    selections_cfg=sel_cfg,
    systematics_cfg=syst_cfg,
    seg_dict=seg_dict,
    systs_to_run=central,
    want_variations=False,
    btag_algo=main_cfg.get("bTagAlgo", "PNet"),
    dnn_payloads=None,
    enable_custom_weights=True,
    split_jet_multiplicity=False,
    skip_validation=True,  # Files come from the passed validation manifest.
)
if not nodes:
    raise RuntimeError("No RDF was prepared")

# These are example scalar features, not a prescribed model architecture.
features = ["pt_mumu", "mu1_pt", "mu2_pt", "m_jj", "delta_eta_jj"]
regions = prepare_region_dataframes(
    nodes["inclusive"],
    mass_regions=[args.region],
    categories=[args.category],
    systs_to_run=central,
    variables=features,
    btag_algo=main_cfg.get("bTagAlgo", "PNet"),
    era=era,
)
# The empty suffix denotes the nominal selection.
rdf, _ = regions[(args.region, args.category, "")]
rdf = rdf.Define("label", str(args.label))

# Keep bookkeeping and evaluation quantities separate from model inputs.
identifiers = ["run", "luminosityBlock", "event", "FullEventId"]
columns = list(dict.fromkeys(
    features + ["m_mumu", "weight__Central", "label"] + identifiers
))
available = {str(name) for name in rdf.GetColumnNames()}
missing = sorted(set(columns) - available)
if missing:
    raise RuntimeError(f"Requested tuple columns are missing: {missing}")

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
if output.exists():
    raise FileExistsError(output)

# Book checks before Snapshot so they can share its event loop.
count = rdf.Count()
sumw = rdf.Sum("weight__Central")
sum_absw = rdf.Define("tuple_abs_weight", "std::abs(weight__Central)").Sum("tuple_abs_weight")
rdf.Snapshot("Events", str(output), columns)

metadata = {
    "manifest": str(Path(args.manifest).resolve()),
    "dataset": dataset,
    "era": era,
    "region": args.region,
    "category": args.category,
    "label": args.label,
    "features": features,
    "columns": columns,
    "events": int(count.GetValue()),
    "sum_weights": float(sumw.GetValue()),
    "sum_abs_weights": float(sum_absw.GetValue()),
    "want_variations": False,
    "enable_custom_weights": True,
    "main_config": main_cfg,
    "selection_config": sel_cfg,
    "systematics_config": syst_cfg,
}
output.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2))
print(f"Wrote {metadata['events']} events to {output}")
```

Run it after sourcing the environment:

```bash
python3 make_train_tuple.py \
  --manifest /path/to/manifests/Run3_2025/VBFHto2Mu_M125_amcatnlo.json \
  --output /path/to/train_tuples/Run3_2025/VBFHto2Mu_M125_amcatnlo.root \
  --label 1 --region Signal_Fit --category VBF
```

For a binary classifier, `1` and `0` commonly represent signal and background.
Define which physical processes belong to each class before producing the
files. The label is a supervised-learning target; it is not an input feature.
Keep the metadata sidecar with each tuple. Record the repository revision and
custom-weight payload versions used for a production as well.

## 4. Understand the weights before training

For MC, `weight__Central` combines the configured central event weight with
generator normalization from the skim reports, followed by the applicable
custom corrections. It is the physical analysis weight, not automatically
the sample weight your training algorithm should consume.

- Use the **full validated dataset's reports** to compute `seg_dict`. If you
  process 5 of 100 ROOT files in a batch, still use the same full-dataset
  normalization as all other batches.
- Do not recompute the denominator after a training-region filter, or
  normalize each batch independently. Doing so changes the physical yield.
- Preserve the signed `weight__Central`. MC can contain negative weights.
  Class balancing, clipping, or an absolute-weight prescription for a chosen
  training method should produce a separate training weight. Do not replace
  the physical weight in the tuple.
- `enable_custom_weights=True` enables the configured DY component, dimuon-pT,
  and jet-multiplicity reweights where applicable. Their payloads must exist.
  Match the choices used by the histogram campaign you want to reproduce;
  its command-line options can override these defaults.
- `enable_custom_weights=False` disables those custom reweights, but retains
  base normalization and the dedicated DY amcatnlo normalization correction.
  It does **not** create unweighted events.

For comparison with nominal analysis histograms, use the same input manifest,
region, category, configuration, and weight options. Compare the selected
entry count, sum of weights, and a few weighted distributions. An event count
and a weighted yield are different quantities; signed weights can cancel.

## 5. Choose features and preserve event identity

`Snapshot` writes only the explicitly requested columns. Save scalar features
needed by the model, the physical weight, the label, and event identifiers.
Other useful quantities, such as `m_mumu`, can be saved for validation without
being passed to the model. Decide separately whether the mass is an intended
feature: in a mass-fit analysis, this affects the training objective and mass
correlations.

The example exports VBF observables after the VBF selection. In a more
inclusive selection, some jet observables can contain sentinel values when
there are too few jets. Inspect the feature definitions in `add_vars.py` and
choose an explicit treatment for unavailable objects. A missing column and a
sentinel value inside an existing column are different cases.

Partition train/validation/test samples by a deterministic event key, keeping
all copies and systematic versions of an event in the same partition. Store
the dataset and era together with the event identifiers. Snapshot row order,
particularly with multithreading, is not an event identifier. `FullEventId`
is useful within a fixed skim production, but its construction also depends
on skim input grouping; do not assume it is unchanged when re-skimming with
different file chunks.

For multiple datasets, keep the output files separate until the training
reader combines them. If you merge them into one tree, first add a dataset ID
and retain its mapping, so sample provenance is not lost.

## 6. Arrays, existing RDFs, and larger productions

For a small selected sample, replace the output action with:

```python
arrays = rdf.AsNumpy(features + ["weight__Central", "label"])
```

This materializes the requested columns in memory. For large samples, write
ROOT tuples and read them in batches during training. If you execute both
`Snapshot` and `AsNumpy`, expect additional event processing; these are two
separate output actions.

You can also pass an existing raw skim RDF as `rdf=raw_rdf` to `prepare_rdf`
instead of supplying `input_files`. Supply the same configurations and MC
normalization. An `additional_cuts` expression is evaluated before derived
observables are built, so it must refer to columns already present in the
raw skim. Apply cuts on newly prepared features to the returned RDF instead.

For batching, slice `root_files` while retaining the full `seg_dict`, and
write distinct output filenames. Do not sum the component nodes returned by
`split_jet_multiplicity=True` together with `inclusive`: they contain
additional views of the same events, not independent samples.

## 7. Systematics and existing DNN outputs

Start with `want_variations=False`, as in the example. This prepares the
nominal graph even if the input skim also stores shifted columns.

Producing variations requires more than switching that flag: pass a
consistent `systs_to_run` mapping, use skims containing the requested shifted
inputs, and provide the required normalization metadata. In particular,
QCD-scale weights use `qcd_scale_seg_dicts`. The histogram producer's
`get_systs_to_run` and related configuration handling show how campaign
variations are selected. Choose explicitly which shifted features and weights
to export; they are not automatically included in a Snapshot column list.

`dnn_payloads=None` leaves optional DNN inference disabled. This is appropriate
when preparing ordinary inputs for a new model. Request existing model
outputs only when they are part of your intended study. When
`prepare_region_dataframes` is given `DNN_NNOutput` in `variables`, supported
sideband regions can trigger mass-remapped DNN evaluation; this differs from
simply filtering the nominal graph. For repeated DNN-enabled productions in
one process, call `common.dnn_application.clear_prediction_registry()` after
all actions using the current predictions have completed.

## 8. Where the shared code lives

| Module | Responsibility |
|---|---|
| `prepare_rdf.py` | Shared RDF preparation for histograms, training, and snapshots |
| `utilities.py` | Configuration, file discovery, reports, chunking, validation, manifests, ROOT setup, and histogram binning |
| `add_vars.py` | Physics observables and region/category definitions |
| `trigger_weights.py` | Event trigger scale factors and variations |
| `jer_split.py` | Split JER variations |
| `apply_custom_weights.py` | Dataset-dependent custom weight corrections |
| `jet_component_splitting.py` | Reco/gen jet matching and component categories |
| `dnn_application.py` | Inference with existing DNN payloads |

Use `prepare_rdf` as the public preparation entry point. Its lower-level
`GetRdfForDataset`, `build_rdf`, and finalization helpers implement successive
parts of the pipeline; calling only the loader/base builder does not reproduce
all final selections and weights. Input TChains are retained by the module so
that lazy RDF actions can still read them.

`utilities.py` loads ROOT or NumPy only inside functions that need them, so
configuration and campaign commands can import it without loading ROOT.
`list_root_files` preserves the input path style; `discover_root_files`
produces absolute paths for validation.

`HistHelper.h` is a legacy header. Current production uses
`analysis/AnalysisTools.h`, initialized by `initialize_root_runtime()`.
Campaign orchestration lives in `campaigns/workflow.py`, invoked by the shell
entry points through `campaigns/submit_campaign.sh`; it manages submission,
completeness checks, merging, and plotting. It is not required to run the
single-dataset training-tuple example above.
