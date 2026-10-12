import os
import sys
import argparse
import zlib
from pathlib import Path
import ROOT
import json

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])

headers_dir = os.path.dirname(os.path.abspath(__file__))
ROOT.gInterpreter.Declare(
    f'#include "{os.path.join(headers_dir, "AnalysisTools.h")}"'
)
import common.utilities as utilities

parser = argparse.ArgumentParser(description="Run the Hmumu skim.")
parser.add_argument("--era", required=True)
parser.add_argument("--input-file",required=True,action="append")
parser.add_argument("--dataset-name", required=True)
parser.add_argument("--output-file", required=True)
parser.add_argument("--report-file",default=None)
parser.add_argument("--n-events", default=-1, type=int,
                    help="Process at most this many input events before selections; -1 processes all events.")
parser.add_argument("--jet-horn-veto", choices=("configured", "with", "without"),
                    default="configured", help="Override the jet horn veto for this skim only.")
parser.add_argument("--want-variations", required=False, action="store_true", help="request for variations from command line")
parser.add_argument("--sync", action="store_true", help="Write a nominal sync skim, event lists and distributions; data sideband only.")
parser.add_argument("--sync-category", choices=("baseline", "VBF", "ggF"), default="baseline")
args = parser.parse_args()
if args.n_events != -1 and args.n_events <= 0:
    parser.error("--n-events must be -1 (all events) or positive")
input_files = [ path for value in args.input_file for path in value.split(",") if path]
if not input_files:
    parser.error("--input-file did not contain any input paths")

## all configurations to load ##
config = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", args.era, "maincfg.yaml"))

dataset_cfg = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", args.era, "samples.yaml"))[args.dataset_name]
sel_config = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", args.era, "selections.yaml"))
from common.jet_horn_policy import configure_horn_veto
configure_horn_veto(sel_config, args.era, args.jet_horn_veto)
trigger_config = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", args.era, "triggers.yaml"))
process_cfg = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", args.era, "process_names.yaml"))
systematics_cfg = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", args.era, "systematics.yaml"))
xs_cfg = utilities.get_config(config["crossSectionsFile"])

## some utilities definitions ##
if args.sync:
    config["want_variations"] = False
nano_version = config.get("nano_version", "v15")
is_data = dataset_cfg.get("is_data", False)
is_signal = dataset_cfg.get("is_signal", False)
process = utilities.process_from_dataset(process_cfg, args.dataset_name) if not is_data else None

want_variations = (config.get("want_variations", False) or args.want_variations) and not is_data and not args.sync

muon_pt_default_suffix = sel_config.get("muon_pt_default_suffix", "")

# columns to save #
cols_to_save = []
# Open all files in one chain so one Snapshot and one report are produced per size-balanced Condor chunk.
input_chain = ROOT.TChain("Events")
for input_file in input_files:
    if input_chain.Add(input_file) == 0:
        raise RuntimeError(f"Could not add input ROOT file: {input_file}")
if args.n_events > 0:
    # RDataFrame.Range requires single-thread execution.
    if ROOT.IsImplicitMTEnabled():
        ROOT.DisableImplicitMT()
df = ROOT.RDataFrame(input_chain)
if args.n_events > 0:
    df = df.Range(args.n_events)
if args.n_events > 0:
    # ROOT's progress helper reports source-size totals even for a Range.
    print(f"[INFO] Test run: processing at most {args.n_events} input events.")
else:
    ROOT.RDF.Experimental.AddProgressBar(ROOT.RDF.AsRNode(df))

# useful definitions #
df = df.Define("period", f"static_cast<int>(Period::{config['era']})")
df = df.Define("is_data", "true" if is_data else "false")
df = df.Define("is_data_int", "1" if is_data else "0")
df = df.Define("is_signal", "true" if is_signal else "false")
dataset_crc = zlib.crc32(args.dataset_name.encode()) & 0xFFFF
input_crc = zlib.crc32("\n".join(input_files).encode()) & 0xFFFF
df = df.Define("FullEventId",f"eventId::encodeFullEventId({dataset_crc}, {input_crc}, rdfentry_)")
# default columns to store: #
cols_to_save.extend(utilities.GetObservablesCols("default", is_data, nano_version))

if not is_data:
    from analysis.mc_splitting import ApplyOrthogonalLumiFilter
    df, ortho_cols = ApplyOrthogonalLumiFilter(df, args.era, seed=12345, keep_tag_column=True)
    cols_to_save.extend(ortho_cols)
    from common.gen_vbf_filter import ApplyGenVBFFilter
    df,cols_to_save = ApplyGenVBFFilter(df,cols_to_save, args.era, args.dataset_name, process)

# define weights #
if not is_data:
    from corrections.general import define_base_weights
    xs_entry = dataset_cfg.get("crossSection", args.dataset_name)
    process_entry = process_cfg[process]
    df,base_weights,json_dict_to_store = define_base_weights(df, config.get("luminosity", ""), xs_entry, xs_cfg,config,dataset_cfg, process_entry,want_variations,systematics_cfg)
    cols_to_save.extend(base_weights)
else:
    from corrections.general import apply_golden_json
    df = apply_golden_json(df, config["lumiFile"])

# apply corrections --> this time also for data (e.g. JEC/ScaRe) #
from corrections.general import apply_corrections
from common.jet_horn_policy import horn_mitigation_enabled
config["apply_jet_horn_mitigation"] = horn_mitigation_enabled(args.era, sel_config)
df = apply_corrections(df, config, dataset_cfg, args.dataset_name, want_variations)

# MET FLAGS
if "MET_flags" in config and 'Flag_goodVertices' in df.GetColumnNames():
    from analysis.other import applyMETFlags
    df = applyMETFlags(df, config, is_data)
elif "MET_flags" in config:
    # FlashSim ntuples have no Flag_* branches: the skim runs without the MET filters.
    print("[WARNING] no Flag_* branches in the input: MET filters not applied")

## muon definitions ##
from analysis.muons import DefineMuonPtAndP4,ProcessMuonVariables,ApplyMuonTriggerMatching,ProcessExtraMuonVariables,ApplyElectronVeto,DefineMuonSelection
# definitions of p4
df = DefineMuonPtAndP4(df,want_variations)

# muon_error_columns = { "Muon_pt_err", "Muon_ptErr", "Muon_bsConstrainedPtErr", "Muon_bsConstrainedChi2"}
# available_muon_columns = {str(column) for column in df.GetColumnNames()}
# cols_to_save.extend(
#     column
#     for column in muon_error_columns
#     if column in available_muon_columns
# )

# trigger application && matching
df, trigger_cols = ApplyMuonTriggerMatching(df, trigger_config, sel_config.get("apply_trg_filter", True), want_variations, systematics_cfg)
cols_to_save.extend(trigger_cols)
# dimuon system definitions
muon_cols_initial = utilities.GetObservablesCols("Muon", is_data, nano_version)
from common.flashsim_patch import patch_flashsim_muon_columns,patch_flashsim_redefine_cols
muon_cols_initial = patch_flashsim_muon_columns(df, muon_cols_initial, args.dataset_name, is_data)
df, new_muon_cols =  ProcessMuonVariables(df,muon_cols_initial,muon_pt_default_suffix,trigger_config,want_variations,sel_config.get("muon_pt_min", 10.0),sel_config.get("dimuon_mass_min", 50.0),sel_config.get("dimuon_mass_max", 200.0),systematics_cfg,apply_trigger_filter=sel_config.get("apply_trg_filter", True))
cols_to_save.extend(new_muon_cols)
# electron veto
df = patch_flashsim_redefine_cols(df, args.dataset_name, is_data)
df = ApplyElectronVeto(df)
# # muon id/iso weights definitions
if not is_data:
    from corrections.mu import apply_muIDIso_weights
    df, mu_weights = apply_muIDIso_weights(df, config, want_variations)
    cols_to_save.extend(mu_weights)
# # extra lepton inclusion
# df, extra_lep_cols = ProcessExtraMuonVariables(df,muon_cols_initial,muon_pt_default_suffix,trigger_config,want_variations,sel_config.get("muon_pt_min", 10.0))
# cols_to_save.extend(extra_lep_cols)
# # muon selection
df,vars_sel = DefineMuonSelection(df,sel_config,want_variations,systematics_cfg)
cols_to_save.extend(vars_sel)

# ## jet definitions ##
from analysis.jets import ProcessAllJetVariables, SelectJetVars,SelectVBFJets
from corrections.btag_wpValues import getBTagWPValues
# define all varied variables for jets
bTagWPDict = getBTagWPValues(config)
jet_cols_initial = utilities.GetObservablesCols("Jet", is_data, nano_version)

df, jet_cols = ProcessAllJetVariables(df,jet_cols_initial,sel_config,config.get("bTagAlgo", "PNet"),bTagWPDict,want_variations,systematics_cfg)
cols_to_save.extend(jet_cols)
# apply jet veto map
from corrections.jetVetoMap import ApplyJetVetoMap
df, jet_veto_map_cols = ApplyJetVetoMap(df, config, muon_pt_default_suffix, False, sel_config.get("define_electron_cleaning", False),(nano_version == "v12"), want_variations,systematics_cfg)
cols_to_save.extend(jet_veto_map_cols)
# define selected jet vars
# Final analysis categories, including their shifted versions, are intentionally
# defined at histogram level by DefineSelections.  The skim stores only
# the nominal/shifted object and selection columns needed to build them.

## additional col to store ##
collections = ["SoftActivityJet"]
if not is_data:
    collections.append('GenJet')
    df = df.Define("GenJet_idx", "CreateIndexes(GenJet_pt.size())")
    cols_to_save.append("GenJet_idx")
    # if "Jet_genJetIdx" in input_columns:
    #     cols_to_save.append("Jet_genJetIdx")
    if "LHE_Vpt" in df.GetColumnNames():
        collections.append("LHE")
    if "LHEScaleWeight" in df.GetColumnNames():
        collections.append("LHEWeight")
for c in collections:
    cols_to_save.extend(utilities.GetObservablesCols(c, is_data, nano_version))
cols_to_save = list(set(cols_to_save))

df, selected_jet_cols = SelectJetVars(df,jet_cols_initial,sel_config,config.get("bTagAlgo", "PNet"),bTagWPDict,want_variations,systematics_cfg)
cols_to_save.extend(selected_jet_cols)
df,vbf_jet_cols = SelectVBFJets(df,want_variations,systematics_cfg)
cols_to_save.extend(vbf_jet_cols)


if args.sync:
    from common.sync_skim import prepare_sync, export_sync
    df, sync_cols = prepare_sync(df, sel_config, is_data, args.sync_category)
    cols_to_save.extend(sync_cols)
    cols_to_save = list(dict.fromkeys(cols_to_save))
    snapshot_options = ROOT.RDF.RSnapshotOptions()
    snapshot_options.fLazy = True
    sync_snapshot = df.Snapshot("Events", args.output_file, utilities.ListToVector(cols_to_save), snapshot_options)
    sync_report = df.Report()
    export_sync(df, sel_config, is_data, args.sync_category, Path(args.output_file).parent, {
        "era": args.era, "dataset": args.dataset_name, "is_data": is_data,
        "input_files": input_files, "input_entries": int(input_chain.GetEntries()),
        "n_events_limit": args.n_events, "nominal_only": True,
    })
cols_to_save = list(dict.fromkeys(cols_to_save))

# The FlashSim Sep2026 ntuples store the collection counters as ULong_t: for
# arrays read straight from the input Snapshot cannot build the leaflist and
# drops them with only a warning.  Redefined as RVec they are written normally.
defined_columns = {str(c) for c in df.GetDefinedColumnNames()}
for col in cols_to_save:
    leaf = input_chain.GetLeaf(col) if col not in defined_columns else None
    count = leaf.GetLeafCount() if leaf else None
    if count and count.GetTypeName() not in ("Int_t", "UInt_t"):
        value_type = "float" if leaf.GetTypeName() == "Float16_t" else leaf.GetTypeName()
        df = df.Redefine(col, f"ROOT::RVec<{value_type}>({col}.begin(), {col}.end())")

## snapshot + report store ##
if args.sync:
    sync_snapshot.GetValue()
    report = sync_report.GetValue()
else:
    df.Snapshot("Events",args.output_file,utilities.ListToVector(cols_to_save))
    report = df.Report().GetValue()
df, report_json = utilities.SaveReport(df, report, verbose=0)
if not is_data:
    for pu_key,pu_dict in json_dict_to_store.items():
        for xs_key,xs_dict in pu_dict.items():
            value_to_extract = xs_dict['value']
            xs_dict['value']=value_to_extract.GetValue()
            if 'value_unsigned' in xs_dict.keys():
                xs_dict['value_unsigned']=xs_dict['value_unsigned'].GetValue()
    report_json.update(json_dict_to_store)

json_file = (
    args.report_file
    if args.report_file is not None
    else os.path.splitext(args.output_file)[0] + "_report.json"
)
with open(json_file, "w") as f:
    json.dump(report_json, f, indent=4)
out_tfile = ROOT.TFile.Open(args.output_file, "UPDATE")
out_tfile.Close()
# df.Report().Print()
