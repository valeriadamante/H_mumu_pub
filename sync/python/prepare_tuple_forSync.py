#!/usr/bin/env python3
"""Sync tuple of one dataset (or of the data streams of one run) from skim_v4.

Loose selection by default: every skim event with 50 < m_mumu < 200 GeV, with
the region/category/muon flags of config/Run3_2024/selections_sync.yaml as
branches. Weights and derived observables come from common.prepare_rdf, as in
the histograms. Next to the tuple: the cutflow, the region/category yields and
<output stem>_metadata.json (datasets, luminosity and cross sections, corrections,
selections and the pT/correction conventions), also stored in the ROOT file as
the TNamed "metadata".
"""
import argparse
import copy
import csv
import datetime
import glob
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import ROOT
import yaml

ROOT.gROOT.SetBatch(True)
MISSING = -1000.0
REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))

os.environ.setdefault("ANALYSIS_PATH", str(REPOSITORY))
from common.add_vars import DefineSelections  # noqa: E402

ERA = "Run3_2024"
CONFIG_DIR = REPOSITORY / "config" / ERA
SYNC_CONFIG_PATH = CONFIG_DIR / "selections_sync.yaml"
with SYNC_CONFIG_PATH.open() as stream:
    SYNC_CONFIG = yaml.safe_load(stream)
REGIONS = SYNC_CONFIG["masses_regions"]
CATEGORIES = SYNC_CONFIG["categories"]
# Systematic variations of the skim (JER/JES jets, muon FSR scale/resolution,
# up/down weights): the sync tuple is nominal only.
VARIATION = re.compile(r"(_(up|down|Up|Down)$|(JER|JES)\w*(up|down)$)")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="ROOT file(s), glob(s), or dataset directory")
    parser.add_argument("-o", "--output", type=Path,
                        help="output ROOT file (required unless --cutflow-only)")
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--mc", action="store_true", help="MC input: also save the event weights")
    parser.add_argument("--dataset", help="dataset name (default: directory of the input files); "
                        "MC needs exactly one dataset for its normalisation")
    parser.add_argument("--selection", default="mass_50_200", help="output filter; use 1 for none")
    parser.add_argument("--no-framework", action="store_true",
                        help="skip common.prepare_rdf: no weights and no framework observables (faster)")
    parser.add_argument("--no-raw-jet-pt", action="store_true",
                        help="do not add the raw jet pT read from the NanoAOD (add_raw_jet_pt.py)")
    parser.add_argument("--nano-workers", type=int, default=24,
                        help="parallel NanoAOD reads for the raw jet pT")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--cutflow-output", type=Path,
                        help="cutflow CSV (default: <output stem>_cutflow.csv)")
    parser.add_argument("--no-cutflow", action="store_true",
                        help="skip report JSON aggregation and cutflow output")
    parser.add_argument("--cutflow-only", type=Path, metavar="TUPLE",
                        help="read an existing tuple and produce only its cutflow CSV")
    parser.add_argument("--tuple-preselection", default="mass_inclusive",
                        help="cut already applied to --cutflow-only tuple; use 1 for none")
    parser.add_argument("--cutflow-key", default="cutflow",
                        help="YAML key holding the cumulative cutflow step list")
    parser.add_argument("--vbf-pair-strategy", choices=("standard", "leading"),
                        default="standard", help=argparse.SUPPRESS)
    return parser.parse_args()


def find_files(items):
    found = []
    for item in items:
        path = Path(item).expanduser()
        matches = (sorted(path.rglob("*.root")) if path.is_dir() else
                   [path] if path.is_file() else [Path(x) for x in sorted(glob.glob(item))])
        found += [str(x.resolve()) for x in matches if x.suffix == ".root"]
    found = list(dict.fromkeys(found))
    if not found:
        raise ValueError("no ROOT input files found")
    return found


def report_for_root(root_file):
    path = Path(root_file)
    if path.stem.startswith("skim_"):
        report = path.with_name(f"report_{path.stem[len('skim_'):]}.json")
    else:
        report = path.with_name(f"{path.stem}_report.json")
    return report if report.is_file() else None


def aggregate_reports(files):
    totals = {}
    order = []
    for root_file in files:
        report_path = report_for_root(root_file)
        if report_path is None:
            raise FileNotFoundError(f"missing report JSON for {root_file}")
        with report_path.open() as stream:
            report = json.load(stream)
        steps = cutflow_report_steps(report)
        if not order:
            order = list(steps)
        elif list(steps) != order:
            raise RuntimeError(f"inconsistent cutflow keys in {report_path}")
        for name, passed in steps.items():
            totals[name] = totals.get(name, 0) + int(passed)
    return [(name, totals[name]) for name in order]


def cutflow_report_steps(report):
    """Event counts of a skim report; the sums of weights (MC) are not steps."""
    return {
        name: value if isinstance(value, (int, float)) else value["pass"]
        for name, value in report.items()
        if isinstance(value, (int, float)) or "pass" in value
    }


def reports_from_inputs(items):
    reports = []
    for item in items:
        path = Path(item).expanduser()
        if path.is_dir():
            reports.extend(sorted(path.rglob("report_*.json")))
        elif path.suffix == ".json" and path.is_file():
            reports.append(path)
    reports = list(dict.fromkeys(path.resolve() for path in reports))
    if not reports:
        raise ValueError("no report JSON files found")
    return reports


def aggregate_report_paths(report_paths):
    totals = {}
    order = []
    for report_path in report_paths:
        with report_path.open() as stream:
            report = json.load(stream)
        steps = cutflow_report_steps(report)
        if not order:
            order = list(steps)
        elif list(steps) != order:
            raise RuntimeError(f"inconsistent cutflow keys in {report_path}")
        for name, passed in steps.items():
            totals[name] = totals.get(name, 0) + int(passed)
    return [(name, totals[name]) for name in order]


def make_cutflow_only(args):
    report_counts = aggregate_report_paths(reports_from_inputs(args.inputs))
    if not args.cutflow_only.is_file():
        raise FileNotFoundError(f"tuple does not exist: {args.cutflow_only}")
    df = ROOT.RDataFrame(args.tree, str(args.cutflow_only))
    available = columns(df)
    configured = {
        name: definition
        for section in ("masses_regions", "muons_selection", "jets_selection", "categories")
        for name, definition in SYNC_CONFIG.get(section, {}).items()
    }
    defining = set()

    def ensure_column(node, name):
        if name in available:
            return node
        if name not in configured:
            raise RuntimeError(f"unknown cutflow selection: {name}")
        if name in defining:
            raise RuntimeError(f"cyclic selection dependency involving {name}")
        defining.add(name)
        expression = configured[name]["expression"].format(
            tot_suff="", mu_suff="", jet_suff=""
        )
        for dependency in configured:
            if dependency != name and re.search(rf"\b{re.escape(dependency)}\b", expression):
                node = ensure_column(node, dependency)
        node = node.Define(name, expression)
        available.add(name)
        defining.remove(name)
        return node

    category_config = SYNC_CONFIG.get("categories", {})
    configured_order = [
        name for name in cutflow_steps(args.cutflow_key)
        if name not in category_config or category_config[name].get("store", False)
    ]
    preselection = args.tuple_preselection
    order = ([preselection] if preselection != "1" else []) + [
        name for name in configured_order if name != preselection
    ]
    count_actions = []
    node = df
    for selection in order:
        node = ensure_column(node, selection)
        node = node.Filter(selection, f"cutflow_{selection}")
        count_actions.append((selection, node.Count()))
    total_action = df.Count()
    yield_actions = [
        ((kind, name), df.Filter(name, f"yield_{kind}_{name}").Count())
        for kind, name in stored_category_region_names()
    ]
    ROOT.RDF.RunGraphs([
        *(action for _, action in count_actions), total_action,
        *(action for _, action in yield_actions),
    ])
    selection_counts = [(name, int(action.GetValue())) for name, action in count_actions]
    path = args.cutflow_output or args.cutflow_only.with_name(
        f"{args.cutflow_only.stem}_cutflow.csv"
    )
    write_cutflow(path, report_counts, selection_counts)
    write_selection_yields(
        yields_path_from_cutflow(path), int(total_action.GetValue()),
        [(key, int(action.GetValue())) for key, action in yield_actions],
    )


def cutflow_steps(key):
    if key not in SYNC_CONFIG:
        raise KeyError(f"{SYNC_CONFIG_PATH} has no cutflow list named {key}")
    return SYNC_CONFIG[key]


def write_cutflow(path, report_counts, selection_counts):
    counts = report_counts + selection_counts
    if not counts:
        raise RuntimeError("empty cutflow")
    initial = counts[0][1]
    previous = initial
    rows = []
    for step, (name, passed) in enumerate(counts):
        rows.append({
            "step": step,
            "selection": name,
            "cumulative_pass": passed,
            "rejected_at_step": previous - passed,
            "relative_efficiency": passed / previous if previous else 0.0,
            "cumulative_efficiency": passed / initial if initial else 0.0,
        })
        previous = passed
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[sync] wrote cutflow {path}")


def stored_category_region_names():
    return [
        ("region" if section == "masses_regions" else "category", name)
        for section in ("masses_regions", "categories")
        for name, definition in SYNC_CONFIG.get(section, {}).items()
        if definition.get("store", False)
    ]


def write_selection_yields(path, total, counts):
    rows = [{
        "kind": kind,
        "selection": name,
        "pass": passed,
        "total_events": total,
        "efficiency": passed / total if total else 0.0,
    } for (kind, name), passed in counts]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[sync] wrote category/region yields {path}")


def yields_path_from_cutflow(path):
    stem = path.stem
    if stem.endswith("_cutflow"):
        stem = stem[:-len("_cutflow")]
    return path.with_name(f"{stem}_categories_regions.csv")


def columns(df):
    return {str(x) for x in df.GetColumnNames()}


def put(df, name, expression):
    return df.Redefine(name, expression) if name in columns(df) else df.Define(name, expression)


def require(df, needed, label):
    missing = sorted(set(needed) - columns(df))
    if missing:
        raise RuntimeError(f"cannot define {label}; missing: {', '.join(missing)}")


VBF_STAGES_HELPER = """
#include <Math/Vector4D.h>
#include <cmath>
namespace sync_vbf {
using P4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
// Stages of FindVBFJets in the order the C++ selection applies them, so a
// cutflow can attribute the VBF loss to the horn cleaning, the jet pT, the
// dijet mass or the pseudorapidity gap. Stage 3 reproduces HasVBF exactly.
template <typename V, typename B>
bool PairStage(const V &pt, const V &eta, const V &phi, const V &mass,
               const B &outside, int stage) {
  for (size_t i = 0; i < pt.size(); ++i) {
    if (!outside[i]) continue;
    for (size_t j = i + 1; j < pt.size(); ++j) {
      if (!outside[j]) continue;
      if (stage == 0) return true;
      if (pt[i] < 35. || pt[j] < 25.) continue;
      if (stage == 1) return true;
      P4 first(pt[i], eta[i], phi[i], mass[i]);
      P4 second(pt[j], eta[j], phi[j], mass[j]);
      if ((first + second).M() < 400.) continue;
      if (stage == 2) return true;
      if (std::abs(eta[i] - eta[j]) < 2.5) continue;
      return true;
    }
  }
  return false;
}
}
"""


def define_vbf_stages(df):
    """Intermediate VBF requirements, only used as cutflow steps."""
    require(df, ["SelectedJet_pt", "SelectedJet_eta", "SelectedJet_phi",
                 "SelectedJet_mass", "SelectedJet_IsOutsideHorn"], "VBF stages")
    ROOT.gInterpreter.Declare(VBF_STAGES_HELPER)
    arguments = ("SelectedJet_pt, SelectedJet_eta, SelectedJet_phi, "
                 "SelectedJet_mass, SelectedJet_IsOutsideHorn")
    stages = {"vbf_jets_ge2": 0, "vbf_pair_pt": 1, "vbf_pair_mjj": 2,
              "vbf_pair_deta": 3}
    df = put(df, "jets_ge2", "N_SelectedJets >= 2")
    for name, stage in stages.items():
        df = put(df, name, f"sync_vbf::PairStage({arguments}, {stage})")
    return df


def define_leading_vbf_pair(df):
    """Choose the pT-leading pair among those passing mjj/deta thresholds."""
    expression = r"""
    ROOT::VecOps::RVec<int> result{-1, -1};
    for (size_t i = 0; i < SelectedJet_pt.size(); ++i) {
        for (size_t j = i + 1; j < SelectedJet_pt.size(); ++j) {
            SyncP4 ji(SelectedJet_pt[i], SelectedJet_eta[i],
                      SelectedJet_phi[i], SelectedJet_mass[i]);
            SyncP4 jj(SelectedJet_pt[j], SelectedJet_eta[j],
                      SelectedJet_phi[j], SelectedJet_mass[j]);
            if ((ji + jj).M() > 400. &&
                std::abs(SelectedJet_eta[i] - SelectedJet_eta[j]) > 2.5) {
                result[0] = static_cast<int>(i);
                result[1] = static_cast<int>(j);
                return result;
            }
        }
    }
    return result;
    """
    df = put(df, "VBFLeadingPairIdx", expression)
    df = put(df, "HasVBF", "VBFLeadingPairIdx[0] >= 0")
    df = put(df, "VBFJetIdx_1", "VBFLeadingPairIdx[0]")
    df = put(df, "VBFJetIdx_2", "VBFLeadingPairIdx[1]")
    return df


def prepare(df, vbf_pair_strategy="standard", framework=None):
    """Sync columns on a skim_v4 node.

    ``framework`` (keyword arguments of common.prepare_rdf.prepare_rdf, without
    the RDF and the selections) builds the weights and the analysis observables
    with the framework; without it only the selections are defined, which is
    what special_sync_skim.py needs on its NanoAOD-level node.
    """
    aliases = {"dimuon_mass": "m_mumu", "dimuon_pt": "pt_mumu",
               "dimuon_eta": "eta_mumu", "dimuon_phi": "phi_mumu"}
    available = columns(df)
    require(df, ["m_mumu"], "dimuon mass")
    missing_kinematics = {"pt_mumu", "eta_mumu", "phi_mumu"} - available
    if missing_kinematics:
        require(
            df,
            ["mu1_pt", "mu1_eta", "mu1_phi", "mu1_mass",
             "mu2_pt", "mu2_eta", "mu2_phi", "mu2_mass"],
            "dimuon kinematics",
        )
        ROOT.gInterpreter.Declare(
            "#include <Math/Vector4D.h>\n"
            "using SyncP4 = ROOT::Math::LorentzVector<"
            "ROOT::Math::PtEtaPhiM4D<double>>;"
        )
        df = put(
            df,
            "sync_dimuon_p4",
            "SyncP4(mu1_pt,mu1_eta,mu1_phi,mu1_mass) + "
            "SyncP4(mu2_pt,mu2_eta,mu2_phi,mu2_mass)",
        )
        fallback = {
            "pt_mumu": "static_cast<float>(sync_dimuon_p4.Pt())",
            "eta_mumu": "static_cast<float>(sync_dimuon_p4.Eta())",
            "phi_mumu": "static_cast<float>(sync_dimuon_p4.Phi())",
        }
        for name in sorted(missing_kinematics):
            df = df.Define(name, fallback[name])
    for target, source in aliases.items():
        df = put(df, target, source)

    df = define_vbf_stages(df)

    selection_config = copy.deepcopy(SYNC_CONFIG)
    if vbf_pair_strategy == "leading":
        ROOT.gInterpreter.Declare(
            "#include <Math/Vector4D.h>\n"
            "using SyncP4 = ROOT::Math::LorentzVector<"
            "ROOT::Math::PtEtaPhiM4D<double>>;"
        )
        df = define_leading_vbf_pair(df)
        selection_config["categories"]["VBF_def"]["expression"] = (
            "HasVBF{jet_suff}"
        )

    if framework is None:
        # Same selection builder as prepare_rdf, without its weights and observables.
        df = DefineSelections(df, selection_config)
    else:
        from common.prepare_rdf import prepare_rdf
        from common.utilities import initialize_root_runtime
        initialize_root_runtime()
        df = prepare_rdf(df, selections_cfg=selection_config, **framework)["inclusive"]
        if "mu1_p4_noCorr" in columns(df):
            # m_mumu_noCorr is defined by the framework; its pT completes the pair.
            df = put(df, "pt_mumu_noCorr", "static_cast<float>((mu1_p4_noCorr + mu2_p4_noCorr).Pt())")

    # Expose every stored SelectedJet quantity also for the first three jets.
    available = columns(df)
    vectors = sorted(x for x in available if x.startswith("SelectedJet_")
                     and "_p4" not in x and not VARIATION.search(x))
    for index, prefix in enumerate(("leadingjet", "subleadingjet", "thirdleadingjet")):
        for source in vectors:
            target = f"{prefix}_{source[len('SelectedJet_') :]}"
            if target not in available:
                df = df.Define(target, f"{source}.size() > {index} ? static_cast<float>({source}.at({index})) : {MISSING}f")
                available.add(target)

    require(df, ["HasVBF", "VBFJetIdx_1", "VBFJetIdx_2", "SelectedJet_pt",
                 "SelectedJet_eta", "SelectedJet_phi", "SelectedJet_mass"], "dijet variables")
    available = columns(df)
    for leg, index_name in ((1, "VBFJetIdx_1"), (2, "VBFJetIdx_2")):
        for source in vectors:
            target = f"vbfjet{leg}_{source[len('SelectedJet_') :]}"
            if target in available:
                continue
            expression = (
                f"HasVBF && {index_name} >= 0 && "
                f"static_cast<size_t>({index_name}) < {source}.size() ? "
                f"static_cast<float>({source}.at({index_name})) : {MISSING}f"
            )
            df = df.Define(target, expression)
            available.add(target)

    ROOT.gInterpreter.Declare("#include <Math/Vector4D.h>\nusing SyncP4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;")
    valid_vbf = "HasVBF && VBFJetIdx_1 >= 0 && VBFJetIdx_2 >= 0 && static_cast<size_t>(VBFJetIdx_1) < SelectedJet_pt.size() && static_cast<size_t>(VBFJetIdx_2) < SelectedJet_pt.size()"
    valid_ls = "SelectedJet_pt.size() >= 2 && SelectedJet_eta.size() >= 2 && SelectedJet_phi.size() >= 2 && SelectedJet_mass.size() >= 2"
    p4 = lambda i: f"SyncP4(SelectedJet_pt.at({i}),SelectedJet_eta.at({i}),SelectedJet_phi.at({i}),SelectedJet_mass.at({i}))"
    df = put(
        df, "sync_dijet_vbf_p4",
        f"{valid_vbf} ? {p4('VBFJetIdx_1')}+{p4('VBFJetIdx_2')} : SyncP4(0.,0.,0.,0.)",
    )
    df = put(
        df, "sync_dijet_ls_p4",
        f"{valid_ls} ? {p4('0')}+{p4('1')} : SyncP4(0.,0.,0.,0.)",
    )
    for name, method in (("mass", "M"), ("pt", "Pt"), ("eta", "Eta"), ("phi", "Phi")):
        df = put(df, f"dijet_vbf_{name}", f"{valid_vbf} ? static_cast<float>(sync_dijet_vbf_p4.{method}()) : {MISSING}f")
        df = put(df, f"dijet_ls_{name}", f"{valid_ls} ? static_cast<float>(sync_dijet_ls_p4.{method}()) : {MISSING}f")
    return df


# Event-level columns of the full sync tuple, when they exist.
EVENT_COLUMNS = (
    "HLT_IsoMu24", "Event_HasTriggerMatching_singleMu", "is_data", "period",
    "OrthogonalEraTag", "nJet", "N_SelectedJets", "JetTagSel", "SelectedJetTagSel",
    "m_mumu", "pt_mumu", "eta_mumu", "phi_mumu", "y_mumu", "dR_mumu",
    "m_mumu_noCorr", "pt_mumu_noCorr",
    "m_jj", "delta_eta_jj", "delta_phi_jj", "pt_vbfj1j2", "m_jj_ls", "delta_eta_jj_ls",
)


def selected_columns(df, is_mc, extended=False):
    """Columns to write. ``extended`` adds the event columns and the jet
    collections (all jets and selected jets) of the loose sync tuple."""
    available = columns(df)
    selection_sections = (
        "masses_regions", "muons_selection", "jets_selection", "categories"
    )
    configured_selections = {
        name
        for section in selection_sections
        for name in SYNC_CONFIG.get(section, {})
    }
    stored_selections = {
        name
        for section in selection_sections
        for name, definition in SYNC_CONFIG.get(section, {}).items()
        if definition.get("store", False)
    }
    exact = {"run", "luminosityBlock", "event", "HasVBF", "N_SelectedJets",
             "VBFJetIdx_1", "VBFJetIdx_2", *stored_selections,
             "dimuon_mass", "dimuon_pt", "dimuon_eta", "dimuon_phi",
             "dijet_ls_mass", "dijet_ls_pt", "dijet_ls_eta", "dijet_ls_phi",
             "dijet_vbf_mass", "dijet_vbf_pt", "dijet_vbf_eta", "dijet_vbf_phi"}
    prefixes = ("mu1_", "mu2_", "leadingjet_", "subleadingjet_",
                "thirdjet_", "thirdleadingjet_", "vbfjet1_", "vbfjet2_")
    result = {
        x for x in available
        if x in exact or (x.startswith(prefixes) and x not in configured_selections)
    }
    if is_mc:
        result |= {x for x in available if x.startswith("weight_") or x.startswith("weight__")}
        result |= {x for x in ("genWeight", "puWeight", "Generator_weight", "LHEWeight_originalXWGTUP",
                               "inv_N_orig") if x in available}
    if extended:
        result |= {x for x in EVENT_COLUMNS if x in available}
        result |= {x for x in available if x.startswith(("Jet_", "SelectedJet_")) and "_p4" not in x}
        result = {x for x in result if not VARIATION.search(x) and "_p4" not in x}
    return sorted(result)


# Point (5) of the sync request: which pT enters which step. Checked against
# analysis/muons.py, analysis/jets.py, corrections/mu.py, common/trigger_weights.py
# and config/Run3_2024/{selections,triggers,systematics}.yaml.
CONVENTIONS = {
    "muon_pt": {
        "mu{i}_pt_raw_noCorr": "NanoAOD pT, no ScaRe, no FSR: Muon_bsConstrainedPt if "
                               "Muon_bsConstrainedChi2 < 30, otherwise Muon_pt",
        "mu{i}_pt_raw_scale": "same choice, ScaRe scale correction only",
        "mu{i}_pt_raw_corr": "same choice, ScaRe scale + resolution smearing (MC)",
        "mu{i}_pt": "nominal: same choice + ScaRe + FSR photon recovery (Muon_pt_FSR_corr); "
                    "mu{i}_eta/phi/mass are those of the same four-vector",
        "mu{i}_eta_raw_noCorr, mu{i}_phi_raw_noCorr": "NanoAOD angles",
        "ordering": "mu1/mu2 ordered by the nominal (ScaRe + FSR) pT",
    },
    "dimuon": {
        "m_mumu, pt_mumu (dimuon_mass, dimuon_pt)": "from the nominal (ScaRe + FSR) muons",
        "m_mumu_noCorr, pt_mumu_noCorr": "from mu{i}_pt_raw_noCorr and the NanoAOD angles",
    },
    "cuts_on_raw_pt": [
        "skim good muons: mu pT_raw_noCorr > 10 GeV, |eta_NanoAOD| < 2.4, mediumId, pfIsoId >= 2",
        "trigger matching: HLT_IsoMu24, offline mu pT_raw_noCorr > 26 GeV, trigger object id 13 "
        "with filter bit 8, DeltaR(raw muon, object) < 0.4",
        "muon_pre_sel_check: pT_raw_noCorr > 10 GeV, |eta_raw_noCorr| < 2.4",
        "sync baseline (base_sel_muons, baseline, VBF, ggF): mu1_rawpt_g26 = mu1_pt_raw_noCorr > 26, "
        "mu2_rawpt_g20 = mu2_pt_raw_noCorr > 20",
    ],
    "cuts_on_corrected_pt": [
        "skim dimuon mass window 50 < m_mumu < 200 GeV",
        "all mass regions",
        "flags mu1_pt_g26, mu2_pt_g26, mu2_pt_g20 (corrected-pT version of the cuts, for comparison only)",
    ],
    "scale_factors": {
        "muon ID/iso/trigger SF inputs": "pT = mu{i}_pt_raw_noCorr, eta = mu{i}_eta (nominal four-vector)",
        "ID/iso SF choice": "mu{i}_pt (corrected) < 30 GeV: TightID + TightPFIso, otherwise "
                            "MediumID + LoosePFIso (systematics.yaml weights.Central)",
        "trigger SF": "IsoMu24 SF of the muon with the higher corrected pT, if it is trigger "
                      "matched, else 1 (weight_TrgSF_singleMu_IsoMu24Central)",
        "pileup": "Pileup_nTrueInt, LUM correctionlib payload (weight_pu)",
    },
    "jets": {
        "Jet_pt, SelectedJet_pt": "JEC reapplied (+ JER smearing in MC)",
        "Jet_pt_nocorr, SelectedJet_pt_nocorr": "NanoAOD jet pT (NanoAOD JEC), before our JEC "
                                                "reapplication and JER",
        "Jet_pt_raw, SelectedJet_pt_raw, *_pt_raw": "raw jet pT (no JEC, no JER) = Jet_pt_nocorr * "
                                                   "(1 - Jet_rawFactor); Jet_rawFactor is not in the skim "
                                                   "and is read from the NanoAOD (add_raw_jet_pt.py)",
        "Jet_passJetIdTight / Jet_passJetIdTightLepVeto": "tight jet ID recomputed from the jet "
                                                          "energy fractions (JetIdNewDef, NanoAOD v15 recipe)",
        "Jet_vetoMap": "jet in the JME veto map (jets in it are removed, no event veto)",
        "Jet_IsInsideHorn": "2.5 <= |eta| < 3.0 and pT <= 50 GeV (removed from the selected jets)",
        "selected jets": "corrected pT > 25 GeV, |eta| < 4.7, tight ID, not in veto map, outside horn, "
                         "DeltaR(jet, mu1/mu2) > 0.4 (nominal muons); ordered by corrected pT",
        "b tagging": "UParTAK4 B score, loose/medium WPs, only |eta| < 2.5; SelectedJetTagSel: "
                     "no medium and fewer than two loose b-tagged selected jets",
        "leadingjet_/subleadingjet_/thirdjet_": "framework definitions, -1000 when absent",
        "VBF pair (HasVBF, VBFJetIdx_1/2, vbfjet1/2_, dijet_vbf_)": "pair of selected jets outside the "
                                                                    "horn with pT > 35/25, m_jj > 400, |Delta eta| > 2.5",
        "dijet_ls_": "two leading selected jets",
    },
    "missing_value": MISSING,
}


def git_state():
    def run(*command):
        return subprocess.run(command, cwd=REPOSITORY, capture_output=True, text=True).stdout.strip()
    return {"commit": run("git", "rev-parse", "HEAD"),
            "modified_files": len(run("git", "status", "--porcelain").splitlines())}


def build_metadata(args, files, datasets, dataset, keep, value_ranges, framework):
    """Points (1)-(5) of the sync request, as one JSON document."""
    from corrections.btag_wpValues import getBTagWPValues
    from sync.python.collect_corrections import build_inventory

    def load(path):
        with path.open() as stream:
            return yaml.safe_load(stream)

    samples = load(CONFIG_DIR / "samples.yaml")
    cross_sections = load(REPOSITORY / "config" / "crossSections13p6TeV.yaml")
    maincfg = load(CONFIG_DIR / "maincfg.yaml")
    skim_selections = load(CONFIG_DIR / "selections.yaml")
    systematics = load(CONFIG_DIR / "systematics.yaml")
    inventory = build_inventory(REPOSITORY, [ERA])

    dataset_info = []
    for name in datasets:
        entry = samples.get(name, {})
        info = {"name": name, "nanoAOD": entry.get("nanoAOD"),
                "n_skim_files": sum(Path(f).parent.name == name for f in files)}
        if entry.get("crossSection"):
            info["cross_section_entry"] = entry["crossSection"]
            info["cross_section"] = cross_sections.get(entry["crossSection"])
        dataset_info.append(info)

    normalisation = {"integrated_luminosity_pb": maincfg["luminosity"],
                     "golden_json": maincfg.get("lumiFile")}
    if args.mc:
        normalisation.update({
            "event_weight": "weight__Central = (" + " ".join(
                systematics["weights"]["Central"]["expression"].split()) + ") * inv_N_orig",
            "inv_N_orig": "1 / sum(genWeight) over the whole dataset (skim report_*.json, "
                          "before any selection)",
            "orthogonal_lumi_split": "2024 MC events are split among the 2024-2026 eras "
                                     "(OrthogonalEraTag); weight_lumi is the luminosity of this era",
            "dy_reweights": "DY only: weight_dy_* multiply weight__Central (common/apply_custom_weights.py)",
            "values_in_tuple": value_ranges,
        })
    return {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "repository": git_state(),
        "era": ERA,
        "is_mc": bool(args.mc),
        "dataset_label": dataset,
        "datasets": dataset_info,                                        # (1)
        "input_files": files,
        "normalisation": normalisation,                                  # (2)
        "corrections": {                                                 # (3)
            "chain": inventory["correction_chain"],
            "payloads": inventory["eras"][ERA],
            "btag_wp_values": getBTagWPValues(maincfg),
        },
        "selections": {                                                  # (4)
            "skim": {
                "muon_pt_min_raw_GeV": skim_selections["muon_pt_min"],
                "dimuon_mass_window_GeV": [skim_selections["dimuon_mass_min"],
                                           skim_selections["dimuon_mass_max"]],
                "trigger_filter": skim_selections["apply_trg_filter"],
                "jet_pt_min_GeV": skim_selections["jet_pt_min"],
                "jet_eta_max": skim_selections["jet_eta_max"],
                "horn_veto": skim_selections["jet_horn_veto_expr"],
                "event_filters": ["golden JSON (data)", "MET filters", "single-muon trigger matching",
                                  "exactly two good muons", "dimuon mass window", "electron veto"],
            },
            "tuple_filter": args.selection,
            "flag_definitions_file": str(SYNC_CONFIG_PATH.relative_to(REPOSITORY)),
            "flags": {section: {name: d["expression"].format(tot_suff="", mu_suff="", jet_suff="")
                                for name, d in SYNC_CONFIG.get(section, {}).items()}
                      for section in ("masses_regions", "muons_selection", "jets_selection", "categories")},
            "flags_stored": [name for _, name in stored_category_region_names()] + [
                name for section in ("muons_selection", "jets_selection")
                for name, d in SYNC_CONFIG.get(section, {}).items() if d.get("store", False)],
        },
        "conventions": CONVENTIONS,                                      # (5)
        "framework_observables": framework is not None,
        "branches": keep,
    }


def write_metadata(output, metadata):
    path = output.with_name(f"{output.stem}_metadata.json")
    text = json.dumps(metadata, indent=2, default=str)
    path.write_text(text)
    handle = ROOT.TFile.Open(str(output), "UPDATE")
    if not handle or handle.IsZombie():
        raise OSError(f"cannot reopen {output} to store the metadata")
    ROOT.TNamed("metadata", text).Write()
    handle.Close()
    print(f"[sync] wrote metadata {path} (also TNamed 'metadata' in the tuple)")


def framework_options(args, files, datasets, dataset):
    """Keyword arguments of prepare_rdf for this input."""
    from common.utilities import get_config, get_segmentation_dict
    options = dict(dataset_name=dataset, era=ERA, is_data=not args.mc, skip_validation=True,
                   systematics_cfg=get_config(str(CONFIG_DIR / "systematics.yaml")))
    if args.mc:
        if len(datasets) != 1:
            raise ValueError(f"MC input must be one dataset, got {datasets}")
        directory = Path(files[0]).parent
        # The normalisation is the whole dataset's, also when only some files are read.
        skims = sorted(directory.glob("skim_*.root"))
        reports = sorted(directory.glob("report_*.json"))
        if len(reports) != len(skims):
            raise RuntimeError(f"{directory}: {len(skims)} skim files but {len(reports)} reports")
        options["seg_dict"] = get_segmentation_dict([str(r) for r in reports])
    return options


def main():
    args = arguments()
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    ROOT.EnableImplicitMT(args.threads)
    if args.cutflow_only:
        make_cutflow_only(args)
        return
    if args.output is None:
        raise ValueError("--output is required unless --cutflow-only is used")
    if args.output.exists() and not (args.overwrite or args.dry_run):
        raise FileExistsError(f"{args.output} exists; pass --overwrite")
    files = find_files(args.inputs)
    datasets = sorted({Path(f).parent.name for f in files})
    dataset = args.dataset or "+".join(datasets)
    report_counts = [] if args.no_cutflow else aggregate_reports(files)
    print(f"[sync] reading {dataset} from {len(files)} file(s)")
    framework = None if args.no_framework else framework_options(args, files, datasets, dataset)
    df = prepare(ROOT.RDataFrame(args.tree, files), args.vbf_pair_strategy, framework)
    keep = selected_columns(df, args.mc, extended=True)
    print(f"[sync] selected {len(keep)} columns" + (" including MC weights" if args.mc else ""))
    if args.dry_run:
        print("\n".join(keep))
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    opts = ROOT.RDF.RSnapshotOptions()
    opts.fMode = "RECREATE"
    opts.fLazy = True
    snapshot = df.Filter(args.selection, "sync selection").Snapshot(
        args.tree, str(args.output), keep, opts
    )
    cut_node = df
    count_actions = []
    for selection in cutflow_steps(args.cutflow_key):
        if selection not in columns(cut_node):
            raise RuntimeError(f"unknown cutflow selection: {selection}")
        cut_node = cut_node.Filter(selection, f"cutflow_{selection}")
        count_actions.append((selection, cut_node.Count()))
    total_action = df.Count()
    yield_actions = [
        ((kind, name), df.Filter(name, f"yield_{kind}_{name}").Count())
        for kind, name in stored_category_region_names()
    ]
    # Normalisation factors as stored per event, to report them in the metadata.
    range_columns = ("weight_lumi", "weight_xs", "inv_N_orig") if args.mc else ()
    range_actions = {name: (df.Min(name), df.Max(name)) for name in range_columns
                     if name in columns(df)}
    ROOT.RDF.RunGraphs([
        snapshot, *(action for _, action in count_actions), total_action,
        *(action for _, action in yield_actions),
        *(action for pair in range_actions.values() for action in pair),
    ])
    print(f"[sync] wrote {args.output}")
    value_ranges = {name: {"min": low.GetValue(), "max": high.GetValue()}
                    for name, (low, high) in range_actions.items()}
    if "inv_N_orig" in value_ranges and value_ranges["inv_N_orig"]["max"]:
        value_ranges["sum_genWeight"] = 1.0 / value_ranges["inv_N_orig"]["max"]
    write_metadata(args.output, build_metadata(args, files, datasets, dataset, keep,
                                               value_ranges, framework))
    if not args.no_cutflow:
        selection_counts = [(name, int(action.GetValue())) for name, action in count_actions]
        cutflow_path = args.cutflow_output or args.output.with_name(
            f"{args.output.stem}_cutflow.csv"
        )
        write_cutflow(cutflow_path, report_counts, selection_counts)
        write_selection_yields(
            yields_path_from_cutflow(cutflow_path), int(total_action.GetValue()),
            [(key, int(action.GetValue())) for key, action in yield_actions],
        )
    if not args.no_raw_jet_pt:
        # A separate process: its NanoAOD readers fork, which is unsafe after ROOT's threads start.
        command = [sys.executable, str(Path(__file__).with_name("add_raw_jet_pt.py")), str(args.output),
                   "--workers", str(args.nano_workers), "--tree", args.tree]
        for name in datasets:
            command += ["--dataset", name]
        subprocess.run(command, check=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
