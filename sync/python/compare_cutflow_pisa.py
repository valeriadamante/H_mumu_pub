#!/usr/bin/env python3
"""Reconcile our VBF cutflow with the Pisa one.

Their chain differs from ours in the cut order and in four definitions: they
keep the horn jets, they take the two pT-leading jets instead of the pair with
the largest m_jj, they have no 35/25 GeV requirement on the VBF pair and they
do not ask the subleading muon for 26 GeV. The script writes our cutflow in
their order and a bridge that switches those four differences on one at a time,
so the final gap is attributed cut by cut.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))

HELPER = """
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
namespace pisa {
using P4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;

// Our goodJet without the horn veto: pT/eta/jet-id preselection, veto map and
// the DeltaR > 0.4 cleaning against the two selected muons.
template <typename V, typename B1, typename B2>
ROOT::VecOps::RVec<int> JetsNoHorn(const V &pt, const V &eta, const V &phi,
                                   const B1 &presel, const B2 &veto_map,
                                   double mu1_eta, double mu1_phi,
                                   double mu2_eta, double mu2_phi) {
  ROOT::VecOps::RVec<int> keep;
  for (size_t i = 0; i < pt.size(); ++i) {
    if (!presel[i] || veto_map[i]) continue;
    const double d1_eta = eta[i] - mu1_eta, d2_eta = eta[i] - mu2_eta;
    const double d1_phi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi[i] - mu1_phi);
    const double d2_phi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi[i] - mu2_phi);
    if (std::sqrt(d1_eta * d1_eta + d1_phi * d1_phi) <= 0.4) continue;
    if (std::sqrt(d2_eta * d2_eta + d2_phi * d2_phi) <= 0.4) continue;
    keep.push_back(static_cast<int>(i));
  }
  // pT-descending, like SelectedJet_sortIdx.
  std::sort(keep.begin(), keep.end(),
            [&pt](int a, int b) { return pt[a] > pt[b]; });
  return keep;
}

// Dijet observable of the two pT-leading jets of an index collection.
template <typename V>
double LeadingPairValue(const ROOT::VecOps::RVec<int> &idx, const V &pt,
                        const V &eta, const V &phi, const V &mass, int what) {
  if (idx.size() < 2) return -1.;
  const int i = idx[0], j = idx[1];
  if (what == 1) return std::abs(eta[i] - eta[j]);
  P4 a(pt[i], eta[i], phi[i], mass[i]);
  P4 b(pt[j], eta[j], phi[j], mass[j]);
  return (a + b).M();
}

// Pair with the largest m_jj passing the thresholds, as FindVBFJets does.
// min_pt1/min_pt2 = 0 removes the 35/25 GeV requirement.
template <typename V, typename B>
bool BestPair(const V &pt, const V &eta, const V &phi, const V &mass,
              const B &usable, double min_pt1, double min_pt2, bool need_mjj,
              bool need_deta) {
  for (size_t i = 0; i < pt.size(); ++i) {
    if (!usable[i]) continue;
    for (size_t j = i + 1; j < pt.size(); ++j) {
      if (!usable[j]) continue;
      if (pt[i] < min_pt1 || pt[j] < min_pt2) continue;
      if (need_mjj || need_deta) {
        P4 a(pt[i], eta[i], phi[i], mass[i]);
        P4 b(pt[j], eta[j], phi[j], mass[j]);
        if (need_mjj && (a + b).M() < 400.) continue;
        if (need_deta && std::abs(eta[i] - eta[j]) < 2.5) continue;
      }
      return true;
    }
  }
  return false;
}

template <typename V>
bool BestPairFromIndex(const ROOT::VecOps::RVec<int> &idx, const V &pt,
                       const V &eta, const V &phi, const V &mass,
                       double min_pt1, double min_pt2, bool need_mjj,
                       bool need_deta) {
  ROOT::VecOps::RVec<bool> usable(pt.size(), false);
  for (int i : idx) usable[i] = true;
  return BestPair(pt, eta, phi, mass, usable, min_pt1, min_pt2, need_mjj,
                  need_deta);
}

// b-tag veto of an index collection, with our working points.
template <typename B1, typename B2>
bool TagVeto(const ROOT::VecOps::RVec<int> &idx, const B1 &loose,
             const B2 &medium) {
  int n_loose = 0, n_medium = 0;
  for (int i : idx) {
    if (medium[i]) ++n_medium;
    if (loose[i]) ++n_loose;
  }
  return n_medium < 1 && n_loose < 2;
}
}
"""

# Their cutflow, as pasted, for the side-by-side table.
PISA = [
    ("pass_BASESELECTION", 91922769), ("pass_NoiseFilters", 91669412),
    ("pass_GoldenJson", 88160908), ("pass_GoodPV", 88160908),
    ("pass_isTriggered", 88160908), ("pass_AtLeastOneTrgMuon", 61708138),
    ("pass_TwoMuons", 4207453), ("pass_TwoOppositeSignMuons", 4202984),
    ("pass_LeadingMuonCut", 4202984), ("pass_TwoJets", 284055),
    ("pass_DijetMass", 38763), ("pass_DijetDeltaEta", 30353),
    ("pass_BJetVeto", 23876), ("pass_EleVeto", 23857),
    ("pass_AdditionalMuonVeto", 23828), ("pass_ZMassAndSideBandMass", 18337),
]

Z_CR_H_SB = ("((m_mumu > 76 && m_mumu < 115) || (m_mumu > 130 && m_mumu < 150))")

# The official VBF_def of config/Run3_2024/selections_sync.yaml.
VBF_DEF = ("HasVBF && VBFJetIdx_1 >= 0 && VBFJetIdx_2 >= 0 && "
           "static_cast<size_t>(VBFJetIdx_1) < SelectedJet_pt.size() && "
           "static_cast<size_t>(VBFJetIdx_2) < SelectedJet_pt.size() && "
           "SelectedJet_pt.at(VBFJetIdx_1) >= 35 && "
           "SelectedJet_pt.at(VBFJetIdx_2) >= 25")


def define(df):
    ROOT.gInterpreter.Declare(HELPER)
    df = df.Define("pisa_jet_idx",
                   "pisa::JetsNoHorn(Jet_pt, Jet_eta, Jet_phi, Jet_preSel, "
                   "Jet_vetoMap, mu1_eta, mu1_phi, mu2_eta, mu2_phi)")
    df = df.Define("horn_jet_idx",
                   "pisa::JetsNoHorn(Jet_pt, Jet_eta, Jet_phi, "
                   "Jet_preSel && Jet_IsOutsideHorn, Jet_vetoMap, "
                   "mu1_eta, mu1_phi, mu2_eta, mu2_phi)")
    jets = "Jet_pt, Jet_eta, Jet_phi, Jet_mass"
    for tag, idx in (("pisa", "pisa_jet_idx"), ("horn", "horn_jet_idx")):
        df = df.Define(f"{tag}_njets", f"(int){idx}.size()")
        df = df.Define(f"{tag}_ls_mjj",
                       f"pisa::LeadingPairValue({idx}, {jets}, 0)")
        df = df.Define(f"{tag}_ls_deta",
                       f"pisa::LeadingPairValue({idx}, {jets}, 1)")
        df = df.Define(f"{tag}_btag_veto",
                       f"pisa::TagVeto({idx}, Jet_btag_loose, Jet_btag_medium)")
        for name, args in (("best_nopt", "0., 0."), ("best_pt", "35., 25.")):
            df = df.Define(
                f"{tag}_{name}_mjj",
                f"pisa::BestPairFromIndex({idx}, {jets}, {args}, true, false)")
            df = df.Define(
                f"{tag}_{name}_full",
                f"pisa::BestPairFromIndex({idx}, {jets}, {args}, true, true)")
    return df


def chain(df, steps):
    node, actions = df, []
    for name, expression in steps:
        node = node.Filter(expression, name)
        actions.append((name, node.Count()))
    return actions


def write(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[pisa] wrote {path}")


def table(counts, initial):
    rows, previous = [], initial
    for step, (name, passed) in enumerate(counts):
        rows.append({"step": step, "selection": name, "cumulative_pass": passed,
                     "rejected_at_step": previous - passed,
                     "relative_efficiency": passed / previous if previous else 0.0,
                     "cumulative_efficiency": passed / initial if initial else 0.0})
        previous = passed
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    ROOT.EnableImplicitMT(args.threads)

    files, reports = [], []
    for item in args.inputs:
        path = Path(item)
        if path.is_file():
            files.append(str(path))
            report = path.with_name(f"report_{path.stem[len('skim_'):]}.json")
            if report.is_file():
                reports.append(report)
            continue
        files += [str(x) for x in sorted(path.rglob("*.root"))]
        reports += sorted(path.rglob("report_*.json"))
    if not files:
        raise SystemExit("no input ROOT files found")
    initial = 0
    for report in reports:
        with report.open() as stream:
            initial += int(json.load(stream)["Initial"])
    print(f"[pisa] {len(files)} files, {initial} input events")

    df = define(ROOT.RDataFrame("Events", files))

    # Our selection re-ordered like theirs, with our definitions.
    ours = chain(df, [
        ("muons_OS", "mu1_charge * mu2_charge < 0"),
        ("mu1_pt>26", "mu1_pt > 26"),
        ("mu2_pt>26", "mu2_pt > 26"),
        ("two_jets", "N_SelectedJets >= 2"),
        ("dijet_mass_and_deta", VBF_DEF),
        ("btag_veto", "SelectedJetTagSel"),
        ("Z_CR_H_SB", Z_CR_H_SB),
    ])

    # Their selection, applied to our skim.
    theirs = chain(df, [
        ("muons_OS", "mu1_charge * mu2_charge < 0"),
        ("mu1_pt>26", "mu1_pt > 26"),
        ("two_jets", "pisa_njets >= 2"),
        ("dijet_mass", "pisa_ls_mjj > 400"),
        ("dijet_deta", "pisa_ls_deta > 2.5"),
        ("btag_veto", "pisa_btag_veto"),
        ("Z_CR_H_SB", Z_CR_H_SB),
    ])

    # Bridge: from their definitions to ours, one difference at a time.
    common = "mu1_charge * mu2_charge < 0 && mu1_pt > 26"
    bridge_steps = [
        ("pisa_like",
         f"{common} && pisa_njets >= 2 && pisa_ls_mjj > 400 && "
         f"pisa_ls_deta > 2.5 && pisa_btag_veto && {Z_CR_H_SB}"),
        ("+horn_veto",
         f"{common} && horn_njets >= 2 && horn_ls_mjj > 400 && "
         f"horn_ls_deta > 2.5 && horn_btag_veto && {Z_CR_H_SB}"),
        ("+best_mjj_pair",
         f"{common} && horn_njets >= 2 && horn_best_nopt_full && "
         f"horn_btag_veto && {Z_CR_H_SB}"),
        ("+vbf_jet_pt_35_25",
         f"{common} && horn_njets >= 2 && horn_best_pt_full && "
         f"horn_btag_veto && {Z_CR_H_SB}"),
        ("+mu2_pt>15", f"{common} && mu2_pt > 15 && horn_njets >= 2 && "
         f"horn_best_pt_full && horn_btag_veto && {Z_CR_H_SB}"),
        ("+mu2_pt>20", f"{common} && mu2_pt > 20 && horn_njets >= 2 && "
         f"horn_best_pt_full && horn_btag_veto && {Z_CR_H_SB}"),
        ("+mu2_pt>26 (ours)",
         f"{common} && mu2_pt > 26 && horn_njets >= 2 && horn_best_pt_full && "
         f"horn_btag_veto && {Z_CR_H_SB}"),
        ("closure: official columns",
         f"{common} && mu2_pt > 26 && {VBF_DEF} && SelectedJetTagSel && "
         f"{Z_CR_H_SB}"),
    ]
    bridge = [(name, df.Filter(expression, name).Count())
              for name, expression in bridge_steps]

    # Is their mass window ours? Fraction kept by Z_CR_H_SB on their sample.
    ROOT.RDF.RunGraphs([action for _, action in
                        [*ours, *theirs, *bridge]])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ours = [(name, int(action.GetValue())) for name, action in ours]
    theirs = [(name, int(action.GetValue())) for name, action in theirs]
    write(args.output_dir / f"cutflow_{args.label}_ours_pisa_order.csv",
          table(ours, initial))
    write(args.output_dir / f"cutflow_{args.label}_pisa_definitions.csv",
          table(theirs, initial))

    rows, previous = [], None
    for name, action in bridge:
        passed = int(action.GetValue())
        rows.append({"variant": name, "events": passed,
                     "delta": 0 if previous is None else passed - previous})
        previous = passed
    write(args.output_dir / f"bridge_{args.label}.csv", rows)
    for row in rows:
        print(f"  {row['variant']:22s} {row['events']:8d} {row['delta']:+8d}")


if __name__ == "__main__":
    main()
