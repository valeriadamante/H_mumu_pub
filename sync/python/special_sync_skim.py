#!/usr/bin/env python3
"""Special sync skim of 2026-09-23: NanoAOD -> compact sync tuple, data only.

Reuses the framework's golden JSON, MET filters, muon ScaRe/FSR and jet
corrections, jet ID and veto map, and redefines only what the special sync
asks for (definitions in config/Run3_2024/selections_sync_special_20260923.yaml):

  * muon pT is always the plain NanoAOD pT (Muon_pt), including sorting,
    trigger matching and the 26/20 GeV thresholds; only the dimuon system
    uses the NanoAOD pT with ScaRe and FSR recovery (Muon_p4_nano_corr_FSR,
    no bsc);
  * good muons: pT>10, |eta|<2.4, mediumId, |dz|<0.1, |dxy|<0.05 and
    newIso = (pfRelIso04_all*pT - pT_gamma)/pT < 0.25, where pT_gamma is the
    FSR photon (Muon_fsrPhotonIdx) accepted by the FSR recovery and inside the
    isolation cone (dR<0.4); exactly 2 good muons;
  * exactly 2 loose muons: as good muons but looseId, |dz|<1, |dxy|<0.5;
    a superset of the good muons, which are counted among them;
  * electron veto on good electrons: pT>10, |eta|<2.5, mvaIso_WP90;
  * good jets cleaned (dR>0.4) against good muons and good electrons;
  * HasVBF: among pT-ordered good jets, the first (pT-leading) pair with
    pT>35/25, m_jj>400 and |deta|>2.5.

--variant filippo_20260925 reproduces Filippo's VBF selection instead
(config/Run3_2024/selections_sync_filippo_20260925.yaml):

  * the dimuon system uses the production muon p4 (bsc + ScaRe + FSR,
    Muon_p4_FSR_corr);
  * trigger matching without the pT > 26 cut on the matched muon;
  * TightLepVeto jet ID, horn veto with |eta| = 2.5 excluded;
  * HasVBF on the two leading good jets (j0, j1) only;
  * b-veto on jets with corrected pT > 20 GeV (BJetVeto20);
  * the electron veto is left to the selection instead of the skim.

One job writes a tuple (mass_inclusive events) and a JSON with the skim
report and the cutflow counts; special_sync_merge.py combines the jobs.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))
os.environ.setdefault("ANALYSIS_PATH", str(REPOSITORY))
ROOT.gInterpreter.Declare(f'#include "{REPOSITORY / "analysis" / "AnalysisTools.h"}"')

import yaml  # noqa: E402
import common.utilities as utilities  # noqa: E402

ERA = "Run3_2024"
VARIANTS = ("special_20260923", "filippo_20260925", "ggfvbf_20261015", "ggfvbf_pisa_20261013")


def selection_config_path(variant):
    if variant == "ggfvbf_20261015":
        return REPOSITORY / "config" / ERA / "selections_sync.yaml"
    if variant == "ggfvbf_pisa_20261013":
        return REPOSITORY / "config" / ERA / "selections_sync_ggfvbf_pisa_20261013.yaml"
    stem = "selections_sync_special_20260923" if variant == "special_20260923" else f"selections_sync_{variant}"
    return REPOSITORY / "config" / ERA / f"{stem}.yaml"

HELPERS = r"""
#ifndef SPECIAL_SYNC_20260923
#define SPECIAL_SYNC_20260923
namespace special_sync {
using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecB = ROOT::VecOps::RVec<bool>;
using P4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using RVecP4 = ROOT::VecOps::RVec<P4>;

// pT of the FSR photon associated with each muon, with the same acceptance
// as fsr_corrected_p4 (corrections/muon_fsr.py) evaluated on the NanoAOD pT.
// maxDR restricts it to photons inside the isolation cone (0.4); a negative
// value keeps the full FSR-recovery acceptance (dR<0.5).
RVecF FsrPhotonPt(const RVecF& mu_pt, const RVecF& mu_eta, const RVecF& mu_phi,
                  const RVecI& fsr_idx, const RVecF& fsr_pt, const RVecF& fsr_eta,
                  const RVecF& fsr_phi, const RVecF& fsr_dROverEt2,
                  const RVecF& fsr_relIso03, const RVecI& fsr_electronIdx,
                  float maxDR) {
  RVecF out(mu_pt.size(), 0.f);
  for (size_t i = 0; i < mu_pt.size(); ++i) {
    const int k = fsr_idx[i];
    if (k < 0 || k >= int(fsr_pt.size())) continue;
    const float dr = ROOT::VecOps::DeltaR(mu_eta[i], fsr_eta[k], mu_phi[i], fsr_phi[k]);
    if (!(dr > 0.0001 && dr < 0.5)) continue;
    if (maxDR > 0 && !(dr < maxDR)) continue;
    if (fsr_electronIdx[k] != -1) continue;
    if (fsr_pt[k] / mu_pt[i] > 0.4) continue;
    if (fsr_dROverEt2[k] > 0.012) continue;
    if (fsr_relIso03[k] / mu_pt[i] > 1.8) continue;
    out[i] = fsr_pt[k];
  }
  return out;
}

RVecF FsrPhotonDR(const RVecF& mu_eta, const RVecF& mu_phi, const RVecI& fsr_idx,
                  const RVecF& fsr_eta, const RVecF& fsr_phi) {
  RVecF out(mu_eta.size(), -1.f);
  for (size_t i = 0; i < mu_eta.size(); ++i) {
    const int k = fsr_idx[i];
    if (k < 0 || k >= int(fsr_eta.size())) continue;
    out[i] = ROOT::VecOps::DeltaR(mu_eta[i], fsr_eta[k], mu_phi[i], fsr_phi[k]);
  }
  return out;
}

// Jets passing preselection, veto map and horn veto, and farther than 0.4
// from every good muon and every good electron.
RVecB CleanJets(const RVecP4& jets, const RVecB& presel, const RVecB& vetoMap,
                const RVecB& outsideHorn, const RVecP4& muons, const RVecP4& electrons) {
  RVecB out(jets.size(), false);
  for (size_t i = 0; i < jets.size(); ++i) {
    bool ok = presel[i] && !vetoMap[i] && outsideHorn[i];
    for (const auto& m : muons)
      if (ok && ROOT::Math::VectorUtil::DeltaR(jets[i], m) <= 0.4) ok = false;
    for (const auto& e : electrons)
      if (ok && ROOT::Math::VectorUtil::DeltaR(jets[i], e) <= 0.4) ok = false;
    out[i] = ok;
  }
  return out;
}

// Jets are pT ordered, so the first pair found in (i, j>i) order is the one
// with the highest leading pT and, for that leading jet, the highest
// subleading pT.
RVecI LeadingVBFPair(const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass) {
  RVecI result{-1000, -1000};
  for (size_t i = 0; i < pt.size(); ++i) {
    if (!(pt[i] > 35.)) continue;
    for (size_t j = i + 1; j < pt.size(); ++j) {
      if (!(pt[j] > 25.)) continue;
      if (!(std::abs(eta[i] - eta[j]) > 2.5)) continue;
      const P4 a(pt[i], eta[i], phi[i], mass[i]);
      const P4 b(pt[j], eta[j], phi[j], mass[j]);
      if (!((a + b).M() > 400.)) continue;
      result[0] = int(i);
      result[1] = int(j);
      return result;
    }
  }
  return result;
}

// Distance of each object from the closest reference object (99 if none).
RVecF MinDeltaR(const RVecP4& objects, const RVecP4& references) {
  RVecF out(objects.size(), 99.f);
  for (size_t i = 0; i < objects.size(); ++i)
    for (const auto& r : references)
      out[i] = std::min(out[i], float(ROOT::Math::VectorUtil::DeltaR(objects[i], r)));
  return out;
}

// Filippo's VBF pair: always the two leading jets.
RVecI TwoLeadingVBFPair(const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass) {
  RVecI result{-1000, -1000};
  if (pt.size() < 2 || !(pt[0] > 35.) || !(pt[1] > 25.)) return result;
  if (!(std::abs(eta[0] - eta[1]) > 2.5)) return result;
  if (!((P4(pt[0], eta[0], phi[0], mass[0]) + P4(pt[1], eta[1], phi[1], mass[1])).M() > 400.)) return result;
  result[0] = 0;
  result[1] = 1;
  return result;
}
}  // namespace special_sync
#endif
"""


def arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset-name", required=True,
                        help="key of config/Run3_2024/samples.yaml: data (e.g. Muon0_Run2024H) or MC")
    parser.add_argument("--mc-lumi-split", action="store_true",
                        help="MC only: keep the Run3_2024 share (110/247) of the events with the production "
                             "algorithm. The hash uses FullEventId (input-file CRC + entry), so the kept events match "
                             "the production skim only for the same job file grouping. Default: all MC events")
    parser.add_argument("--input-file", required=True, action="append")
    parser.add_argument("--output-file", required=True, type=Path)
    parser.add_argument("--report-file", type=Path)
    parser.add_argument("--n-events", type=int, default=-1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--variant", choices=VARIANTS, default="special_20260923")
    return parser.parse_args()


def load(path):
    with open(path) as stream:
        return yaml.safe_load(stream)


def define_muons(df, muon_columns, trigger_config, filippo, inclusive=False, raw_pt_min=20.0, pisa=False):
    """Good/loose muons with NanoAOD pT, trigger matching and dimuon system."""
    new_cols = []

    def track(node, name, expression):
        if name in new_cols:
            return node
        new_cols.append(name)
        return node.Define(name, expression)

    df = df.Define("Muon_p4_nano", "GetP4(Muon_pt, Muon_eta, Muon_phi, Muon_mass)")
    df = df.Define("Muon_p4_nano_corr", "GetP4(Muon_pt_nano_corr, Muon_eta, Muon_phi, Muon_mass)")
    fsr_args = ("Muon_pt, Muon_eta, Muon_phi, Muon_fsrPhotonIdx, FsrPhoton_pt, FsrPhoton_eta, "
                "FsrPhoton_phi, FsrPhoton_dROverEt2, FsrPhoton_relIso03, FsrPhoton_electronIdx")
    df = df.Define("Muon_fsrPhotonPtIso", f"special_sync::FsrPhotonPt({fsr_args}, 0.4f)")
    df = df.Define("Muon_fsrPhotonPtRecovered", f"special_sync::FsrPhotonPt({fsr_args}, -1.f)")
    df = df.Define("Muon_fsrPhotonDR", "special_sync::FsrPhotonDR(Muon_eta, Muon_phi, "
                   "Muon_fsrPhotonIdx, FsrPhoton_eta, FsrPhoton_phi)")
    df = df.Define("Muon_newIso", "(Muon_pfRelIso04_all * Muon_pt - Muon_fsrPhotonPtIso) / Muon_pt")

    kinematic = "Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_newIso < 0.25"
    df = df.Define("good_muons", f"{kinematic} && Muon_mediumId && abs(Muon_dz) < 0.1 && abs(Muon_dxy) < 0.05")
    df = df.Define("loose_muons", f"{kinematic} && Muon_looseId && abs(Muon_dz) < 1 && abs(Muon_dxy) < 0.5")
    df = track(df, "n_goodMuons", "int(ROOT::VecOps::Sum(good_muons))")
    df = track(df, "n_looseMuons", "int(ROOT::VecOps::Sum(loose_muons))")

    # Trigger matching with the NanoAOD pT and direction.
    df = df.Define("TrigObj_idx", "CreateIndexes(TrigObj_pt.size())")
    df = df.Define("TrigObj_mass", "RVecF(TrigObj_pt.size(), 0.f)")
    df = df.Define("TrigObj_p4", "GetP4(TrigObj_pt, TrigObj_eta, TrigObj_phi, TrigObj_mass, TrigObj_idx)")
    trigger_filters = []
    for path, config in trigger_config.items():
        path_name = config["path"][0]
        leg = config["legs"][0]
        df = df.Define(f"TrigObj_passOnlineCut_{path}", leg["online_obj"]["cut"])
        offline_cut = leg["offline_obj"]["cut"].format(obj="Muon", pt="pt")
        # Filippo matches without a pT cut on the matched muon; the canonical
        # matching (pT > 26) is kept as a diagnostic column.
        df = df.Define(f"Muon_passOfflineCut_{path}", "Muon_pt > 0" if filippo else offline_cut)
        df = df.Define(f"Muon_TriggerMatchingIdx_{path}",
                       f"FindMatching(Muon_passOfflineCut_{path}, TrigObj_passOnlineCut_{path}, "
                       f"Muon_p4_nano, TrigObj_p4, 0.4)")
        if filippo:
            df = df.Define(f"Muon_TriggerMatchingIdx_{path}_pt26",
                           f"FindMatching({offline_cut}, TrigObj_passOnlineCut_{path}, "
                           f"Muon_p4_nano, TrigObj_p4, 0.4)")
        df = track(df, f"Event_HasTriggerMatching_{path}",
                   f"{path_name} && Any(Muon_TriggerMatchingIdx_{path} > -1)")
        new_cols.append(path_name)
        trigger_filters.append(f"Event_HasTriggerMatching_{path}")
    if not inclusive:
        df = df.Filter(" || ".join(trigger_filters), "Trigger matching for " + "__".join(trigger_config))

    index_mask = f"Muon_pt > {raw_pt_min}" if inclusive else "good_muons"
    df = df.Define("sync_preselection_muons", index_mask)
    df = df.Define("good_idx", "ROOT::VecOps::Nonzero(sync_preselection_muons)")
    df = df.Define("sorted_idx_raw", "Reverse(Take(good_idx, Argsort(Take(Muon_pt, good_idx))))")
    if pisa:
        # Pisa candidate: exactly two good and exactly two loose muons.  Their pair
        # replaces the raw-pT pairing, whatever its pT (good muons start at 10 GeV).
        df = track(df, "pisa_candidate", "n_goodMuons == 2 && n_looseMuons == 2")
        df = df.Define("pisa_idx", "ROOT::VecOps::Nonzero(good_muons)")
        df = df.Define("sorted_idx_pisa", "Reverse(Take(pisa_idx, Argsort(Take(Muon_pt, pisa_idx))))")
        df = df.Define("sorted_idx", "pisa_candidate ? sorted_idx_pisa : sorted_idx_raw")
    else:
        df = df.Define("sorted_idx", "sorted_idx_raw")
    df = track(df, "mu1_idx", "sorted_idx.size() > 0 ? int(sorted_idx[0]) : -1")
    df = track(df, "mu2_idx", "sorted_idx.size() > 1 ? int(sorted_idx[1]) : -1")
    df = df.Filter("sorted_idx.size() >= 2" if inclusive else "sorted_idx.size() == 2",
                   f"At least 2 raw-pT>{raw_pt_min:g} muons" if inclusive else "Exactly 2 muons")
    if not inclusive:
        df = df.Filter("n_looseMuons == 2", "Exactly 2 loose muons")

    stage_columns = {
        "pt_raw_noCorr": "Muon_pt_raw_noCorr", "pt_raw_corr": "Muon_pt_raw_corr",
        "pt_raw_scale": "Muon_pt_raw_scale", "pt_FSR_noCorr": "Muon_pt_FSR_noCorr",
        "pt_FSR_corr": "Muon_pt_FSR_corr", "pt_FSR_scale": "Muon_pt_FSR_scale",
        "pt_err": "Muon_pt_err", "pt_nano": "Muon_pt", "pt_nano_corr": "Muon_pt_nano_corr",
        "pt_nano_corr_FSR": "Muon_pt_nano_corr_FSR", "newIso": "Muon_newIso",
        "fsrPhoton_ptIso": "Muon_fsrPhotonPtIso", "fsrPhoton_ptRecovered": "Muon_fsrPhotonPtRecovered",
        "fsrPhoton_dR": "Muon_fsrPhotonDR",
    }
    for i in (1, 2):
        idx = f"mu{i}_idx"
        df = df.Define(f"mu{i}_p4", f"Muon_p4_nano.at({idx})")
        df = df.Define(f"mu{i}_p4_nano_corr", f"Muon_p4_nano_corr.at({idx})")
        df = df.Define(f"mu{i}_p4_nano_corr_FSR", f"Muon_p4_nano_corr_FSR.at({idx})")
        df = df.Define(f"mu{i}_p4_FSR_corr", f"Muon_p4_FSR_corr.at({idx})")
        df = track(df, f"mu{i}_pt", f"Muon_pt[{idx}]")
        df = track(df, f"mu{i}_eta", f"Muon_eta[{idx}]")
        df = track(df, f"mu{i}_phi", f"Muon_phi[{idx}]")
        df = track(df, f"mu{i}_mass", f"Muon_mass[{idx}]")
        for column in muon_columns:
            df = track(df, f"mu{i}_{column.split('_', 1)[1]}", f"{column}[{idx}]")
        for stage, column in stage_columns.items():
            df = track(df, f"mu{i}_{stage}", f"{column}[{idx}]")
        for stage in ("raw_noCorr", "raw_corr", "raw_scale", "FSR_noCorr", "FSR_corr", "FSR_scale"):
            df = track(df, f"mu{i}_eta_{stage}", f"Muon_p4_{stage}[{idx}].Eta()")
            df = track(df, f"mu{i}_phi_{stage}", f"Muon_p4_{stage}[{idx}].Phi()")
        gen_match = (f"Muon_genPartFlav[{idx}] == 1 || Muon_genPartFlav[{idx}] == 15"
                     if "Muon_genPartFlav" in muon_columns else "false")
        df = track(df, f"mu{i}_GenMatched", gen_match)
        for path in trigger_config:
            df = track(df, f"mu{i}_HasTriggerMatching_{path}", f"Muon_TriggerMatchingIdx_{path}[{idx}] >= 0")
            if filippo:
                df = track(df, f"mu{i}_HasTriggerMatching_{path}_pt26",
                           f"Muon_TriggerMatchingIdx_{path}_pt26[{idx}] >= 0")

    pair_filters = [f"(Event_HasTriggerMatching_{path} && (mu1_HasTriggerMatching_{path} || "
                    f"mu2_HasTriggerMatching_{path}))" for path in trigger_config]
    if not inclusive:
        df = df.Filter(" || ".join(pair_filters), "Selected dimuon trigger matching")

    # Dimuon system: NanoAOD pT with ScaRe and FSR recovery (23/09), or the
    # production p4 with bsc + ScaRe + FSR (Filippo).
    df = df.Define("dimuon_p4_nano_corr_FSR", "mu1_p4_nano_corr_FSR + mu2_p4_nano_corr_FSR")
    df = df.Define("dimuon_p4_FSR_corr", "mu1_p4_FSR_corr + mu2_p4_FSR_corr")
    dimuon = "dimuon_p4_FSR_corr" if filippo else "dimuon_p4_nano_corr_FSR"
    df = track(df, "m_mumu", f"{dimuon}.M()")
    df = track(df, "pt_mumu", f"{dimuon}.Pt()")
    df = track(df, "eta_mumu", f"{dimuon}.Eta()")
    df = track(df, "phi_mumu", f"{dimuon}.Phi()")
    df = track(df, "m_mumu_nano_corr_FSR", "dimuon_p4_nano_corr_FSR.M()")
    df = track(df, "m_mumu_nano", "(mu1_p4 + mu2_p4).M()")
    df = track(df, "m_mumu_nano_corr", "(mu1_p4_nano_corr + mu2_p4_nano_corr).M()")
    df = track(df, "m_mumu_FSR_corr", "(mu1_p4_FSR_corr + mu2_p4_FSR_corr).M()")
    df = track(df, "pt_mumu_nano", "(mu1_p4 + mu2_p4).Pt()")
    if inclusive:
        if pisa:
            # Raw OR corrected mass, so that Pisa's corrected-mass window is a subset.
            df = df.Define("sync_preselection", "(m_mumu_nano > 50 && m_mumu_nano < 200) || "
                           "(m_mumu_FSR_corr > 50 && m_mumu_FSR_corr < 200)")
        else:
            df = df.Define("sync_preselection", "m_mumu_nano > 50 && m_mumu_nano < 200")
        df = df.Filter("sync_preselection", "50 < dimuon mass < 200 GeV")
    else:
        df = df.Filter("m_mumu > 50 && m_mumu < 200", "dimuon mass cut")
    return df, new_cols


def define_electrons(df, apply_veto):
    df = df.Define("Electron_p4", "GetP4(Electron_pt, Electron_eta, Electron_phi, Electron_mass)")
    df = df.Define("good_electrons", "Electron_pt > 10 && abs(Electron_eta) < 2.5 && Electron_mvaIso_WP90")
    df = df.Define("n_goodElectrons", "int(ROOT::VecOps::Sum(good_electrons))")
    df = df.Define("GoodElectron_p4", "Electron_p4[good_electrons]")
    if apply_veto:
        df = df.Filter("n_goodElectrons == 0", "No extra electrons")
        return df, ["n_goodElectrons"]
    # Alternative electron vetoes, to find the one Filippo uses.
    df = df.Define("Electron_minDRGoodMuon", "special_sync::MinDeltaR(Electron_p4, GoodMuon_p4)")
    df = df.Define("n_goodElectrons_WP80", "int(Sum(Electron_pt > 10 && abs(Electron_eta) < 2.5 && Electron_mvaIso_WP80))")
    df = df.Define("n_goodElectrons_pt20", "int(Sum(good_electrons && Electron_pt > 20))")
    df = df.Define("n_goodElectrons_dRmu04", "int(Sum(good_electrons && Electron_minDRGoodMuon > 0.4))")
    return df, ["n_goodElectrons", "n_goodElectrons_WP80", "n_goodElectrons_pt20", "n_goodElectrons_dRmu04"]


def main():
    args = arguments()
    if args.threads > 1 and args.n_events < 0:
        ROOT.EnableImplicitMT(args.threads)
    ROOT.gInterpreter.Declare(HELPERS)
    config_dir = REPOSITORY / "config" / ERA
    config = utilities.get_config(str(config_dir / "maincfg.yaml"))
    config["want_variations"] = False
    samples = utilities.get_config(str(config_dir / "samples.yaml"))
    if args.dataset_name not in samples:
        raise SystemExit(f"{args.dataset_name} is not in {config_dir / 'samples.yaml'}")
    dataset_cfg = samples[args.dataset_name]
    is_data = bool(dataset_cfg.get("is_data", False))
    inclusive = args.variant.startswith("ggfvbf")
    pisa = args.variant == "ggfvbf_pisa_20261013"
    raw_pt_min = 10.0 if args.dataset_name in {
        "TTH_Hto2Mu", "ZH_Hto2Mu", "WplusH_Hto2Mu",
    } else 20.0
    sel_config = utilities.get_config(str(config_dir / "selections.yaml"))
    from common.jet_horn_policy import configure_horn_veto, horn_mitigation_enabled
    configure_horn_veto(sel_config, ERA, "configured")
    filippo = args.variant in ("filippo_20260925", "ggfvbf_pisa_20261013")
    if filippo:
        # Filippo's horn: 2.5 < |eta| < 3.0, |eta| = 2.5 excluded.
        sel_config["jet_horn_veto_expr"] = ("(abs(v_ops::eta(Jet_p4)) > 2.5 && abs(v_ops::eta(Jet_p4)) < 3.0 "
                                            "&& v_ops::pt(Jet_p4) <= 50.0)")
    special_config = selection_config_path(args.variant)
    trigger_config = utilities.get_config(str(config_dir / "triggers.yaml"))
    syst_cfg = utilities.get_config(str(config_dir / "systematics.yaml"))
    nano_version = config.get("nano_version", "v15")
    input_files = [path for value in args.input_file for path in value.split(",") if path]

    chain = ROOT.TChain("Events")
    for path in input_files:
        if chain.Add(path) == 0:
            raise RuntimeError(f"could not add input file {path}")
    df = ROOT.RDataFrame(chain)
    if args.n_events > 0:
        df = df.Range(args.n_events)
    df = df.Define("period", f"static_cast<int>(Period::{config['era']})")
    df = df.Define("is_data", "true" if is_data else "false").Define("is_data_int", "1" if is_data else "0")
    df = df.Define("is_signal", "true" if dataset_cfg.get("is_signal", False) else "false")
    df = df.Define("sample", f'std::string("{args.dataset_name}")')
    df = df.Define("sample_id", str(abs(__import__("zlib").crc32(args.dataset_name.encode()))))
    df = df.Define("MET_pt", "PuppiMET_pt")
    df = df.Define("MET_phi", "PuppiMET_phi")
    df = df.Define("MET_collection", 'std::string("PuppiMET")')

    from corrections.general import apply_golden_json, apply_corrections
    mc_cols = []
    if is_data:
        if inclusive:
            from corrections.lumi import load_lumi_map, _build_cpp_function
            lumi_map = load_lumi_map(config["lumiFile"])
            function_name = "passesLumiMask_ggfvbf_sync_20261015"
            ROOT.gInterpreter.Declare(_build_cpp_function(function_name, lumi_map))
            df = df.Define("pass_golden_json", f"{function_name}(run, luminosityBlock)")
        else:
            df = apply_golden_json(df, config["lumiFile"])
    else:
        mc_cols = ["genWeight"]
        if inclusive:
            process_cfg = utilities.get_config(str(config_dir / "process_names.yaml"))
            xs_cfg = utilities.get_config(config["crossSectionsFile"])
            process = utilities.process_from_dataset(process_cfg, args.dataset_name)
            from corrections.general import define_base_weights
            df, base_weights, _ = define_base_weights(
                df, config.get("luminosity", ""), dataset_cfg.get("crossSection", args.dataset_name),
                xs_cfg, config, dataset_cfg, process_cfg[process], False, syst_cfg,
            )
            mc_cols += base_weights
        if args.mc_lumi_split:
            import zlib
            from analysis.mc_splitting import ApplyOrthogonalLumiFilter
            # FullEventId as in analysis/skim.py; the job's file list enters the hash
            dataset_crc = zlib.crc32(args.dataset_name.encode()) & 0xFFFF
            input_crc = zlib.crc32("\n".join(input_files).encode()) & 0xFFFF
            df = df.Define("FullEventId", f"eventId::encodeFullEventId({dataset_crc}, {input_crc}, rdfentry_)")
            df, _ = ApplyOrthogonalLumiFilter(df, ERA, seed=12345, keep_tag_column=False)
    config["apply_jet_horn_mitigation"] = horn_mitigation_enabled(ERA, sel_config)
    df = apply_corrections(df, config, dataset_cfg, args.dataset_name, False)
    if "MET_flags" in config and not inclusive:
        from analysis.other import applyMETFlags
        df = applyMETFlags(df, config, is_data)
    if pisa:
        # Pisa applies the MET filters as a skim filter: keep the decision as a flag.
        from analysis.other import applyBadMETfilter
        if config.get("badMET_flag_runs"):
            df = applyBadMETfilter(df, config["badMET_flag_runs"], is_data)
        df = df.Define("pass_met_filters", " && ".join(config.get("MET_flags", [])) or "true")
        if not is_data:
            df = df.Define("pass_golden_json", "true")

    from analysis.muons import DefineMuonPtAndP4
    df = DefineMuonPtAndP4(df, False)
    muon_columns = utilities.GetObservablesCols("Muon", is_data, nano_version)
    df, cols = define_muons(df, muon_columns, trigger_config, filippo, inclusive, raw_pt_min, pisa)
    if not is_data and inclusive:
        from corrections.mu import apply_muIDIso_weights
        from common.trigger_weights import AddTriggerWeightsAndErrors
        if pisa:
            # As Pisa (and sync_mc_weights.py): scale factors at the bsc+ScaRe+FSR pT and eta.
            for i in (1, 2):
                df = df.Redefine(f"mu{i}_pt", f"mu{i}_pt_FSR_corr").Redefine(f"mu{i}_eta", f"mu{i}_eta_FSR_corr")
        df, mu_weights = apply_muIDIso_weights(df, config, False)
        df = AddTriggerWeightsAndErrors(df, WantErrors=False)
        cols += mu_weights + ["weight_TrgSF_singleMu_IsoMu24Central"]
        df = df.Define("weight_generator", "genWeight")
        df = df.Define("weight_pileup", "weight_pu")
        df = df.Define("weight_trigger", "weight_TrgSF_singleMu_IsoMu24Central")
        df = df.Define(
            "weight_muon_idiso",
            "((mu1_pt < 30) ? weight_mu1_TightPFIso_TightID_Central : weight_mu1_LoosePFIso_MediumID_Central) * "
            "((mu1_pt < 30) ? weight_mu1_TightID_Trk_Central : weight_mu1_MediumID_Trk_Central) * "
            "((mu2_pt < 30) ? weight_mu2_TightPFIso_TightID_Central : weight_mu2_LoosePFIso_MediumID_Central) * "
            "((mu2_pt < 30) ? weight_mu2_TightID_Trk_Central : weight_mu2_MediumID_Trk_Central)",
        )
        cols += ["weight_generator", "weight_pileup", "weight_trigger", "weight_muon_idiso"]
        if pisa:
            for i in (1, 2):
                df = df.Redefine(f"mu{i}_pt", f"Muon_pt[mu{i}_idx]").Redefine(f"mu{i}_eta", f"Muon_eta[mu{i}_idx]")
    df = df.Define("GoodMuon_p4", "Muon_p4_nano[good_muons]")
    df, electron_cols = define_electrons(df, not filippo and not inclusive)
    if inclusive and not pisa:  # the Pisa selection yaml defines electron_veto itself
        df = df.Define("electron_veto", "n_goodElectrons == 0")
        electron_cols.append("electron_veto")
    cols += electron_cols

    from analysis.jets import ProcessAllJetVariables, SelectJetVars
    from corrections.btag_wpValues import getBTagWPValues
    from corrections.jetVetoMap import ApplyJetVetoMap
    btag_wps = getBTagWPValues(config)
    btag_algo = config.get("bTagAlgo", "PNet")
    jet_columns = utilities.GetObservablesCols("Jet", is_data, nano_version)
    df, _ = ProcessAllJetVariables(df, jet_columns, sel_config, btag_algo, btag_wps, False, syst_cfg)
    df, _ = ApplyJetVetoMap(df, config, "nano", False, False, nano_version == "v12", False, syst_cfg)
    # Good-jet definition with the lepton cleaning of the special sync;
    # SelectJetVars keeps predefined columns and builds SelectedJet_* on them.
    pt_min = sel_config.get("jet_pt_min", 25.0)
    eta_max = sel_config.get("jet_eta_max", 4.7)
    jet_id = "Jet_passJetIdTightLepVeto" if filippo else "Jet_passJetIdTight"
    df = df.Define("Jet_preSel", f"v_ops::pt(Jet_p4) > {pt_min} && abs(v_ops::eta(Jet_p4)) < {eta_max} && {jet_id}")
    df = df.Define("Jet_NoOverlapWithMuons", "special_sync::CleanJets(Jet_p4, Jet_preSel, Jet_vetoMap, "
                   "Jet_IsOutsideHorn, GoodMuon_p4, GoodElectron_p4)")
    df = df.Define("goodJet", "Jet_NoOverlapWithMuons")
    df, _ = SelectJetVars(df, jet_columns, sel_config, btag_algo, btag_wps, False, syst_cfg)
    if inclusive:
        df = df.Define("Jet_pt_raw", "Jet_pt_nocorr * (1.f - Jet_rawFactor)")
        df = df.Define("SelectedJet_pt_raw", "Take(Jet_pt_raw, SelectedJet_idx)")
        cols += ["Jet_pt_raw", "SelectedJet_pt_raw", "Jet_rawFactor", "SelectedJet_idx"]
    if inclusive and args.dataset_name in {"ZH_Hto2Mu", "WplusH_Hto2Mu"}:
        available = {str(name) for name in df.GetColumnNames()}
        fatjet_cols = sorted(name for name in available if name.startswith("FatJet_"))
        if "FatJet_pt" in available and "FatJet_rawFactor" in available:
            df = df.Define("FatJet_pt_beforeJERC", "FatJet_pt * (1.f - FatJet_rawFactor)")
            fatjet_cols.append("FatJet_pt_beforeJERC")
        # NanoAOD normally stores the collection counter already.  A few input
        # variants do not, so only synthesize it when it is actually missing.
        if "nFatJet" not in available:
            df = df.Define("nFatJet", "static_cast<int>(FatJet_pt.size())")
        cols += fatjet_cols + ["nFatJet"]
    pair = "TwoLeadingVBFPair" if filippo else "LeadingVBFPair"
    df = df.Define("VBFPairIdx", f"special_sync::{pair}(SelectedJet_pt, SelectedJet_eta, "
                   "SelectedJet_phi, SelectedJet_mass)")
    if filippo:
        # b-veto on the good-jet definition with corrected pT > 20 GeV.
        df = df.Define("BVetoJet_preSel", f"v_ops::pt(Jet_p4) > 20 && abs(v_ops::eta(Jet_p4)) < 2.5 && {jet_id}")
        df = df.Define("BVetoJet", "special_sync::CleanJets(Jet_p4, BVetoJet_preSel, Jet_vetoMap, "
                       "Jet_IsOutsideHorn, GoodMuon_p4, GoodElectron_p4)")
        df = df.Define("n_BVetoJets_loose", "int(Sum(BVetoJet && Jet_btag_loose))")
        df = df.Define("n_BVetoJets_medium", "int(Sum(BVetoJet && Jet_btag_medium))")
        df = df.Define("BJetVeto20", "n_BVetoJets_medium < 1 && n_BVetoJets_loose < 2")
        # Steps of Filippo's cutflow on (j0, j1).
        lead = "SyncP4(SelectedJet_pt[0], SelectedJet_eta[0], SelectedJet_phi[0], SelectedJet_mass[0])"
        sub = "SyncP4(SelectedJet_pt[1], SelectedJet_eta[1], SelectedJet_phi[1], SelectedJet_mass[1])"
        ROOT.gInterpreter.Declare("#include <Math/Vector4D.h>\nusing SyncP4 = "
                                  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;")
        df = df.Define("filippo_LeadingJetPt", "N_SelectedJets >= 2 && SelectedJet_pt[0] > 35 && SelectedJet_pt[1] > 25")
        df = df.Define("filippo_DijetMass", f"filippo_LeadingJetPt && ({lead} + {sub}).M() > 400")
        df = df.Define("filippo_DijetDeltaEta", "filippo_DijetMass && abs(SelectedJet_eta[0] - SelectedJet_eta[1]) > 2.5")
        cols += ["n_BVetoJets_loose", "n_BVetoJets_medium", "BJetVeto20", "filippo_LeadingJetPt",
                 "filippo_DijetMass", "filippo_DijetDeltaEta", "SelectedJetTagSel"]
    df = df.Define("HasVBF", "VBFPairIdx[0] >= 0")
    df = df.Define("VBFJetIdx_1", "VBFPairIdx[0]")
    df = df.Define("VBFJetIdx_2", "VBFPairIdx[1]")

    # Compact tuple with the column set of the previous sync tuples.
    import sync.python.prepare_tuple_forSync as sync_tuple
    special = load(special_config)
    sync_tuple.SYNC_CONFIG = special
    sync_tuple.SYNC_CONFIG_PATH = special_config
    df = sync_tuple.prepare(df, "standard")
    report = df.Report()
    keep = sync_tuple.selected_columns(df, False, extended=inclusive)
    extra = cols + mc_cols + ["sample", "sample_id", "MET_pt", "MET_phi", "MET_collection",
                    "pass_golden_json", "pass_met_filters", "sync_preselection",
                    "n_goodMuons", "n_looseMuons", "pisa_candidate",
                    "SelectedJetTagSel", "JetTagSel", "exactly_2_loose_muons", "bothMuons_ML",
                    "electron_veto", "m_mumu_nano",
                    "pt_mumu_nano", "m_mumu_nano_corr", "m_mumu_FSR_corr"]
    available = sync_tuple.columns(df)
    if inclusive:
        extra += sorted(name for name in available if name.startswith("Flag_"))
    keep = sorted(c for c in set(keep) | {c for c in extra if c in available} if "_p4" not in c)

    cutflows = {}
    count_actions = []
    for key in (k for k in special if k.startswith("cutflow")):
        node = df
        cutflows[key] = []
        for step in special[key]:
            node = node.Filter(step, f"{key}_{step}")
            action = node.Count()
            cutflows[key].append((step, action))
            count_actions.append(action)
    yields = [((kind, name), df.Filter(name).Count()) for kind, name in sync_tuple.stored_category_region_names()]
    total = df.Count()

    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    options = ROOT.RDF.RSnapshotOptions()
    options.fLazy = True
    options.fMode = "RECREATE"
    output_selection = "sync_preselection" if inclusive else "mass_inclusive"
    snapshot = df.Filter(output_selection, "sync output selection").Snapshot(
        "Events", str(args.output_file), utilities.ListToVector(keep), options)
    ROOT.RDF.RunGraphs([snapshot, total, *count_actions, *(a for _, a in yields)])

    _, report_json = utilities.SaveReport(df, report.GetValue())
    summary = {
        "dataset": args.dataset_name, "input_files": input_files,
        "variant": args.variant, "selection_config": str(special_config),
        "raw_muon_pt_min_GeV": raw_pt_min,
        "met_collection": "PuppiMET",
        "report": report_json,
        "cutflows": {key: [[step, int(a.GetValue())] for step, a in steps] for key, steps in cutflows.items()},
        "yields": [[kind, name, int(a.GetValue())] for (kind, name), a in yields],
        "events_after_skim": int(total.GetValue()),
    }
    report_path = args.report_file or args.output_file.with_name(args.output_file.stem + "_report.json")
    report_path.write_text(json.dumps(summary, indent=2))
    print(f"[special-sync] wrote {args.output_file} and {report_path}")


if __name__ == "__main__":
    main()
