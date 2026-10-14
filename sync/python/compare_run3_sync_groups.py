#!/usr/bin/env python3
"""Compare the Run 3 H->mumu exchange tuples of two groups, region by region and channel by channel.

Every group publishes its own ntuple schema in /eos/cms/store/group/phys_higgs/cmshmm/run3_sync.
A ``treat_<group>_ntuple`` function projects one group onto a common description
(mass, ggH/VBF channel, weights, flags and the common variables); the comparison itself only
sees that description.  A new group is added by writing one GroupSpec and its treat function.

Layout of the output:

    <output>/<GROUP1>-<GROUP2>/<sample>/<region>/<channel>/SUMMARY.md + plots

  python3 sync/python/compare_run3_sync_groups.py --output $HMUMU_WWW/skim_v4/studies/sync/groups_comparison_13Oct
  python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 PKU --sample Data_Run2024H ...
  python3 sync/python/compare_run3_sync_groups.py --list
"""
from __future__ import annotations

import os
import argparse
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

BASE = Path("/eos/cms/store/group/phys_higgs/cmshmm/run3_sync")
CHANNELS = ("ggH", "VBF")
REGIONS = {
    "Z_sideband": (lambda m: (m > 70) & (m < 110), "70 < m(mumu) < 110 GeV"),
    "H_sideband": (lambda m: ((m > 110) & (m < 115)) | ((m > 135) & (m < 150)),
                   "110 < m(mumu) < 115 or 135 < m(mumu) < 150 GeV"),
    "Signal_Fit": (lambda m: (m > 115) & (m < 135), "115 < m(mumu) < 135 GeV"),
}
BLINDED = {("Data_Run2024H", "Signal_Fit")}   # no distributions of data in the signal window
MC_MODES = ("genWeights", "all_common_weights")
RESIDUAL_TOLERANCE = 1e-3                        # relative, on events selected by both groups
UNDEFINED = -9999.0                              # a jet that does not exist
HELPERS = {"abs": np.abs, "sqrt": np.sqrt, "pi": math.pi, "cos": np.cos, "sin": np.sin, "sinh": np.sinh,
           "arcsinh": np.arcsinh, "arctan2": np.arctan2, "hypot": np.hypot,
           "dphi": lambda a, b: (a - b + math.pi) % (2 * math.pi) - math.pi}


# ----------------------------------------------------------------------------------------------
# Common description of the groups
# ----------------------------------------------------------------------------------------------

@dataclass
class GroupSpec:
    """How one native ntuple is projected onto the common schema."""
    name: str
    files: dict            # sample -> path relative to BASE
    tree: str
    lumi_branch: str
    mass: str              # corrected dimuon mass
    channels: dict         # channel -> boolean expression, defined on the whole tuple
    weights: dict          # mode -> expression
    flags: list            # boolean branches shown as efficiencies
    notes: list = field(default_factory=list)
    variables: dict = field(default_factory=dict)   # filled from VARIABLES


@dataclass
class Variable:
    name: str
    title: str
    bins: np.ndarray
    expr: dict             # group -> expression, a missing group means "not available"
    min_jets: int = 0      # jets needed for the value to be defined


def V(name, title, bins, TAMU=None, Marmara=None, PKU=None, NNU=None, min_jets=0, Pisa=None):
    return Variable(name, title, np.asarray(bins), {"TAMU": TAMU, "Marmara": Marmara, "PKU": PKU, "NNU": NNU, "Pisa": Pisa,
                                                      "Purdue": PURDUE_EXPR.get(name)}, min_jets)


# Purdue branches (variable name -> expression); jets are the two leading selected jets (`*_nominal`)
PURDUE_EXPR = {
    "m_mumu": "dimuon_mass", "m_mumu_uncorrected": "dimuon_mass_raw", "pt_mumu": "dimuon_pt", "eta_mumu": "dimuon_eta",
    "phi_mumu": "dimuon_phi", "mu1_pt": "mu1_pt", "mu2_pt": "mu2_pt", "mu1_pt_nano": "mu1_pt_nano", "mu2_pt_nano": "mu2_pt_nano",
    "mu1_eta": "mu1_eta", "mu2_eta": "mu2_eta", "mu1_phi": "mu1_phi", "mu2_phi": "mu2_phi",
    "met_pt": "PuppiMET_pt", "met_phi": "PuppiMET_phi", "njets": "njets_nominal",
    "n_bjets_loose": "nBtagLoose_nominal", "n_bjets_medium": "nBtagMedium_nominal",
    "jet1_pt": "jet1_pt_nominal", "jet1_eta": "jet1_eta_nominal", "jet1_phi": "jet1_phi_nominal",
    "jet2_pt": "jet2_pt_nominal", "jet2_eta": "jet2_eta_nominal", "jet2_phi": "jet2_phi_nominal",
    "mjj": "jj_mass_nominal", "pt_jj": "jj_pt_nominal", "eta_jj": "jj_eta_nominal", "deta_jj": "abs(jj_dEta_nominal)",
    "dphi_jj": "abs(jj_dPhi_nominal)", "soft_activity_njets5": "nsoftjets5_nominal",
}


PHI = np.linspace(-math.pi, math.pi, 33)
ETA = np.linspace(-2.5, 2.5, 51)
JETETA = np.linspace(-4.7, 4.7, 48)
VARIABLES = [
    V("m_mumu", r"$m_{\mu\mu}$ [GeV]", np.linspace(70, 150, 81), "m_mumu_FSR_corr", "mm_mass_fsr", "diMufsr_kit_BSC_mass", "dimuon_fsr_mass", Pisa="Dimuon_mass"),
    V("m_mumu_uncorrected", r"uncorrected $m_{\mu\mu}$ [GeV]", np.linspace(60, 160, 101), "m_mumu_nano", "mm_mass_uncorr", None, "dimuon_raw_mass", Pisa="Dimuon_rawMass"),
    V("pt_mumu", r"$p_T^{\mu\mu}$ [GeV]", np.linspace(0, 250, 51), "hypot((mu1_pt_FSR_corr*cos(mu1_phi_FSR_corr)+mu2_pt_FSR_corr*cos(mu2_phi_FSR_corr)),(mu1_pt_FSR_corr*sin(mu1_phi_FSR_corr)+mu2_pt_FSR_corr*sin(mu2_phi_FSR_corr)))", "mm_pt_fsr", "diMufsr_kit_BSC_pt", "dimuon_fsr_pt", Pisa="Dimuon_pt"),
    V("eta_mumu", r"$\eta_{\mu\mu}$", np.linspace(-5, 5, 51), "arcsinh((mu1_pt_FSR_corr*sinh(mu1_eta_FSR_corr)+mu2_pt_FSR_corr*sinh(mu2_eta_FSR_corr))/hypot((mu1_pt_FSR_corr*cos(mu1_phi_FSR_corr)+mu2_pt_FSR_corr*cos(mu2_phi_FSR_corr)),(mu1_pt_FSR_corr*sin(mu1_phi_FSR_corr)+mu2_pt_FSR_corr*sin(mu2_phi_FSR_corr))))", "mm_eta_corr", "diMufsr_kit_BSC_eta", "dimuon_fsr_eta"),
    V("phi_mumu", r"$\phi_{\mu\mu}$", PHI, "arctan2((mu1_pt_FSR_corr*sin(mu1_phi_FSR_corr)+mu2_pt_FSR_corr*sin(mu2_phi_FSR_corr)),(mu1_pt_FSR_corr*cos(mu1_phi_FSR_corr)+mu2_pt_FSR_corr*cos(mu2_phi_FSR_corr)))", "mm_phi_corr", "diMufsr_kit_BSC_phi", "dimuon_fsr_phi"),
    V("mu1_pt", r"leading $\mu$ corrected $p_T$ [GeV]", np.linspace(0, 200, 51), "mu1_pt_FSR_corr", "mu_pt_corr[:,0]", "Mu1_fsr_kit_BSC_pt", "mu1_pt_fsr", Pisa="LeadingMuon_ptToUse"),
    V("mu2_pt", r"subleading $\mu$ corrected $p_T$ [GeV]", np.linspace(0, 120, 49), "mu2_pt_FSR_corr", "mu_pt_corr[:,1]", "Mu2_fsr_kit_BSC_pt", "mu2_pt_fsr", Pisa="SubLeadingMuon_ptToUse"),
    V("mu1_pt_nano", r"leading $\mu$ NanoAOD $p_T$ [GeV]", np.linspace(0, 200, 51), "mu1_pt", "mu_pt_nano[:,0]", None, "mu1_pt_raw", Pisa="LeadingMuon_pt"),
    V("mu2_pt_nano", r"subleading $\mu$ NanoAOD $p_T$ [GeV]", np.linspace(0, 120, 49), "mu2_pt", "mu_pt_nano[:,1]", None, "mu2_pt_raw", Pisa="SubLeadingMuon_pt"),
    V("mu1_eta", r"leading $\mu$ $\eta$", ETA, "mu1_eta", "mu_eta[:,0]", "mu1_fromH_eta", "mu1_eta_fsr", Pisa="LeadingMuon_eta"),
    V("mu2_eta", r"subleading $\mu$ $\eta$", ETA, "mu2_eta", "mu_eta[:,1]", "mu2_fromH_eta", "mu2_eta_fsr", Pisa="SubLeadingMuon_eta"),
    V("mu1_phi", r"leading $\mu$ $\phi$", PHI, "mu1_phi", "mu_phi[:,0]", "mu1_fromH_phi", "mu1_phi_fsr", Pisa="LeadingMuon_phi"),
    V("mu2_phi", r"subleading $\mu$ $\phi$", PHI, "mu2_phi", "mu_phi[:,1]", "mu2_fromH_phi", "mu2_phi_fsr", Pisa="SubLeadingMuon_phi"),
    V("met_pt", r"MET $p_T$ [GeV]", np.linspace(0, 250, 51), "MET_pt", "puppimet_pt", "met_pt", "MET_pt", Pisa="PuppiMET_pt"),
    V("met_phi", r"MET $\phi$", PHI, "MET_phi", "puppimet_phi", "met_phi", "MET_phi", Pisa="PuppiMET_phi"),
    V("njets", "selected jet multiplicity", np.arange(-0.5, 10.5, 1), "N_SelectedJets", "nJetsAna", "njets", "nGoodJets", Pisa="nSelectedJet"),
    V("n_bjets_loose", "loose b-tagged jets", np.arange(-0.5, 5.5, 1), "n_BVetoJets_loose", "nLooseB", "nbjets_loose", "nBtagLoose", Pisa="nLooseBTagJet"),
    V("n_bjets_medium", "medium b-tagged jets", np.arange(-0.5, 5.5, 1), "n_BVetoJets_medium", "nMediumB", "nbjets_medium", "nBtagMedium", Pisa="nMediumBTagJet"),
    V("jet1_pt", r"leading jet $p_T$ [GeV]", np.linspace(0, 300, 61), "leadingjet_pt", "byrank(jet_pt, jet_passAnalysis, jet_pt, 0)", "jet1_pt", "jet1pt", 1, Pisa="LeadingJet_ptToUse"),
    V("jet1_eta", r"leading jet $\eta$", JETETA, "leadingjet_eta", "byrank(jet_eta, jet_passAnalysis, jet_pt, 0)", "jet1_eta", "jet1eta", 1, Pisa="LeadingJet_eta"),
    V("jet1_phi", r"leading jet $\phi$", PHI, "leadingjet_phi", "byrank(jet_phi, jet_passAnalysis, jet_pt, 0)", "jet1_phi", None, 1, Pisa="LeadingJet_phi"),
    V("jet2_pt", r"subleading jet $p_T$ [GeV]", np.linspace(0, 200, 51), "subleadingjet_pt", "byrank(jet_pt, jet_passAnalysis, jet_pt, 1)", "jet2_pt", "jet2pt", 2, Pisa="SubLeadingJet_ptToUse"),
    V("jet2_eta", r"subleading jet $\eta$", JETETA, "subleadingjet_eta", "byrank(jet_eta, jet_passAnalysis, jet_pt, 1)", "jet2_eta", "jet2eta", 2, Pisa="SubLeadingJet_eta"),
    V("jet2_phi", r"subleading jet $\phi$", PHI, "subleadingjet_phi", "byrank(jet_phi, jet_passAnalysis, jet_pt, 1)", "jet2_phi", None, 2, Pisa="SubLeadingJet_phi"),
    V("mjj", r"$m_{jj}$ [GeV]", np.linspace(0, 2000, 51), "dijet_ls_mass", "jj_mass", "dijet_mass", "dijetmass", 2, Pisa="Dijet_mass"),
    V("pt_jj", r"$p_T^{jj}$ [GeV]", np.linspace(0, 400, 51), "dijet_ls_pt", "jj_pt", "dijet_pt", "dijetpt", 2, Pisa="Dijet_pt"),
    V("eta_jj", r"$\eta_{jj}$", np.linspace(-8, 8, 51), "dijet_ls_eta", "jj_eta", "dijet_eta", "dijeteta", 2, Pisa="Dijet_eta"),
    V("deta_jj", r"$|\Delta\eta_{jj}|$", np.linspace(0, 10, 51), "abs(leadingjet_eta-subleadingjet_eta)", "abs(jj_deta)", "abs(delta_eta_jj)", "absDetajj", 2, Pisa="abs(Dijet_DeltaEta)"),
    V("dphi_jj", r"$|\Delta\phi_{jj}|$", np.linspace(0, math.pi, 33), "abs(dphi(leadingjet_phi,subleadingjet_phi))", None, "abs(delta_phi_jj)", "abs(dijetdphi)", 2),
    V("soft_activity_njets5", "soft-activity jets, $p_T>5$ GeV", np.arange(-0.5, 15.5, 1), None, None, "SoftActivityJetNjets5", "SoftActivityJetNjets5_tmva"),
]


def tamu_files():
    """Exchange area by default; HMUMU_TAMU_STAGING=<dir with TAMU-<dataset>/sync_tuple.root> reads the unpublished merged tuples."""
    staging = os.environ.get("HMUMU_TAMU_STAGING")
    names = {"Data_Run2024H": "Run2024H", "ggH": "GluGluHto2Mu", "VBFH": "VBFHto2Mu_M125_powheg", "ttH": "TTH_Hto2Mu",
             "ZH": "ZH_Hto2Mu", "WplusH": "WplusH_Hto2Mu"}
    if staging:
        return {sample: f"{staging}/TAMU-{name}/sync_tuple.root" for sample, name in names.items()}
    return {sample: f"TAMU/TAMU-{name}.root" for sample, name in names.items()}


def treat_TAMU_ntuple():
    return GroupSpec(
        "TAMU", tamu_files(),
        "Events", "luminosityBlock", "m_mumu_FSR_corr",
        {"ggH": "ggF", "VBF": "VBF"},
        {"genWeights": "weight_generator",
         "all_common_weights": "weight_generator*weight_pileup*weight_muon_idiso*weight_trigger"},
        ["pass_golden_json", "pass_met_filters", "HLT_IsoMu24", "Event_HasTriggerMatching_singleMu", "muons_trigger",
         "muons_OS", "mu_pt_sel", "bothMuons_ML", "exactly_2_loose_muons", "electron_veto", "jets_btag_clean_sel",
         "pisa_candidate", "baseline", "VBF", "ggF"],
        ["Channels are the stored Pisa-style `ggF` / `VBF` flags (baseline, then VBF = two leading jets pass 35/25 GeV, "
         "m_jj > 400, |deta| > 2.5; ggF = baseline and not VBF). Pisa's exchange tuples use 30/30 GeV jets, so jet-related "
         "differences are expected (decision 2026-10-14: TAMU keeps 35/25).",
         "The dimuon is the pair of the single good muon pair (exactly 2 good and 2 loose muons) whenever it exists "
         "(`pisa_candidate`); otherwise the two highest-pT NanoAOD muons."])


def treat_Marmara_ntuple():
    return GroupSpec(
        "Marmara", {"Data_Run2024H": "Marmara/tree_DATA_2024_loose_larger_stat_Data_File.root", "ggH": "Marmara/tree_GGH_2024_loose_LAST.root",
                    "VBFH": "Marmara/tree_VBF_2024_loose_LAST.root"},
        "HmmTree", "lumi", "mm_mass_fsr",
        {"ggH": "sel_full & cat_isGGH[:,0]", "VBF": "sel_full & cat_isVBF[:,0]"},
        {"genWeights": "w_gen", "all_common_weights": "w_total"},
        ["HLT_IsoMu24", "trigMatched", "passMuonEta", "passMuonID", "passMuonIso", "passMuonPt", "passMuonSel",
         "passJetVetoMap", "passBVeto", "passEleVeto", "passMETFilters", "sel_full"],
        ["`cat_isGGH`/`cat_isVBF` only classify the jets (they are true whenever the b-veto passes), so the channel "
         "is `sel_full && cat_*[0]`.  Index 0 of the two-element category arrays is used; the meaning of index 1 "
         "should be confirmed by Marmara.",
         "Jets are the entries of the jet arrays with `jet_passAnalysis` (`nJetsAna` of them), ranked by `jet_pt` because the arrays are not pT-ordered (92% of events have the largest-pT jet first); `jj_*` is their leading pair. Which jet-energy corrections `jet_pt` carries (it differs from the other groups' corrected pT by a median of a few GeV) has to be confirmed with Marmara.",
         "Data: `tree_DATA_2024_loose_larger_stat_Data_File.root` is a **partial** sample (about 0.8M events, i.e. a few input files), so coverage against the full-statistics tuples is low by construction and only the events in common are meaningful; `data_Run2024H_Muon*.root` are NanoAOD inputs.",
         "The `GluGluH-...root` / `VBFH-...root` files are NanoAOD too; the `tree_*_loose_LAST.root` tuples are used."])


def treat_PKU_ntuple():
    return GroupSpec(
        "PKU", {"Data_Run2024H": "PKU/Muon_24H.root", "ggH": "PKU/GluGluH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root",
                "VBFH": "PKU/VBFH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root",
                "ttH": "PKU/TTH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root",
                "ZH": "PKU/ZH-Hto2Mu_Par-M-125_powhegMINLO_RunIII2024Summer24.root",
                "WplusH": "PKU/WplusH-Hto2Mu_Par-M-125_powhegMINLO_RunIII2024Summer24.root"},
        "ntuple", "lumi", "diMufsr_kit_BSC_mass",
        {"ggH": "GGH", "VBF": "VBF"},
        {"genWeights": "weight", "all_common_weights": "eventWeight"},
        ["trg_single_mu24", "VBF", "GGH", "TTH", "VH"],
        ["PKU stores one exclusive category per event (`VBF`, `GGH`, `TTH`, `VH`); the tuple only contains events "
         "that passed the PKU analysis selection.  A ttH or VH event is therefore never in the ggH or VBF channel."])


def treat_NNU_ntuple():
    base = ("passLumiJSON & passTrigger & passTriggerMatch & passOS & passCandidateGoodMuons & passElectronVeto & "
            "passBTagVeto & passMETFilters & passJetVetoMap & passLeadPt & passPV & passExactly2GoodMuons")
    return GroupSpec(
        "NNU", {"Data_Run2024H": "NNU/NNU-dataH_2024_muon01.root", "VBFH": "NNU/NNU-VBFH-Hto2Mu_M125_RunIII2024Summer24.root"},
        "Events", "luminosityBlock", "dimuon_fsr_mass",
        {"VBF": "passAnalysisSelection", "ggH": f"({base}) & (category != 3)"},
        {"genWeights": "genWeight", "all_common_weights": "weightBaseNoNorm"},
        ["passLumiJSON", "passTrigger", "passTriggerMatch", "passOS", "passCandidateGoodMuons", "passElectronVeto",
         "passBTagVeto", "passMETFilters", "passJetVetoMap", "passLeadPt", "passExactly2GoodMuons",
         "passVBFKinematics", "passAnalysisSelection"],
        ["`category` is assigned to every event (code 3 = VBF kinematics, 0-2 = other); only `passAnalysisSelection` "
         "(the 11 common flags, exactly two good muons and category 3) is a selection.  The ggH channel is "
         "the same common selection with category 0-2; this is our reading of the tuple and should be confirmed by NNU."])


def treat_Pisa_ntuple():
    """Pisa's exchange tuples (`PISA/PISA-<dataset>/Events.root`).  Pisa runs a VBF cutflow only: `first_failed_step` is the index of the
    first failed `pass_*` cut (0 = all cuts pass); step 18 is the mass window (Z+sidebands in data, Higgs mass in MC), so events with
    first_failed_step 0 or >= 18 pass every other cut and our mass region is applied afterwards."""
    return GroupSpec(
        "Pisa", {"Data_Run2024H": ["PISA/PISA-Muon0_Run2024H-MINIv6NANOv15-v1/Events.root", "PISA/PISA-Muon1_Run2024H-MINIv6NANOv15-v2/Events.root"],
                 "ggH": "PISA/PISA-GluGluH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15-v2/Events.root",
                 "VBFH": "PISA/PISA-VBFH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15-v2/Events.root"},
        "Events", "luminosityBlock", "Dimuon_mass",
        {"ggH": "first_failed_step < -1", "VBF": "(first_failed_step == 0) | (first_failed_step >= 18)"},
        {"genWeights": "Generator_weight", "all_common_weights": "HiggsMass_cumulativeWeight"},
        ["pass_NoiseFilters", "pass_GoldenJson", "pass_isTriggered", "pass_AtLeastOneTrgMuon", "pass_TwoMuons", "pass_TwoOppositeSignMuons",
         "pass_LeadingMuonCut", "pass_EleVeto", "pass_AdditionalMuonVeto", "pass_RawDimuonMassWindow", "pass_TwoJets", "pass_LeadingJetPt",
         "pass_DijetMass", "pass_DijetDeltaEta", "pass_BJetVeto"],
        ["Pisa publishes only the **VBF cutflow**: there is no ggH category (the ggH cells are empty).  The tuple holds the events that "
         "reach the raw-mass window (about 1.9M per data dataset), with `first_failed_step` telling where each one stops; the VBF channel is "
         "every cut but the mass window passing.  Data in the signal window are blinded and absent.",
         "Variable mapping checked against the TAMU Pisa-like tuple on common events: `Dimuon_mass` is the FSR-recovered, bs-constrained, ScaRe mass "
         "(0.4% of events differ from TAMU `m_mumu_FSR_corr`); `LeadingMuon_pt` is the NanoAOD pT and `LeadingMuon_ptToUse` the corrected one; "
         "`LeadingJet_ptToUse` is the corrected jet pT (`LeadingJet_pt` is uncorrected).",
         "MET is PuppiMET.  There is no dimuon eta/phi and no dphi_jj.  `all_common_weights` is `HiggsMass_cumulativeWeight` (MC only)."])


def treat_Purdue_ntuple():
    """Purdue exchange tuples (`Purdue/Purdue-<dataset>/Purdue-<dataset>_<n>.root`, several parts per sample, tree `Events`)."""
    def parts(dataset, n):
        """The n part files of a dataset (ordered; n is the number published on 2026-10-14)."""
        return [f"Purdue/Purdue-{dataset}/Purdue-{dataset}_{i}.root" for i in range(n)]
    return GroupSpec(
        "Purdue", {"Data_Run2024H": parts("Muon0_Run2024H-MINIv6NANOv15-v1", 4) + parts("Muon1_Run2024H-MINIv6NANOv15-v2", 4),
                   "ggH": parts("GluGluH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8", 2),
                   "VBFH": parts("VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8", 4),
                   "ttH": parts("TTH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8", 3),
                   "ZH": parts("ZH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8", 2),
                   "WplusH": parts("WplusH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8", 1)},
        "Events", "luminosityBlock", "dimuon_mass",
        {"ggH": "(pass_analysis_preselection == 1) & (cat_ggh == 1)", "VBF": "(pass_analysis_preselection == 1) & (cat_vbf == 1)"},
        {"genWeights": "separate_wgt_genWeight", "all_common_weights": "wgt_nominal"},
        ["pass_analysis_preselection", "pass_HLT", "pass_trigger_match", "pass_MET_filters", "pass_PV", "pass_nmuons", "pass_OS",
         "pass_electron_veto", "pass_jet_veto_map", "pass_bveto", "pass_vbf_cuts", "cat_vbf", "cat_ggh"],
        ["Channels: `pass_analysis_preselection` and `cat_ggh` / `cat_vbf` (the category flags alone also fire on events that fail the preselection). "
         "VBF jets have pT > 35 / 25 GeV like TAMU's (minimum leading-jet pT 35.01 in `cat_vbf`).",
         "`dimuon_mass` equals `dimuon_mass_corr` (corrected mass); `dimuon_mass_raw` is the uncorrected one. Jets are the `*_nominal` leading jets.",
         "`all_common_weights` is `wgt_nominal` (1 in data)."])


TREATERS = {"TAMU": treat_TAMU_ntuple, "Marmara": treat_Marmara_ntuple,
            "PKU": treat_PKU_ntuple, "NNU": treat_NNU_ntuple, "Pisa": treat_Pisa_ntuple, "Purdue": treat_Purdue_ntuple}


def treat(group):
    spec = TREATERS[group]()
    spec.variables = {v.name: v.expr[group] for v in VARIABLES if v.expr.get(group)}
    return spec


# ----------------------------------------------------------------------------------------------
# Reading one native ntuple
# ----------------------------------------------------------------------------------------------

def preprocess(expr):
    """`name[:,i]` -> first(name, i), `name[mask][:,i]` -> first(sel(name, mask), i): padded access to a jagged branch."""
    expr = re.sub(r"(\w+)\[(\w+)\]\[:,(\d)\]", r"first(sel(\1,\2),\3)", expr)
    return re.sub(r"(\w+)\[:,(\d)\]", r"first(\1,\2)", expr)


def names(expr):
    return [x for x in dict.fromkeys(re.findall(r"[A-Za-z_]\w*", preprocess(expr))) if x not in HELPERS and x not in ("first", "sel", "byrank")]


def first_element(array, index):
    """The index-th entry of a jagged branch; missing entries are NaN (False for booleans)."""
    padded = ak.pad_none(array, index + 1, axis=1)[:, index]
    return ak.fill_none(padded, False if "bool" in str(ak.type(array)) else np.nan)


def byrank(array, mask, key, rank):
    """The rank-th entry (0 = largest `key`) of `array` among the entries where `mask` is true; for branches that are not pT-ordered."""
    mask = ak.values_astype(mask, bool)
    return first_element(array[mask][ak.argsort(key[mask], axis=1, ascending=False)], rank)


def evaluate(arrays, expr):
    scope = {n: arrays[n] for n in names(expr)}
    scope.update(HELPERS)
    scope["first"] = first_element
    scope["sel"] = lambda a, mask: a[ak.values_astype(mask, bool)]
    scope["byrank"] = byrank
    value = eval(preprocess(expr), {"__builtins__": {}}, scope)  # expressions are literals of this file
    value = np.asarray(ak.to_numpy(value))
    return value[:, 0] if value.ndim > 1 else value


def pack_keys(run, lumi, event, is_data):
    run, lumi, event = (np.asarray(x).astype(np.uint64) for x in (run, lumi, event))
    if is_data:
        return (run - np.uint64(370000)) << np.uint64(48) | lumi << np.uint64(34) | event
    return np.uint64(1) << np.uint64(63) | run << np.uint64(59) | lumi << np.uint64(39) | event


@dataclass
class Extract:
    group: str
    sample: str
    entries: int
    full_keys: np.ndarray                  # every event of the file, sorted and unique
    keys: np.ndarray                       # events with 70 < m < 150 in at least one channel, sorted and unique
    mass: np.ndarray
    channel: dict                          # channel -> bool
    values: dict                           # variable -> float32
    flags: dict                            # flag -> bool
    weights: dict                          # mode -> float64


def sample_paths(spec, sample):
    """One path or a list of paths (e.g. Muon0 + Muon1); absolute paths are used as they are."""
    files = spec.files[sample]
    return [BASE / f for f in ([files] if isinstance(files, str) else files)]


def needed_branches(spec, sample, tree_keys):
    exprs = [spec.mass, *spec.channels.values(), *spec.variables.values()]
    if sample != "Data_Run2024H":
        exprs += list(spec.weights.values())
    wanted = ["run", spec.lumi_branch, "event", *spec.flags]
    for expr in exprs:
        wanted += names(expr)
    return [b for b in dict.fromkeys(wanted) if b in tree_keys]


def extract(group, sample, cache_dir):
    spec = treat(group)
    paths = sample_paths(spec, sample)
    cache = cache_dir / f"{group}_{sample}.npz"
    signature = {"path": str(paths[0]) if len(paths) == 1 else [str(p) for p in paths],
                 "size": paths[0].stat().st_size if len(paths) == 1 else [p.stat().st_size for p in paths],
                 "mtime": paths[0].stat().st_mtime if len(paths) == 1 else [p.stat().st_mtime for p in paths],
                 "variables": spec.variables, "channels": spec.channels, "weights": spec.weights,
                 "flags": spec.flags, "version": 2}
    if cache.exists():
        stored = np.load(cache, allow_pickle=True)
        if json.loads(str(stored["signature"])) == signature:
            return Extract(group, sample, int(stored["entries"]), stored["full_keys"], stored["keys"], stored["mass"],
                           {k: stored[f"ch_{k}"] for k in CHANNELS}, {k[2:]: stored[k] for k in stored.files if k.startswith("v_")},
                           {k[2:]: stored[k] for k in stored.files if k.startswith("f_")},
                           {k[2:]: stored[k] for k in stored.files if k.startswith("w_")})
    trees = [uproot.open(p)[spec.tree] for p in paths]
    present = set(trees[0].keys())
    is_data = sample == "Data_Run2024H"
    modes = () if is_data else MC_MODES
    usable_vars = {k: e for k, e in spec.variables.items() if set(names(e)) <= present}
    usable_flags = [f for f in spec.flags if f in present]
    full, parts = [], []
    entries = sum(t.num_entries for t in trees)
    print(f"[read] {group} {sample}: {entries:,} entries", flush=True)
    chunks = (chunk for tree in trees for chunk in tree.iterate(needed_branches(spec, sample, present), step_size="200 MB", library="ak"))
    for chunk in chunks:
        keys = pack_keys(chunk["run"], chunk[spec.lumi_branch], chunk["event"], is_data)
        mass = np.asarray(evaluate(chunk, spec.mass), dtype=float)
        channel = {c: evaluate(chunk, e).astype(bool) for c, e in spec.channels.items()}
        keep = (mass > 70) & (mass < 150) & (channel["ggH"] | channel["VBF"])
        full.append(keys)
        values = {k: evaluate(chunk, e).astype(np.float32)[keep] for k, e in usable_vars.items()}
        parts.append({"keys": keys[keep], "mass": mass[keep], "channel": {c: v[keep] for c, v in channel.items()},
                      "values": values,
                      "flags": {f: np.asarray(ak.to_numpy(ak.fill_none(chunk[f], False))).reshape(len(keys), -1)[:, 0].astype(bool)[keep]
                                for f in usable_flags},
                      "weights": {m: evaluate(chunk, spec.weights[m]).astype(float)[keep] for m in modes}})
    full_keys = np.unique(np.concatenate(full))
    cat = lambda f: np.concatenate([f(p) for p in parts])
    keys = cat(lambda p: p["keys"])
    order = np.argsort(keys, kind="stable")
    keys = keys[order]
    first = np.concatenate([[True], keys[1:] != keys[:-1]])          # first occurrence of every event
    pick = lambda a: a[order][first]
    out = Extract(group, sample, entries, full_keys, keys[first], pick(cat(lambda p: p["mass"])),
                  {c: pick(cat(lambda p: p["channel"][c])) for c in CHANNELS},
                  {k: pick(cat(lambda p: p["values"][k])) for k in usable_vars},
                  {k: pick(cat(lambda p: p["flags"][k])) for k in usable_flags},
                  {m: pick(cat(lambda p: p["weights"][m])) for m in modes})
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez(cache, signature=json.dumps(signature), entries=out.entries, full_keys=out.full_keys, keys=out.keys,
             mass=out.mass, **{f"ch_{k}": v for k, v in out.channel.items()},
             **{f"v_{k}": v for k, v in out.values.items()}, **{f"f_{k}": v for k, v in out.flags.items()},
             **{f"w_{k}": v for k, v in out.weights.items()})
    print(f"[read] {group} {sample}: kept {len(out.keys):,} events in 70-150 GeV and ggH|VBF", flush=True)
    return out


# ----------------------------------------------------------------------------------------------
# Comparison
# ----------------------------------------------------------------------------------------------

def selection(e, region, channel):
    return REGIONS[region][0](e.mass) & e.channel[channel]


def destination_labels(e):
    """For every kept event the first 'channel / region' it belongs to (index into the label list)."""
    names_ = [f"{channel} / {region}" for channel in CHANNELS for region in REGIONS] + ["outside 70-150 GeV or ggH/VBF"]
    code = np.full(len(e.keys), len(names_) - 1)
    for i, (channel, region) in reversed(list(enumerate((c, r) for c in CHANNELS for r in REGIONS))):
        code[selection(e, region, channel)] = i
    return code, names_


def undefined_jets(values, name, njets, min_jets):
    if min_jets and njets is not None:
        return np.where(njets >= min_jets, values, UNDEFINED)
    return values


def usable(e, variable):
    values = e.values[variable.name].astype(float)
    return undefined_jets(values, variable.name, e.values.get("njets"), variable.min_jets)


def histogram(x, w, bins):
    valid = np.isfinite(x) & np.isfinite(w) & (x >= bins[0]) & (x <= bins[-1])
    x, w = x[valid], w[valid]
    h = np.histogram(np.clip(x, bins[0], np.nextafter(bins[-1], bins[0])), bins, weights=w)[0]
    mean = float(np.average(x, weights=w)) if w.sum() else math.nan
    std = float(math.sqrt(np.average((x - mean) ** 2, weights=w))) if w.sum() else math.nan
    return h, int(valid.sum()), mean, std


def distance(a, b):
    return float(0.5 * np.abs(a / a.sum() - b / b.sum()).sum()) if a.sum() and b.sum() else math.nan


def draw(ha, hb, bins, title, labels, heading, out):
    width = np.diff(bins)
    ya = ha / ha.sum() / width if ha.sum() else ha
    yb = hb / hb.sum() / width if hb.sum() else hb
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": .05})
    ax.stairs(ya, bins, label=labels[0], linewidth=1.7)
    ax.stairs(yb, bins, label=labels[1], linewidth=1.7)
    ax.set_ylabel("unit-normalised density")
    ax.set_title(heading, fontsize=10)
    ax.legend()
    # a bin empty in either histogram has no meaningful ratio: leave it undefined (NaN) instead of a point at zero
    ratio = np.divide(ya, yb, out=np.full_like(ya, np.nan), where=(ya > 0) & (yb > 0))
    rax.axhline(1, color="black", linewidth=.8)
    centres = 0.5 * (bins[1:] + bins[:-1])
    rax.errorbar(centres, ratio, xerr=width / 2, fmt="o", color="C0", markersize=3, linestyle="none", elinewidth=1)
    finite = ratio[np.isfinite(ratio)]
    spread = float(np.max(np.abs(finite - 1))) if finite.size else 0.
    half = min(max(1.2 * spread, 0.02), 0.5)      # zoom on the differences, never wider than 0.5-1.5
    rax.set_ylim(1 - half, 1 + half)
    rax.set_ylabel(f"{labels[0]} /\n{labels[1]}")
    rax.set_xlabel(title)
    fig.savefig(out.with_suffix(".png"), dpi=130, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def migration(selected_keys, other, other_selected_keys):
    """Where events selected by one group but not by the other go in the other tuple."""
    code, label_names = destination_labels(other)
    lost = selected_keys[~np.isin(selected_keys, other_selected_keys) & np.isin(selected_keys, other.full_keys)]
    outside = len(label_names) - 1
    if not len(other.keys):
        labels = np.full(len(lost), outside)
    else:
        position = np.minimum(np.searchsorted(other.keys, lost), len(other.keys) - 1)
        labels = np.where(other.keys[position] == lost, code[position], outside)
    counts = {label_names[i]: int(n) for i, n in zip(*np.unique(labels, return_counts=True))}
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))


def compare_cell(a, b, region, channel):
    """Everything known about one (region, channel) of one sample for the pair (a, b)."""
    ga, gb = a.group, b.group
    sa, sb = selection(a, region, channel), selection(b, region, channel)
    ka, kb = a.keys[sa], b.keys[sb]
    both = np.isin(ka, kb)
    a_in_b, b_in_a = np.isin(ka, b.full_keys), np.isin(kb, a.full_keys)
    cell = {"selected": {ga: int(sa.sum()), gb: int(sb.sum())}, "both": int(both.sum()),
            f"{ga}_only_absent_from_{gb}": int((~a_in_b).sum()), f"{ga}_only_present_in_{gb}": int((a_in_b & ~both).sum()),
            f"{gb}_only_absent_from_{ga}": int((~b_in_a).sum()), f"{gb}_only_present_in_{ga}": int((b_in_a & ~np.isin(kb, ka)).sum()),
            "migration": {ga: migration(ka, b, kb), gb: migration(kb, a, ka)}}
    ia, ib = np.searchsorted(a.keys, ka[both]), np.searchsorted(b.keys, ka[both])
    residuals = {}
    for v in VARIABLES:
        if v.name not in a.values or v.name not in b.values:
            continue
        xa, xb = usable(a, v)[ia], usable(b, v)[ib]
        valid = np.isfinite(xa) & np.isfinite(xb) & (xa > UNDEFINED / 2) & (xb > UNDEFINED / 2)
        if not valid.any():
            continue
        d = xa[valid] - xb[valid]
        scale = np.maximum(np.maximum(np.abs(xa[valid]), np.abs(xb[valid])), 1.0)
        residuals[v.name] = {"events": int(valid.sum()), "fraction_different": float((np.abs(d) > RESIDUAL_TOLERANCE * scale).mean()),
                             "median_difference": float(np.median(d)), "p99_abs_difference": float(np.quantile(np.abs(d), .99))}
    cell["residuals"] = residuals
    cell["flags"] = {g: {f: {"passed": int(e.flags[f][s].sum()), "total": int(s.sum())} for f in e.flags}
                     for g, e, s in ((ga, a, sa), (gb, b, sb))}
    return cell, (sa, sb)


def compare_distributions(a, b, region, channel, sample, selections, leaf):
    """Shapes, weighted sums and plots of every common variable."""
    sa, sb = selections
    modes = ("unweighted",) if sample == "Data_Run2024H" else MC_MODES
    result = {}
    for mode in modes:
        wa = np.ones(sa.sum()) if mode == "unweighted" else a.weights[mode][sa]
        wb = np.ones(sb.sum()) if mode == "unweighted" else b.weights[mode][sb]
        entry = {"sum_weights": {a.group: float(np.nansum(wa)), b.group: float(np.nansum(wb))}, "variables": {}}
        folder = leaf / mode
        folder.mkdir(parents=True, exist_ok=True)
        for v in VARIABLES:
            if v.name not in a.values or v.name not in b.values:
                continue
            ha, na, ma, sda = histogram(usable(a, v)[sa], wa, v.bins)
            hb, nb, mb, sdb = histogram(usable(b, v)[sb], wb, v.bins)
            entry["variables"][v.name] = {"title": v.title, "entries": {a.group: na, b.group: nb},
                                          "mean": {a.group: ma, b.group: mb}, "std": {a.group: sda, b.group: sdb},
                                          "shape_distance": distance(ha, hb)}
            draw(ha, hb, v.bins, v.title, (f"{a.group} ({na:,})", f"{b.group} ({nb:,})"), f"{sample} - {region} - {channel} - {mode}", folder / v.name)
        result[mode] = entry
    return result


# ----------------------------------------------------------------------------------------------
# Markdown
# ----------------------------------------------------------------------------------------------

def fmt(x, digits=4):
    return "n/a" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{digits}g}"


def explain(r, region, channel, sample, ga, gb):
    cell = r["cells"][region][channel]
    s = cell["selected"]
    lines = []
    a_total, b_total = s[ga], s[gb]
    both = cell["both"]
    if not a_total or not b_total:
        return [f"- One group has no event in this cell ({ga}: {a_total:,}, {gb}: {b_total:,}), so nothing can be compared."]
    lines.append(f"- {ga} selects {a_total:,} events and {gb} selects {b_total:,} ({a_total / b_total:.3f} ratio); "
                 f"{both:,} events are selected by both ({both / a_total:.1%} of {ga}, {both / b_total:.1%} of {gb}).")
    for x, y in ((ga, gb), (gb, ga)):
        missing, present = cell[f"{x}_only_absent_from_{y}"], cell[f"{x}_only_present_in_{y}"]
        if missing + present:
            lines.append(f"- {x} selects {missing + present:,} events that {y} does not: {missing:,} are not in the {y} tuple at all "
                         f"(input coverage or preselection difference) and {present:,} are in the {y} tuple but fail this region/channel "
                         "(selection or definition difference).")
        top = list(cell["migration"][x].items())[:3]
        if top:
            lines.append(f"  - where those {present:,} {x}-only events go in {y}: " + "; ".join(f"{n:,} to *{label}*" for label, n in top) + ".")
    modes = r["distributions"].get(region, {}).get(channel, {})
    reference = "unweighted" if sample == "Data_Run2024H" else "genWeights"
    ranking = sorted(((v["shape_distance"], k) for k, v in modes.get(reference, {}).get("variables", {}).items()
                      if math.isfinite(v["shape_distance"])), reverse=True)
    if ranking:
        top = ", ".join(f"`{k}` ({d:.3f})" for d, k in ranking[:3])
        lines.append(f"- The three variables whose unit-normalised shapes differ most are {top} (total-variation distance; 0 = identical).")
        leaders = {k for _, k in ranking[:3]}
        if leaders & {"m_mumu", "m_mumu_uncorrected", "mu1_pt", "mu2_pt", "mu1_pt_nano", "mu2_pt_nano", "pt_mumu"}:
            lines.append("  - Muon/dimuon quantities lead: look at muon pairing, pT thresholds, momentum calibration (ScaRe/bsc) and FSR recovery.")
        if leaders & {"njets", "jet1_pt", "jet2_pt", "mjj", "deta_jj", "dphi_jj", "pt_jj", "eta_jj", "n_bjets_loose", "n_bjets_medium"}:
            lines.append("  - Jet quantities lead: look at JEC/JER, jet ID, veto maps, lepton cleaning, pT thresholds, b-veto and the choice of the dijet pair.")
        if leaders & {"met_pt", "met_phi"}:
            lines.append("  - MET leads: the groups probably use different MET collections or corrections (TAMU/Marmara use PuppiMET).")
    residuals = cell["residuals"]
    bad = sorted(((v["fraction_different"], k) for k, v in residuals.items() if v["fraction_different"] > .01), reverse=True)
    if both:
        if bad:
            lines.append("- On events selected by both groups, these variables differ by more than 0.1% in over 1% of events: "
                         + ", ".join(f"`{k}` ({f:.1%})" for f, k in bad[:5])
                         + ". The same event is therefore reconstructed differently, not just selected differently.")
        elif residuals:
            lines.append("- On events selected by both groups every common variable agrees to better than 0.1% in more than 99% of the events.")
    return lines


def leaf_markdown(r, region, channel, sample, path):
    ga, gb = r["groups"]
    cell = r["cells"][region][channel]
    out = [f"# {ga} vs {gb}: {sample}, {region}, {channel}", "",
           f"Region: {REGIONS[region][1]}. Channel: `{channel}` as defined by each group (below).", "",
           "## Summary of the differences", "", *explain(r, region, channel, sample, ga, gb), "",
           "## Native-to-common treatment", ""]
    for g in (ga, gb):
        spec = r["treatment"][g]
        out += [f"- **{g}**: mass `{spec['mass']}`, channel `{spec['channels'][channel]}`"]
    out += ["", "## Event accounting", "", f"| | {ga} | {gb} |", "|---|---:|---:|",
            f"| events selected | {cell['selected'][ga]:,} | {cell['selected'][gb]:,} |",
            f"| selected by both | {cell['both']:,} | {cell['both']:,} |",
            f"| selected only by this group, absent from the other tuple | {cell[f'{ga}_only_absent_from_{gb}']:,} | {cell[f'{gb}_only_absent_from_{ga}']:,} |",
            f"| selected only by this group, present in the other tuple | {cell[f'{ga}_only_present_in_{gb}']:,} | {cell[f'{gb}_only_present_in_{ga}']:,} |", ""]
    for g, other in ((ga, gb), (gb, ga)):
        if cell["migration"][g]:
            out += [f"Where the {g}-only events that are present in the {other} tuple go in {other}:", "",
                    "| destination in " + other + " | events |", "|---|---:|"]
            out += [f"| {label} | {n:,} |" for label, n in list(cell["migration"][g].items())[:8]] + [""]
    out += ["## Event-by-event agreement on events selected by both", "",
            "| variable | events | fraction differing by > 0.1% | median difference | 99th percentile of the absolute difference |",
            "|---|---:|---:|---:|---:|"]
    out += [f"| `{k}` | {v['events']:,} | {v['fraction_different']:.3%} | {fmt(v['median_difference'])} | {fmt(v['p99_abs_difference'])} |"
            for k, v in cell["residuals"].items()]
    modes = r["distributions"].get(region, {}).get(channel)
    if modes is None:
        out += ["", "## Distributions", "", "Data are blinded in the signal window: only event accounting is reported.", ""]
    else:
        for mode, entry in modes.items():
            out += ["", f"## Distributions: {mode}", "",
                    f"Sum of weights: {ga} {fmt(entry['sum_weights'][ga], 6)}, {gb} {fmt(entry['sum_weights'][gb], 6)} "
                    f"(ratio {fmt(entry['sum_weights'][ga] / entry['sum_weights'][gb] if entry['sum_weights'][gb] else None, 5)}). "
                    "The groups normalise their weights differently, so the ratio of the sums is not a yield comparison. Shapes are unit-normalised; the total-variation distance is 0 for identical shapes and 1 for disjoint ones.", "",
                    f"| variable | {ga} mean | {gb} mean | {ga} std | {gb} std | shape distance | plot |", "|---|---:|---:|---:|---:|---:|---|"]
            for k, v in sorted(entry["variables"].items(), key=lambda kv: -(kv[1]["shape_distance"] if math.isfinite(kv[1]["shape_distance"]) else -1)):
                out.append(f"| `{k}` | {fmt(v['mean'][ga])} | {fmt(v['mean'][gb])} | {fmt(v['std'][ga])} | {fmt(v['std'][gb])} | "
                           f"{fmt(v['shape_distance'])} | [png]({mode}/{k}.png) / [pdf]({mode}/{k}.pdf) |")
            out += [""] + [f"![{k}]({mode}/{k}.png)" for k in list(entry["variables"])[:6]] + [""]
    out += ["## Selection flags among the selected events", "", "| group | flag | passing / selected |", "|---|---|---:|"]
    for g in (ga, gb):
        out += [f"| {g} | `{f}` | {v['passed']:,} / {v['total']:,} |" for f, v in cell["flags"][g].items()]
    out += ["", "## Limits of this diagnosis", "",
            "- Flags with similar names are not assumed to have identical definitions.",
            "- Every group's corrected quantities are used as published; shapes are unit-normalised, so overall normalisation differences do not show up in the plots.",
            "- A definitive cause requires the object-level (muon index, first-failed cut) comparison on the events selected by only one group.", ""]
    for g in (ga, gb):
        out += [f"- {g}: {note}" for note in r["treatment"][g]["notes"]]
    path.write_text("\n".join(out) + "\n")


def compare_pair(ga, gb, sample, output, cache_dir, regions, channels):
    a, b = extract(ga, sample, cache_dir), extract(gb, sample, cache_dir)
    root = output / f"{ga}-{gb}" / sample
    root.mkdir(parents=True, exist_ok=True)
    r = {"groups": [ga, gb], "sample": sample, "files": [", ".join(str(p) for p in sample_paths(treat(g), sample)) for g in (ga, gb)],
         "treatment": {g: dict(vars(treat(g))) for g in (ga, gb)},
         "entries": [a.entries, b.entries], "unique_events": [len(a.full_keys), len(b.full_keys)],
         "common_events": int(np.isin(a.full_keys, b.full_keys).sum()), "cells": {}, "distributions": {}}
    for region in regions:
        r["cells"][region], r["distributions"][region] = {}, {}
        region_dir = root / region
        for channel in channels:
            leaf = region_dir / channel
            leaf.mkdir(parents=True, exist_ok=True)
            cell, selections = compare_cell(a, b, region, channel)
            r["cells"][region][channel] = cell
            if (sample, region) not in BLINDED:
                r["distributions"][region][channel] = compare_distributions(a, b, region, channel, sample, selections, leaf)
            else:
                r["distributions"][region].pop(channel, None)
            leaf_markdown(r, region, channel, sample, leaf / "SUMMARY.md")
        write_region_summary(r, region, region_dir / "SUMMARY.md", channels)
    (root / "summary.json").write_text(json.dumps(r, indent=2, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)) + "\n")
    write_sample_summary(r, root / "SUMMARY.md", regions, channels)
    return r


def write_region_summary(r, region, path, channels):
    ga, gb = r["groups"]
    out = [f"# {ga} vs {gb}: {r['sample']}, {region}", "", f"Region: {REGIONS[region][1]}.", "",
           f"| channel | {ga} | {gb} | both | {ga} only | {gb} only | report |", "|---|---:|---:|---:|---:|---:|---|"]
    for channel in channels:
        c = r["cells"][region][channel]
        out.append(f"| {channel} | {c['selected'][ga]:,} | {c['selected'][gb]:,} | {c['both']:,} | "
                   f"{c['selected'][ga] - c['both']:,} | {c['selected'][gb] - c['both']:,} | [{channel}]({channel}/SUMMARY.md) |")
    path.write_text("\n".join(out) + "\n")


def write_sample_summary(r, path, regions, channels):
    ga, gb = r["groups"]
    ea, eb = r["unique_events"]
    c = r["common_events"]
    out = [f"# {ga} vs {gb}: {r['sample']}", "", f"- {ga}: `{r['files'][0]}`", f"- {gb}: `{r['files'][1]}`", "",
           "## Input coverage", "",
           f"{ga} has **{ea:,}** unique events, {gb} has **{eb:,}**, and **{c:,}** are in both "
           f"({c / ea:.1%} of {ga}, {c / eb:.1%} of {gb}).  A low overlap is a difference of input coverage or loose preselection, "
           "not of reconstruction; it is kept separate from the selection migration reported per region and channel.", "",
           "## Regions and channels", "",
           f"| region | channel | {ga} | {gb} | both | report |", "|---|---|---:|---:|---:|---|"]
    for region in regions:
        for channel in channels:
            x = r["cells"][region][channel]
            out.append(f"| {region} | {channel} | {x['selected'][ga]:,} | {x['selected'][gb]:,} | {x['both']:,} | [report]({region}/{channel}/SUMMARY.md) |")
    out += [""] + [f"- {g}: {n}" for g in (ga, gb) for n in r["treatment"][g]["notes"]]
    path.write_text("\n".join(out) + "\n")


def write_root_summary(output):
    rows = []
    for summary in sorted(output.glob("*/*/summary.json")):
        r = json.loads(summary.read_text())
        ga, gb = r["groups"]
        ea, eb = r["unique_events"]
        rows.append(f"| [{ga}-{gb}]({ga}-{gb}/) | [{r['sample']}]({ga}-{gb}/{r['sample']}/SUMMARY.md) | {ea:,} | {eb:,} | "
                    f"{r['common_events']:,} | {r['common_events'] / ea:.1%} | {r['common_events'] / eb:.1%} |")
    text = ["# Run 3 H->mumu synchronisation between groups", "",
            "Layout: `GROUP-GROUP/<sample>/<region>/<channel>/` with plots of every variable available in both tuples and a "
            "`SUMMARY.md` that lists the differences and explains them.  Regions: `Z_sideband`, `H_sideband`, `Signal_Fit`; "
            "channels: `ggH` and `VBF`, each as defined by the group (see the `Native-to-common treatment` of every report).", "",
            "Raw tuple counts are not yields: the groups use different loose preselections.  Selected events are therefore always "
            "counted per region and channel, and events absent from one tuple are separated from events present but classified differently.", "",
            "| pair | sample | unique events A | unique events B | common | common / A | common / B |", "|---|---|---:|---:|---:|---:|---:|", *rows, "",
            "Regenerate with `sync/python/compare_run3_sync_groups.py`; `--group1/--group2` restrict the pair and a new group is added with a `treat_<group>_ntuple` function.", ""]
    (output / "SUMMARY.md").write_text("\n".join(text))


# ----------------------------------------------------------------------------------------------
# Pair summary: quantitative and qualitative reading of one GROUP1-GROUP2 folder
# ----------------------------------------------------------------------------------------------

FAMILIES = {
    "dimuon mass": {"m_mumu", "m_mumu_uncorrected"},
    "dimuon kinematics": {"pt_mumu", "eta_mumu", "phi_mumu"},
    "muon": {"mu1_pt", "mu2_pt", "mu1_pt_nano", "mu2_pt_nano", "mu1_eta", "mu2_eta", "mu1_phi", "mu2_phi"},
    "MET": {"met_pt", "met_phi"},
    "jet multiplicity and b-tag": {"njets", "n_bjets_loose", "n_bjets_medium", "soft_activity_njets5"},
    "jet kinematics": {"jet1_pt", "jet1_eta", "jet1_phi", "jet2_pt", "jet2_eta", "jet2_phi"},
    "dijet system": {"mjj", "pt_jj", "eta_jj", "deta_jj", "dphi_jj"},
}
FAMILY_CAUSES = {
    "dimuon mass": "muon momentum calibration (Rochester/ScaRe/bsc), FSR recovery, or a different muon pair chosen for the dimuon",
    "dimuon kinematics": "different muon pair or different muon corrections; a different pairing also moves pT, eta and phi of the dimuon",
    "muon": "muon pT/eta thresholds, muon ID and isolation working points, momentum corrections, or the choice of the two muons",
    "MET": "different MET collection (PuppiMET or PFMET) or different MET corrections",
    "jet multiplicity and b-tag": "jet ID, pT/eta thresholds, lepton-jet cleaning, pileup jet ID, JEC/JER, or the b-tag working point",
    "jet kinematics": "JEC/JER versions, jet pT thresholds, jet veto maps, lepton cleaning, or which jets are called leading and subleading",
    "dijet system": "the choice of the dijet pair and the jet thresholds; the ggH/VBF split can also sculpt these shapes",
}
SYNC, MOSTLY, NOT = 0.99, 0.90, 0.0


def family_of(name):
    return next((f for f, members in FAMILIES.items() if name in members), "other")


def verdict(overlap):
    return "in sync" if overlap >= SYNC else "mostly in sync" if overlap >= MOSTLY else "NOT in sync"


def cell_numbers(r, region, channel):
    ga, gb = r["groups"]
    c = r["cells"][region][channel]
    sa, sb, both = c["selected"][ga], c["selected"][gb], c["both"]
    only_a_absent, only_a_present = c[f"{ga}_only_absent_from_{gb}"], c[f"{ga}_only_present_in_{gb}"]
    only_b_absent, only_b_present = c[f"{gb}_only_absent_from_{ga}"], c[f"{gb}_only_present_in_{ga}"]
    modes = r["distributions"].get(region, {}).get(channel, {})
    reference = "unweighted" if r["sample"] == "Data_Run2024H" else "genWeights"
    variables = modes.get(reference, {}).get("variables", {})
    distances = {k: v["shape_distance"] for k, v in variables.items() if math.isfinite(v["shape_distance"])}
    worst = sorted(distances.items(), key=lambda kv: -kv[1])
    residual = {k: v["fraction_different"] for k, v in c["residuals"].items()}
    return {"a": sa, "b": sb, "both": both, "overlap": both / max(min(sa, sb), 1), "ratio": sa / sb if sb else math.nan,
            "only_a": sa - both, "only_b": sb - both, "only_a_absent": only_a_absent, "only_a_present": only_a_present,
            "only_b_absent": only_b_absent, "only_b_present": only_b_present, "distances": distances, "worst": worst,
            "residual": residual, "weights": {m: e["sum_weights"] for m, e in modes.items()},
            "migration": c["migration"]}


def sample_reading(r, ga, gb):
    """Qualitative sentences for one sample: what differs and the likely causes."""
    sample = r["sample"]
    lines = []
    ea, eb = r["unique_events"]
    com = r["common_events"]
    lines.append(f"- **Input coverage**: {ga} {ea:,} events, {gb} {eb:,}, {com:,} in both ({com / ea:.1%} of {ga}, {com / eb:.1%} of {gb}).")
    if min(com / ea, com / eb) < 0.95:
        lines.append("  - The tuples do not contain the same events. Candidate causes: different loose preselection (muon pT or mass window, "
                     "trigger requirement applied or only stored), different input datasets or runs, or a different golden-JSON handling.")
    cells = {(rg, ch): cell_numbers(r, rg, ch) for rg in r["cells"] for ch in r["cells"][rg]}
    for (rg, ch), n in cells.items():
        if not n["a"] or not n["b"]:
            continue
        absent = n["only_a_absent"] + n["only_b_absent"]
        present = n["only_a_present"] + n["only_b_present"]
        one_sided = absent + present
        if verdict(n["overlap"]) == "in sync" or one_sided == 0:
            continue
        kind = ("mostly input coverage (events missing from the other tuple)" if absent > 2 * present else
                "mostly classification (events present in both tuples but put in a different region/channel)" if present > 2 * absent
                else "a mix of input coverage and classification")
        top = []
        for g, other in ((ga, gb), (gb, ga)):
            for label, count in list(n["migration"][g].items())[:1]:
                top.append(f"{g}-only events land in *{label}* of {other} ({count:,})")
        lines.append(f"- **{rg} / {ch}**: {verdict(n['overlap'])} ({n['overlap']:.1%} of the smaller selection is common); "
                     f"the {one_sided:,} one-sided events are {kind}" + (": " + "; ".join(top) if top else "") + ".")
    shapes = {}
    for (rg, ch), n in cells.items():
        for var, d in n["distances"].items():
            if d >= 0.05:
                shapes.setdefault(family_of(var), []).append((d, var, rg, ch))
    for family, items in sorted(shapes.items(), key=lambda kv: -max(i[0] for i in kv[1])):
        items.sort(reverse=True)
        d, var, rg, ch = items[0]
        lines.append(f"- **Shape difference in {family}** (total-variation distance > 0.05 in {len(items)} cells; worst `{var}` = {d:.3f} in {rg}/{ch}). "
                     f"Possible causes: {FAMILY_CAUSES.get(family, 'not classified')}.")
    resid = {}
    for (rg, ch), n in cells.items():
        for var, f in n["residual"].items():
            if f > 0.01:
                resid.setdefault(family_of(var), []).append((f, var, rg, ch))
    for family, items in sorted(resid.items(), key=lambda kv: -max(i[0] for i in kv[1])):
        items.sort(reverse=True)
        f, var, rg, ch = items[0]
        lines.append(f"- **Same event, different value in {family}**: on events selected by both groups `{var}` differs by more than 0.1% in {f:.1%} "
                     f"of the events ({rg}/{ch}). This is a reconstruction or correction difference, not a selection one.")
    for obj, pt, angles in (("jets", ("jet1_pt", "jet2_pt"), ("jet1_eta", "jet2_eta")), ("muons", ("mu1_pt", "mu2_pt"), ("mu1_eta", "mu2_eta"))):
        worst = {}
        for (rg, ch), n in cells.items():
            if n["both"] and all(k in n["residual"] for k in pt + angles):
                worst[(rg, ch)] = (min(n["residual"][k] for k in pt), max(n["residual"][k] for k in angles))
        if worst and all(p_ > 0.5 and a_ < 0.2 for p_, a_ in worst.values()):
            lines.append(f"- **Same {obj}, different pT scale**: on common events the {obj} directions (eta) agree while their pT differs in more than 50% of the events. "
                         "The groups are using the same objects but a different energy/momentum correction (jet energy corrections and their version, or muon momentum calibration), "
                         "so the pT-dependent selections and all derived masses move together.")
    if not resid and any(n["both"] for n in cells.values()):
        lines.append("- On the events selected by both groups every common variable agrees to better than 0.1% in more than 99% of the events: the reconstruction is the same.")
    ratios = [n["weights"].get("genWeights", {}) for n in cells.values() if n["weights"]]
    for rg, ch in (("Signal_Fit", "VBF"), ("Signal_Fit", "ggH")):
        n = cells.get((rg, ch))
        if n and sample != "Data_Run2024H" and n["weights"].get("all_common_weights"):
            w = n["weights"]
            if all(w[m].get(gb) for m in w):
                lines.append(f"- **Weights in {rg}/{ch}**: sum of weights ratio {ga}/{gb} = "
                             f"{w['genWeights'][ga] / w['genWeights'][gb]:.3f} with generator weights only and "
                             f"{w['all_common_weights'][ga] / w['all_common_weights'][gb]:.3f} with all weights. A large change between the two points to the "
                             "scale factors (pileup, muon ID/iso/trigger) or to the luminosity/cross-section normalisation, which each group applies differently.")
                break
    return lines


def write_pair_summary(output, ga, gb):
    folder = output / f"{ga}-{gb}"
    reports = [json.loads(p.read_text()) for p in sorted(folder.glob("*/summary.json"))]
    if not reports:
        return
    out = [f"# {ga} vs {gb}: summary of the differences", "",
           "Generated by `compare_run3_sync_groups.py` from the `summary.json` of every sample. Quantitative tables first, then a qualitative reading "
           "with the likely causes. The causes are hypotheses from the pattern of the differences and from the groups' published definitions: "
           "a proof needs the object-level comparison (muon indices, first failed cut) of the one-sided events.", "",
           "Thresholds: *in sync* means at least 99% of the smaller selection is common, *mostly in sync* 90-99%. A shape distance is the total-variation "
           "distance of the unit-normalised distributions (0 identical, 1 disjoint); above 0.05 is a visible difference.", ""]
    if "TAMU" in (ga, gb):
        out += ["> **Provisional**: this pair uses the TAMU exchange files published before the Pisa-like regeneration "
                "(variant `ggfvbf_20261015`). It must be regenerated with the new `ggfvbf_pisa_20261013` tuples.", ""]
    out += ["## Quantitative overview", "",
            f"| sample | region | channel | {ga} | {gb} | both | {ga}/{gb} | common / smaller | status | {ga} only (absent / present) | {gb} only (absent / present) | largest shape distance | largest event-by-event mismatch |",
            "|---|---|---|---:|---:|---:|---:|---:|---|---|---|---|---|"]
    for r in reports:
        for rg in r["cells"]:
            for ch in r["cells"][rg]:
                n = cell_numbers(r, rg, ch)
                worst = f"`{n['worst'][0][0]}` {n['worst'][0][1]:.3f}" if n["worst"] else "n/a (blinded)" if r["sample"] == "Data_Run2024H" and rg == "Signal_Fit" else "n/a"
                mism = max(n["residual"].items(), key=lambda kv: kv[1], default=None)
                mism = f"`{mism[0]}` {mism[1]:.1%}" if mism else "n/a"
                out.append(f"| {r['sample']} | {rg} | {ch} | {n['a']:,} | {n['b']:,} | {n['both']:,} | {n['ratio']:.3f} | {n['overlap']:.1%} | {verdict(n['overlap']) if n['a'] and n['b'] else 'empty'} | "
                           f"{n['only_a_absent']:,} / {n['only_a_present']:,} | {n['only_b_absent']:,} / {n['only_b_present']:,} | {worst} | {mism} |")
    out += ["", "## Qualitative reading, sample by sample", ""]
    for r in reports:
        out += [f"### {r['sample']}", "", *sample_reading(r, ga, gb), "", f"Details: [{r['sample']}/SUMMARY.md]({r['sample']}/SUMMARY.md)", ""]
    out += ["## Definitions that can explain differences", ""]
    for g in (ga, gb):
        out += [f"- **{g}**: {note}" for note in reports[0]["treatment"][g]["notes"]]
    out += ["", "## How to read the cells", "",
            "- *absent*: the event is not in the other group's tuple (input coverage or preselection). *present*: it is in the other tuple but is not selected in this region/channel (definition or selection).",
            "- Data in `Signal_Fit` are blinded, so only counts are given there.", ""]
    (folder / "SUMMARY.md").write_text("\n".join(out))


# ----------------------------------------------------------------------------------------------
# All groups on the same plot
# ----------------------------------------------------------------------------------------------

def draw_all(hists, bins, title, heading, out):
    """hists: ordered {group: (histogram, entries)}; the ratio is to the first group."""
    width = np.diff(bins)
    densities = {g: (h / h.sum() / width if h.sum() else h, n) for g, (h, n) in hists.items()}
    reference = next(iter(densities))
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": .05})
    centres = 0.5 * (bins[1:] + bins[:-1])
    rax.axhline(1, color="black", linewidth=.8)
    spread = 0.
    for i, (g, (y, n)) in enumerate(densities.items()):
        ax.stairs(y, bins, label=f"{g} ({n:,})", linewidth=1.7, color=f"C{i}")
        if g == reference:
            continue
        yr = densities[reference][0]
        ratio = np.divide(y, yr, out=np.full_like(y, np.nan), where=(y > 0) & (yr > 0))
        finite = ratio[np.isfinite(ratio)]
        if finite.size:
            spread = max(spread, float(np.max(np.abs(finite - 1))))
        rax.errorbar(centres + 0.0, ratio, xerr=width / 2, fmt="o", color=f"C{i}", markersize=3, linestyle="none", elinewidth=1)
    half = min(max(1.2 * spread, 0.02), 0.5)
    rax.set_ylim(1 - half, 1 + half)
    ax.set_ylabel("unit-normalised density")
    ax.set_title(heading, fontsize=10)
    ax.legend()
    rax.set_ylabel(f"group / {reference}")
    rax.set_xlabel(title)
    fig.savefig(out.with_suffix(".png"), dpi=130, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def compare_all_groups(sample, output, cache_dir, regions, channels, groups):
    """Every group that provides the sample on one plot per variable, region, channel and weight mode."""
    present = [g for g in groups if sample in treat(g).files]
    if len(present) < 2:
        return None
    extracts = {g: extract(g, sample, cache_dir) for g in present}
    root = output / "ALL_GROUPS" / sample
    root.mkdir(parents=True, exist_ok=True)
    info = {"sample": sample, "groups": present, "cells": {}}
    out = [f"# All groups together: {sample}", "",
           f"Groups: {', '.join(present)}. Every plot overlays the unit-normalised distribution of each group; the ratio panel is to **{present[0]}**. "
           "Every group uses its own definition of the channel and of the corrected mass (see the pair reports for the details).", ""]
    for region in regions:
        for channel in channels:
            selections = {g: selection(e, region, channel) for g, e in extracts.items()}
            keys = {g: extracts[g].keys[s] for g, s in selections.items()}
            common = keys[present[0]]
            for g in present[1:]:
                common = np.intersect1d(common, keys[g])
            row = {"selected": {g: int(s.sum()) for g, s in selections.items()}, "common_to_all": int(common.size), "distances": {}}
            info["cells"].setdefault(region, {})[channel] = row
            out += [f"## {region} / {channel}", "", "| group | selected | common to all groups |", "|---|---:|---:|"]
            out += [f"| {g} | {row['selected'][g]:,} | {common.size:,} |" for g in present] + [""]
            if (sample, region) in BLINDED:
                out += ["Data are blinded in the signal window: counts only.", ""]
                continue
            modes = ("unweighted",) if sample == "Data_Run2024H" else MC_MODES
            for mode in modes:
                folder = root / region / channel / mode
                folder.mkdir(parents=True, exist_ok=True)
                table = []
                for v in VARIABLES:
                    hists = {}
                    for g, e in extracts.items():
                        if v.name not in e.values:
                            continue
                        w = np.ones(selections[g].sum()) if mode == "unweighted" else e.weights[mode][selections[g]]
                        h, n, _, _ = histogram(usable(e, v)[selections[g]], w, v.bins)
                        hists[g] = (h, n)
                    if len(hists) < 2:
                        continue
                    reference = next(iter(hists))
                    d = {g: distance(hists[reference][0], h) for g, (h, _) in hists.items() if g != reference}
                    table.append((v.name, reference, d))
                    draw_all(hists, v.bins, v.title, f"{sample} - {region} - {channel} - {mode}", folder / v.name)
                row["distances"][mode] = {name: d for name, _, d in table}
                if mode == ("unweighted" if sample == "Data_Run2024H" else "genWeights") and table:
                    others = [g for g in present if g != present[0]]
                    out += [f"Shape distance to {present[0]} (0 identical, 1 disjoint), {mode}:", "",
                            "| variable | " + " | ".join(others) + " | plot |", "|---|" + "---:|" * len(others) + "---|"]
                    for name, _, d in sorted(table, key=lambda t: -max([x for x in t[2].values() if math.isfinite(x)] or [0])):
                        out.append(f"| `{name}` | " + " | ".join(fmt(d.get(g)) for g in others) + f" | [png]({region}/{channel}/{mode}/{name}.png) |")
                    out.append("")
    (root / "SUMMARY.md").write_text("\n".join(out) + "\n")
    (root / "overlay.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", type=Path)
    p.add_argument("--cache-dir", type=Path, default=Path("/tmp/vdamante/run3_sync_extracts"))
    p.add_argument("--group1", choices=TREATERS)
    p.add_argument("--group2", choices=TREATERS)
    p.add_argument("--sample", action="append", help="restrict to this sample (repeatable)")
    p.add_argument("--region", action="append", choices=REGIONS)
    p.add_argument("--channel", action="append", choices=CHANNELS)
    p.add_argument("--resume", action="store_true", help="keep pair/sample directories that already have a summary.json")
    p.add_argument("--summaries-only", action="store_true", help="rewrite the pair summaries and the root index from the existing summary.json files")
    p.add_argument("--no-overlay", action="store_true", help="skip the all-groups overlay plots")
    p.add_argument("--list", action="store_true", help="list the groups, samples and pairs and exit")
    a = p.parse_args()
    if bool(a.group1) != bool(a.group2):
        p.error("--group1 and --group2 must be given together")
    if a.group1 and a.group1 == a.group2:
        p.error("--group1 and --group2 must differ")
    groups = list(TREATERS)
    pairs = [(a.group1, a.group2)] if a.group1 else [(x, y) for i, x in enumerate(groups) for y in groups[i + 1:]]
    work = [(x, y, s) for x, y in pairs for s in sorted(set(treat(x).files) & set(treat(y).files)) if not a.sample or s in a.sample]
    if a.list:
        for x, y, s in work:
            print(f"{x}-{y}/{s}")
        return
    if not a.output:
        p.error("--output is required")
    a.output.mkdir(parents=True, exist_ok=True)
    if a.summaries_only:
        work = []
    for x, y, s in work:
        if a.resume and (a.output / f"{x}-{y}" / s / "summary.json").exists():
            print(f"[resume] {x}-{y} {s}")
            continue
        print(f"[sync] {x}-{y} {s}", flush=True)
        compare_pair(x, y, s, a.output, a.cache_dir, a.region or list(REGIONS), a.channel or list(CHANNELS))
    for x, y in pairs:
        write_pair_summary(a.output, x, y)
    if not a.group1 and not a.summaries_only and not a.no_overlay:
        for s in sorted({s for g in groups for s in treat(g).files if not a.sample or s in a.sample}):
            print(f"[overlay] {s}", flush=True)
            compare_all_groups(s, a.output, a.cache_dir, a.region or list(REGIONS), a.channel or list(CHANNELS), groups)
    write_root_summary(a.output)


if __name__ == "__main__":
    main()
