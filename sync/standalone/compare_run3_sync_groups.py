#!/usr/bin/env python3
"""Compare the Run 3 H->mumu exchange tuples of two groups, region by region and channel by channel.

Every group publishes its own ntuple schema in /eos/cms/store/group/phys_higgs/cmshmm/run3_sync.
A ``treat_<group>_ntuple`` function projects one group onto a common description
(mass, ggH/VBF channel, weights, flags and the common variables); the comparison itself only
sees that description.  A new group is added by writing one GroupSpec and its treat function.

Layout of the output:

    <output>/<GROUP1>-<GROUP2>/<sample>/<region>/<channel>/SUMMARY.md + plots

  export HMUMU_SYNC_VARIANT=update_after_Hmm_15Oct     # the tuples published after the Hmm meeting of 15 Oct 2026
  python3 compare_run3_sync_groups.py --group1 TAMU --group2 Pisa --output pairs --cache-dir cache --no-overlay
  python3 compare_run3_sync_groups.py --group1 TAMU --group2 Iowa --sample Data_Run2024H --output pairs
  python3 compare_run3_sync_groups.py --list

Requires numpy, awkward>=2, uproot>=5, matplotlib. HMUMU_SYNC_BASE sets the folder with one sub-folder per group.
For yield tables and bar charts only, use the simpler sync_compare.py.
"""
from __future__ import annotations

import os
import argparse
import datetime
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

BASE = Path(os.environ.get("HMUMU_SYNC_BASE", "/eos/cms/store/group/phys_higgs/cmshmm/run3_sync"))   # directory with one folder per group
CHANNELS = ("ggH", "VBF", "inclusive")
BASE_CHANNELS = ("ggH", "VBF")      # the groups' own categories; inclusive = ggH | VBF
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
HELPERS = {"abs": np.abs, "sqrt": np.sqrt, "pi": math.pi, "cos": np.cos, "sin": np.sin, "sinh": np.sinh, "cosh": np.cosh,
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
                                                      "Purdue": PURDUE_EXPR.get(name), "MIT": MIT_EXPR.get(name),
                                                      "Caltech": CALTECH_EXPR.get(name), "Iowa": IOWA_EXPR.get(name)}, min_jets)


# Caltech branches (`update_after_Hmm_15Oct` only, raw-level): `*_bscPt`/`diMuon_bsc*` are the beam-spot-constrained values, `*_noCorr` the
# uncorrected ones (the update variant compares the raw mass and pT); the jets are the two leading ones (`leading_jet_*`, `diJet_*`)
CALTECH_EXPR = {
    "m_mumu": "diMuon_bscMass_noCorr", "m_mumu_uncorrected": "diMuon_bscMass_noCorr", "pt_mumu": "diMuon_bscPt_noCorr",
    "eta_mumu": "diMuon_eta", "phi_mumu": "diMuon_phi", "mu1_pt": "mu1_bscPt", "mu2_pt": "mu2_bscPt",
    "mu1_pt_nano": "mu1_bscPt_noCorr", "mu2_pt_nano": "mu2_bscPt_noCorr", "mu1_eta": "mu1_eta", "mu2_eta": "mu2_eta",
    "mu1_phi": "mu1_phi", "mu2_phi": "mu2_phi", "met_phi": "PuppiMET_phi", "njets": "n_jet",
    "jet1_pt": "leading_jet_pt", "jet1_eta": "leading_jet_eta", "jet1_phi": "leading_jet_phi",
    "jet2_pt": "subleading_jet_pt", "jet2_eta": "subleading_jet_eta", "jet2_phi": "subleading_jet_phi",
    "mjj": "diJet_mass", "pt_jj": "diJet_pt", "eta_jj": "diJet_eta", "deta_jj": "abs(delta_eta_diJet)", "dphi_jj": "abs(delta_phi_diJet)",
}


# Iowa branches (`update_after_Hmm_15Oct`, tree `SyncEvents`, raw-level): the muon branches are two-entry arrays (the analysis pair, in storage order),
# so leading/subleading are ranked by `muonRawPt`; the jets are the entries with `jetPassQuality` (the selected jets, `nSelectedJets` of them)
_IOWA_MU = "muonRawPt>0"
def _iowa_mu(branch, rank):
    return f"byrank({branch}, {_IOWA_MU}, muonRawPt, {rank})"
def _iowa_jet(branch, rank):
    return f"byrank({branch}, jetPassQuality, jetCorrectedPt, {rank})"
IOWA_EXPR = {
    "m_mumu": "massRaw", "m_mumu_uncorrected": "massRaw", "pt_mumu": "ptRaw",
    "mu1_pt": _iowa_mu("muonPostFSRPt", 0), "mu2_pt": _iowa_mu("muonPostFSRPt", 1),
    "mu1_pt_nano": _iowa_mu("muonRawPt", 0), "mu2_pt_nano": _iowa_mu("muonRawPt", 1),
    "mu1_eta": _iowa_mu("muonEta", 0), "mu2_eta": _iowa_mu("muonEta", 1),
    "mu1_phi": _iowa_mu("muonPhi", 0), "mu2_phi": _iowa_mu("muonPhi", 1),
    "met_pt": "metPt", "met_phi": "metPhi", "njets": "nSelectedJets", "n_bjets_loose": "nBJetsLoose", "n_bjets_medium": "nBJetsMedium",
    "jet1_pt": _iowa_jet("jetCorrectedPt", 0), "jet2_pt": _iowa_jet("jetCorrectedPt", 1),
    "jet1_eta": _iowa_jet("jetEta", 0), "jet2_eta": _iowa_jet("jetEta", 1),
    "jet1_phi": _iowa_jet("jetPhi", 0), "jet2_phi": _iowa_jet("jetPhi", 1),
    "mjj": "leadingDijetMass", "pt_jj": "leadingDijetPt", "deta_jj": "abs(leadingDijetDEta)",
}


# MIT branches (variable name -> expression); the jets are the VBF dijet (-1 when there is no dijet)
MIT_EXPR = {
    "m_mumu": "HiggsCandCorrMass", "mu1_pt": "Muon1_pt", "mu2_pt": "Muon2_pt", "mu1_eta": "Muon1_eta", "mu2_eta": "Muon2_eta",
    "mu1_phi": "Muon1_phi", "mu2_phi": "Muon2_phi", "met_pt": "PuppiMET_pt", "met_phi": "PuppiMET_phi", "njets": "nGoodJetsAll",
    "n_bjets_medium": "nBMJets", "jet1_pt": "jetVBF1_Pt", "jet2_pt": "jetVBF2_Pt", "jet1_eta": "jetVBF1_Eta", "jet2_eta": "jetVBF2_Eta",
    "mjj": "Mjj", "deta_jj": "abs(jetVBF1_Eta - jetVBF2_Eta)",
}


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
    V("mjj", r"$m_{jj}$ [GeV]", np.linspace(0, 2000, 51), "dijet_ls_mass", "jj_mass[:,0]", "dijet_mass", "dijetmass", 2, Pisa="Dijet_mass"),
    V("pt_jj", r"$p_T^{jj}$ [GeV]", np.linspace(0, 400, 51), "dijet_ls_pt", "jj_pt", "dijet_pt", "dijetpt", 2, Pisa="Dijet_pt"),
    V("eta_jj", r"$\eta_{jj}$", np.linspace(-8, 8, 51), "dijet_ls_eta", "jj_eta", "dijet_eta", "dijeteta", 2, Pisa="Dijet_eta"),
    V("deta_jj", r"$|\Delta\eta_{jj}|$", np.linspace(0, 10, 51), "abs(leadingjet_eta-subleadingjet_eta)", "abs(jj_deta[:,0])", "abs(delta_eta_jj)", "absDetajj", 2, Pisa="abs(Dijet_DeltaEta)"),
    V("dphi_jj", r"$|\Delta\phi_{jj}|$", np.linspace(0, math.pi, 33), "abs(dphi(leadingjet_phi,subleadingjet_phi))", None, "abs(delta_phi_jj)", "abs(dijetdphi)", 2),
    V("soft_activity_njets5", "soft-activity jets, $p_T>5$ GeV", np.arange(-0.5, 15.5, 1), None, None, "SoftActivityJetNjets5", "SoftActivityJetNjets5_tmva"),
]


def _expr(group, name):
    return next(v.expr.get(group) for v in VARIABLES if v.name == name)


def _set(name, **per_group):
    next(v for v in VARIABLES if v.name == name).expr.update(per_group)


# MIT: `Muon*_rawpt` is taken as the uncorrected (NanoAOD) muon pT, the uncorrected dimuon quantities are built from it
_set("mu1_pt_nano", MIT="Muon1_rawpt")
_set("mu2_pt_nano", MIT="Muon2_rawpt")
_set("m_mumu_uncorrected", MIT="sqrt(2*Muon1_rawpt*Muon2_rawpt*(cosh(Muon1_eta-Muon2_eta)-cos(Muon1_phi-Muon2_phi)))")

_ALL = ("TAMU", "Marmara", "PKU", "NNU", "Pisa", "Purdue", "MIT", "Caltech", "Iowa")


def _derived(name, title, bins, build, min_jets=0, native=None):
    """A variable built group by group from the expressions of other variables (None if one of them is missing)."""
    expr = {g: (native or {}).get(g) or build(g) for g in _ALL}
    VARIABLES.append(Variable(name, title, np.asarray(bins), expr, min_jets))


def _pair(fn, *names):
    def build(g):
        parts = [_expr(g, n) for n in names]
        return None if any(p is None for p in parts) else fn(*[f"({p})" for p in parts])
    return build


_derived("pt_mumu_uncorrected", r"uncorrected $p_T^{\mu\mu}$ [GeV]", np.linspace(0, 250, 51),
         lambda g: None if None in (_expr(g, "mu1_pt_nano"), _expr(g, "mu2_pt_nano"), _expr(g, "mu1_phi"), _expr(g, "mu2_phi")) else
         f"sqrt(({_expr(g, 'mu1_pt_nano')})**2+({_expr(g, 'mu2_pt_nano')})**2+2*({_expr(g, 'mu1_pt_nano')})*({_expr(g, 'mu2_pt_nano')})"
         f"*cos(({_expr(g, 'mu1_phi')})-({_expr(g, 'mu2_phi')})))",
         native={"TAMU": "pt_mumu_nano", "Marmara": "mm_pt_uncorr", "NNU": "dimuon_raw_pt", "Pisa": "Dimuon_rawPt", "Purdue": "dimuon_pt_raw",
                 "Caltech": "diMuon_bscPt_noCorr", "Iowa": "ptRaw"})
_derived("mu1_charge", "leading $\mu$ charge", np.arange(-1.5, 2.5, 1), lambda g: None,
         native={"TAMU": "mu1_charge", "Marmara": "mu_charge[:,0]", "PKU": "mu1_fromH_charge", "NNU": "mu1_charge", "Pisa": "LeadingMuon_charge",
                 "Purdue": "mu1_charge", "MIT": "Muon1_charge", "Iowa": _iowa_mu("muonCharge", 0)})
_derived("mu2_charge", "subleading $\mu$ charge", np.arange(-1.5, 2.5, 1), lambda g: None,
         native={"TAMU": "mu2_charge", "Marmara": "mu_charge[:,1]", "PKU": "mu2_fromH_charge", "NNU": "mu2_charge", "Pisa": "SubLeadingMuon_charge",
                 "Purdue": "mu2_charge", "MIT": "Muon2_charge", "Iowa": _iowa_mu("muonCharge", 1)})
_derived("deta_mumu", r"$|\Delta\eta_{\mu\mu}|$", np.linspace(0, 5, 51), _pair(lambda a, b: f"abs({a}-{b})", "mu1_eta", "mu2_eta"))
_derived("dphi_mumu", r"$|\Delta\phi_{\mu\mu}|$", np.linspace(0, math.pi, 33), _pair(lambda a, b: f"abs(dphi({a},{b}))", "mu1_phi", "mu2_phi"))
_derived("dR_mumu", r"$\Delta R_{\mu\mu}$", np.linspace(0, 6, 61),
         lambda g: None if None in (_expr(g, "mu1_eta"), _expr(g, "mu2_eta"), _expr(g, "mu1_phi"), _expr(g, "mu2_phi")) else
         f"hypot(({_expr(g, 'mu1_eta')})-({_expr(g, 'mu2_eta')}),dphi({_expr(g, 'mu1_phi')},{_expr(g, 'mu2_phi')}))")
_derived("mu1_pt_over_mmumu", r"leading $\mu$ $p_T/m_{\mu\mu}$", np.linspace(0, 3, 61), _pair(lambda a, b: f"{a}/{b}", "mu1_pt", "m_mumu"))
_derived("mu2_pt_over_mmumu", r"subleading $\mu$ $p_T/m_{\mu\mu}$", np.linspace(0, 2, 41), _pair(lambda a, b: f"{a}/{b}", "mu2_pt", "m_mumu"))
_derived("zeppenfeld", r"Zeppenfeld $\eta_{\mu\mu}-(\eta_{j1}+\eta_{j2})/2$", np.linspace(-6, 6, 61),
         lambda g: None if None in (_expr(g, "eta_mumu"), _expr(g, "jet1_eta"), _expr(g, "jet2_eta")) else
         f"({_expr(g, 'eta_mumu')})-(({_expr(g, 'jet1_eta')})+({_expr(g, 'jet2_eta')}))/2", min_jets=2)


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
         "baseline", "VBF", "ggF"],
        ["Channels are the stored `ggF` / `VBF` flags (baseline, then VBF = two leading jets pass 35/25 GeV, m_jj > 400, |deta| > 2.5; "
         "ggF = baseline and not VBF).",
         "The dimuon is formed by the two candidate muons stored as `mu1` and `mu2`."])


def treat_Marmara_ntuple():
    return GroupSpec(
        "Marmara", {"Data_Run2024H": "Marmara/tree_DATA_2024_loose_for_ALL_2024_H_Data_Samples.root", "ggH": "Marmara/tree_GGH_2024_loose_LAST.root",
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
         "Data: `tree_DATA_2024_loose_for_ALL_2024_H_Data_Samples.root` (4.5M events, full Run2024H, replaces the earlier partial `larger_stat` file); `data_Run2024H_Muon*.root` are NanoAOD inputs.",
         "The `GluGluH-...root` / `VBFH-...root` files are NanoAOD too; the `tree_*_loose_LAST.root` tuples are used."])


def treat_Caltech_ntuple():
    """Caltech exchange tuples (`Caltech/update_after_Hmm_15Oct/{ggH,VBF,Muon0_2024H,Muon1_2024H}_skim.root`, tree `tree_output`)."""
    d = f"Caltech/{UPDATE}/"
    return GroupSpec(
        "Caltech", {"Data_Run2024H": [d + "Muon0_2024H_skim.root", d + "Muon1_2024H_skim.root"], "ggH": d + "ggH_skim.root", "VBFH": d + "VBF_skim.root"},
        "tree_output", "luminosityBlock", "diMuon_bscMass_noCorr",
        {"ggH": "(is_ggH_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto",
         "VBF": "(is_VBF_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto"},
        {"genWeights": "gen_weight", "all_common_weights": "gen_weight*pileup_weight*mu_trig_eff_SF*mu_ID_eff_SF*mu_ISO_eff_SF"},
        ["event_pass_muon_trigger_match", "event_pass_jet_veto", "is_ggH_category", "is_VBF_category"],
        [f"Tuple of `{UPDATE}` (raw-level): regions and m_mumu on `diMuon_bscMass_noCorr`, pT on `diMuon_bscPt_noCorr`, muon pT corrected = `mu*_bscPt`, "
         "NanoAOD-like = `mu*_bscPt_noCorr`; jets are `leading_jet_*`/`subleading_jet_*`/`diJet_*`.",
         "Channels: `is_ggH_category`/`is_VBF_category` (mutually exclusive) together with `event_pass_muon_trigger_match` and `event_pass_jet_veto`; "
         "this definition of the baseline is our assumption and has to be confirmed by Caltech. There is no b-tag flag/count in the tuple, no MET pT and no charge.",
         "Data have `gen_weight=-999` (unused); no events in the 115-135 GeV window, i.e. blinded.",
         "Weights: Caltech's `weight` = gen_weight * pileup_weight * xsec[fb] * lumi[fb-1] * mu_trig_eff_SF * mu_ID_eff_SF * mu_ISO_eff_SF / sum(gen_weights); "
         "checked: it equals our un-normalised product times a constant per sample (VBF 4.177e-06, ggH 1.365e-05). `all_common_weights` is therefore the product without "
         "the normalisation (like the other groups' un-normalised weights); `weight` is only comparable after dividing by that constant."])


def iowa_files():
    """Iowa publishes per-job files: data = Muon0 + Muon1, VBFH = VBF (no ggH sample)."""
    d = BASE / "Iowa" / UPDATE
    rel = lambda sub: sorted(str(p.relative_to(BASE)) for p in (d / sub).glob("*.root"))
    return {"Data_Run2024H": rel("Muon0") + rel("Muon1"), "VBFH": rel("VBF")}


def treat_Iowa_ntuple():
    """Iowa exchange tuples (`Iowa/update_after_Hmm_15Oct/{Muon0,Muon1,VBF}/*_sync_round2.root`, tree `SyncEvents`)."""
    base = ("passGolden & passFilters & passHLTAndObjectMatch & passAnalysisMuonSelection & passExactlyTwoAnalysisMuons & passElectronVeto & "
            "passAdditionalLeptonVeto & passJetVetoMap & passBJetVeto & passVertex & analysisPairPassesRawWindow & oppositeSign")
    return GroupSpec(
        "Iowa", iowa_files() if (BASE / "Iowa" / UPDATE).exists() else {},
        "SyncEvents", "luminosityBlock", "massRaw",
        {"ggH": f"{base} & ~passVBFCategoryApprox", "VBF": f"{base} & passVBFCategoryApprox"},
        {"genWeights": "weightGenerator", "all_common_weights": "weightGenerator*weightPileup*weightMuonIDISO*weightTrigger"},
        ["passGolden", "passFilters", "passHLTAndObjectMatch", "passAnalysisMuonSelection", "passExactlyTwoAnalysisMuons", "passElectronVeto",
         "passAdditionalLeptonVeto", "passJetVetoMap", "passBJetVeto", "passVertex", "passVBFCategoryApprox"],
        [f"Tuple of `{UPDATE}` (raw-level): regions and m_mumu on `massRaw`, pT on `ptRaw`; muon pT corrected = `muonPostFSRPt`, NanoAOD-like = `muonRawPt` "
         "(not confirmed that this is the NanoAOD pT; `muonBasePt` is also stored).",
         "Baseline = the stored pass flags (golden, filters, HLT + object match, analysis muon pair, exactly two analysis muons, electron/additional-lepton veto, "
         "jet veto map, b-veto, vertex, raw mass window, OS); VBF = baseline and `passVBFCategoryApprox`; ggH = baseline and not VBF. "
         "Both are our reading of the tuple and have to be confirmed by Iowa. There is no ggH sample.",
         "The muon arrays are in storage order; leading/subleading are ranked by `muonRawPt`. Jets are the entries with `jetPassQuality` ranked by "
         "`jetCorrectedPt`; `leadingDijet*` is the leading pair. No dijet phi/eta and no soft activity.",
         "Data blinding is on the raw mass (`massRaw`); ~20 events per file have the post-FSR mass in 115-135 GeV."])


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
        {"ggH": "(first_failed_step >= 13) & (first_failed_step <= 16) & (nMediumBTagJet == 0) & (nLooseBTagJet < 2)",
         "VBF": "(first_failed_step == 0) | (first_failed_step >= 18)"},
        {"genWeights": "Generator_weight", "all_common_weights": "Generator_weight*PuSF*MuIdIsoSF*MuTrgSF"},
        ["pass_NoiseFilters", "pass_GoldenJson", "pass_isTriggered", "pass_AtLeastOneTrgMuon", "pass_TwoMuons", "pass_TwoOppositeSignMuons",
         "pass_LeadingMuonCut", "pass_EleVeto", "pass_AdditionalMuonVeto", "pass_RawDimuonMassWindow", "pass_TwoJets", "pass_LeadingJetPt",
         "pass_DijetMass", "pass_DijetDeltaEta", "pass_BJetVeto"],
        ["Pisa publishes only the **VBF cutflow**, with no ggH category.  The tuple holds the events that reach the raw-mass window "
         "(about 1.9M per data dataset), with `first_failed_step` telling where each one stops; the VBF channel is every cut but the mass "
         "window passing.  Our reconstruction of a ggH channel (to be confirmed by Pisa): events that pass every cut before the jets "
         "(`first_failed_step` 13-16, i.e. fail the VBF jet requirements) and the b-veto recomputed from the jet counts (no medium b jet and fewer "
         "than two loose ones; the published `pass_BJetVeto` is cumulative and not evaluated for those events), the analogue of TAMU's "
         "`ggF = baseline and not VBF`.  Data in the signal window are blinded and absent.",
         "Variable mapping checked against the TAMU tuple on common events: `Dimuon_mass` is the FSR-recovered, bs-constrained, ScaRe mass "
         "(0.4% of events differ from TAMU `m_mumu_FSR_corr`); `LeadingMuon_pt` is the NanoAOD pT and `LeadingMuon_ptToUse` the corrected one; "
         "`LeadingJet_ptToUse` is the corrected jet pT (`LeadingJet_pt` is uncorrected).",
         "MET is PuppiMET.  There is no dimuon eta/phi and no dphi_jj.  `all_common_weights` is `Generator_weight*PuSF*MuIdIsoSF*MuTrgSF` (MC only): the published `HiggsMass_cumulativeWeight` is the same product but is zero for every event that does not pass the whole VBF selection, so it cannot be used outside the VBF channel and the mass window."])


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


def treat_MIT_ntuple():
    """MIT exchange tuples (`MIT/VBF-GGH/snapshot_mc_<id>_2024_ggHcat_synch.root`, tree `events`): ggH + VBF selection with the
    b-veto relaxed.  mc_10 = VBFH, mc_11 = ggH (GluGluHto2Mu, powheg), mc_-51 / mc_-52 = Run2024H Muon0 / Muon1 (info.txt of the group)."""
    d = "MIT/VBF-GGH/snapshot_mc_{}_2024_ggHcat_synch.root"
    return GroupSpec(
        "MIT", {"Data_Run2024H": [d.format("-51"), d.format("-52")], "ggH": d.format("11"), "VBFH": d.format("10")},
        "events", "luminosityBlock", "HiggsCandCorrMass",
        {"ggH": "(nBMJets == 0) & ~((Mjj > 400) & (abs(jetVBF1_Eta - jetVBF2_Eta) > 2.5))",
         "VBF": "(nBMJets == 0) & (Mjj > 400) & (abs(jetVBF1_Eta - jetVBF2_Eta) > 2.5)"},
        {"genWeights": "w", "all_common_weights": "w_allSF"},
        [],
        ["The tuple has no category flag: it stores the ggH + VBF selection with the b-veto relaxed (`nBMJets` up to 3) and the extra mu-e veto on. "
         "Our reading of the channels: `nBMJets == 0` and VBF = `Mjj > 400` and |deta(jetVBF1, jetVBF2)| > 2.5, ggH = the rest; to be confirmed by MIT.",
         "Jets: the VBF dijet (`jetVBF1/2_*`, pT > 30 GeV, -1 when absent) and `nGoodJetsAll`; `Muon1/2_rawpt` is not mapped. Region flags `isZ/isH/isHSB` are not used: the regions are rebuilt from `HiggsCandCorrMass`.",
         "Weights: `w` (cross-section, luminosity and generator weight, without scale factors) and `w_allSF` (with all scale factors); both 1 in data."])


TREATERS = {"TAMU": treat_TAMU_ntuple, "Marmara": treat_Marmara_ntuple,
            "PKU": treat_PKU_ntuple, "NNU": treat_NNU_ntuple, "Pisa": treat_Pisa_ntuple, "Purdue": treat_Purdue_ntuple, "MIT": treat_MIT_ntuple, "Caltech": treat_Caltech_ntuple, "Iowa": treat_Iowa_ntuple}


def add_inclusive_channel(spec):
    """inclusive: ggH | VBF of the group (the leading-jet VBF definition is the standard one: there are no further derived categories)."""
    spec.channels["inclusive"] = f"(({spec.channels['ggH']}) | ({spec.channels['VBF']}))"


# Update after the Hmm meeting of 15 Oct 2026: TAMU, Marmara and Pisa published new tuples in `<group>/update_after_Hmm_15Oct/`
# (HMUMU_SYNC_VARIANT=update_after_Hmm_15Oct).  They are raw-level: the mass regions and m_mumu use the uncorrected (NanoAOD) dimuon mass,
# which is the only mass Pisa stores.  Pisa lost the corrected mass/pT, the MET and Dijet_eta, so those variables are not compared.
VARIANT = os.environ.get("HMUMU_SYNC_VARIANT", "")
UPDATE = "update_after_Hmm_15Oct"
LUMI_2024H = 5440.0       # pb-1, Run2024H


def apply_update_variant(spec):
    group = spec.name
    if group == "TAMU":
        spec.files = {s: f.replace("TAMU/", f"TAMU/{UPDATE}/") for s, f in spec.files.items()}
        # VBFH: the main production and its _ext1 extension (Pisa's file contains both)
        spec.files["VBFH"] = [spec.files["VBFH"], spec.files["VBFH"].replace(".root", "_ext1.root")]
        spec.mass = "m_mumu_nano"
        # the tuples published on 19 Oct have weight_lumi = 109,948 pb-1 (all 2024); the data are Run2024H only: patch to LUMI_2024H on reading.
        # Tuples produced with the 2024H luminosity already have weight_lumi = LUMI_2024H and the factor is 1.
        spec.weights = {"genWeights": f"weight_generator*({LUMI_2024H}/weight_lumi)",
                        "all_common_weights": f"weight_generator*({LUMI_2024H}/weight_lumi)*weight_pileup*weight_muon_idiso*weight_trigger"}
        spec.notes.append(f"Tuple of `{UPDATE}`; regions and m_mumu on the raw (NanoAOD) dimuon mass `m_mumu_nano`. MC weights are scaled to the Run2024H "
                          f"luminosity ({LUMI_2024H:g} pb-1) on reading (`weight_generator * {LUMI_2024H:g} / weight_lumi`).")
    elif group == "Caltech":
        # Caltech `weight` = gen*PU*SF*xsec*lumi/sum(gen) with the lumi of all 2024 (109.82 fb-1): scaled to Run2024H
        f = "(5440.0/109820.0)"
        spec.weights = {"genWeights": f"weight/(pileup_weight*mu_trig_eff_SF*mu_ID_eff_SF*mu_ISO_eff_SF)*{f}", "all_common_weights": f"weight*{f}"}
        spec.notes.append("MC weights: Caltech's `weight` (normalised to all 2024, 109.82 fb-1) times 5.44/109.82, i.e. Run2024H; `genWeights` = `weight` without the PU and muon SF factors. "
                          "`VBF_skim.root` has 6.83M entries, i.e. it should contain the main production and `ext1`.")
    elif group == "Marmara":
        spec.files = {"Data_Run2024H": f"Marmara/{UPDATE}/tree_DATA_2024_loose_round2.root", "ggH": f"Marmara/{UPDATE}/tree_GGH_2024_loose_round2_upd2_nobtag_weight.root", "VBFH": f"Marmara/{UPDATE}/tree_VBF_2024_loose_round2_upd2_nobtag_weight.root"}
        spec.mass = "mm_mass_raw"
        # `w_total` includes the b-tag SF (it is exactly the product of its factors): the b-tag weight is left out here
        spec.weights = {"genWeights": "w_gen", "all_common_weights": "w_gen*w_norm*w_pu*w_muID_tot"}
        spec.notes.append(f"Tuples of `{UPDATE}` (Data_Run2024H, ggH and VBFH); regions and m_mumu on `mm_mass_raw`.")
    elif group == "NNU":
        d = f"NNU/{UPDATE}/"
        spec.files = {"Data_Run2024H": d + "NNU-2024Hmuon01_update.root", "ggH": d + "NNU-GluGluH-Hto2Mu_M125_RunIII2024Summer24_update.root",
                      "VBFH": d + "NNU-VBFH-Hto2Mu_M125_RunIII2024Summer24_update.root"}
        spec.mass = "dimuon_raw_mass"
        base = ("passLumiJSON & passTrigger & passPV & passMETFilters & passTriggerMatch & passRawPreselection & passMuonSelection & "
                "passThirdLeptonVeto & passBTagVeto")
        # flaginfo.txt (21 Oct): ggH = baseline (no jet selection) with passVBFKinematics == 0; VBF = passAnalysisSelection
        spec.channels = {"VBF": "passAnalysisSelection", "ggH": f"({base}) & (passVBFKinematics==0)"}
        spec.weights = {"genWeights": "genWeight", "all_common_weights": "weightBaseNoNorm"}
        spec.flags = ["passLumiJSON", "passTrigger", "passPV", "passMETFilters", "passTriggerMatch", "passRawPreselection", "passMuonSelection",
                      "passThirdLeptonVeto", "passJetSelection", "passBTagVeto", "passVBFKinematics", "passAnalysisSelection"]
        spec.notes = [f"Tuple of `{UPDATE}` (20 Oct): raw-level, regions and m_mumu on `dimuon_raw_mass`, pT on `dimuon_raw_pt`. VBF = `passAnalysisSelection`; "
                      "ggH = the baseline of `flaginfo.txt` (no jet selection) with `passVBFKinematics == 0`. Data are blinded on the raw mass; NNU asks to also exclude 115-135 on `dimuon_fsr_mass`, "
                      "which is not needed here because the regions use the raw mass."]
    elif group == "MIT":
        d = f"MIT/VBF-GGH/{UPDATE}/snapshot_mc_{{}}_2024_ggHcat_synch.root"
        spec.files = {"Data_Run2024H": [d.format("-51"), d.format("-52")], "ggH": d.format("11"), "VBFH": d.format("10")}
        spec.mass = MIT_RAW_MASS
        spec.channels = {"ggH": "isGGH", "VBF": "~isGGH"}      # MIT (21 Oct): the tuple is ggH + VBF with the flag isGGH; VBF = !isGGH
        spec.flags = ["isGGH"]
        # w and w_allSF are per fb-1 (MIT): the weight for the era is w * lumi [fb-1]; here the Run2024H luminosity
        spec.weights = {"genWeights": f"w*{LUMI_2024H / 1000:g}", "all_common_weights": f"w_allSF*{LUMI_2024H / 1000:g}"}
        spec.notes = ["Tuple of `update_after_Hmm_15Oct` (21 Oct, 25 branches, tree `events`): ggH + VBF selection with the b-veto relaxed and the extra mu-e veto on; "
                      "ggH = `isGGH`, VBF = `!isGGH` (as told by MIT). There are no jet or b-tag branches except `nGoodJetsAll`, so only the muon and mass variables are compared. "
                      "Regions and m_mumu on the raw mass rebuilt from `Muon1/2_rawpt`, eta and phi; `w` and `w_allSF` are per fb-1 (`lumiIntegrated` = 109.82 fb-1 for all 2024): the weights used here are `w` and `w_allSF` times the Run2024H luminosity, 5.44 fb-1. "
                      "The data are not blinded in 115-135 GeV."]
    elif group == "PKU":
        d = f"PKU/{UPDATE}/"
        spec.files = {k: d + v.split("/", 1)[1] for k, v in spec.files.items()}
        spec.mass = "H_mass"
        # PKU `weight` has the same normalisation constant as Caltech's (ggH: 1.3719e-5 vs 1.3653e-5, i.e. all 2024, inferred): scaled to Run2024H
        spec.weights = {"genWeights": "weight*(5440.0/109820.0)", "all_common_weights": "eventWeight*(5440.0/109820.0)"}
        spec.channels = {"ggH": "GGH", "VBF": "is_vbf_leading"}      # PKU (21 Oct): `is_vbf_leading` is the VBF selection to use (`is_vbf_maxmjj` is the max-m_jj one)
        spec.flags = ["trg_single_mu24", "VBF", "GGH", "is_vbf_leading", "is_vbf_maxmjj"]
        spec.notes.append(f"Tuple of `{UPDATE}` (21 Oct): regions and m_mumu on the raw dimuon mass `H_mass`, pT on `H_pt`. Channels: `GGH` and `is_vbf_leading` (the VBF "
                          "selection on the two leading jets, as PKU asked; it coincides event by event with the `VBF` flag; `is_vbf_maxmjj` is the max-m_jj selection). The data are not blinded on `H_mass` (415 events in 115-135 GeV "
                          "in the data file, 1,417 on `diMufsr_kit_BSC_mass`).")
    elif group == "Purdue":
        d = f"Purdue/{UPDATE}/Purdue-"
        spec.files = {"Data_Run2024H": [d + "Muon0_Run2024H-MINIv6NANOv15-v1.root", d + "Muon1_Run2024H-MINIv6NANOv15-v2.root"],
                      "ggH": d + "GluGluH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8.root",
                      "VBFH": d + "VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8.root",
                      "ttH": d + "TTH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8.root",
                      "ZH": d + "ZH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8.root",
                      "WplusH": d + "WplusH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8.root"}
        spec.mass = "dimuon_mass_raw_agreed"
        spec.channels = {"ggH": "(pass_agreed_selection == 1) & (cat_ggh_agreed == 1)", "VBF": "(pass_agreed_selection == 1) & (cat_vbf_agreed == 1)"}
        # wgt_nominal = gen * norm * xsec * lumi * PU * NNLOPS * muID * muIso * muTrig, with lumi = separate_wgt_lumi = 109,820 pb-1 (all 2024): scaled to Run2024H on reading
        spec.weights = {"genWeights": "separate_wgt_genWeight", "all_common_weights": f"wgt_nominal*({LUMI_2024H}/separate_wgt_lumi)"}
        spec.flags = ["pass_agreed_selection", "pass_HLT", "pass_trigger_match_agreed", "pass_MET_filters", "pass_PV", "pass_nmuons_agreed", "pass_OS_agreed",
                      "pass_electron_veto_agreed", "pass_jet_veto_map", "pass_bveto_agreed", "pass_vbf_agreed", "cat_vbf_agreed", "cat_ggh_agreed"]
        spec.notes = [f"Tuple of `{UPDATE}` (20 Oct, one file per sample, tree `Events`, `sync_info` string inside): raw-level, regions and m_mumu on "
                      "`dimuon_mass_raw_agreed`, pT on `dimuon_pt_raw`. Channels: `pass_agreed_selection` and `cat_ggh_agreed` / `cat_vbf_agreed` "
                      "(the category flags alone also fire on events that fail the selection). Data are blinded on the raw mass; the corrected "
                      "`dimuon_mass` still has events in 115-135 GeV.",
                      "Muon pT = `mu*_pt` (corrected), NanoAOD-like = `mu*_pt_nano`. Jets are the `*_nominal` leading jets; the VBF cuts of the agreed "
                      "selection (mjj > 400, |deta| > 2.5, pT 35/30) are applied on the leading and max-mjj pairs."]
    elif group == "Pisa":
        d = f"PISA/{UPDATE}/"
        spec.files = {"Data_Run2024H": [d + "PISA-Muon0_Run2024H-MINIv6NANOv15-v1/Events.root", d + "PISA-Muon1_Run2024H-MINIv6NANOv15-v2/Events.root"],
                      "ggH": d + "PISA-GluGluH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15/Events.root",
                      "VBFH": d + "PISA-VBFH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15/Events.root"}
        spec.mass = "Dimuon_rawMass"
        # the pass_* flags are cumulative; `first_failed_step` is shifted by one in data (extra BlindingSR step), so the flags are used
        spec.channels = {"VBF": "pass_BJetVeto",
                         "ggH": "pass_EleVeto & ~pass_DijetDeltaEta & (nMediumBTagJet == 0) & (nLooseBTagJet < 2)"}
        spec.weights = {"genWeights": "Generator_weight", "all_common_weights": "Generator_weight*PuSF*MuIdSF*MuIsoSF*MuTrgSF"}
        spec.flags = ["pass_BASESELECTION", "pass_NoiseFilters", "pass_GoldenJson", "pass_GoodPV", "pass_isTriggered", "pass_AtLeastOneTrgMuon",
                      "pass_TwoMuons", "pass_TwoOppositeSignMuons", "pass_RawDimuonMassWindow", "pass_LeadingMuonCut", "pass_EleVeto",
                      "pass_TwoJets", "pass_LeadingJetPt", "pass_DijetMass", "pass_DijetDeltaEta", "pass_BJetVeto"]
        spec.notes = [f"Tuple of `{UPDATE}`: raw-level, no corrected mass/pT, no MET, no `Dijet_eta`; mass = `Dimuon_rawMass` (50-200 GeV window). "
                      "VBF = `pass_BJetVeto` (cumulative flags), ggH = passes up to the electron veto, fails the VBF jet requirements, b-veto "
                      "recomputed from the jet counts. Jets: `*_jesPt`."]
    return spec


PISA_UPDATE_VARIABLES = {
    "m_mumu": "Dimuon_rawMass", "pt_mumu": "Dimuon_rawPt", "mu1_pt": None, "mu2_pt": None, "met_pt": None, "met_phi": None,
    "jet1_pt": "LeadingJet_jesPt", "jet2_pt": "SubLeadingJet_jesPt", "eta_jj": None}
MIT_RAW_MASS = "sqrt(2*Muon1_rawpt*Muon2_rawpt*(cosh(Muon1_eta-Muon2_eta)-cos(Muon1_phi-Muon2_phi)))"
RAW_MASS = {"MIT": MIT_RAW_MASS, "PKU": "H_mass", "TAMU": "m_mumu_nano", "Marmara": "mm_mass_raw", "NNU": "dimuon_raw_mass", "Purdue": "dimuon_mass_raw_agreed"}
# Pisa only publishes the raw dimuon pT in this variant, so `pt_mumu` is the raw one for every group
RAW_PT = {"PKU": "H_pt", "TAMU": "pt_mumu_nano", "Marmara": "mm_pt_uncorr", "Pisa": "Dimuon_rawPt", "NNU": "dimuon_raw_pt", "Purdue": "dimuon_pt_raw"}


def treat(group):
    if VARIANT == UPDATE:
        # Pisa defines the jet1_* variables only with at least 2 jets, TAMU with at least 1: compare them on events with >= 2 jets
        for variable in VARIABLES:
            if variable.name.startswith("jet1_"):
                variable.min_jets = 2
    spec = TREATERS[group]()
    spec.variables = {v.name: v.expr[group] for v in VARIABLES if v.expr.get(group)}
    if VARIANT == UPDATE and group == "Caltech":
        apply_update_variant(spec)       # only the Run2024H weights
    if VARIANT == UPDATE and group in ("TAMU", "Marmara", "Pisa", "NNU", "Purdue", "PKU", "MIT"):
        apply_update_variant(spec)
        if group == "Pisa":
            spec.variables = {k: v for k, v in spec.variables.items() if PISA_UPDATE_VARIABLES.get(k, v) is not None}
            spec.variables.update({k: v for k, v in PISA_UPDATE_VARIABLES.items() if v is not None})
        else:
            spec.variables["m_mumu"] = RAW_MASS[group]
        if group in RAW_PT:
            spec.variables["pt_mumu"] = RAW_PT[group]
        else:
            spec.variables.pop("pt_mumu", None)
    spec.channels = dict(spec.channels)
    add_inclusive_channel(spec)
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
    array, mask, key = (ak.from_regular(x, axis=1) for x in (array, mask, key))   # fixed-size lists would be flattened by the mask
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
    print(f"[files] {group} {sample}: tree `{spec.tree}`, lumi branch `{spec.lumi_branch}`, mass `{spec.mass}`", flush=True)
    for p in paths:
        st = p.stat()
        print(f"[files]   {p}  ({st.st_size / 1e6:,.0f} MB, modified {datetime.datetime.fromtimestamp(st.st_mtime):%Y-%m-%d %H:%M})", flush=True)
    for c, e in spec.channels.items():
        print(f"[files]   channel {c}: {e}", flush=True)
    for m, e in spec.weights.items():
        print(f"[files]   weight {m}: {e}", flush=True)
    signature = {"path": str(paths[0]) if len(paths) == 1 else [str(p) for p in paths],
                 "size": paths[0].stat().st_size if len(paths) == 1 else [p.stat().st_size for p in paths],
                 "mtime": paths[0].stat().st_mtime if len(paths) == 1 else [p.stat().st_mtime for p in paths],
                 "variables": spec.variables, "channels": spec.channels, "weights": spec.weights,
                 "flags": spec.flags, "version": 2}
    if cache.exists():
        stored = np.load(cache, allow_pickle=True)
        old = json.loads(str(stored["signature"]))
        old["channels"] = {k: v for k, v in old.get("channels", {}).items() if not k.endswith("onlyLeading")}   # caches written before these channels were removed
        if old == json.loads(json.dumps(signature)):
            print(f"[cache] {group} {sample}: using {cache}", flush=True)
            return Extract(group, sample, int(stored["entries"]), stored["full_keys"], stored["keys"], stored["mass"],
                           {k: stored[f"ch_{k}"] for k in CHANNELS}, {k[2:]: stored[k] for k in stored.files if k.startswith("v_")},
                           {k[2:]: stored[k] for k in stored.files if k.startswith("f_")},
                           {k[2:]: stored[k] for k in stored.files if k.startswith("w_")})
    print(f"[cache] {group} {sample}: not found or outdated ({cache}), reading the ntuple", flush=True)
    trees = [uproot.open(p)[spec.tree] for p in paths]
    present = set(trees[0].keys())
    is_data = sample == "Data_Run2024H"
    modes = () if is_data else MC_MODES
    usable_vars = {k: e for k, e in spec.variables.items() if set(names(e)) <= present}
    usable_flags = [f for f in spec.flags if f in present]
    full, parts = [], []
    entries = sum(t.num_entries for t in trees)
    print(f"[read] {group} {sample}: {entries:,} entries", flush=True)
    missing = [k for k in spec.variables if k not in usable_vars]
    print(f"[read] {group} {sample}: {len(usable_vars)} variables, {len(usable_flags)}/{len(spec.flags)} flags found"
          + (f"; not in the tuple: {', '.join(missing)}" if missing else ""), flush=True)
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
    print(f"[read] {group} {sample}: {len(out.full_keys):,} unique (run, lumi, event) of {entries:,} rows ({entries - len(out.full_keys):,} duplicates counted once); "
          + ", ".join(f"{c}: {int(np.count_nonzero(v)):,}" for c, v in out.channel.items()), flush=True)
    return out


# ----------------------------------------------------------------------------------------------
# Comparison
# ----------------------------------------------------------------------------------------------

def selection(e, region, channel):
    return REGIONS[region][0](e.mass) & e.channel[channel]


def destination_labels(e):
    """For every kept event the first 'channel / region' it belongs to (index into the label list)."""
    names_ = [f"{channel} / {region}" for channel in BASE_CHANNELS for region in REGIONS] + ["outside 70-150 GeV or ggH/VBF"]
    code = np.full(len(e.keys), len(names_) - 1)
    for i, (channel, region) in reversed(list(enumerate((c, r) for c in BASE_CHANNELS for r in REGIONS))):
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


def draw(ha, hb, bins, title, labels, heading, out, na=None, nb=None):
    width = np.diff(bins)
    # top panel: events per bin (sum of weights for MC); the ratio panel compares the unit-normalised shapes
    ya = ha / ha.sum() / width if ha.sum() else ha
    yb = hb / hb.sum() / width if hb.sum() else hb
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": .05})
    # each group's (weighted) shape is scaled to its number of selected events, so that groups with incompatible weight
    # normalisations (raw generator weight, cross-section x luminosity, ...) are drawn on a comparable scale
    ca = ha / ha.sum() * na if (na and ha.sum()) else ha
    cb = hb / hb.sum() * nb if (nb and hb.sum()) else hb
    ax.stairs(ca, bins, label=labels[0], linewidth=1.7)
    ax.stairs(cb, bins, label=labels[1], linewidth=1.7)
    ax.set_ylabel("events / bin" if heading.endswith("unweighted") else "weighted events / bin")
    ax.set_title(heading, fontsize=10)
    ax.legend()
    # a bin empty in either histogram has no meaningful ratio: leave it undefined (NaN) instead of a point at zero
    ratio = np.divide(ya, yb, out=np.full_like(ya, np.nan), where=(ya > 0) & (yb > 0))
    rax.axhline(1, color="black", linewidth=.8)
    centres = 0.5 * (bins[1:] + bins[:-1])
    rax.errorbar(centres, ratio, xerr=width / 2, fmt="o", color="C0", markersize=3, linestyle="none", elinewidth=1)
    finite = ratio[np.isfinite(ratio)]
    spread = float(np.max(np.abs(finite - 1))) if finite.size else 0.
    ratio_axis(rax, spread)
    rax.set_ylabel(f"{labels[0]} /\n{labels[1]}\n(shapes)")
    rax.set_xlabel(title)
    fig.savefig(out.with_suffix(".png"), dpi=130, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def ratio_axis(rax, spread):
    """Ratio panel: 5% (green) and 10% (yellow) bands, range 0.8-1.2 extended only if the points exceed it."""
    rax.axhspan(0.90, 1.10, color="gold", alpha=0.25, linewidth=0, zorder=0)
    rax.axhspan(0.95, 1.05, color="limegreen", alpha=0.35, linewidth=0, zorder=0)
    half = min(max(0.2, 1.1 * spread), 1.0)
    rax.set_ylim(1 - half, 1 + half)


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


def draw_flags(flags, ga, gb, heading, out):
    """Fraction of the selected events passing each boolean flag, one panel per group (the flags are the groups' own)."""
    panels = [(g, flags[g]) for g in (ga, gb) if flags.get(g)]
    if not panels or all(v["total"] == 0 for _, d in panels for v in d.values()):
        return
    height = max(len(d) for _, d in panels) * 0.32 + 1.6
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), height), squeeze=False)
    for ax, (g, d) in zip(axes[0], panels):
        names_ = list(d)
        frac = [v["passed"] / v["total"] if v["total"] else math.nan for v in d.values()]
        ax.barh(range(len(names_)), frac, color="C0" if g == ga else "C1")
        ax.set_yticks(range(len(names_)))
        ax.set_yticklabels(names_, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("fraction of selected events passing")
        ax.set_title(f"{g}", fontsize=10)
        for i, f in enumerate(frac):
            if math.isfinite(f):
                ax.text(min(f + 0.01, 0.98), i, f"{f:.4f}", va="center", fontsize=7)
    fig.suptitle(heading, fontsize=10)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=120, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


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
            draw(ha, hb, v.bins, v.title, (f"{a.group} ({na:,})", f"{b.group} ({nb:,})"), f"{sample} - {region} - {channel} - {mode}", folder / v.name, na, nb)
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


def weights_and_branches_markdown(r, region, channel, sample, ga, gb):
    """Which weight each group's expression uses in each mode, and which native branch/expression is behind every plotted variable."""
    cell_ = lambda text: f"`{text}`".replace("|", "\\|") if text else "-"
    out = ["## Weights used", ""]
    if sample == "Data_Run2024H":
        out += ["Data: every event has weight 1 for every group (the MC weights below are not used).", ""]
    else:
        out += [f"| mode | {ga} | {gb} |", "|---|---|---|"]
        out += [f"| {m} | {cell_(r['treatment'][ga]['weights'].get(m))} | {cell_(r['treatment'][gb]['weights'].get(m))} |" for m in MC_MODES]
        out += ["", "`genWeights` is the generator weight only. `all_common_weights` is each group's own product of generator weight and scale factors, as written in the table: "
                "groups that publish the factors separately are used without cross section, luminosity or sum-of-weights normalisation; groups that publish only a combined weight use it as published. "
                "The normalisation therefore differs between groups, and the ratio of the sums of weights is not a yield comparison.", ""]
    modes = r["distributions"].get(region, {}).get(channel)
    names = list(next(iter(modes.values()))["variables"]) if modes else []
    if not names:
        return out
    va, vb = r["treatment"][ga]["variables"], r["treatment"][gb]["variables"]
    out += ["## Branches plotted", "", f"| variable | {ga} branch / expression | {gb} branch / expression |", "|---|---|---|"]
    out += [f"| `{k}` | {cell_(va.get(k))} | {cell_(vb.get(k))} |" for k in names]
    out += [""]
    return out


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
    out += ["", *weights_and_branches_markdown(r, region, channel, sample, ga, gb)]
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
                    "The groups normalise their weights differently, so the ratio of the sums is not a yield comparison. The plots show events per bin (for MC arbitrary units, A.U.: the weighted shape of each group scaled to its number of selected events, because the weight normalisations differ); the ratio panel and the total-variation distance use the unit-normalised shapes; the distance is 0 for identical shapes and 1 for disjoint ones.", "",
                    f"| variable | {ga} mean | {gb} mean | {ga} std | {gb} std | shape distance | plot |", "|---|---:|---:|---:|---:|---:|---|"]
            for k, v in sorted(entry["variables"].items(), key=lambda kv: -(kv[1]["shape_distance"] if math.isfinite(kv[1]["shape_distance"]) else -1)):
                out.append(f"| `{k}` | {fmt(v['mean'][ga])} | {fmt(v['mean'][gb])} | {fmt(v['std'][ga])} | {fmt(v['std'][gb])} | "
                           f"{fmt(v['shape_distance'])} | [png]({mode}/{k}.png) / [pdf]({mode}/{k}.pdf) |")
            out += [""] + [f"![{k}]({mode}/{k}.png)" for k in list(entry["variables"])[:6]] + [""]
    out += ["## Selection flags among the selected events", "", "| group | flag | passing / selected |", "|---|---|---:|"]
    for g in (ga, gb):
        out += [f"| {g} | `{f}` | {v['passed']:,} / {v['total']:,} |" for f, v in cell["flags"][g].items()]
    if any(cell["flags"][g] for g in (ga, gb)) and cell["selected"][ga] + cell["selected"][gb]:
        out += ["", "![flags](flags.png)"]
    out += ["", "## Limits of this diagnosis", "",
            "- Flags with similar names are not assumed to have identical definitions.",
            "- Every group's corrected quantities are used as published; the top panel shows events per bin (MC: A.U., weighted shape scaled to the selected events), but the ratio panel and the shape distances are unit-normalised, so overall normalisation differences do not enter them.",
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
            draw_flags(cell["flags"], ga, gb, f"{sample} - {region} - {channel} - flags", leaf / "flags")
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
        h_g = hists[g][0]
        ax.stairs(h_g / h_g.sum() * n if h_g.sum() else h_g, bins, label=f"{g} ({n:,})", linewidth=1.7, color=f"C{i}")
        if g == reference:
            continue
        yr = densities[reference][0]
        ratio = np.divide(y, yr, out=np.full_like(y, np.nan), where=(y > 0) & (yr > 0))
        finite = ratio[np.isfinite(ratio)]
        if finite.size:
            spread = max(spread, float(np.max(np.abs(finite - 1))))
        rax.errorbar(centres + 0.0, ratio, xerr=width / 2, fmt="o", color=f"C{i}", markersize=3, linestyle="none", elinewidth=1)
    ratio_axis(rax, spread)
    ax.set_ylabel("events / bin" if heading.endswith("unweighted") else "weighted events / bin")
    ax.set_title(heading, fontsize=10)
    ax.legend()
    rax.set_ylabel(f"group / {reference}\n(shapes)")
    rax.set_xlabel(title)
    fig.savefig(out.with_suffix(".png"), dpi=130, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def compare_all_groups(sample, output, cache_dir, regions, channels, groups, folder_name="ALL_GROUPS"):
    """Every group that provides the sample on one plot per variable, region, channel and weight mode."""
    present = [g for g in groups if sample in treat(g).files]
    if len(present) < 2:
        return None
    extracts = {g: extract(g, sample, cache_dir) for g in present}
    root = output / folder_name / sample
    root.mkdir(parents=True, exist_ok=True)
    info = {"sample": sample, "groups": present, "cells": {}}
    out = [f"# {'All groups together' if folder_name == 'ALL_GROUPS' else folder_name.replace('_VS_OTHERS', ' against all the other groups')}: {sample}", "",
           f"Groups: {', '.join(present)}. Every plot overlays the events per bin of each group (MC: A.U., weighted shape scaled to the selected events); the ratio panel compares the unit-normalised shapes and is to **{present[0]}**. "
           "Every group uses its own definition of the channel and of the corrected mass (see the pair reports for the details). "
           "The plots in `<mode>/` have the ratio panel to the first group; `<mode>/ratio_to_<GROUP>/` has the same plots with the ratio to each other group.", ""]
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
                    # the same plot with the ratio panel to each of the other groups
                    for other in hists:
                        if other == reference:
                            continue
                        sub = folder / f"ratio_to_{other}"
                        sub.mkdir(exist_ok=True)
                        ordered = {other: hists[other], **{g: h for g, h in hists.items() if g != other}}
                        draw_all(ordered, v.bins, v.title, f"{sample} - {region} - {channel} - {mode}", sub / v.name)
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
    p.add_argument("--cache-dir", type=Path, default=Path("run3_sync_extracts"))
    p.add_argument("--group1", choices=TREATERS)
    p.add_argument("--group2", choices=TREATERS)
    p.add_argument("--sample", action="append", help="restrict to this sample (repeatable)")
    p.add_argument("--region", action="append", choices=REGIONS)
    p.add_argument("--channel", action="append", choices=CHANNELS)
    p.add_argument("--resume", action="store_true", help="keep pair/sample directories that already have a summary.json")
    p.add_argument("--summaries-only", action="store_true", help="rewrite the pair summaries and the root index from the existing summary.json files")
    p.add_argument("--reference", choices=TREATERS, help="one group against all the others: only the overlay plots, with the ratio panels "
                   "to this group, written in <OUTPUT>/<REFERENCE>_VS_OTHERS/ (no pairwise comparison)")
    p.add_argument("--no-overlay", action="store_true", help="skip the all-groups overlay plots")
    p.add_argument("--list", action="store_true", help="list the groups, samples and pairs and exit")
    a = p.parse_args()
    if bool(a.group1) != bool(a.group2):
        p.error("--group1 and --group2 must be given together")
    if a.group1 and a.group1 == a.group2:
        p.error("--group1 and --group2 must differ")
    groups = list(TREATERS)
    if a.reference:
        if not a.output:
            p.error("--output is required")
        a.output.mkdir(parents=True, exist_ok=True)
        ordered = [a.reference] + [g for g in groups if g != a.reference]
        for s in sorted({s for g in groups for s in treat(g).files if not a.sample or s in a.sample}):
            print(f"[reference {a.reference}] {s}", flush=True)
            compare_all_groups(s, a.output, a.cache_dir, a.region or list(REGIONS), a.channel or list(CHANNELS), ordered,
                               f"{a.reference}_VS_OTHERS")
        return
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
