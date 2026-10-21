#!/usr/bin/env python3
"""Standalone tool to compare the Run 3 H->mumu synchronisation ntuples of several groups.

It needs only numpy, awkward, uproot and matplotlib (see requirements.txt); nothing else from the analysis framework.

Two steps (full instructions in HOWTO.txt):

  1. yields : read the ntuples (events with the same run, lumi, event are counted once) of every group and write tables
       yields_<sample>.csv   events and sum of weights per group, mass region and channel
       overlap_<sample>.csv  events present in both groups of every pair (matched on run, lumi, event)
       mass_<sample>.csv     m(mumu) distribution of every group (selected ggH or VBF events, 70-150 GeV)
  2. plots  : draw bar charts of the yields, the ratio to every group (or to --reference), and the overlap matrices (common events / events of the row group) from those tables.

  python3 sync_compare.py yields --round after_15Oct --base /eos/cms/store/group/phys_higgs/cmshmm/run3_sync --outdir out
  python3 sync_compare.py plots  --round after_15Oct --outdir out --reference TAMU
  python3 sync_compare.py markdown --round after_15Oct --outdir out      # tables.md with the yield and overlap tables

Rounds (ROUNDS below): before_15Oct = tuples published 12-15 Oct 2026 directly in the group folders (TAMU, Pisa, Marmara, PKU,
Purdue, MIT, NNU, Caltech); after_15Oct = tuples of <group>/update_after_Hmm_15Oct (TAMU, Pisa, Marmara, Caltech, NNU, Iowa, Purdue, PKU).
Add --samples ttH ZH WplusH for the other signal samples (TAMU, PKU, Purdue).

A group is one entry of a round in ROUNDS below: the files of every sample, the tree name, the dimuon-mass branch used for the
regions, one boolean expression per channel (ggH, VBF) and the weight expressions.  A new group is added by adding one entry.
Expressions are Python expressions on branch names (with abs, sqrt, hypot, cos, sin, cosh, sinh); `name[:,0]` is the first entry of a jagged branch.
"""
import argparse
import csv
import glob
import re
import sys
from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

# ----------------------------------------------------------------------------------------------
# Description of the groups (update_after_Hmm_15Oct tuples; paths are relative to --base)
# ----------------------------------------------------------------------------------------------
ROUNDS = {
    # tuples published 12-15 Oct 2026, directly in the group folder (corrected dimuon mass; Caltech: first small skim of 15 Oct)
    'before_15Oct': {'TAMU': {'tree': 'Events',
          'lumi': 'luminosityBlock',
          'mass': 'm_mumu_FSR_corr',
          'files': {'Data_Run2024H': ['TAMU/TAMU-Run2024H.root'],
                    'ggH': ['TAMU/TAMU-GluGluHto2Mu.root'],
                    'VBFH': ['TAMU/TAMU-VBFHto2Mu_M125_powheg.root'],
                    'ttH': ['TAMU/TAMU-TTH_Hto2Mu.root'],
                    'ZH': ['TAMU/TAMU-ZH_Hto2Mu.root'],
                    'WplusH': ['TAMU/TAMU-WplusH_Hto2Mu.root']},
          'channels': {'ggH': 'ggF', 'VBF': 'VBF'},
          'weights': {'gen': 'weight_generator', 'all': 'weight_generator*weight_pileup*weight_muon_idiso*weight_trigger'}},
 'Pisa': {'tree': 'Events',
          'lumi': 'luminosityBlock',
          'mass': 'Dimuon_mass',
          'files': {'Data_Run2024H': ['PISA/PISA-Muon0_Run2024H-MINIv6NANOv15-v1/Events.root',
                                      'PISA/PISA-Muon1_Run2024H-MINIv6NANOv15-v2/Events.root'],
                    'ggH': ['PISA/PISA-GluGluH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15-v2/Events.root'],
                    'VBFH': ['PISA/PISA-VBFH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15-v2/Events.root']},
          'channels': {'ggH': '(first_failed_step >= 13) & (first_failed_step <= 16) & (nMediumBTagJet == 0) & (nLooseBTagJet < 2)',
                       'VBF': '(first_failed_step == 0) | (first_failed_step >= 18)'},
          'weights': {'gen': 'Generator_weight', 'all': 'Generator_weight*PuSF*MuIdIsoSF*MuTrgSF'}},
 'Marmara': {'tree': 'HmmTree',
             'lumi': 'lumi',
             'mass': 'mm_mass_fsr',
             'files': {'Data_Run2024H': ['Marmara/tree_DATA_2024_loose_for_ALL_2024_H_Data_Samples.root'],
                       'ggH': ['Marmara/tree_GGH_2024_loose_LAST.root'],
                       'VBFH': ['Marmara/tree_VBF_2024_loose_LAST.root']},
             'channels': {'ggH': 'sel_full & cat_isGGH[:,0]', 'VBF': 'sel_full & cat_isVBF[:,0]'},
             'weights': {'gen': 'w_gen', 'all': 'w_total'}},
 'PKU': {'tree': 'ntuple',
         'lumi': 'lumi',
         'mass': 'diMufsr_kit_BSC_mass',
         'files': {'Data_Run2024H': ['PKU/Muon_24H.root'],
                   'ggH': ['PKU/GluGluH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root'],
                   'VBFH': ['PKU/VBFH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root'],
                   'ttH': ['PKU/TTH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root'],
                   'ZH': ['PKU/ZH-Hto2Mu_Par-M-125_powhegMINLO_RunIII2024Summer24.root'],
                   'WplusH': ['PKU/WplusH-Hto2Mu_Par-M-125_powhegMINLO_RunIII2024Summer24.root']},
         'channels': {'ggH': 'GGH', 'VBF': 'VBF'},
         'weights': {'gen': 'weight', 'all': 'eventWeight'}},
 'Purdue': {'tree': 'Events',
            'lumi': 'luminosityBlock',
            'mass': 'dimuon_mass',
            'files': {'Data_Run2024H': ['Purdue/Purdue-Muon0_Run2024H-MINIv6NANOv15-v1/*.root',
                                        'Purdue/Purdue-Muon1_Run2024H-MINIv6NANOv15-v2/*.root'],
                      'ggH': ['Purdue/Purdue-GluGluH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/*.root'],
                      'VBFH': ['Purdue/Purdue-VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/*.root'],
                      'ttH': ['Purdue/Purdue-TTH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/*.root'],
                      'ZH': ['Purdue/Purdue-ZH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8/*.root'],
                      'WplusH': ['Purdue/Purdue-WplusH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8/*.root']},
            'channels': {'ggH': '(pass_analysis_preselection == 1) & (cat_ggh == 1)',
                         'VBF': '(pass_analysis_preselection == 1) & (cat_vbf == 1)'},
            'weights': {'gen': 'separate_wgt_genWeight', 'all': 'wgt_nominal'}},
 'MIT': {'tree': 'events',
         'lumi': 'luminosityBlock',
         'mass': 'HiggsCandCorrMass',
         'files': {'Data_Run2024H': ['MIT/VBF-GGH/snapshot_mc_-51_2024_ggHcat_synch.root',
                                     'MIT/VBF-GGH/snapshot_mc_-52_2024_ggHcat_synch.root'],
                   'ggH': ['MIT/VBF-GGH/snapshot_mc_11_2024_ggHcat_synch.root'],
                   'VBFH': ['MIT/VBF-GGH/snapshot_mc_10_2024_ggHcat_synch.root']},
         'channels': {'ggH': '(nBMJets == 0) & ~((Mjj > 400) & (abs(jetVBF1_Eta - jetVBF2_Eta) > 2.5))',
                      'VBF': '(nBMJets == 0) & (Mjj > 400) & (abs(jetVBF1_Eta - jetVBF2_Eta) > 2.5)'},
         'weights': {'gen': 'w', 'all': 'w_allSF'}},
 'NNU': {'tree': 'Events',
         'lumi': 'luminosityBlock',
         'mass': 'dimuon_fsr_mass',
         'files': {'Data_Run2024H': ['NNU/NNU-dataH_2024_muon01.root'], 'VBFH': ['NNU/NNU-VBFH-Hto2Mu_M125_RunIII2024Summer24.root']},
         'channels': {'VBF': 'passAnalysisSelection',
                      'ggH': '(passLumiJSON & passTrigger & passTriggerMatch & passOS & passCandidateGoodMuons & passElectronVeto & '
                             'passBTagVeto & passMETFilters & passJetVetoMap & passLeadPt & passPV & passExactly2GoodMuons) & (category != '
                             '3)'},
         'weights': {'gen': 'genWeight', 'all': 'weightBaseNoNorm'}},
 'Caltech': {'tree': 'tree_output',
             'lumi': 'luminosityBlock',
             'mass': 'diMuon_bscMass_noCorr',
             'files': {'Data_Run2024H': ['Caltech/Muon0_2024H_skim.root', 'Caltech/Muon1_2024H_skim.root'],
                       'ggH': ['Caltech/ggH_skim.root'],
                       'VBFH': ['Caltech/VBF_skim.root']},
             'channels': {'ggH': '(is_ggH_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto',
                          'VBF': '(is_VBF_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto'},
             'weights': {'gen': 'gen_weight', 'all': 'gen_weight*pileup_weight*mu_trig_eff_SF*mu_ID_eff_SF*mu_ISO_eff_SF'}}},
    # tuples of <group>/update_after_Hmm_15Oct (raw dimuon mass)
    'after_15Oct': {'TAMU': {'tree': 'Events',
          'lumi': 'luminosityBlock',
          'mass': 'm_mumu_nano',
          'files': {'Data_Run2024H': ['TAMU/update_after_Hmm_15Oct/TAMU-Run2024H.root'],
                    'ggH': ['TAMU/update_after_Hmm_15Oct/TAMU-GluGluHto2Mu.root'],
                    'VBFH': ['TAMU/update_after_Hmm_15Oct/TAMU-VBFHto2Mu_M125_powheg.root',
                             'TAMU/update_after_Hmm_15Oct/TAMU-VBFHto2Mu_M125_powheg_ext1.root'],
                    'ttH': ['TAMU/update_after_Hmm_15Oct/TAMU-TTH_Hto2Mu.root'],
                    'ZH': ['TAMU/update_after_Hmm_15Oct/TAMU-ZH_Hto2Mu.root'],
                    'WplusH': ['TAMU/update_after_Hmm_15Oct/TAMU-WplusH_Hto2Mu.root']},
          'channels': {'ggH': 'ggF', 'VBF': 'VBF'},
          'weights': {'gen': 'weight_generator*(5440.0/weight_lumi)', 'all': 'weight_generator*(5440.0/weight_lumi)*weight_pileup*weight_muon_idiso*weight_trigger'}},   # published 19 Oct files: weight_lumi = 109,948 pb-1; scaled to Run2024H
 'Pisa': {'tree': 'Events',
          'lumi': 'luminosityBlock',
          'mass': 'Dimuon_rawMass',
          'files': {'Data_Run2024H': ['PISA/update_after_Hmm_15Oct/PISA-Muon0_Run2024H-MINIv6NANOv15-v1/Events.root',
                                      'PISA/update_after_Hmm_15Oct/PISA-Muon1_Run2024H-MINIv6NANOv15-v2/Events.root'],
                    'ggH': ['PISA/update_after_Hmm_15Oct/PISA-GluGluH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15/Events.root'],
                    'VBFH': ['PISA/update_after_Hmm_15Oct/PISA-VBFH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15/Events.root']},
          'channels': {'VBF': 'pass_BJetVeto', 'ggH': 'pass_EleVeto & ~pass_DijetDeltaEta & (nMediumBTagJet == 0) & (nLooseBTagJet < 2)'},
          'weights': {'gen': 'Generator_weight', 'all': 'Generator_weight*PuSF*MuIdSF*MuIsoSF*MuTrgSF'}},
 'Marmara': {'tree': 'HmmTree',
             'lumi': 'lumi',
             'mass': 'mm_mass_raw',
             'files': {'Data_Run2024H': ['Marmara/update_after_Hmm_15Oct/tree_DATA_2024_loose_round2.root'],
                       'ggH': ['Marmara/update_after_Hmm_15Oct/tree_GGH_2024_loose_round2_upd2_nobtag_weight.root'],
                       'VBFH': ['Marmara/update_after_Hmm_15Oct/tree_VBF_2024_loose_round2_upd2_nobtag_weight.root']},
             'channels': {'ggH': 'sel_full & cat_isGGH[:,0]', 'VBF': 'sel_full & cat_isVBF[:,0]'},
             'weights': {'gen': 'w_gen', 'all': 'w_gen*w_norm*w_pu*w_muID_tot'}},
 'Caltech': {'tree': 'tree_output',
             'lumi': 'luminosityBlock',
             'mass': 'diMuon_bscMass_noCorr',
             'files': {'Data_Run2024H': ['Caltech/update_after_Hmm_15Oct/Muon0_2024H_skim.root',
                                         'Caltech/update_after_Hmm_15Oct/Muon1_2024H_skim.root'],
                       'ggH': ['Caltech/update_after_Hmm_15Oct/ggH_skim.root'],
                       'VBFH': ['Caltech/update_after_Hmm_15Oct/VBF_skim.root']},
             'channels': {'ggH': '(is_ggH_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto',
                          'VBF': '(is_VBF_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto'},
             'weights': {'gen': 'weight/(pileup_weight*mu_trig_eff_SF*mu_ID_eff_SF*mu_ISO_eff_SF)*(5.44/109.82)', 'all': 'weight*(5.44/109.82)'}},   # Caltech `weight` is normalised to all 2024 (109.82 fb-1): scaled to Run2024H
 'NNU': {'tree': 'Events',
         'lumi': 'luminosityBlock',
         'mass': 'dimuon_raw_mass',
         'files': {'Data_Run2024H': ['NNU/update_after_Hmm_15Oct/NNU-2024Hmuon01_update.root'],
                   'ggH': ['NNU/update_after_Hmm_15Oct/NNU-GluGluH-Hto2Mu_M125_RunIII2024Summer24_update.root'],
                   'VBFH': ['NNU/update_after_Hmm_15Oct/NNU-VBFH-Hto2Mu_M125_RunIII2024Summer24_update.root']},
         'channels': {'VBF': 'passAnalysisSelection',
                      'ggH': '(passLumiJSON & passTrigger & passPV & passMETFilters & passTriggerMatch & passRawPreselection & '
                             'passMuonSelection & passThirdLeptonVeto & passBTagVeto) & (passAnalysisSelection==0)'},
         'weights': {'gen': 'genWeight', 'all': 'weightBaseNoNorm'}},
 'MIT': {'tree': 'events', 'lumi': 'luminosityBlock',
            # no raw-mass branch: built from the raw muon pT; the tuple has no jet or b-tag branches (ggH+VBF selection with relaxed b-veto)
            'mass': 'sqrt(2*Muon1_rawpt*Muon2_rawpt*(cosh(Muon1_eta-Muon2_eta)-cos(Muon1_phi-Muon2_phi)))',
            'files': {'Data_Run2024H': ['MIT/VBF-GGH/update_after_Hmm_15Oct/snapshot_mc_-51_2024_ggHcat_synch.root',
                                        'MIT/VBF-GGH/update_after_Hmm_15Oct/snapshot_mc_-52_2024_ggHcat_synch.root'],
                      'ggH': ['MIT/VBF-GGH/update_after_Hmm_15Oct/snapshot_mc_11_2024_ggHcat_synch.root'],
                      'VBFH': ['MIT/VBF-GGH/update_after_Hmm_15Oct/snapshot_mc_10_2024_ggHcat_synch.root']},
            'channels': {'ggH': 'isGGH', 'VBF': '~isGGH'},     # MIT 21 Oct: ggH = isGGH, VBF = !isGGH
            'weights': {'gen': 'w*5.44', 'all': 'w_allSF*5.44'}},   # MIT weights are per fb-1: times the Run2024H luminosity (5.44 fb-1)
 'PKU': {'tree': 'ntuple',
         'lumi': 'lumi',
         'mass': 'H_mass',
         'files': {'Data_Run2024H': ['PKU/update_after_Hmm_15Oct/Muon_24H.root'],
                   'ggH': ['PKU/update_after_Hmm_15Oct/GluGluH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root'],
                   'VBFH': ['PKU/update_after_Hmm_15Oct/VBFH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root'],
                   'ttH': ['PKU/update_after_Hmm_15Oct/TTH-Hto2Mu_Par-M-125_powheg_RunIII2024Summer24.root'],
                   'ZH': ['PKU/update_after_Hmm_15Oct/ZH-Hto2Mu_Par-M-125_powhegMINLO_RunIII2024Summer24.root'],
                   'WplusH': ['PKU/update_after_Hmm_15Oct/WplusH-Hto2Mu_Par-M-125_powhegMINLO_RunIII2024Summer24.root']},
         'channels': {'ggH': 'GGH', 'VBF': 'is_vbf_leading'},
         'weights': {'gen': 'weight*(5.44/109.82)', 'all': 'eventWeight*(5.44/109.82)'}},   # PKU `weight` has the same normalisation as Caltech (all 2024, inferred): scaled to Run2024H
 'Purdue': {'tree': 'Events',
            'lumi': 'luminosityBlock',
            'mass': 'dimuon_mass_raw_agreed',
            'files': {'Data_Run2024H': ['Purdue/update_after_Hmm_15Oct/Purdue-Muon0_Run2024H-MINIv6NANOv15-v1.root',
                                        'Purdue/update_after_Hmm_15Oct/Purdue-Muon1_Run2024H-MINIv6NANOv15-v2.root'],
                      'ggH': ['Purdue/update_after_Hmm_15Oct/Purdue-GluGluH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8.root'],
                      'VBFH': ['Purdue/update_after_Hmm_15Oct/Purdue-VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8.root'],
                      'ttH': ['Purdue/update_after_Hmm_15Oct/Purdue-TTH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8.root'],
                      'ZH': ['Purdue/update_after_Hmm_15Oct/Purdue-ZH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8.root'],
                      'WplusH': ['Purdue/update_after_Hmm_15Oct/Purdue-WplusH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8.root']},
            'channels': {'ggH': '(pass_agreed_selection == 1) & (cat_ggh_agreed == 1)', 'VBF': '(pass_agreed_selection == 1) & (cat_vbf_agreed == 1)'},
            'weights': {'gen': 'separate_wgt_genWeight', 'all': 'wgt_nominal*(5440.0/separate_wgt_lumi)'}},   # wgt_nominal has lumi = 109,820 pb-1 (all 2024): scaled to Run2024H
 'Iowa': {'tree': 'SyncEvents',
          'lumi': 'luminosityBlock',
          'mass': 'massRaw',
          'files': {'Data_Run2024H': ['Iowa/update_fullScaRe_21Oct/Muon0/*.root', 'Iowa/update_fullScaRe_21Oct/Muon1/*.root'],
                    'VBFH': ['Iowa/update_fullScaRe_21Oct/VBF/*.root']},
          'channels': {'ggH': 'passGolden & passFilters & passHLTAndObjectMatch & passAnalysisMuonSelection & passExactlyTwoAnalysisMuons '
                              '& passElectronVeto & passAdditionalLeptonVeto & passJetVetoMap & passBJetVeto & passVertex & '
                              'analysisPairPassesRawWindow & oppositeSign & ~passVBFCategoryApprox',
                       'VBF': 'passGolden & passFilters & passHLTAndObjectMatch & passAnalysisMuonSelection & passExactlyTwoAnalysisMuons '
                              '& passElectronVeto & passAdditionalLeptonVeto & passJetVetoMap & passBJetVeto & passVertex & '
                              'analysisPairPassesRawWindow & oppositeSign & passVBFCategoryApprox'},
          'weights': {'gen': 'weightGenerator', 'all': 'weightGenerator*weightPileup*weightMuonIDISO*weightTrigger'}}},
}
ROUND = "after_15Oct"                  # selected with --round
GROUPS = ROUNDS[ROUND]

SAMPLES = ["Data_Run2024H", "ggH", "VBFH"]            # default
ALL_SAMPLES = SAMPLES + ["ttH", "ZH", "WplusH"]       # available for some groups (--samples)
CHANNELS = ["ggH", "VBF", "inclusive"]    # inclusive = ggH or VBF (a group that only has an inclusive selection defines it directly)
# mass regions on the group's raw dimuon mass; the data are blinded in the signal window
REGIONS = {"Z_sideband": (70, 110), "H_sideband": None, "Signal_Fit": (115, 135)}
H_SIDEBAND = ((110, 115), (135, 150))
MASS_RANGE = (70, 150)


def defined(group, channel):
    """Whether a group provides a channel (inclusive exists whenever ggH and VBF exist)."""
    channels = GROUPS[group]["channels"]
    return channel in channels or (channel == "inclusive" and "ggH" in channels and "VBF" in channels)


def in_region(mass, region):
    if region == "H_sideband":
        return ((mass > 110) & (mass < 115)) | ((mass > 135) & (mass < 150))
    low, high = REGIONS[region]
    return (mass > low) & (mass < high)


# ----------------------------------------------------------------------------------------------
# Reading the ntuples
# ----------------------------------------------------------------------------------------------

def first(array):
    """First entry of a jagged branch, False/NaN when missing."""
    padded = ak.pad_none(array, 1, axis=1)[:, 0]
    return ak.fill_none(padded, False if "bool" in str(ak.type(array)) else np.nan)


def evaluate(arrays, expression):
    """Evaluate an expression on the branches of a chunk and return a flat numpy array."""
    expression = re.sub(r"(\w+)\[:,0\]", r"first(\1)", expression)
    scope = {name: arrays[name] for name in branch_names(expression) if name in arrays.fields}
    scope["first"] = first
    scope.update(abs=np.abs, sqrt=np.sqrt, hypot=np.hypot, cos=np.cos, sin=np.sin, cosh=np.cosh, sinh=np.sinh)
    return np.asarray(ak.to_numpy(eval(expression, {"__builtins__": {}}, scope)))


def branch_names(expression):
    expression = re.sub(r"(\w+)\[:,0\]", r"first(\1)", expression)
    return [n for n in dict.fromkeys(re.findall(r"[A-Za-z_]\w*", expression)) if n != "first"]


def sample_files(base, spec, sample):
    paths = []
    for pattern in spec["files"][sample]:
        paths += sorted(glob.glob(str(Path(base) / pattern)))
    if not paths:
        raise FileNotFoundError(f"no file for {sample}: {spec['files'][sample]} under {base}")
    return paths


def pack_keys(run, lumi, event, is_data):
    """One integer per event (run, lumi, event), used to match events between groups."""
    run, lumi, event = (np.asarray(x).astype(np.uint64) for x in (run, lumi, event))
    if is_data:
        return (run - np.uint64(370000)) << np.uint64(48) | lumi << np.uint64(34) | event
    return np.uint64(1) << np.uint64(63) | run << np.uint64(59) | lumi << np.uint64(39) | event


def read_sample(base, spec, sample, step_size=500_000):
    """Selected events of one group and sample: key, mass, channel flags and weights (70 < m < 150 GeV only)."""
    is_data = sample.startswith("Data")
    expressions = [spec["mass"], *spec["channels"].values()] + ([] if is_data else list(spec["weights"].values()))
    wanted = ["run", spec["lumi"], "event"]
    for expression in expressions:
        wanted += branch_names(expression)
    out = {k: [] for k in ("key", "mass", "ggH", "VBF", "inclusive", "w_gen", "w_all")}
    entries = 0
    for path in sample_files(base, spec, sample):
        tree = uproot.open(path)[spec["tree"]]
        names = [n for n in dict.fromkeys(wanted) if n in tree.keys()]
        for chunk in tree.iterate(names, step_size=step_size, library="ak"):
            entries += len(chunk)
            mass = evaluate(chunk, spec["mass"]).astype(np.float64)
            keep = (mass > MASS_RANGE[0]) & (mass < MASS_RANGE[1])
            channel = {c: evaluate(chunk, spec["channels"][c]).astype(bool) for c in CHANNELS if c in spec["channels"]}
            for c in CHANNELS:
                channel.setdefault(c, np.zeros(len(chunk), bool))
            if "inclusive" not in spec["channels"]:
                channel["inclusive"] = channel["ggH"] | channel["VBF"]
            keep &= channel["inclusive"]
            out["key"].append(pack_keys(chunk["run"], chunk[spec["lumi"]], chunk["event"], is_data)[keep])
            out["mass"].append(mass[keep])
            for c in CHANNELS:
                out[c].append(channel[c][keep])
            for mode, name in (("w_gen", "gen"), ("w_all", "all")):
                w = np.ones(len(chunk)) if is_data else evaluate(chunk, spec["weights"][name]).astype(np.float64)
                out[mode].append(w[keep])
    result = {k: np.concatenate(v) for k, v in out.items()}
    result["entries"] = np.array(entries)
    return result


def load_or_read(base, group, sample, outdir):
    """Read the ntuple once; the selected events are cached in <outdir>/cache/<group>_<sample>.npz."""
    cache = Path(outdir) / "cache" / ROUND / f"{group}_{sample}.npz"
    if cache.exists():
        d = dict(np.load(cache))
        if "inclusive" not in d:          # caches written before the inclusive channel existed
            d["inclusive"] = d["ggH"] | d["VBF"]
        return d
    print(f"[read] {ROUND} {group} {sample}", flush=True)
    result = read_sample(base, GROUPS[group], sample)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **result)
    return result


def deduplicate(d):
    """Count every (run, lumi, event) once per group and sample (some files contain events twice).
    `duplicates` is the number of rows dropped; for VBFH (main + ext1) a few identical identifiers of distinct events (<0.01%) are dropped too."""
    _, index = np.unique(d["key"], return_index=True)
    index = np.sort(index)
    out = {k: v[index] for k, v in d.items() if k != "entries"}
    out["entries"] = d["entries"]
    out["duplicates"] = np.array(len(d["key"]) - len(index))
    return out


# ----------------------------------------------------------------------------------------------
# Step 1: tables
# ----------------------------------------------------------------------------------------------

def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def cells(sample):
    """(region, channel) pairs of a sample; the data signal window is blinded."""
    return [(r, c) for r in REGIONS for c in CHANNELS if not (sample.startswith("Data") and r == "Signal_Fit")]


def number_of(rows, group, region, channel):
    """Events of a cell, None when the group does not provide the channel."""
    return next((int(r["events"]) for r in rows if (r["group"], r["region"], r["channel"]) == (group, region, channel)), None)


def step_yields(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    groups = args.groups or list(GROUPS)
    for sample in args.samples or SAMPLES:
        data = {g: deduplicate(load_or_read(args.base, g, sample, outdir)) for g in groups if sample in GROUPS[g]["files"]}
        yields, overlap, mass_rows = [], [], []
        for g, d in data.items():
            for region, channel in cells(sample):
                if not defined(g, channel):
                    continue
                sel = d[channel] & in_region(d["mass"], region)
                yields.append([g, sample, region, channel, int(sel.sum()), f"{d['w_gen'][sel].sum():.6g}", f"{d['w_all'][sel].sum():.6g}", int(d["duplicates"])])
            for channel in CHANNELS:
                if not defined(g, channel):
                    continue
                sel = d[channel]
                counts, edges = np.histogram(d["mass"][sel], bins=np.arange(MASS_RANGE[0], MASS_RANGE[1] + 1, 1.0))
                mass_rows += [[g, channel, edges[i], edges[i + 1], int(counts[i])] for i in range(len(counts))]
        names = list(data)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                for region, channel in cells(sample):
                    if not (defined(a, channel) and defined(b, channel)):
                        continue
                    ka = np.unique(data[a]["key"][data[a][channel] & in_region(data[a]["mass"], region)])
                    kb = np.unique(data[b]["key"][data[b][channel] & in_region(data[b]["mass"], region)])
                    common = len(np.intersect1d(ka, kb, assume_unique=True))
                    overlap.append([a, b, sample, region, channel, len(ka), len(kb), common,
                                    f"{common / max(len(ka), 1):.5f}", f"{common / max(len(kb), 1):.5f}"])
        write_csv(outdir / f"yields_{sample}.csv", ["group", "sample", "region", "channel", "events", "sum_gen_weight", "sum_all_weights", "duplicate_rows_dropped"], yields)
        write_csv(outdir / f"overlap_{sample}.csv", ["groupA", "groupB", "sample", "region", "channel", "nA", "nB", "common", "common_over_A", "common_over_B"], overlap)
        write_csv(outdir / f"mass_{sample}.csv", ["group", "channel", "low", "high", "events"], mass_rows)
        print(f"[tables] {sample}: {len(yields)} yield rows, {len(overlap)} overlap rows", flush=True)


# ----------------------------------------------------------------------------------------------
# Step 2: plots
# ----------------------------------------------------------------------------------------------

def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def colours(groups):
    return {g: plt.cm.tab10(i % 10) for i, g in enumerate(groups)}


def plot_bars(outdir, sample, rows):
    groups = list(dict.fromkeys(r["group"] for r in rows))
    labels = [f"{r}\n{c}" for r, c in cells(sample)]
    colour = colours(groups)
    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.8 / len(groups)
    for i, g in enumerate(groups):
        y = [number_of(rows, g, reg, ch) or 0 for reg, ch in cells(sample)]
        ax.bar(np.arange(len(labels)) + (i - (len(groups) - 1) / 2) * width, y, width, label=g, color=colour[g])
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("selected events")
    ax.set_title(f"{sample}: events per mass region and channel")
    ax.legend(ncol=len(groups), fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"yield_bars_{sample}.png", dpi=130)
    plt.close(fig)


def plot_ratio(outdir, sample, rows, reference):
    groups = [g for g in dict.fromkeys(r["group"] for r in rows) if g != reference]
    if reference not in {r["group"] for r in rows}:
        return
    labels = [f"{r}\n{c}" for r, c in cells(sample)]
    colour = colours(list(dict.fromkeys(r["group"] for r in rows)))
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.axhspan(0.99, 1.01, color="0.85", label="±1%")
    ax.axhspan(0.95, 1.05, color="0.93", zorder=0, label="±5%")
    for i, g in enumerate(groups):
        ratios = []
        for reg, ch in cells(sample):
            n_g, n_ref = number_of(rows, g, reg, ch), number_of(rows, reference, reg, ch)
            ratios.append(n_g / n_ref if n_g is not None and n_ref else np.nan)
        x = np.arange(len(labels)) + (i - (len(groups) - 1) / 2) * 0.12
        ax.plot(x, ratios, "o", color=colour[g], label=g)
    ax.axhline(1, color="k", lw=0.8)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 3)
    ax.set_yticks([0.1, 0.2, 0.5, 0.8, 0.9, 1, 1.1, 1.25, 2, 3])
    ax.set_yticklabels(["0.1", "0.2", "0.5", "0.8", "0.9", "1", "1.1", "1.25", "2", "3"])
    ax.set_ylabel(f"events / {reference}")
    ax.set_title(f"{sample}: yield agreement with {reference} (missing points: channel not provided)")
    ax.legend(ncol=len(groups) + 2, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"yield_ratio_{sample}_vs_{reference}.png", dpi=130)
    plt.close(fig)


def plot_overlap(outdir, sample, rows):
    groups = sorted({r["groupA"] for r in rows} | {r["groupB"] for r in rows})
    if len(groups) < 2:
        return
    fig, axes = plt.subplots(1, len(CHANNELS), figsize=(6 * len(CHANNELS), 5.2))
    for ax, channel in zip(axes, CHANNELS):
        matrix = np.full((len(groups), len(groups)), np.nan)
        for a_i, a in enumerate(groups):
            for b_i, b in enumerate(groups):
                if a == b:
                    continue
                nA = nB = common = 0
                for r in rows:
                    if r["channel"] == channel and {r["groupA"], r["groupB"]} == {a, b}:
                        nA += int(r["nA"]) if r["groupA"] == a else int(r["nB"])
                        nB += int(r["nB"]) if r["groupA"] == a else int(r["nA"])
                        common += int(r["common"])
                if nA and nB:
                    matrix[a_i, b_i] = common / nA     # fraction of the row group's events that the column group also has
        image = ax.imshow(matrix, vmin=0.5, vmax=1.0, cmap="viridis")
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=45, ha="right")
        ax.set_yticks(range(len(groups)))
        ax.set_yticklabels(groups)
        for i in range(len(groups)):
            for j in range(len(groups)):
                if not np.isnan(matrix[i, j]):
                    ax.text(j, i, f"{100 * matrix[i, j]:.1f}", ha="center", va="center", color="w", fontsize=8)
        ax.set_title(f"{sample} {channel}: common / events of the row group [%]")
        fig.colorbar(image, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(outdir / f"overlap_matrix_{sample}.png", dpi=130)
    plt.close(fig)


def plot_mass(outdir, sample, rows):
    fig, axes = plt.subplots(1, len(CHANNELS), figsize=(6 * len(CHANNELS), 4.5))
    groups = list(dict.fromkeys(r["group"] for r in rows))
    colour = colours(groups)
    for ax, channel in zip(axes, CHANNELS):
        for g in groups:
            sel = [r for r in rows if r["group"] == g and r["channel"] == channel]
            if not sel:
                continue
            edges = np.array([float(r["low"]) for r in sel] + [float(sel[-1]["high"])])
            counts = np.array([float(r["events"]) for r in sel])
            ax.stairs(counts, edges, label=g, color=colour[g])
        ax.set_xlabel("m(mumu) used for the regions [GeV]")
        ax.set_ylabel("events / GeV")
        ax.set_title(f"{sample} {channel}")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"mass_{sample}.png", dpi=130)
    plt.close(fig)


def step_plots(args):
    outdir = Path(args.outdir)
    for sample in args.samples or SAMPLES:
        if (outdir / f"yields_{sample}.csv").exists():
            rows = read_csv(outdir / f"yields_{sample}.csv")
            plot_bars(outdir, sample, rows)
            groups = list(dict.fromkeys(r["group"] for r in rows))
            for reference in (groups if args.reference == "all" else [args.reference]):
                plot_ratio(outdir, sample, rows, reference)
            plot_overlap(outdir, sample, read_csv(outdir / f"overlap_{sample}.csv"))
            plot_mass(outdir, sample, read_csv(outdir / f"mass_{sample}.csv"))
            print(f"[plots] {sample}", flush=True)


def md_table(header, rows):
    return ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] + ["| " + " | ".join(r) + " |" for r in rows] + [""]


def step_markdown(args):
    """tables.md: the yield tables and the overlap tables (common / events of the row group) of every sample, from the CSV files."""
    outdir = Path(args.outdir)
    lines = [f"# Yields and event overlap, round {ROUND}\n",
             "Events with the same (run, lumi, event) are counted once per group and sample; regions are on the mass branch of each group "
             "(see HOWTO.txt); data Signal_Fit is blinded.\n"]
    for sample in args.samples or SAMPLES:
        if not (outdir / f"yields_{sample}.csv").exists():
            continue
        rows = read_csv(outdir / f"yields_{sample}.csv")
        groups = list(dict.fromkeys(r["group"] for r in rows))
        lines += [f"## {sample}: selected events per region and channel\n"]
        cell_list = list(dict.fromkeys((r["region"], r["channel"]) for r in rows))
        lines += md_table(["region / channel"] + groups,
                          [[f"{reg} / {ch}"] + [("n/a" if number_of(rows, g, reg, ch) is None else f"{number_of(rows, g, reg, ch):,}") for g in groups]
                           for reg, ch in cell_list])
        lines += ["Duplicate rows dropped: " + ", ".join(f"{g} {int(next(r['duplicate_rows_dropped'] for r in rows if r['group'] == g)):,}" for g in groups) + "\n"]
        overlap = read_csv(outdir / f"overlap_{sample}.csv")
        names = sorted({r["groupA"] for r in overlap} | {r["groupB"] for r in overlap})
        for channel in CHANNELS:
            lines += [f"### {sample} {channel}: common events / events of the row group [%] (all regions summed)\n"]
            table = []
            for a in names:
                row = [a]
                for b in names:
                    n_a = common = 0
                    for r in overlap:
                        if a != b and r["channel"] == channel and {r["groupA"], r["groupB"]} == {a, b}:
                            n_a += int(r["nA"]) if r["groupA"] == a else int(r["nB"])
                            common += int(r["common"])
                    row.append("-" if a == b or not n_a else f"{100 * common / n_a:.1f}")
                table.append(row)
            lines += md_table(["row \\ column"] + names, table)
    (outdir / "tables.md").write_text("\n".join(lines))
    print(f"[markdown] {outdir / 'tables.md'}", flush=True)


def main():
    global ROUND, GROUPS
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="step", required=True)
    for name in ("yields", "plots", "markdown"):
        p = sub.add_parser(name)
        p.add_argument("--round", default="after_15Oct", choices=list(ROUNDS), help="which set of published ntuples")
        p.add_argument("--outdir", default=None, help="where tables, plots and the cache are written (default sync_out_<round>)")
        p.add_argument("--samples", nargs="+", choices=ALL_SAMPLES)
        p.add_argument("--groups", nargs="+", help="groups of the round (default: all of them)")
        if name == "yields":
            p.add_argument("--base", default="/eos/cms/store/group/phys_higgs/cmshmm/run3_sync", help="directory with one folder per group")
        elif name == "plots":
            p.add_argument("--reference", default="all", help="group the yield ratios are computed against (default: one plot per group)")
    args = parser.parse_args()
    ROUND, GROUPS = args.round, ROUNDS[args.round]
    args.outdir = args.outdir or f"sync_out_{ROUND}"
    unknown = [g for g in args.groups or [] if g not in GROUPS]
    if unknown:
        parser.error(f"groups not in round {ROUND}: {unknown}; available: {list(GROUPS)}")
    {"yields": step_yields, "plots": step_plots, "markdown": step_markdown}[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
