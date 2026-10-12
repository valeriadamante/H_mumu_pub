from pathlib import Path
import os
import ROOT
from .general import pog_folder_names,period_names


##### ID + trigger ####

MediumMu_SF_Sources_dict = {
    # reco SF
    "NUM_TrackerMuons_DEN_genTracks": "Reco",  # --> used in Run 2 - NOT FOR RUN 3!! https://muon-wiki.docs.cern.ch/guidelines/corrections/#medium-pt-30-gev-pt-200-gev "No correction is recommended for Run 3 data. Data/MC scale factors are expected to be 1 and therefore no correction is needed/provided."
    # ID SF - with tracker muons - RECOMMENDED
    "NUM_LooseID_DEN_TrackerMuons": "LooseID_Trk",
    "NUM_MediumID_DEN_genTracks": "MediumID_genTrk",
    "NUM_MediumPromptID_DEN_TrackerMuons": "MediumPromptID_Trk",
    "NUM_MediumID_DEN_TrackerMuons": "MediumID_Trk",
    "NUM_TightID_DEN_TrackerMuons": "TightID_Trk",
    "NUM_SoftID_DEN_TrackerMuons": "SoftID_Trk",
    "NUM_HighPtID_DEN_TrackerMuons": "HighPtID_Trk",
    "NUM_TrkHighPtID_DEN_TrackerMuons": "TrkHighPtID_Trk",
    # Iso SF
    "NUM_LoosePFIso_DEN_LooseID": "LoosePFIso_LooseID",  # loose ID, loose iso
    "NUM_LoosePFIso_DEN_MediumID": "LoosePFIso_MediumID",  # medium ID, loose  iso
    "NUM_LoosePFIso_DEN_MediumPromptID": "LoosePFIso_MediumPromptID",  # medium prompt ID, loose iso
    "NUM_LoosePFIso_DEN_TightID": "LoosePFIso_TightID",  # tight ID, tight PF Iso
    "NUM_LooseRelTkIso_DEN_HighPtID": "LooseRelTkIso_HighPtID",  # highPtID, loose tkRelIso
    "NUM_LooseRelTkIso_DEN_TrkHighPtID": "LooseRelTkIso_TrkHighPtID",  # trkHighPtID, loose tkRelIso
    "NUM_LooseRelIso_DEN_MediumID": "LooseRelIso_MediumID",
    "NUM_TightPFIso_DEN_MediumID": "TightPFIso_MediumID",  # medium ID, tight PF Iso
    "NUM_TightPFIso_DEN_MediumPromptID": "TightPFIso_MediumPromptID",  # medium prompt ID, tight PF Iso
    "NUM_TightRelTkIso_DEN_HighPtID": "TightRelIso_HighPtID",  # highPtID, tight tkRelIso
    "NUM_TightRelTkIso_DEN_TrkHighPtID": "TightRelIso_TrkHighPtID",  # trkHighPtID, tight tkRelIso
    "NUM_TightPFIso_DEN_TightID": "TightPFIso_TightID",  # tight ID, tight PF Iso
    "NUM_LooseMiniIso_DEN_LooseID": "LooseMiniIso_LooseID",
    "NUM_LooseMiniIso_DEN_MediumID": "LooseMiniIso_MediumID",
    "NUM_MediumMiniIso_DEN_MediumID": "MediumMiniIso_MediumID",
    "NUM_TightMiniIso_DEN_MediumID": "TightMiniIso_MediumID",
    "NUM_TightRelIso_DEN_MediumPromptID": "MediumRelIso",  # medium ID, tight iso --> old, Run 2 only
    "NUM_TightRelIso_DEN_TightIDandIPCut": "TightRelIso",  # tight ID, tight iso --> old, Run 2 only
    "NUM_TightRelTkIso_DEN_TrkHighPtIDandIPCut": "HighPtIdRelTkIso",  # highPtID, tight tkRelIso --> old, Run 2 only
    # Trigger
    # "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight": "TightIso24",  # trg --> FOR ALL PT RANGE!!
    # "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium": "MediumIso24",  # trg for medium muons
    "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight": "TightIso27",  # trg --> FOR ALL PT RANGE!!  --> old, Run 2 only
    "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": "Mu50",  # trg --> FOR ALL PT RANGE!! --> old, Run 2 only
    "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": "Mu50_tkMu50",  # trg --> FOR ALL PT RANGE!! --> old, Run 2 only
    "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight": "TightIso24OrTightIsoTk24",  # trg --> FOR ALL PT RANGE!! --> old, Run 2 only
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium": "IsoMu24_CutBasedIdMedium_and_PFIsoMedium",
    "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight": "IsoMu24_CutBasedIdTight_and_PFIsoTight",
    "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": "IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_CutBasedIdGlobalHighPt_and_TkIsoLoose",
    "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdMedium_and_PFIsoMedium": "IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_CutBasedIdMedium_and_PFIsoMedium",
    "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight": "IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_CutBasedIdTight_and_PFIsoTight",
    "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTrkHighPt_and_TkIsoLoose": "IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_CutBasedIdTrkHighPt_and_TkIsoLoose",
    "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": "Mu50_or_CascadeMu100_or_HighPtTkMu100_CutBasedIdGlobalHighPt_and_TkIsoLoose",
    "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTrkHighPt_and_TkIsoLoose": "Mu50_or_CascadeMu100_or_HighPtTkMu100_CutBasedIdTrkHighPt_and_TkIsoLoose",
}

MediumMuReco_SF_sources = {
    "2016preVFP_UL": ["NUM_TrackerMuons_DEN_genTracks"],
    "2016postVFP_UL": ["NUM_TrackerMuons_DEN_genTracks"],
    "2017_UL": ["NUM_TrackerMuons_DEN_genTracks"],
    "2018_UL": ["NUM_TrackerMuons_DEN_genTracks"],
    "2022_Summer22": [],  # 2022 has no genTrack to TrackerMuons key
    "2022_Summer22EE": [],
    "2023_Summer23": [],
    "2023_Summer23BPix": [],
    "2024_Winter24": [],
    "2024_Summer24": [],
    "2025_Summer24": [],
    "2025_Winter25": [],
}
MediumMuIDIso_SF_Sources = {
    "2016preVFP_UL": [
        "NUM_TightID_DEN_TrackerMuons",
        "NUM_TightRelIso_DEN_TightIDandIPCut",
    ],
    "2016postVFP_UL": [
        "NUM_TightID_DEN_TrackerMuons",
        "NUM_TightRelIso_DEN_TightIDandIPCut",
    ],
    "2017_UL": [
        "NUM_TightID_DEN_TrackerMuons",
        "NUM_TightRelIso_DEN_TightIDandIPCut",
    ],
    "2018_UL": [
        "NUM_TightID_DEN_TrackerMuons",
        "NUM_TightRelIso_DEN_TightIDandIPCut",
    ],
    "2022_Summer22": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2022_Summer22EE": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2023_Summer23": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2023_Summer23BPix": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2024_Summer24": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        # "NUM_LooseMiniIso_DEN_LooseID",
        # "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2024_Winter24": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2025_Summer24": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
    "2025_Winter25": [
        "NUM_LooseID_DEN_TrackerMuons",
        # "NUM_MediumID_DEN_genTracks",
        # "NUM_MediumPromptID_DEN_TrackerMuons",
        "NUM_MediumID_DEN_TrackerMuons",
        "NUM_TightID_DEN_TrackerMuons",
        # "NUM_SoftID_DEN_TrackerMuons",
        # "NUM_HighPtID_DEN_TrackerMuons",
        # "NUM_TrkHighPtID_DEN_TrackerMuons",
        "NUM_LoosePFIso_DEN_LooseID",
        "NUM_LoosePFIso_DEN_MediumID",
        # "NUM_LoosePFIso_DEN_MediumPromptID",
        "NUM_LoosePFIso_DEN_TightID",
        # "NUM_LooseRelTkIso_DEN_HighPtID",
        # "NUM_LooseRelTkIso_DEN_TrkHighPtID",
        # "NUM_LooseRelIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumID",
        # "NUM_TightPFIso_DEN_MediumPromptID",
        # "NUM_TightRelTkIso_DEN_HighPtID",
        # "NUM_TightRelTkIso_DEN_TrkHighPtID",
        "NUM_TightPFIso_DEN_TightID",
        "NUM_LooseMiniIso_DEN_LooseID",
        "NUM_LooseMiniIso_DEN_MediumID",
        # "NUM_MediumMiniIso_DEN_MediumID",
        # "NUM_TightMiniIso_DEN_MediumID",
    ],
}

MediumMuTrg_SF_Sources = {
    "2022_Summer22": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2022_Summer22EE": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2023_Summer23": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2023_Summer23BPix": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2024_Summer24": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2024_Winter24": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2025_Summer24": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
    "2025_Winter25": [
    "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"],
}
pt_mins = {
    "2026": 10,
    "2025": 10,
    "2024": 10,
    "2023BPix": 15,
    "2023": 15,
    "2022EE": 15,
    "2022": 15,
}
#!/usr/bin/env python3
import os
import ROOT
import correctionlib

correctionlib.register_pyroot_binding()



def apply_muIDIso_weights(df, config, want_variations):
    era = config.get("era")
    requested_SFs = config.get("requested_SFs", [])
    year = era.split("_")[1]
    pt_min = pt_mins.get(year, 15)
    pt_min_trg = 26.
    # Setup JSON paths using environment and shared functions
    period_unc = period_names[era]
    muIDEff_JsonPath = (
        "/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/{}/latest/muon_Z.json.gz"
    )
    if era == "Run3_2025":
        # Not on cvmfs yet: the file ships with the repository (relative
        # to ANALYSIS_PATH through the join below).
        muIDEff_JsonPath = (
            "corrections/data/MUO/SF/Run3-25Prompt-Summer24-NanoAODv15/muon_Z.json.gz"
        )
    jsonFile_path = os.path.join(
        os.environ["ANALYSIS_PATH"],
        muIDEff_JsonPath.format(pog_folder_names["MUO"][period_unc]),
    )
    # Global registration for access in C++ context
    ROOT.gROOT.ProcessLine(
        f'auto cset = correction::CorrectionSet::from_file("{jsonFile_path}");'
    )
    headers_dir = os.path.dirname(os.path.abspath(__file__))
    header_path = os.path.join(headers_dir, "mu.h")
    ROOT.gInterpreter.Declare(f'#include "{header_path}"')

    # Gather available keys for this era
    available_sources = (
        MediumMuIDIso_SF_Sources.get(period_unc, [])
        + MediumMuReco_SF_sources.get(period_unc, [])
    )
    available_trg_sf_sources = MediumMuTrg_SF_Sources.get(period_unc, [])
    # Dict mapping python scale terminology to correctionlib JSON parameters
    scale_map = {"Central": "nominal", "up": "systup", "down": "systdown"}

    SF_branches = []


    # Define standard input column strings matching standard NanoAOD layout
    for leg_idx in [1,2]:
        p4_pt = f"mu{leg_idx}_pt_raw_noCorr" # no corr?
        p4_eta = f"mu{leg_idx}_eta"
        pfRelIso04_all = f"mu{leg_idx}_pfRelIso04_all"
        tightId = f"mu{leg_idx}_tightId"
        tkRelIso = f"mu{leg_idx}_tkRelIso"
        # highPtId = f"mu{leg_idx}_highPtId"
        mediumId = f"mu{leg_idx}_mediumId"
        looseId = f"mu{leg_idx}_looseId"
        gen_kind = f"mu{leg_idx}_genPartFlav"
        trg_matching = f"mu{leg_idx}_HasTriggerMatching_singleMu"
        trg_path = "HLT_IsoMu24"

        # genMatch_bool = f"{gen_kind} == 1 || {gen_kind} == 15" # for MC matching to status==1 muons: 1 = prompt muon (including gamma*->mu mu), 15 = muon from prompt tau, 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched

        for source in available_sources + available_trg_sf_sources:
            short_name = MediumMu_SF_Sources_dict.get(source)

            # Skip execution entirely if a target list was given and this key isn't in it
            if requested_SFs and (short_name not in requested_SFs):
                continue
            pt_min_sf = pt_min_trg if source in available_trg_sf_sources else pt_min
            # Loop through all scales to ALWAYS define them inside the RDataFrame
            for scale, cset_syst_string in scale_map.items():
                branch_name = f"weight_mu{leg_idx}_{short_name}_{scale}"
                # print(branch_name)
                # Direct definition inside RDataFrame using your new C++ signature
                df = df.Define(
                    branch_name,
                    f"""({p4_pt} > {pt_min_sf} && std::abs({p4_eta}) < 2.4)
                        ? static_cast<float>(::correction::getMuonSF_simple(
                            "{source}", "{cset_syst_string}",
                            {p4_pt}, {p4_eta}, {pfRelIso04_all},
                            {tightId}, {tkRelIso}, {mediumId}, {looseId},{trg_matching},{trg_path}, event))
                        : 1.0f""",
                )

                # if source == "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium":
                # df.Display({branch_name}).Print()

                # Control what goes to the downstream storing list:
                # Always keep Central; only add up/down to track if want_variations is True
                if scale == "Central" or want_variations:
                    SF_branches.append(branch_name)

    return df, SF_branches
