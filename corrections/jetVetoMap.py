from common.jet_variation_suffixes import jet_variation_suffixes
#!/usr/bin/env python3
import os
import ROOT
import correctionlib

from .general import pog_folder_names, period_names

correctionlib.register_pyroot_binding()
def _column_names(df):
    return {str(col) for col in df.GetColumnNames()}

def InitializeVetoMap(config):
    JME_vetoMap_JsonPath = (
        "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{}/latest/jetvetomaps.json.gz"
    )
    era = config.get("era")
    period_unc = period_names[era]

    jetvetomap_names = {
        "2022_Summer22": "Summer22_23Sep2023_RunCD_V1",  # period: , entryname
        "2022_Summer22EE": "Summer22EE_23Sep2023_RunEFG_V1",
        "2023_Summer23BPix": "Summer23BPixPrompt23_RunD_V1",
        "2023_Summer23": "Summer23Prompt23_RunC_V1",
        "2024_Summer24": "Summer24Prompt24_RunBCDEFGHI_V1",
        "2025_Winter25": "Winter25Prompt25_RunCDEFG_V1",
        "2025_Summer24": "Summer24Prompt25_RunCDEFG_V1",
        "2026_Summer24": "Summer24Prompt26_RunBCD_V1", # tmp patch as there is no Jet veto maps for 2026
    }
    entry_name = jetvetomap_names[period_unc]
    JME_vetoMap_JsonFile =JME_vetoMap_JsonPath.format(pog_folder_names["JERC"][period_unc])

    headers_dir = os.path.dirname(os.path.abspath(__file__))
    header_path = os.path.join(headers_dir, "jetVetoMap.h")
    ROOT.gInterpreter.Declare(f'#include "{header_path}"')
    ROOT.gInterpreter.ProcessLine(
        f'::correction::JetVetoMapProvider::Initialize("{JME_vetoMap_JsonFile}", "{entry_name}")'
    )


def ApplyJetVetoMap(df, config, muon_default_suffix, apply_filter, defineElectronCleaning, isV12, want_variations,syst_cfg):
    InitializeVetoMap(config)
    new_cols = []

    def track(df, name, expr):
        if name not in new_cols:
            new_cols.append(name)
        if name in _column_names(df): return df
        return df.Define(name, expr)

    cols = _column_names(df)
    syst_suffixes = jet_variation_suffixes(df, want_variations, syst_cfg)

    for suff in syst_suffixes:
        p4_branch = f"Jet_p4{suff}"
        function_for_jetId = (
            f"JetIdNewDef::RedefineJet_passJetIdTight_v12({p4_branch}, Jet_neHEF, Jet_neEmEF, Jet_jetId)"
            if isV12
            else f"JetIdNewDef::RedefineJet_passJetIdTight_v13({p4_branch}, Jet_neHEF, Jet_neEmEF, Jet_chHEF, Jet_chMultiplicity, Jet_neMultiplicity )"
        )
        df = track(df,f"Jet_passJetIdTight{suff}", function_for_jetId)
        df = track(df,f"Jet_passJetIdTightLepVeto{suff}", f"JetIdNewDef::Redefine_Jet_passJetIdTightLepVeto({p4_branch}, Jet_passJetIdTight, Jet_muEF, Jet_chEmEF)")
        df = track(df,f"Jet_isInsideVetoRegion{suff}",f"""::correction::JetVetoMapProvider::getGlobal().GetJetVetoMapValues({p4_branch})""")
        df = track(df,f"Jet_vetoMapLooseRegion_presel{suff}",f"Jet_pt{suff} > 15 && ( Jet_passJetIdTightLepVeto ) && (Jet_chEmEF + Jet_neEmEF < 0.9) && Jet_isInsideVetoRegion")
        df = track(df,f"Jet_vetoMap{suff}",f"RemoveOverlaps({p4_branch}, Jet_vetoMapLooseRegion_presel{suff}, Muon_p4_{muon_default_suffix}[Muon_isPFcand], 0.2)")


        if defineElectronCleaning:
            df = track(df,f"Jet_vetoMapEle{suff}",f" RemoveOverlaps({p4_branch}, Jet_vetoMap, Electron_p4[Electron_isPFcand], 0.2)")

        if apply_filter:
            return df.Filter(f"{p4_branch}[{jet_veto_map_string}].size()==0", "Jet Veto Map filter")
    return df,new_cols
