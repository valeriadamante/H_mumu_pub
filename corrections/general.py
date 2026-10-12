import os
import sys
import ROOT
import json

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])

import common.utilities as utilities


pog_folder_names = {
    "BTV": {
        "2016postVFP_UL": "Run2-2016postVFP-UL-NanoAODv9",
        "2016preVFP_UL": "Run2-2016preVFP-UL-NanoAODv9",
        "2017_UL": "Run2-2017-UL-NanoAODv9",
        "2018_UL": "Run2-2018-UL-NanoAODv9",
        "JER": "Run3-22CDJun23-Summer22-NanoAODv11",
        "2022_Summer22": "Run3-22CDSep23-Summer22-NanoAODv12",
        "2022_Summer22EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        "2023_Summer23": "Run3-23CSep23-Summer23-NanoAODv12",
        "2023_Summer23BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2024_Summer24": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2025_Summer24": "Run3-25Prompt-Summer24-NanoAODv15",
        # "2025_Winter25": "",
        # "2026_Summer24":"Run3-26Prompt-Summer24-NanoAODv15",
        "2026_Summer24":"Run3-25Prompt-Summer24-NanoAODv15", # no 25 for BTV
    },
    "JERC": {
        "2018_UL": "Run2-2018-UL-NanoAODv9",
        "2017_UL": "Run2-2017-UL-NanoAODv9",
        "2016preVFP_UL": "Run2-2016preVFP-UL-NanoAODv9",
        "2016postVFP_UL": "Run2-2016postVFP-UL-NanoAODv9",
        "JER": "JER-Smearing",
        "2022_Summer22": "Run3-22CDSep23-Summer22-NanoAODv12",  #
        "2022_Prompt": "Run3-22Prompt-Winter22-NanoAODv12",
        "2022_Summer22EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        "2023_Summer23BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2023_Summer23": "Run3-23CSep23-Summer23-NanoAODv12",
        "2024_Winter24": "Run3-24Prompt-Winter24-NanoAODv14",
        "2024_Summer24": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",  # https://cms-jerc.web.cern.ch/Recommendations/#2024
        "2025_Summer24": "Run3-25Prompt-Summer24-NanoAODv15",
        # "2025_Winter25": "Run3-25Prompt-Winter25-NanoAODv15",
        "2026_Summer24":"Run3-26Prompt-Summer24-NanoAODv15",
    },
    "EGM": {
        "2016postVFP_UL": "Run2-2016postVFP-UL-NanoAODv15",
        # "2016postVFP_UL":"Run2-2016postVFP-UL-NanoAODv9",
        "2016preVFP_UL": "Run2-2016preVFP-UL-NanoAODv15",
        # "2016preVFP_UL":"Run2-2016preVFP-UL-NanoAODv9",
        "2017_UL": "Run2-2017-UL-NanoAODv15",
        # "2017_UL":"Run2-2017-UL-NanoAODv9",
        "2018_UL": "Run2-2018-UL-NanoAODv15",
        # "2018_UL":"Run2-2018-UL-NanoAODv9",
        "2022_Summer22": "Run3-22CDSep23-Summer22-NanoAODv12",
        "2022_Summer22EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        "2023_Summer23": "Run3-23CSep23-Summer23-NanoAODv12",
        "2023_Summer23BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2024_Summer24": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2025_Summer24": "Run3-25Prompt-Summer24-NanoAODv15",
        # "2025_Winter25": "",
        "2026_Summer24":"Run3-26Prompt-Summer24-NanoAODv15"
    },
    "LUM": {
        "2016postVFP_UL": "Run2-2016postVFP-UL-NanoAODv9",
        "2016preVFP_UL": "Run2-2016preVFP-UL-NanoAODv9",
        "2017_UL": "Run2-2017-UL-NanoAODv9",
        "2018_UL": "Run2-2018-UL-NanoAODv9",
        "2022_Summer22": "Run3-22CDSep23-Summer22-NanoAODv12",
        "2022_Summer22EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        "2023_Summer23": "Run3-23CSep23-Summer23-NanoAODv12",
        "2023_Summer23BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2024_Summer24": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2025_Summer24": "Run3-25Prompt-Summer24-NanoAODv15",
        # "2025_Winter25": "",
        # "2026_Summer24": "Run3-26Prompt-Summer24-NanoAODv15"
        "2026_Summer24": "Run3-25Prompt-Summer24-NanoAODv15"
    },
    "MUO": {
        "2016postVFP_UL": "Run2-2016postVFP-UL-NanoAODv9",
        "2016preVFP_UL": "Run2-2016preVFP-UL-NanoAODv9",
        "2017_UL": "Run2-2017-UL-NanoAODv9",
        "2018_UL": "Run2-2018-UL-NanoAODv9",
        "2022_Summer22": "Run3-22CDSep23-Summer22-NanoAODv12",
        "2022_Summer22EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        "2023_Summer23": "Run3-23CSep23-Summer23-NanoAODv12",
        "2023_Summer23BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2024_Summer24": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2025_Summer24": "Run3-25Prompt-Summer24-NanoAODv15",
        # "2025_Winter25": "",
        "2026_Summer24": "Run3-25Prompt-Summer24-NanoAODv15" #  PATCH
        # "2026_Summer24": "Run3-26Prompt-Summer24-NanoAODv15"
    },
    "TAU": {
        "2016postVFP_UL": "Run2-2016postVFP-UL-NanoAODv15",
        # "2016postVFP_UL":"Run2-2016postVFP-UL-NanoAODv9",
        "2016preVFP_UL": "Run2-2016preVFP-UL-NanoAODv15",
        # "2016preVFP_UL":"Run2-2016preVFP-UL-NanoAODv9",
        "2017_UL": "Run2-2017-UL-NanoAODv15",
        # "2017_UL":"Run2-2017-UL-NanoAODv9",
        "2018_UL": "Run2-2018-UL-NanoAODv15",
        # "2018_UL":"Run2-2018-UL-NanoAODv9",
        "2022_Summer22": "Run3-22CDSep23-Summer22-NanoAODv12",
        "2022_Summer22EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        "2023_Summer23": "Run3-23CSep23-Summer23-NanoAODv12",
        "2023_Summer23BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2024_Summer24": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2025_Summer24": "",
        # "2025_Winter25": "",
        "2026_Summer24":""
    },
}

period_names = {
    "Run2_2016_HIPM": "2016preVFP_UL",
    "Run2_2016": "2016postVFP_UL",
    "Run2_2017": "2017_UL",
    "Run2_2018": "2018_UL",
    "Run3_2022": "2022_Summer22",
    "Run3_2022EE": "2022_Summer22EE",
    "Run3_2023": "2023_Summer23",
    "Run3_2023BPix": "2023_Summer23BPix",
    "Run3_2024": "2024_Summer24",  # 2024_Winter24
    "Run3_2025": "2025_Summer24",  # "2025_Winter25" is also a valid entry, but has files only only for JME
    "Run3_2026": "2026_Summer24",  # TEMPORARY PATCH TO CHECK THAT WORKS FOR 2026 DATA
}

periods = {
    "2026_Summer24": "2026",
    "2025_Winter25": "2025",
    "2025_Summer24": "2025",
    "2024_Summer24": "2024",
    "2023_Summer23BPix": "2023",
    "2023_Summer23": "2023",
    "2022_Summer22EE": "2022",
    "2022_Summer22": "2022",
    "2018_UL": "2018",
    "2017_UL": "2017",
    "2016preVFP_UL": "2016",
    "2016postVFP_UL": "2016",
}

def define_base_weights(df, lumi, xs_entry, xs_cfg,config,dataset_cfg, process_entry, want_variations_from_skim=False, systematics_cfg=None):
    want_variations = config.get("want_variations", False) or want_variations_from_skim
    base_weights_to_store = []
    json_dict_to_store = {}

    lumi_weight_name = "weight_lumi"
    df = df.Define(lumi_weight_name, f"float({lumi})")
    base_weights_to_store.append(lumi_weight_name)

    df = df.Define("weight_gen_sign", "std::copysign<float>(1.f, genWeight)")
    base_weights_to_store.append("genWeight")
    base_weights_to_store.append("weight_gen_sign")
    if 'gen' not in json_dict_to_store.keys():
        json_dict_to_store['gen'] = {}
    json_dict_to_store['gen'][f"total"] = {}
    json_dict_to_store['gen'][f"total"]['selection']= "return true;"
    json_dict_to_store['gen'][f"total"]['value']= df.Sum("genWeight")

    from .pu import apply_pu_weights
    df,pu_branches,json_dict_to_store = apply_pu_weights(df, json_dict_to_store, config, "Pileup_nTrueInt",want_variations)
    base_weights_to_store.extend(pu_branches)

    weight_xs_name = "weight_xs"
    processor_file = process_entry.get('processors','')
    if processor_file:
        # print(f"considering processor file {processor_file} for {process_entry}")
        from .dy_cross_section import apply_dy_cross_section
        processor_cfg = utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], processor_file))
        df,json_dict_xs = apply_dy_cross_section(df, weight_xs_name, xs_cfg, processor_cfg,json_dict_to_store)
    else:
        xs_value = str(xs_cfg[xs_entry]['crossSec'])
        df = df.Define("weight_xs", xs_value)
    base_weights_to_store.append(weight_xs_name)

    if systematics_cfg and systematics_cfg.get("qcd_scale"):
        from .qcd_scale import define_qcd_scale_sum_columns
        df, json_dict_to_store = define_qcd_scale_sum_columns(
            df,
            systematics_cfg["qcd_scale"],
            json_dict_to_store,
        )

    ### to add --> 2024/2025/2026 --> weight 0 or 1 to decide what MC use for what prod


    return df,base_weights_to_store,json_dict_to_store

def apply_golden_json(df, lumiFile_path):
    from .lumi import apply_lumi_filter
    df = apply_lumi_filter(df, lumiFile_path)
    return df

def apply_corrections(df, config, dataset_cfg, dataset_name, want_variations_from_skim=False):
    is_data = dataset_cfg.get("is_data", False)
    want_variations = config.get("want_variations", False) or want_variations_from_skim
    # Compatibility patch for missing branches in Flashsim samples.
    from common.flashsim_patch import apply_flashsim_patch
    df = apply_flashsim_patch(df, dataset_name, is_data)
    # muons ScaRe
    from .muon_scare import apply_muon_scare
    df = apply_muon_scare(df, config, dataset_cfg,want_variations)
    # muons FSR
    from .muon_fsr import apply_muon_fsr
    df = apply_muon_fsr(df, is_data, want_variations)
    # JEC / JER / JES
    from .jets import apply_jet_corrections
    df = apply_jet_corrections(df, config, dataset_cfg, dataset_name, want_variations)
    return df

