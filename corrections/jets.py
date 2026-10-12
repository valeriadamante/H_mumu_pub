import os
import ROOT
from .general import pog_folder_names, period_names

# Jet correction configuration and helper utilities for Hmumu_analysis.

jsonPath_btag = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/{}/btagging.json.gz"
jet_jsonPath = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{}/latest/jet_jerc.json.gz"
jetsmear_jsonFile = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/JER-Smearing/latest/jer_smear.json.gz"
jet_algorithm = "AK4PFPuppi"
uncSources_minimal = ["Total"]
unc_sources_regrouped = ["Total", "Regrouped_RelativeBal", "Regrouped_HF", "Regrouped_BBEC1", "Regrouped_EC2", "Regrouped_Absolute", "Regrouped_FlavorQCD", "Regrouped_BBEC1_{year}", "Regrouped_Absolute_{year}", "Regrouped_EC2_{year}", "Regrouped_HF_{year}", "RelativeSample_{year}"]

unc_source_enum = {
    "Central": "Central",
    "JER": "JER",
    "JESTotal": "Total",
    "Total": "Total",
    "RelativeBal": "RelativeBal",
    "HF": "HF",
    "BBEC1": "BBEC1",
    "EC2": "EC2",
    "Absolute": "Absolute",
    "FlavorQCD": "FlavorQCD",
    "BBEC1_year": "BBEC1_year",
    "Absolute_year": "Absolute_year",
    "EC2_year": "EC2_year",
    "HF_year": "HF_year",
    "RelativeSample_year": "RelativeSample_year",
}

def _jet_uncertainty_year(period):
    return {
        "2022_Summer22EE": "2022EE",
        "2023_Summer23BPix": "2023BPix",
    }.get(period, period.split("_")[0])


def _jet_uncertainty_enum(source, period):
    source = source.format(year=_jet_uncertainty_year(period))
    if source.startswith("Regrouped_"):
        source = source[len("Regrouped_"):]
    suffix = "_" + _jet_uncertainty_year(period)
    if source.endswith(suffix):
        source = source[:-len(suffix)] + "_year"
    return _debug_map_get(unc_source_enum, source, "unc_source_enum")


jer_tag_map = {
    "2022_Summer22": "Summer22_22Sep2023_JRV2_MC",
    "2022_Prompt": "JR_Winter22Run3_V2_MC",
    "2022_Summer22EE": "Summer22EE_22Sep2023_JRV2_MC",
    "2023_Summer23BPix": "Summer23BPixPrompt23_RunD_JRV3_MC",
    "2023_Summer23": "Summer23Prompt23_RunCv1234_JRV3_MC",
    "2024_Summer24": "Summer24Prompt24_JRV2_MC",
    "2025_Summer24": "Summer24Prompt25_JRV2_MC",
    # "2025_Winter25": "Summer24Prompt25_JRV2_MC",
    "2026_Summer24": "Summer24Prompt26_RunBD_JRV1_MC"
}

jec_tag_map_mc = {
    "2022_Prompt": ["Winter22Run3_V4_MC"],
    "2022_Summer22": ["Summer22_22Sep2023_V4_MC"],
    "2022_Summer22EE": ["Summer22EE_22Sep2023_V4_MC"],
    "2023_Summer23BPix": ["Summer23BPixPrompt23_V4_MC"],
    "2023_Summer23": ["Summer23Prompt23_V4_MC"],
    "2024_Summer24": ["Summer24Prompt24_V5_MC"],
    "2025_Summer24": ["Summer24Prompt25_V3_MC"],
    # "2025_Winter25": ["Winter25Prompt25_V3_MC"],
    "2026_Summer24": ["Summer24Prompt26_V1_MC"],
}

jec_tag_map_data = {
    "2022_Prompt": ["Winter22Run3_Run{}_V4_DATA"],
    "2022_Summer22": ["Summer22_22Sep2023_V4_DATA"],
    "2022_Summer22EE": ["Summer22EE_22Sep2023_Run{}_V4_DATA", "Summer22EE_22Sep2023_V4_DATA"],
    "2023_Summer23BPix": ["Summer23BPixPrompt23_Run{}_V4_DATA", "Summer23BPixPrompt23_V4_DATA"],
    "2023_Summer23": ["Summer23Prompt23_Run{}_V4_DATA", "Summer23Prompt23_V4_DATA"],
    "2024_Summer24": ["Summer24Prompt24_V5_DATA"],
    "2025_Summer24": ["Summer24Prompt25_V3_DATA"],
    # "2025_Winter25": ["Winter25Prompt25_Run{}_V3_DATA", "Winter25Prompt25_V3_DATA"],
    "2026_Summer24": ["Summer24Prompt26_Run{}_V1_DATA", "Summer24Prompt26_V1_DATA"], # tmp patch as there is no JEC tag for 2026 right now
}

run_versions = {
    "2022_Summer22": [],
    "2023_Summer23BPix": [],
    "2022_Prompt": [],
    "2023_Summer23": ["v123", "v4"],
    "2022_Summer22EE": [],
    "2024_Winter24": [],
    "2024_Summer24": [],
    "2025_Summer24": [],
    "2025_Winter25": [],
    "2026_Summer24": [],

}

run_letters = {
    "2022_Summer22": ["CD"],
    "2023_Summer23BPix": ["D"],
    "2022_Prompt": ["C", "D"],
    "2023_Summer23": ["C"],
    "2022_Summer22EE": ["E", "F", "G"],
    "2024_Winter24": ["BCD", "E", "F", "G", "H"],
    "2024_Summer24": ["CDEReprocessing", "FGHIPrompt"],
    "2025_Winter25": ["C", "D", "E", "F","G"],
    "2025_Summer24": ["C", "D", "E", "F", "G"],
    "2026_Summer24": ["A", "B", "C", "D"],
}

_jet_correction_state = {
    "initialized": False,
    "period": None,
    "is_data": False,
    "use_regrouped": True,
    "sample_name": None,
    "uncSources_toUse": unc_sources_regrouped,
}


def _jets_debug(message):
    if os.environ.get("JETS_DEBUG", "0") != "0":
        print(f"[jets.py][DEBUG] {message}", flush=True)


def _debug_map_get(mapping, key, map_name):
    try:
        value = mapping[key]
    except KeyError:
        print(
            f"[jets.py][MAP LOOKUP FAILED] {map_name}[{key!r}]; "
            f"available keys={list(mapping)}",
            flush=True,
        )
        raise
    _jets_debug(f"{map_name}[{key!r}] -> {value!r}")
    return value


def _format_data_jec_tags(period, sample_name, jec_tag_array):
    sample_letter = ""
    sample_version = ""
    if sample_name[-1].isalpha():
        sample_letter = sample_name[-1]
    elif sample_name[-1].isnumeric():
        tokens = sample_name.split("_")
        sample_version = tokens[-1]
        sample_letter = tokens[-2][-1]
    if period == "2025_Summer24" or period=="2026_Summer24": # tmp patch as it does not exist for 2026
        sample_version = ""

    if sample_letter not in run_letters[period]:
        matches = [let for let in run_letters[period] if sample_letter in let]
        if len(matches) != 1:
            raise RuntimeError(
                f"ambiguous deduction of sample letter for {sample_name}: got letter options {matches}"
            )
        sample_letter = matches[0]

    version_list = run_versions[period]
    if version_list and sample_version not in version_list:
        matches = [v for v in version_list if sample_version[1:] in v]
        if len(matches) != 1:
            raise RuntimeError(
                f"ambiguous deduction of sample version for {sample_name}: got version options {matches}"
            )
        sample_version = matches[0]

    if not sample_letter and not sample_version:
        raise RuntimeError(
            f"sample name {sample_name} doesn't follow expected pattern base_letter_version"
        )

    letters = sample_letter + sample_version
    return [tag.format(letters) for tag in jec_tag_array]


def initialize_jet_corrections(
    period,
    is_data,
    sample_name,
    use_regrouped=False,
):

    if _jet_correction_state["initialized"]:
        same_state = (
            _jet_correction_state["period"] == period
            and _jet_correction_state["is_data"] == is_data
            and _jet_correction_state["sample_name"] == sample_name
            and _jet_correction_state["use_regrouped"] == use_regrouped
        )
        if same_state:
            return
        raise RuntimeError(
            "Jet corrections already initialized with a different configuration"
        )

    _jet_correction_state["period"] = period
    _jet_correction_state["is_data"] = is_data
    _jet_correction_state["use_regrouped"] = use_regrouped
    _jet_correction_state["sample_name"] = sample_name
    year_unc = _jet_uncertainty_year(period)
    _jet_correction_state["uncSources_toUse"] = ["JER"] + [
        source.format(year=year_unc)
        for source in (unc_sources_regrouped if use_regrouped else uncSources_minimal)
    ]

    jec_jsonFile = jet_jsonPath.format(pog_folder_names["JERC"][period])
    jer_jsonFile = jet_jsonPath.format(pog_folder_names["JERC"][period])
    print(f"jec_jsonFile. = {jec_jsonFile}, jer_jsonFile={jer_jsonFile}")
    year = period.split("_")[0]
    jec_year = _jet_uncertainty_year(period)

    jec_tag_map = jec_tag_map_data if is_data else jec_tag_map_mc
    jec_tag_array = _debug_map_get(jec_tag_map, period, "jec_tag_map")
    if is_data:
        jec_tag_array = _format_data_jec_tags(period, sample_name, jec_tag_array)
    print(jec_tag_array)
    jec_tag = jec_tag_array[0]
    other_jec_tag = jec_tag_array[1] if len(jec_tag_array) > 1 else jec_tag_array[0]
    jer_tag = _debug_map_get(jer_tag_map, period, "jer_tag_map")
    if "C" in jec_tag_array and period == "2026_Summer24":
        jer_tag = "Summer24Prompt26_RunC_JRV1_MC"
    print(jer_tag)
    _jets_debug(
        f"period={period}, is_data={is_data}, sample={sample_name}, "
        f"regrouped={use_regrouped}, period={period}, "
        f"jec_jsonFile={jec_jsonFile}, jer_jsonFile={jer_jsonFile}, "
        f"jec_tag={jec_tag}, other_jec_tag={other_jec_tag}, jer_tag={jer_tag}, "
        f"algorithm={jet_algorithm}, nominal_year={year}, jec_year={jec_year}"
    )

    headers_dir = os.path.dirname(os.path.abspath(__file__))
    header_path = os.path.join(headers_dir, "jets.h")
    ROOT.gInterpreter.Declare(f'#include "{header_path}"')

    is_data_str = "true" if is_data else "false"
    regrouped_str = "true" if use_regrouped else "false"
    apply_compound = "true"

    ROOT.gInterpreter.ProcessLine(
        f"""::correction::JetCorrectionProvider::Initialize(\
            \"{jec_jsonFile}\",\
            \"{jer_jsonFile}\",\
            \"{jetsmear_jsonFile}\",\
            \"{jec_tag}\",\
            \"{other_jec_tag}\",\
            \"{jer_tag}\",\
            \"{jet_algorithm}\",\
            \"{year}\",\
            \"{jec_year}\",\
            {is_data_str},\
            {regrouped_str},\
            {apply_compound})"""
    )

    _jet_correction_state["initialized"] = True

# ============================================================
# JET CORRECTIONS ONLY (JEC/JER PROVIDER LAYER)
# ============================================================

def define_jet_p4_variations(
    df,
    want_variations,
    apply_JER,
    apply_JES,
    apply_horn_mitigation=True,
):
    if not _jet_correction_state["initialized"]:
        raise RuntimeError("Jet corrections are not initialized")

    is_data = _jet_correction_state["is_data"]
    period = _jet_correction_state["period"]

    apply_jer = "true" if apply_JER and not is_data else "false"
    horn_mitigation = "true" if apply_horn_mitigation else "false"
    reapply_jec = "true"
    require_run_number = "true" if is_data or period == "2023_Summer23BPix" else "false"

    wantPhi = (
        "true"
        if ((period in ["2023_Summer23BPix"] and is_data) or (period in ["2024_Summer24", "2025_Summer24", "2026_Summer24"])) # tmp patch as it does not exist for 2026
        else "false"
    )


    df = df.Define("seed_jersmearing","(run << 20) + (luminosityBlock << 10) + event + (Jet_eta.size() > 0 ? Jet_eta[0]/0.01 : 0)")
    # ===========================
    # CORE SHIFTED MAP
    # ===========================
    if not is_data:
        df = df.Define(
            "Jet_p4_shifted_map",
            f"""::correction::JetCorrectionProvider::getGlobal().getShiftedP4(
                Jet_pt, Jet_eta, Jet_phi, Jet_mass,
                Jet_rawFactor, Jet_area,
                Rho_fixedGridRhoFastjetAll,
                event,
                seed_jersmearing,
                {apply_jer},
                {reapply_jec},
                {require_run_number},
                run,
                {wantPhi},
                GenJet_pt, GenJet_eta, GenJet_phi, Jet_genJetIdx,
                {horn_mitigation}
            )"""
        )
    else:
        df = df.Define(
            "Jet_p4_shifted_map",
            f"""::correction::JetCorrectionProvider::getGlobal().getShiftedP4(
                Jet_pt, Jet_eta, Jet_phi, Jet_mass,
                Jet_rawFactor, Jet_area,
                Rho_fixedGridRhoFastjetAll,
                event,
                seed_jersmearing,
                {apply_jer},
                {reapply_jec},
                {require_run_number},
                run,
                {wantPhi},
                {{}}, {{}}, {{}}, {{}}, {horn_mitigation}
            )"""
        )

    # ===========================
    # helper: extract p4
    # ===========================
    def p4_from(unc_source, unc_scale):
        src = _jet_uncertainty_enum(unc_source, period)
        expression = (
            "::correction::JetCorrectionProvider::getP4FromMap("
            "Jet_p4_shifted_map, "
            f"::correction::JetCorrectionProvider::UncSource::{src}, "
            f"::correction::UncScale::{unc_scale})"
        )
        _jets_debug(f"RDataFrame map lookup {unc_source}:{unc_scale} -> {expression}")
        return expression

    # ===========================
    # nominal only p4 (central)
    # ===========================
    df = df.Define("Jet_p4", p4_from("Central", "Central"))

    # ONLY VARIATIONS (NO DERIVED FEATURES HERE)
    if not is_data and want_variations:
        sources = []
        sources.append("JER")

        sources += [
            s.format(year=_jet_uncertainty_year(period))
            for s in _jet_correction_state["uncSources_toUse"] if s != "JER"
        ]
        for unc in sources:
            # JES belongs to the output suffix, not the provider source name.
            suffix = unc if unc == "JER" else f"JES{unc}"
            for scale in ["up", "down"]:
                df = df.Define(
                    f"Jet_p4_{suffix}{scale}",
                    p4_from(unc, scale)
                )

    return df



def get_jet_correction_period(config, is_data):
    era = config.get("era")
    period = period_names[era]

    mode = get_2025_mc_jerc_mode(config, is_data)
    if mode in ("jec2024_jer2025", "2024"):
        return period_names["Run3_2024"]

    return period





def apply_jet_corrections(df, config, dataset_cfg, dataset_name, want_variations):
    # Extract configuration parameters
    is_data = dataset_cfg.get("is_data", False)

    era = config.get("era")
    period = period_names[era]

    apply_JER = config.get("apply_JER", True)
    apply_JES = config.get("apply_JES", True)
    print("apply_JER? ", apply_JER)
    print("apply_JES? ", apply_JES)
    use_regrouped = config.get("use_regrouped", True)

    # Apply jet corrections
    initialize_jet_corrections(
        period,
        is_data,
        dataset_name,
        use_regrouped
    )
    df = define_jet_p4_variations(
        df,
        want_variations,
        apply_JER,
        apply_JES,
        config.get("apply_jet_horn_mitigation", True),
    )

    return df

__all__ = [
    "initialize_jet_corrections",
    "define_jet_p4_variations",
    "get_2025_mc_jerc_mode",
    "get_jet_correction_period",
    "get_jet_correction_periods",
    "apply_jet_corrections",
]
