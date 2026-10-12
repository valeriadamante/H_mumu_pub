
#!/usr/bin/env python3
from pathlib import Path
import os
from .general import pog_folder_names,period_names

__all__ = ["apply_pu_weights", "_resolve_pu_file"]

_initialized = set()

# minimal golden json mapping (used to pick the correction key inside the PU JSON)
golden_json_dict = {
    "2026_Summer24": "Collisions25_goldenJSON", # TMP PATCH AS IT IS NOT AVAILABLE FOR 2026

    "2025_Winter25": "Collisions25_goldenJSON",
    "2025_Summer24": "Collisions25_goldenJSON",
    "2024_Summer24": "Collisions24_BCDEFGHI_goldenJSON",
    "2023_Summer23BPix": "Collisions2023_369803_370790_eraD_GoldenJson",
    "2023_Summer23": "Collisions2023_366403_369802_eraBC_GoldenJson",
    "2022_Summer22EE": "Collisions2022_359022_362760_eraEFG_GoldenJson",
    "2022_Summer22": "Collisions2022_355100_357900_eraBCD_GoldenJson",
    "2018_UL": "Collisions18_UltraLegacy_goldenJSON",
    "2017_UL": "Collisions17_UltraLegacy_goldenJSON",
    "2016preVFP_UL": "Collisions16_UltraLegacy_goldenJSON",
    "2016postVFP_UL": "Collisions16_UltraLegacy_goldenJSON",
}

def apply_pu_weights(
    df,
    json_dict_to_store,
    config,
    pileup_column,
    want_variations,
):


    import ROOT

    json_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/{folder}/latest/puWeights{suffix}.json.gz"
    era = config.get("era")
    period_unc = period_names[era]
    folder_name = pog_folder_names["LUM"][period_unc]
    suffix = "_BCDEFGHI" if period_unc == "2024_Summer24" else ""  # tmp patch ? 
    if period_unc=="2025_Summer24" or period_unc=="2026_Summer24": # TMP PATCH AS IT IS NOT AVAILABLE FOR 2026
        suffix="_2025pp_Golden_Summer24_25ns_69200ub"
    pu_path = json_path.format(
        folder=folder_name, suffix=suffix
    )
    import correctionlib

    cs = correctionlib.CorrectionSet.from_file(str(pu_path))

    json_key = golden_json_dict.get(period_unc)
    corr = cs[json_key]
    # sample pileup values (0..max_pu) to create lookup arrays
    max_pu = 200
    central_arr = []
    up_arr = []
    down_arr = []
    for i in range(max_pu + 1):
        # correctionlib.Correction.evaluate expects separate positional inputs
        # (one per input in the correction), not a single list argument.
        central_arr.append(float(corr.evaluate(float(i), "nominal")))
        up_arr.append(float(corr.evaluate(float(i), "up")))
        down_arr.append(float(corr.evaluate(float(i), "down")))

    # use stable, single symbol names for the embedded arrays/accessors
    func_c = "pu_weight_central"
    func_u = "pu_weight_up"
    func_d = "pu_weight_down"

    resolved_path = str(Path(pu_path).resolve())
    if resolved_path not in _initialized:
        # create C++ source: static arrays and accessor functions
        cpp_lines = [
            "#include <vector>",
            "#include <cmath>",
            "using namespace std;",
            f"static const vector<double> pu_central = {{ {', '.join(map(str, central_arr))} }};",
            f"static const vector<double> pu_up = {{ {', '.join(map(str, up_arr))} }};",
            f"static const vector<double> pu_down = {{ {', '.join(map(str, down_arr))} }};",
            "static inline double _get_from_vec(const vector<double>& v, double x) {",
            "  int idx = (int)floor(x);",
            "  if (idx < 0) return v.front();",
            "  if (idx >= (int)v.size()) return v.back();",
            "  return v[idx];",
            "}",
            f"static inline double {func_c}(double x) {{ return _get_from_vec(pu_central, x); }}",
            f"static inline double {func_u}(double x) {{ return _get_from_vec(pu_up, x); }}",
            f"static inline double {func_d}(double x) {{ return _get_from_vec(pu_down, x); }}",
        ]

        cpp_code = "\n".join(cpp_lines)
        ROOT.gInterpreter.Declare(cpp_code)

        # record initialization to avoid redeclaring the same arrays
        _initialized.add(resolved_path)

    branches = []
    scales = ["Central"]
    if want_variations:
        scales += ["up", "down"]

    # json_dict_to_store = {"pu":{},"pu_up":{},"pu_down":{}}
    for scale in scales:
        branch_name = f"weight_pu"
        if scale != "Central": branch_name += f"_{scale}"
        if scale == "Central":
            expr = f"{func_c}({pileup_column})"
        elif scale == "up":
            expr = f"{func_u}({pileup_column})"
        else:
            expr = f"{func_d}({pileup_column})"
        df = df.Define(branch_name, expr)
        branches.append(branch_name)
    json_dict_to_store['pu'] = {}
    json_dict_to_store['pu'][f"total"] = {}
    json_dict_to_store['pu'][f"total"]['selection']= "return true;"
    json_dict_to_store['pu'][f"total"]['value']= df.Define("signed_weight_pu","weight_pu*weight_gen_sign").Sum("signed_weight_pu")
    json_dict_to_store['pu'][f"total"]['value_unsigned']= df.Sum("weight_pu")
    if want_variations:
        for scale in ["up", "down"]:
            json_dict_to_store[f"pu_{scale}"] = {}
            json_dict_to_store[f"pu_{scale}"][f"total"] = {}
            json_dict_to_store[f"pu_{scale}"][f"total"]['selection']= "return true;"
            json_dict_to_store[f"pu_{scale}"][f"total"]['value_unsigned']= df.Sum(f"weight_pu_{scale}")
            json_dict_to_store[f"pu_{scale}"][f"total"]['value']= df.Define(f"signed_weight_pu_{scale}",f"weight_pu_{scale}*weight_gen_sign").Sum(f"signed_weight_pu_{scale}")
    return df, branches, json_dict_to_store

