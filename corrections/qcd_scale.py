"""QCD renormalization/factorization scale variations from NanoAOD LHE weights."""

import re

import ROOT


_QCD_SCALE_HELPER = r"""
#include <ROOT/RVec.hxx>

namespace qcd_scale {
inline float weightAt(
    const ROOT::VecOps::RVec<float>& weights,
    const unsigned int index
) {
    return index < weights.size() ? static_cast<float>(weights[index]) : 1.f;
}
}  // namespace qcd_scale
"""

if not ROOT.gInterpreter.Declare(_QCD_SCALE_HELPER):
    raise RuntimeError("Failed to declare the QCD scale weight helper")


def get_qcd_scale_points(config):
    """Return configured six-point scale variations."""
    return config.get(
        "points",
        [
            {"name": "muR0p5_muF0p5", "index": 0},
            {"name": "muR0p5_muF1", "index": 1},
            {"name": "muR1_muF0p5", "index": 3},
            {"name": "muR1_muF2", "index": 5},
            {"name": "muR2_muF1", "index": 7},
            {"name": "muR2_muF2", "index": 8},
        ],
    )


def _df_for_selection(df, selection):
    selection = str(selection or "").strip()
    if not selection or selection.lower() == "return true;":
        return df
    return df.Filter(selection)


def _safe_column_fragment(value):
    return re.sub(r"[^A-Za-z0-9_]", "_", str(value))


def define_qcd_scale_sum_columns(df, config, json_dict):
    """Book exact selected-event denominator sums for every scale variation."""
    config = config or {}
    if not config.get("enabled", True):
        return df, json_dict

    branch = config.get("branch", "LHEScaleWeight")
    available_columns = {str(column) for column in df.GetColumnNames()}
    if branch not in available_columns:
        print(f"[WARNING] QCD scale branch '{branch}' is missing.")
        return df, json_dict

    gen_entries = dict(json_dict.get("gen", {}))
    pu_entries = dict(json_dict.get("pu", {}))

    for point in get_qcd_scale_points(config):
        name = point["name"]
        index = int(point["index"])
        weight_column = f"weight_qcd_scale_sum__{_safe_column_fragment(name)}"
        df = df.Define(
            weight_column,
            f"genWeight * qcd_scale::weightAt({branch}, {index}u)",
        )

        gen_node = f"gen_qcdScale_{name}"
        json_dict[gen_node] = {}
        for entry_name, entry in gen_entries.items():
            selection = entry.get("selection", "return true;")
            selected_df = _df_for_selection(df, selection)
            json_dict[gen_node][entry_name] = {
                "selection": selection,
                "value": selected_df.Sum(weight_column),
            }

        if "weight_pu" not in available_columns:
            continue

        pu_node = f"pu_qcdScale_{name}"
        json_dict[pu_node] = {}
        for entry_name, entry in pu_entries.items():
            selection = entry.get("selection", "return true;")
            selected_df = _df_for_selection(df, selection)
            safe_name = _safe_column_fragment(f"{name}_{entry_name}")
            signed_column = f"signed_weight_pu_qcd_scale_sum__{safe_name}"
            unsigned_column = f"unsigned_weight_pu_qcd_scale_sum__{safe_name}"
            json_dict[pu_node][entry_name] = {
                "selection": selection,
                "value": selected_df.Define(
                    signed_column,
                    f"weight_pu * weight_gen_sign * qcd_scale::weightAt({branch}, {index}u)",
                ).Sum(signed_column),
                "value_unsigned": selected_df.Define(
                    unsigned_column,
                    f"weight_pu * qcd_scale::weightAt({branch}, {index}u)",
                ).Sum(unsigned_column),
            }
    return df, json_dict


__all__ = ["define_qcd_scale_sum_columns", "get_qcd_scale_points"]
