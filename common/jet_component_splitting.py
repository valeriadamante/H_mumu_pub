"""Reco-jet and gen-matching components used by DY control fits."""

from __future__ import annotations

from copy import deepcopy


GGF_COMPONENT_VARIABLES = {
    "ggF_0J_Hard": "m_mumu",
    "ggF_1J_Hard": "eta_signed_vs_pt_leadingjet",
    "ggF_1J_PU": "eta_signed_vs_pt_leadingjet",
    "ggF_2J_Hard": "eta_signed_vs_pt_subleadingjet",
    "ggF_2J_PU1": "eta_signed_vs_pt_subleadingjet",
    "ggF_2J_PU2": "eta_signed_vs_pt_subleadingjet",
}

VBF_COMPONENTS = ("VBF_Hard", "VBF_PU1", "VBF_PU2")
DY_JET_COMPONENTS = (*GGF_COMPONENT_VARIABLES, *VBF_COMPONENTS)
DY_COMPONENT_FILE_LABELS = {
    "ggF_0J_Hard": "DY_0J",
    "ggF_1J_Hard": "DY_1J_Hard",
    "ggF_1J_PU": "DY_1J_PU",
    "ggF_2J_Hard": "DY_2J_Hard",
    "ggF_2J_PU1": "DY_2J_PU1",
    "ggF_2J_PU2": "DY_2J_PU2",
    "VBF_Hard": "DY_2J_Hard",
    "VBF_PU1": "DY_2J_PU1",
    "VBF_PU2": "DY_2J_PU2",
}
VBF_ETA_REGIONS = ("incl", "CC", "CF", "FF")
PU_HARD_COMPONENT_STYLES = (
    ("_0J", "0J"),
    ("_1J_Hard", "1J Hard"),
    ("_1J_PU", "1J PU"),
    ("_2J_Hard", "2J Hard"),
    ("_2J_PU1", "2J PU1"),
    ("_2J_PU2", "2J PU2"),
    ("_VBF_Hard", "VBF Hard"),
    ("_VBF_PU1", "VBF PU1"),
    ("_VBF_PU2", "VBF PU2"),
)


def component_output_directory(mass_region, category, eta_region=None):
    """Return the public ROOT layout for PU/hard jet components.

    PU/hard splitting and the optional VBF eta splitting are independent:
    ggF is always flat, while VBF is nested only when an eta region was
    explicitly requested.
    """
    path = f"{mass_region}_{category}"
    return f"{path}/{eta_region}" if eta_region is not None else path


def pu_hard_component_style(sample_name, styles=None):
    """Return family, component label and configured color for components."""
    for suffix, component_label in PU_HARD_COMPONENT_STYLES:
        if sample_name.endswith(suffix):
            process_prefix = sample_name.removesuffix(suffix).rstrip("_")
            configured_families = [
                family
                for family in (styles or {})
                if family != "default" and process_prefix.startswith(family)
            ]
            family = (
                max(configured_families, key=len)
                if configured_families
                else (
                    "DY" if process_prefix.startswith("DY")
                    else "EWK" if process_prefix.startswith("EWK")
                    else process_prefix
                )
            )
            family_styles = (styles or {}).get(
                family, (styles or {}).get("default", {})
            )
            color = family_styles.get(component_label)
            return family, component_label, color
    return None


def jet_components_enabled_for_dataset(
    allowed_groups, dataset_name, process_name, *, is_signal=False
):
    """Match either process names or dataset-campaign group aliases."""
    allowed = set(allowed_groups)
    if process_name in allowed:
        return True
    if "DY_amcatnlo" in allowed and dataset_name in {
        "DYto2L_M_50_amcatnloFXFX",
        "DYto2Mu_M_50_amcatnloFXFX",
        "DYto2Tau_M_50_amcatnloFXFX",
        "DYto2E_M_50_amcatnloFXFX",
    }:
        return True
    if (
        "DY_amcatnlo_105_160" in allowed
        and dataset_name.startswith("DYto2Mu_MLL_105to160_amcatnloFXFX")
    ):
        return True
    if "EWK" in allowed and dataset_name == "EWK_2L2J_madgraph_herwig":
        return True
    if (
        "EWK_105_160" in allowed
        and dataset_name.startswith("EWK_2Mu2J_MLL_105to160_")
    ):
        return True
    return "signals" in allowed and is_signal


def vbf_eta_region_expressions(base_expression):
    """Split a VBF selection at |eta|=2.5 using the selected VBF pair."""
    eta1 = "abs(SelectedJet_eta{jet_suff}[VBFJetIdx_1{jet_suff}])"
    eta2 = "abs(SelectedJet_eta{jet_suff}[VBFJetIdx_2{jet_suff}])"
    central1 = f"{eta1} < 2.5"
    central2 = f"{eta2} < 2.5"
    return {
        "incl": base_expression,
        "CC": f"({base_expression}) && {central1} && {central2}",
        "CF": (
            f"({base_expression}) && (({central1} && !({central2})) || "
            f"(!({central1}) && {central2}))"
        ),
        "FF": f"({base_expression}) && !({central1}) && !({central2})",
    }


def expanded_jet_component_categories(
    include_vbf_eta_regions=False, requested_categories=None
):
    """Internal staging categories needed for split component ROOT files."""
    requested_sequence = tuple(requested_categories or ("ggF", "VBF"))
    requested = set(requested_sequence)
    eta_regions = VBF_ETA_REGIONS if include_vbf_eta_regions else ("incl",)
    # Categories unrelated to the ggF/VBF split (for example ``baseline``)
    # remain ordinary inclusive outputs alongside the component staging dirs.
    categories = [
        category
        for category in requested_sequence
        if category not in {"ggF", "VBF"}
    ]
    if "ggF" in requested:
        categories.append("DY_inclusive_ggF")
        categories.extend(GGF_COMPONENT_VARIABLES)
    if "VBF" in requested:
        categories.extend(f"DY_inclusive_VBF_{region}" for region in eta_regions)
        for component in VBF_COMPONENTS:
            categories.extend(f"{component}_{region}" for region in eta_regions)
    return tuple(categories)


def add_vbf_eta_region_categories(selection_config):
    """Add inclusive/CC/CF/FF VBF categories independently of DY splitting."""
    config = deepcopy(selection_config)
    categories = config.setdefault("categories", {})
    for region, expression in vbf_eta_region_expressions(
        "VBF{tot_suff}"
    ).items():
        categories[f"VBF_eta_{region}"] = {
            "expression": expression,
            "store": True,
        }
    return config


def add_jet_component_categories(
    selection_config, include_vbf_eta_regions=False, requested_categories=None
):
    """Return a config with mutually exclusive reco/gen-matching categories."""
    config = deepcopy(selection_config)
    categories = config.setdefault("categories", {})
    requested = set(requested_categories or ("ggF", "VBF"))
    eta_regions = VBF_ETA_REGIONS if include_vbf_eta_regions else ("incl",)
    components = {
        "ggF_0J_Hard": "ggF{tot_suff} && N_SelectedJets{jet_suff} == 0",
        "ggF_1J_Hard": (
            "ggF{tot_suff} && N_SelectedJets{jet_suff} == 1 "
            "&& N_PU_FirstTwoJets{jet_suff} == 0"
        ),
        "ggF_1J_PU": (
            "ggF{tot_suff} && N_SelectedJets{jet_suff} == 1 "
            "&& N_PU_FirstTwoJets{jet_suff} == 1"
        ),
        "ggF_2J_Hard": (
            "ggF{tot_suff} && N_SelectedJets{jet_suff} >= 2 "
            "&& N_PU_FirstTwoJets{jet_suff} == 0"
        ),
        "ggF_2J_PU1": (
            "ggF{tot_suff} && N_SelectedJets{jet_suff} >= 2 "
            "&& N_PU_FirstTwoJets{jet_suff} == 1"
        ),
        "ggF_2J_PU2": (
            "ggF{tot_suff} && N_SelectedJets{jet_suff} >= 2 "
            "&& N_PU_FirstTwoJets{jet_suff} == 2"
        ),
        "VBF_Hard": "VBF{tot_suff} && N_PU_VBFJets{jet_suff} == 0",
        "VBF_PU1": "VBF{tot_suff} && N_PU_VBFJets{jet_suff} == 1",
        "VBF_PU2": "VBF{tot_suff} && N_PU_VBFJets{jet_suff} == 2",
    }
    staging = {}
    if "ggF" in requested:
        staging["DY_inclusive_ggF"] = "ggF{tot_suff}"
    if "VBF" in requested:
        staging.update({
            f"DY_inclusive_VBF_{region}": expression
            for region, expression in vbf_eta_region_expressions(
                "VBF{tot_suff}"
            ).items()
            if region in eta_regions
        })
    for name, expression in components.items():
        if name in VBF_COMPONENTS:
            if "VBF" not in requested:
                continue
            staging.update(
                {
                    f"{name}_{region}": region_expression
                    for region, region_expression in vbf_eta_region_expressions(
                        expression
                    ).items()
                    if region in eta_regions
                }
            )
        else:
            if "ggF" not in requested:
                continue
            staging[name] = expression
    for name, expression in staging.items():
        categories[name] = {"expression": expression, "store": True}
    return config


def define_jet_gen_matching(df, selection_suffixes):
    """Define PU-jet counts for the first two reco jets and the VBF pair.

    A reconstructed jet is classified as hard when ``Jet_genJetIdx >= 0`` and
    as PU when no generator-level jet is matched (index < 0).
    """
    columns = {str(column) for column in df.GetColumnNames()}
    for suffix in selection_suffixes:
        selected_index = f"SelectedJet_idx{suffix}"
        selected_gen_index = f"SelectedJet_genJetIdx{suffix}"
        n_selected = f"N_SelectedJets{suffix}"
        if selected_index not in columns:
            raise RuntimeError(
                f"Jet component splitting requires column {selected_index}"
            )
        if selected_gen_index not in columns:
            if "Jet_genJetIdx" not in columns:
                raise RuntimeError(
                    "Jet component splitting requires SelectedJet_genJetIdx "
                    "or the raw MC column Jet_genJetIdx. Regenerate this DY "
                    "skim with the current analysis/jets.py."
                )
            df = df.Define(
                selected_gen_index,
                f"ROOT::VecOps::Take(Jet_genJetIdx, {selected_index})",
            )
            columns.add(selected_gen_index)

        first_two_count = f"N_PU_FirstTwoJets{suffix}"
        if first_two_count not in columns:
            df = df.Define(
                first_two_count,
                (
                    f"int(({n_selected} > 0 && {selected_gen_index}[0] < 0) + "
                    f"({n_selected} > 1 && {selected_gen_index}[1] < 0))"
                ),
            )
            columns.add(first_two_count)

        # Generic, process-independent 0/1/2J flags. Here J denotes the
        # number of the first two selected reco jets without a GenJet match.
        for count in (0, 1, 2):
            flag = f"RecoGenJetMatch_{count}J{suffix}"
            if flag not in columns:
                df = df.Define(flag, f"{first_two_count} == {count}")
                columns.add(flag)

        vbf_count = f"N_PU_VBFJets{suffix}"
        if vbf_count not in columns:
            df = df.Define(
                vbf_count,
                (
                    f"HasVBF{suffix} ? int("
                    f"({selected_gen_index}[VBFJetIdx_1{suffix}] < 0) + "
                    f"({selected_gen_index}[VBFJetIdx_2{suffix}] < 0)) : -1"
                ),
            )
            columns.add(vbf_count)
    return df


def variable_for_component(category, requested_variables):
    """Keep requested observables and add the component-specific 2D one."""
    variables = list(requested_variables)
    component_variable = None
    if category in GGF_COMPONENT_VARIABLES:
        component_variable = GGF_COMPONENT_VARIABLES[category]
    for component, variable in GGF_COMPONENT_VARIABLES.items():
        if category.startswith(f"{component}_"):
            component_variable = variable
            break
    if component_variable is not None and component_variable not in variables:
        variables.append(component_variable)
    return tuple(variables)
