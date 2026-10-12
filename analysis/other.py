import ROOT
import os

def applyMETFlags(df, config, is_data):
    MET_flags = config.get("MET_flags", [])
    badMET_flag_runs = config.get("badMET_flag_runs", [])
    if badMET_flag_runs:
        df = applyBadMETfilter(df, badMET_flag_runs, is_data)
    MET_flags_string = " && ".join(MET_flags)
    return df.Filter(MET_flags_string, "MET filters")

def applyBadMETfilter(df, badMET_flag_runs, is_data):
    if not is_data:
        return df
    else:
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
        df = df.Define(
            f"Flag_badMET_calib",
            f""" !( PuppiMET_pt>100 &&
                                                Any(Jet_pt > 50
                                                && Jet_eta >= -0.5 && Jet_eta <= -0.1
                                                && Jet_phi >= -2.1 && Jet_phi <= -1.8
                                                && abs(PuppiMET_phi - Jet_phi) > 2.9
                                                && (Jet_neEmEF > 0.9 || Jet_chEmEF > 0.9)
                                                ) )""",
        )

        df = df.Redefine(
            f"Flag_ecalBadCalibFilter",
            f" ( run >= {badMET_flag_runs[0]} && run <= {badMET_flag_runs[1]} ) ? Flag_badMET_calib : Flag_ecalBadCalibFilter",
        )
        return df


import ROOT

def _column_names(df):
    return {str(col) for col in df.GetColumnNames()}


def DefineCategories(df, sel_config, is_data, want_variations, syst_cfg):
    categories_dict = {}
    sections_to_load = ["masses_regions", "categories"]
    for section in sections_to_load:
        categories_dict.update(sel_config.get(section, {}))
    vars_to_store = []
    tot_suffixes = [""]
    mu_suffixes = [""]
    jet_suffixes = [""]
    if want_variations and syst_cfg:
        scales = syst_cfg.get('scales', ['up', 'down'])
        for tot_name, syst_subdict in syst_cfg.get("systematics", {}).items():
            if tot_name == "Central":
                continue
            for scale in scales:
                tot_suffixes.append(f"_{tot_name}{scale.capitalize()}")
                mu_suffixes.append(syst_subdict.get('muon_suffix', '').format(scale=scale))
                jet_suffixes.append(syst_subdict.get('jet_suffix', '').format(scale=scale))

    defined_columns = set(_column_names(df))
    for cat_name, cat_content in categories_dict.items():
        if isinstance(cat_content, dict):
            base_expression = cat_content.get("expression", "")
            should_store = cat_content.get("store", True)
        else:
            base_expression = cat_content
            should_store = True
        if not base_expression:
            print("no base expression" )
            continue
        for tot_suff, mu_suff, jet_suff in zip(tot_suffixes, mu_suffixes, jet_suffixes):
            formatted_expr = base_expression.format(tot_suff=tot_suff,mu_suff=mu_suff,jet_suff=jet_suff)
            final_column_name = f"{cat_name}{tot_suff}"
            # print(f"defining {final_column_name} with expression: {formatted_expr}")
            df = df.Define(final_column_name, formatted_expr)
            if should_store:
                vars_to_store.append(final_column_name)
    return df, vars_to_store