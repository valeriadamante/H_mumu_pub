import ROOT
import sys
import os

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])

from common.utilities import _column_names, _define_if_missing


def _declare_muon_helpers():
    ROOT.gInterpreter.Declare(
        """
        #ifndef NEW_SKIM_MUON_ANALYSIS_HELPERS
        #define NEW_SKIM_MUON_ANALYSIS_HELPERS
        using RVecF = ROOT::VecOps::RVec<float>;

        template <typename P4>
        ROOT::VecOps::RVec<P4> Muon_p4_sel(const ROOT::VecOps::RVec<P4>& nano,
                                         const ROOT::VecOps::RVec<P4>& bsc,
                                         const RVecF& chi2) {
            ROOT::VecOps::RVec<P4> out(nano.size());
            for (size_t i = 0; i < out.size(); ++i) {
                out[i] = chi2[i] < 30 ? bsc[i] : nano[i];
            }
            return out;
        }

        RVecF Muon_pt_err_sel(const RVecF& Muon_nano_pt_err, const RVecF& Muon_bsc_pt_err, const RVecF& Muon_bsc_chi2) {
            RVecF out(Muon_nano_pt_err.size());
            for (size_t i = 0; i < out.size(); ++i) {
                out[i] = (Muon_bsc_chi2[i] < 30) ? Muon_bsc_pt_err[i] : Muon_nano_pt_err[i];
            }
            return out;
        }

        RVecF Muon_pt_sel(const RVecF& Muon_nano_pt, const RVecF& Muon_bsc_pt, const RVecF& Muon_bsc_chi2) {
            RVecF out(Muon_nano_pt.size());
            for (size_t i = 0; i < out.size(); ++i) {
                out[i] = (Muon_bsc_chi2[i] < 30) ? Muon_bsc_pt[i] : Muon_nano_pt[i];
            }
            return out;
        }
        #endif
        """
    )

def GetPtConfigurations(want_variations):
    configs = {
        "Muon_pt_raw_noCorr": ["Muon_pt", "Muon_bsConstrainedPt"],
        "Muon_pt_raw_corr": ["Muon_pt_nano_corr", "Muon_pt_bsc_corr"],
        "Muon_pt_raw_scale": ["Muon_pt_nano_scale", "Muon_pt_bsc_scale"],
        "Muon_pt_FSR_noCorr": ["Muon_pt_nano_FSR", "Muon_pt_bsc_FSR"],
        "Muon_pt_FSR_corr": ["Muon_pt_nano_corr_FSR", "Muon_pt_bsc_corr_FSR"],
        "Muon_pt_FSR_scale": ["Muon_pt_nano_scale_FSR", "Muon_pt_bsc_scale_FSR"]
    }
    if want_variations:
        configs.update({
            "Muon_pt_raw_scale_up": ["Muon_pt_nano_scale_up", "Muon_pt_bsc_scale_up"],
            "Muon_pt_raw_scale_down": ["Muon_pt_nano_scale_down", "Muon_pt_bsc_scale_down"],
            "Muon_pt_raw_res_up": ["Muon_pt_nano_res_up", "Muon_pt_bsc_res_up"],
            "Muon_pt_raw_res_down": ["Muon_pt_nano_res_down", "Muon_pt_bsc_res_down"],
            "Muon_pt_FSR_scale_up": ["Muon_pt_nano_FSR_scale_up", "Muon_pt_bsc_FSR_scale_up"],
            "Muon_pt_FSR_scale_down": ["Muon_pt_nano_FSR_scale_down", "Muon_pt_bsc_FSR_scale_down"],
            "Muon_pt_FSR_res_up": ["Muon_pt_nano_FSR_res_up", "Muon_pt_bsc_FSR_res_up"],
            "Muon_pt_FSR_res_down": ["Muon_pt_nano_FSR_res_down", "Muon_pt_bsc_FSR_res_down"],
        })
    err_configs = {
        "Muon_pt_err": ["Muon_ptErr", "Muon_bsConstrainedPtErr"],
    }
    return configs,err_configs

def DefineMuonPtAndP4(df, want_variations):
    _declare_muon_helpers()
    configs,err_configs = GetPtConfigurations(want_variations)
    cols = _column_names(df)
    for name_pt, (nano, bsc) in configs.items():
        df = _define_if_missing(df,name_pt,f"Muon_pt_sel({nano}, {bsc}, Muon_bsConstrainedChi2)")
        if "_FSR_" in name_pt:
            # Keep the photon-added direction and mass as well as its pT.
            p4_expr = f"Muon_p4_sel({nano.replace('Muon_pt_', 'Muon_p4_')}, {bsc.replace('Muon_pt_', 'Muon_p4_')}, Muon_bsConstrainedChi2)"
        else:
            p4_expr = f"GetP4({name_pt}, Muon_eta, Muon_phi, Muon_mass)"
        df = _define_if_missing(df, name_pt.replace("pt", "p4"), p4_expr)
    for name_err, (nano_err, bsc_err) in err_configs.items():
         df = _define_if_missing(df,name_err,f"Muon_pt_err_sel({nano_err}, {bsc_err}, Muon_bsConstrainedChi2)")
    return df

def ApplyMuonTriggerMatching(df, trigger_config, apply_filter, want_variations, syst_cfg):
    cols = _column_names(df)
    if "TrigObj_pt" in cols:
        df = _define_if_missing(df, "TrigObj_idx", "CreateIndexes(TrigObj_pt.size())")
        df = _define_if_missing(df, "TrigObj_mass", "RVecF(TrigObj_pt.size(), 0.f)")
        df = _define_if_missing(df, "TrigObj_p4","GetP4(TrigObj_pt, TrigObj_eta, TrigObj_phi, TrigObj_mass, TrigObj_idx)")

    syst_suffixes = [""]
    if want_variations:
        scales = syst_cfg.get("scales", ["up", "down"])
        syst_suffixes.extend(
            syst_cfg["systematics"][syst]["muon_suffix"].format(scale=scale)
            for syst in ("MuonScale", "MuonRes")
            for scale in scales
        )

    filters = []
    cols_to_save = []
    for path, config in trigger_config.items():
        path_name = config["path"][0]
        cols_to_save.append(path_name)
        leg = config["legs"][0]
        online = leg["online_obj"]["cut"]
        df = _define_if_missing(df, f"TrigObj_passOnlineCut_{path}", online)
        cols_to_save.append(f"TrigObj_passOnlineCut_{path}")

        for suff in syst_suffixes:
            # Trigger matching uses the uncorrected pT (no ScaRe, no FSR) for
            # the nominal and for every muon variation, so the shifted
            # matching columns are identical to the nominal ones.
            pt = "pt_raw_noCorr"
            muon_p4 = "Muon_p4_raw_noCorr"
            offline_col = f"Muon_passOfflineCut_{path}{suff}"
            matching_col = f"Muon_TriggerMatchingIdx_{path}{suff}"
            evt = f"Event_HasTriggerMatching_{path}{suff}"

            offline = leg["offline_obj"]["cut"].format(obj="Muon", pt=pt)
            df = _define_if_missing(df, offline_col, offline)
            df = _define_if_missing(
                df,
                matching_col,
                f"FindMatching({offline_col}, TrigObj_passOnlineCut_{path}, "
                f"{muon_p4}, TrigObj_p4, 0.4)",
            )
            df = _define_if_missing(df, evt, f"{path_name} && Any({matching_col} > -1)")

            filters.append(evt)
            cols_to_save.extend([offline_col, matching_col, evt])
    if apply_filter:
        df = df.Filter(" || ".join(filters), "Trigger matching for " + "__".join(trigger_config.keys()))
    return df, cols_to_save

def ProcessMuonVariables(df,muon_columns,default_suffix,trigger_config,want_variations,pt_min,lower_mass_cut,upper_mass_cut,syst_cfg,apply_trigger_filter=True):
    cols = _column_names(df)
    selection_pt = [f"Muon_pt_{default_suffix}"]
    syst_suffixes = [""]
    if want_variations:
        scales = syst_cfg.get('scales',['up','down'])
        syst_suffixes.extend([syst_cfg['systematics']['MuonScale']['muon_suffix'].format(scale=scale) for scale in scales])
        syst_suffixes.extend([syst_cfg['systematics']['MuonRes']['muon_suffix'].format(scale=scale) for scale in scales])
        selection_pt += ["Muon_pt{syst}" for syst in syst_suffixes]
    new_cols = []

    def track(df, name, expr):
        if name in _column_names(df): return df
        if "p4" not in name: new_cols.append(name)
        return df.Define(name, expr)

    event_filters = []
    mass_filters = []
    pair_trigger_filters = []
    pt_branches,err_pt_branches = GetPtConfigurations(want_variations)
    for suff in syst_suffixes:
        is_nominal = (suff == "")
        pt=f"Muon_pt{suff}"
        p4=f"Muon_p4{suff}"
        if is_nominal:
            pt = "Muon_pt_"+default_suffix
            p4 = "Muon_p4_"+default_suffix
        df = df.Define(f"good_muons{suff}",f"Muon_pt_raw_noCorr > {pt_min} && abs(Muon_eta) < 2.4 && Muon_mediumId && Muon_pfIsoId >= 2")
        df = df.Define(f"good_idx{suff}",f"ROOT::VecOps::Nonzero(good_muons{suff})")
        df = df.Define(f"sorted_idx{suff}",f"Reverse(Take(good_idx{suff}, Argsort(Take({pt}, good_idx{suff}))))")
        df = track(df, f"mu1_idx{suff}", f"sorted_idx{suff}.size()>0 ? (int)sorted_idx{suff}[0] : -1")
        df = track(df, f"mu2_idx{suff}", f"sorted_idx{suff}.size()>1 ? (int)sorted_idx{suff}[1] : -1")
        event_filters.append(f"sorted_idx{suff}.size() == 2")
        for i in [1, 2]:
            idx = f"mu{i}_idx{suff}"
            df = track(df, f"mu{i}_p4{suff}", f"{idx}>=0 ? {p4}.at({idx}) : ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(0,0,0,0)")
            df = track(df, f"mu{i}_pt{suff}", f"{idx}>=0 ? {pt}[{idx}] : -999.f")
            df = track(df, f"mu{i}_phi{suff}", f"{idx}>=0 ? {p4}[{idx}].Phi() : -999.f")
            df = track(df, f"mu{i}_eta{suff}", f"{idx}>=0 ? {p4}[{idx}].Eta() : -999.f")
            df = track(df, f"mu{i}_mass{suff}", f"{idx}>=0 ? {p4}[{idx}].M() : -999.f")
            for muon_col in muon_columns:
                suffix_clean = "_".join(c for c in muon_col.split("_")[1:])
                df = track(df, f"mu{i}_{suffix_clean}{suff}", f"{idx}>=0 ? {muon_col}[{idx}] : -999.f")
            # Use NanoAOD's existing reco-to-gen association, including nonprompt muons.
            gen_match = "false"
            if "Muon_genPartIdx" in cols:
                gen_idx = f"Muon_genPartIdx[{idx}]"
                gen_match = f"{idx} >= 0 && {gen_idx} >= 0"
                if "GenPart_pdgId" in cols:
                    gen_match += f" && {gen_idx} < int(GenPart_pdgId.size()) && abs(GenPart_pdgId[{gen_idx}]) == 13"
            elif "Muon_genPartFlav" in cols:
                gen_match = f"{idx} >= 0 && Muon_genPartFlav[{idx}] > 0"
            df = track(df, f"mu{i}_GenMatched{suff}", gen_match)
            for path in trigger_config.keys():
                df = track(df, f"mu{i}_HasTriggerMatching_{path}{suff}", f"{idx} >= 0 ? (Muon_TriggerMatchingIdx_{path}{suff}[{idx}] >= 0) : false")

        for path in trigger_config:
            pair_trigger_filters.append(
                f"(Event_HasTriggerMatching_{path}{suff} && "
                f"(mu1_HasTriggerMatching_{path}{suff} || mu2_HasTriggerMatching_{path}{suff}))"
            )

        p4 = f"(mu1_p4{suff} + mu2_p4{suff})"
        df = track(df, f"m_mumu{suff}", f"{p4}.M()")
        mass_filters.append(f"m_mumu{suff} > {lower_mass_cut} && m_mumu{suff} < {upper_mass_cut}")

    df = df.Filter(" && ".join(event_filters), "Exactly 2 muons")
    if apply_trigger_filter and pair_trigger_filters:
        df = df.Filter(" || ".join(pair_trigger_filters), "Selected dimuon trigger matching")
    df = df.Filter(" && ".join(mass_filters), "dimuon mass cut")

    # These observables are evaluated for the nominally selected muons. The
    # FSR scale/res branches used to build shifted categories are already
    # stored above with their own shifted muon indices.

    idx1 = "mu1_idx"
    idx2 = "mu2_idx"
    # Store every correction stage for the same nominal muon indices.
    # Shifted selections already have their own eta/phi branches above.
    for pt_branch in pt_branches:
        stage = pt_branch.removeprefix("Muon_pt_")
        p4_branch = pt_branch.replace("Muon_pt_", "Muon_p4_", 1)
        for i in (1, 2):
            idx = f"mu{i}_idx"
            for coordinate, accessor in (("eta", "Eta"), ("phi", "Phi")):
                df = track(df, f"mu{i}_{coordinate}_{stage}",
                           f"{idx}>=0 ? {p4_branch}[{idx}].{accessor}() : -999.f")
    for centr_br in list(pt_branches.keys()) + list(err_pt_branches.keys()):
        suffix_br = "_".join(c for c in centr_br.split("_")[1:])
        df = track(df, f"mu1_{suffix_br}", f"{idx1}>=0 ? {centr_br}[{idx1}] : -999.f")
        df = track(df, f"mu2_{suffix_br}", f"{idx2}>=0 ? {centr_br}[{idx2}] : -999.f")
    return df, new_cols


def ProcessExtraMuonVariables(df,muon_columns,default_suffix,trigger_config,want_variations,pt_min):
    pt_branches,err_pt_branches = GetPtConfigurations(want_variations)
    cols = _column_names(df)
    new_cols = []
    pt_list = [f"Muon_pt_{default_suffix}"]
    def track(df, name, expr):
        if name in _column_names(df): return df
        new_cols.append(name)
        return df.Define(name, expr)
    for pt in pt_list:
        if pt not in cols:
            print(f"pt column {pt} not found, skipping extra muon variables")
            continue
        df = _define_if_missing(df, "Muon_idx", f"CreateIndexes({pt}.size())")
        df = df.Define( "extra_good_muons", f"{pt} > {pt_min} && abs(Muon_eta) < 2.4 && Muon_looseId && Muon_pfIsoId >= 2 && " f"Muon_idx != mu1_idx && " f"Muon_idx != mu2_idx")
        df = df.Define("extra_good_idx", "ROOT::VecOps::Nonzero(extra_good_muons)")
        df = df.Define("extra_sorted_idx", f"Reverse(Take(extra_good_idx, Argsort(Take({pt}, extra_good_idx))))")
        df = track(df, "extraMuon_idx", "extra_sorted_idx")
        df = track(df, "n_extraMuon", "int(extraMuon_idx.size())")
        df = track(df, "extraMuon_pt", f"Take({pt}, extraMuon_idx)")
        df = track(df, "extraMuon_eta", "Take(Muon_eta, extraMuon_idx)")
        df = track(df, "extraMuon_phi", "Take(Muon_phi, extraMuon_idx)")
        df = track(df, "extraMuon_charge", "Take(Muon_charge, extraMuon_idx)")
        for col in muon_columns+list(pt_branches.keys()) + list(err_pt_branches.keys()):
            suffix_clean = "_".join(c for c in col.split("_")[1:])
            if f"extraMuon_{suffix_clean}" not in _column_names(df):
                df = track(df, f"extraMuon_{suffix_clean}", f"Take({col}, extraMuon_idx)")
    return df, new_cols

def ApplyElectronVeto(df):
    df = df.Define("Electron_p4", "GetP4(Electron_pt, Electron_eta, Electron_phi, Electron_mass)")
    df = _define_if_missing(df, "veto_electrons","Electron_pt > 20 && abs(Electron_eta) < 2.5 && Electron_mvaIso_WP90")
    return df.Filter("ROOT::VecOps::Nonzero(veto_electrons).size() == 0", "No extra electrons")

def DefineMuonSelection(df,sel_config,want_variations,syst_cfg):
    sel_dict = sel_config.get("muons_selection", {})
    vars_to_store = []
    syst_suffixes = [""]
    print(sel_dict)
    if want_variations:
        scales = syst_cfg.get('scales',['up','down'])
        syst_suffixes.extend([syst_cfg['systematics']['MuonScale']['muon_suffix'].format(scale=scale) for scale in scales])
        syst_suffixes.extend([syst_cfg['systematics']['MuonRes']['muon_suffix'].format(scale=scale) for scale in scales])
    for suff in syst_suffixes:
        for sel_name, sel_subdict in sel_dict.items():
            sel_str = sel_subdict["expression"]
            full_name = f"{sel_name}{suff}"
            df = df.Define(full_name,sel_str.format(mu_suff=suff))
            if sel_subdict.get("store", False):
                vars_to_store.append(full_name)
    return df, vars_to_store
