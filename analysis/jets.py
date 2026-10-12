from common.jet_variation_suffixes import jet_variation_suffixes

import ROOT
import sys

def _column_names(df):
    return {str(col) for col in df.GetColumnNames()}


def _vec_pt(df, p4, out,new_cols):
    new_cols.append(out)
    if out not in df.GetColumnNames():
        return df.Define(out, f"""
            ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
            for (size_t i = 0; i < {p4}.size(); ++i)
                v[i] = {p4}[i].Pt();
            return v;
        """)
    return df.Redefine(out, f"""
        ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
        for (size_t i = 0; i < {p4}.size(); ++i)
            v[i] = {p4}[i].Pt();
        return v;
    """)

def _vec_mass(df, p4, out,new_cols):
    new_cols.append(out)
    if out not in df.GetColumnNames():
        return df.Define(out, f"""
            ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
            for (size_t i = 0; i < {p4}.size(); ++i)
                v[i] = {p4}[i].M();
            return v;
        """)
    return df.Redefine(out, f"""
        ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
        for (size_t i = 0; i < {p4}.size(); ++i)
            v[i] = {p4}[i].M();
        return v;
    """)

def _vec_eta(df, p4, out,new_cols):
    new_cols.append(out)
    if out not in df.GetColumnNames():
        return df.Define(out, f"""
            ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
            for (size_t i = 0; i < {p4}.size(); ++i)
                v[i] = {p4}[i].Eta();
            return v;
        """)
    return df.Redefine(out, f"""
        ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
        for (size_t i = 0; i < {p4}.size(); ++i)
            v[i] = {p4}[i].Eta();
        return v;
    """)

def _vec_phi(df, p4, out,new_cols):
    new_cols.append(out)
    if out not in df.GetColumnNames():
        return df.Define(out, f"""
            ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
            for (size_t i = 0; i < {p4}.size(); ++i)
                v[i] = {p4}[i].Phi();
            return v;
        """)
    return df.Redefine(out, f"""
        ROOT::VecOps::RVec<float> v({p4}.size(), 0.);
        for (size_t i = 0; i < {p4}.size(); ++i)
            v[i] = {p4}[i].Phi();
        return v;
    """)

def _define_if_missing(df, name, expression):
    if name in _column_names(df): return df
    return df.Define(name, expression)



def ProcessAllJetVariables(df,jet_columns,config,bTagAlgo,bTagDict,want_variations,syst_cfg):
    cols = _column_names(df)
    pt_min = config.get("jet_pt_min", 25.0)
    eta_max = config.get("jet_eta_max", 4.7)
    horn_expr = config.get("jet_horn_veto_expr", "false")
    loose_wp = bTagDict["L"]
    medium_wp = bTagDict["M"]

    ### to be fixed here

    cols = _column_names(df)
    syst_suffixes = jet_variation_suffixes(df, want_variations, syst_cfg)

    jet_extra = {}
    for col in jet_columns:
        if (
            any(x in col.lower() for x in ["pt", "eta", "phi", "mass", "idx"])
            and col != "Jet_genJetIdx"
        ):
            continue
        jet_extra[col.split("Jet_")[-1]] = col

    new_cols = []

    def track(df, name, expr):
        if name not in new_cols:
            new_cols.append(name)
        # Some legacy skim campaigns persisted these derived ordering columns
        # only in a subset of their files.  Reading such files in one TChain
        # makes ROOT bind the branch from the first file and then fail when it
        # reaches a file produced with the other schema.  Recompute the jet
        # ordering from the stable p4/goodJet inputs even when a stored copy is
        # present, so mixed-schema skim campaigns remain readable.
        if name.startswith("SelectedJet_sortIdx") and name in _column_names(df):
            return df.Redefine(name, expr)
        if name in _column_names(df): return df
        return df.Define(name, expr)

    for suff in syst_suffixes:
        p4_branch = f"Jet_p4{suff}"
        df = track(df, f"Jet_idx{suff}", f"CreateIndexes(Jet_p4{suff}.size())")
        horn = horn_expr.replace("Jet_p4", p4_branch)
        df = track(df, f"Jet_IsInsideHorn{suff}", horn)
        df = track(df, f"Jet_IsOutsideHorn{suff}", f"!{horn}")
        if suff=="":
            df = track(df, f"Jet_pt_nocorr", f"Jet_pt")
            df = track(df, f"Jet_mass_nocorr", f"Jet_mass")

        df = _vec_pt(df, p4_branch, f"Jet_pt{suff}",new_cols)
        df = _vec_mass(df, p4_branch, f"Jet_mass{suff}",new_cols)
        df = _vec_eta(df, p4_branch, f"Jet_eta{suff}",new_cols)
        df = _vec_phi(df, p4_branch, f"Jet_phi{suff}",new_cols)

        for b, col in jet_extra.items():
            if col in cols:
                df = track(df, col, f"{col}")
                df = track(
                    df,
                    f"Jet_btag_loose{suff}",
                    f"Jet_btag{bTagAlgo}B >= {loose_wp} && abs(Jet_eta{suff}) < 2.5"
                )

                df = track(
                    df,
                    f"Jet_btag_medium{suff}",
                    f"Jet_btag{bTagAlgo}B >= {medium_wp} && abs(Jet_eta{suff}) < 2.5"
                )

        df = track(
            df,
            f"JetTagSel{suff}",
            f"Jet_p4{suff}[Jet_btag_medium{suff}].size() < 1 && "
            f"Jet_p4{suff}[Jet_btag_loose{suff}].size() < 2"
        )
    return df, new_cols

def SelectJetVars(df,jet_columns,config,bTagAlgo,bTagDict,want_variations,syst_cfg):
    cols = _column_names(df)
    pt_min = config.get("jet_pt_min", 25.0)
    eta_max = config.get("jet_eta_max", 4.7)
    horn_expr = config.get("jet_horn_veto_expr", "false")
    loose_wp = bTagDict["L"]
    medium_wp = bTagDict["M"]

    ### to be fixed here


    cols = _column_names(df)
    syst_suffixes = jet_variation_suffixes(df, want_variations, syst_cfg)

    jet_extra = {}
    for col in jet_columns:
        if (
            any(x in col.lower() for x in ["pt", "eta", "phi", "mass", "idx"])
            and col != "Jet_genJetIdx"
        ):
            continue
        jet_extra[col.split("Jet_")[-1]] = col

    new_cols = []

    def track(df, name, expr):
        if name not in new_cols:
            new_cols.append(name)
        if name in _column_names(df): return df
        return df.Define(name, expr)

    for suff in syst_suffixes:
        p4_branch = f"Jet_p4{suff}"
        df = track(df,f"Jet_preSel{suff}",f"v_ops::pt({p4_branch}) > {pt_min} && abs(v_ops::eta({p4_branch})) < {eta_max} && Jet_passJetIdTight")
        # -----------------------------------------------------
        # MUON CLEANING
        # -----------------------------------------------------
        overlap_expr = f"""
        ROOT::VecOps::RVec<bool> out;
        out.reserve({p4_branch}.size());
        for (size_t i = 0; i < {p4_branch}.size(); ++i) {{
            bool ok = Jet_preSel{suff}[i] && !Jet_vetoMap{suff}[i]
                      && Jet_IsOutsideHorn{suff}[i]
                      && ROOT::Math::VectorUtil::DeltaR({p4_branch}[i], mu1_p4) > 0.4
                      && ROOT::Math::VectorUtil::DeltaR({p4_branch}[i], mu2_p4) > 0.4;
            out.push_back(ok);
        }}
        return out;
        """

        df = track(df, f"Jet_NoOverlapWithMuons{suff}", overlap_expr)
        df = track(df, f"goodJet{suff}", f"Jet_NoOverlapWithMuons{suff}")

        # Sort the selected jets by decreasing pT.  The ordering is recomputed
        # for every jet-energy variation and then applied consistently to all
        # SelectedJet_* properties, so their elements remain aligned.
        df = track(
            df,
            f"SelectedJet_sortIdx{suff}",
            f"Reverse(Argsort(v_ops::pt({p4_branch}[goodJet{suff}])))",
        )
        df = track(
            df,
            f"SelectedJet_idx{suff}",
            f"Take(Jet_idx{suff}[goodJet{suff}], SelectedJet_sortIdx{suff})",
        )
        df = track(df, f"N_SelectedJets{suff}", f"(int)SelectedJet_idx{suff}.size()")
        df = track(
            df,
            f"SelectedJet_IsInsideHorn{suff}",
            f"Take(Jet_IsInsideHorn{suff}[goodJet{suff}], SelectedJet_sortIdx{suff})",
        )
        df = track(
            df,
            f"SelectedJet_IsOutsideHorn{suff}",
            f"Take(Jet_IsOutsideHorn{suff}[goodJet{suff}], SelectedJet_sortIdx{suff})",
        )

        if suff=="":
            df = track(
                df,
                "SelectedJet_pt_nocorr",
                f"Take(Jet_pt_nocorr[goodJet{suff}], SelectedJet_sortIdx{suff})",
            )
            df = track(
                df,
                "SelectedJet_mass_nocorr",
                f"Take(Jet_mass_nocorr[goodJet{suff}], SelectedJet_sortIdx{suff})",
            )

        df = df.Define(
            f"Selected{p4_branch}",
            f"Take({p4_branch}[goodJet{suff}], SelectedJet_sortIdx{suff})",
        )

        df = _vec_pt(df, f"Selected{p4_branch}", f"SelectedJet_pt{suff}",new_cols)
        df = _vec_mass(df, f"Selected{p4_branch}", f"SelectedJet_mass{suff}",new_cols)
        df = _vec_eta(df, f"Selected{p4_branch}", f"SelectedJet_eta{suff}",new_cols)
        df = _vec_phi(df, f"Selected{p4_branch}", f"SelectedJet_phi{suff}",new_cols)

        for b, col in jet_extra.items():
            if col in cols:
                df = track(
                    df,
                    f"SelectedJet_{b}{suff}",
                    f"Take({col}[goodJet{suff}], SelectedJet_sortIdx{suff})",
                )
        df = track(
            df,
            f"SelectedJet_btag_loose{suff}",
            f"SelectedJet_btag{bTagAlgo}B{suff} >= {loose_wp} && abs(SelectedJet_eta{suff}) < 2.5"
        )

        df = track(
            df,
            f"SelectedJet_btag_medium{suff}",
            f"SelectedJet_btag{bTagAlgo}B{suff} >= {medium_wp} && abs(SelectedJet_eta{suff}) < 2.5"
        )

        df = track(
            df,
            f"SelectedJetTagSel{suff}",
            f"SelectedJet_p4{suff}[SelectedJet_btag_medium{suff}].size() < 1 && "
            f"SelectedJet_p4{suff}[SelectedJet_btag_loose{suff}].size() < 2"
        )
    return df, new_cols

def SelectVBFJets(df,want_variations,syst_cfg):

    syst_suffixes = jet_variation_suffixes(df, want_variations, syst_cfg)

    new_cols = []

    def track(df, name, expr):
        if name not in new_cols:
            new_cols.append(name)
        if name in _column_names(df): return df
        return df.Define(name, expr)

    for suff in syst_suffixes:
        df = df.Define(
            f"VBFJetCand{suff}",
            f"FindVBFJets(SelectedJet_p4{suff}, SelectedJet_IsOutsideHorn{suff})"
        )

        df = track(df, f"HasVBF{suff}", f"static_cast<bool>(VBFJetCand{suff}.isVBF)")

        df = track(
            df,
            f"VBFJetIdx_1{suff}",
            f"HasVBF{suff} ? int(VBFJetCand{suff}.leg_index[0]) : -1000"
        )

        df = track(
            df,
            f"VBFJetIdx_2{suff}",
            f"HasVBF{suff} ? int(VBFJetCand{suff}.leg_index[1]) : -1000"
        )

    return df, new_cols
