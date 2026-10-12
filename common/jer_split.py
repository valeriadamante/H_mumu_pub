"""Auxiliary functions for the CMS JER Case-I histogram split."""

JER_SPLIT_COMPONENTS = {
    "JEReta0pt0": "aeta < 1.93",
    "JEReta1pt0": "aeta >= 1.93 && aeta < 2.5",
    "JEReta2pt0": "aeta >= 2.5 && aeta < 3.0 && pt < 50.0",
    "JEReta2pt1": "aeta >= 2.5 && aeta < 3.0 && pt >= 50.0",
    "JEReta3pt0": "aeta >= 3.0 && aeta < 5.0 && pt < 50.0",
    "JEReta3pt1": "aeta >= 3.0 && aeta < 5.0 && pt >= 50.0",
}


def define_split_jer_collections(rdf, systs_to_run):
    """Build six hybrid collections from the nominal and full-JER jets.

    A jet receives the full JER up/down shift only in its nominal pt/|eta|
    component. All other jets retain nominal pt and mass. Jets are matched by
    their original index and re-sorted after the component shift.
    """
    columns = {str(column) for column in rdf.GetColumnNames()}
    nominal_properties = sorted(
        column
        for column in columns
        if column.startswith("SelectedJet_")
        and not any(token in column for token in ("_JER", "_JES"))
        and column not in ("SelectedJet_pt_nocorr", "SelectedJet_mass_nocorr")
        and column != "SelectedJet_p4"
    )
    for info in systs_to_run.values():
        component = info.get("jer_component")
        target_suffix = info.get("jet_suffix", "")
        source_suffix = info.get("source_jet_suffix", "")
        if not component or not target_suffix:
            continue
        condition = JER_SPLIT_COMPONENTS[component]
        required = {
            "SelectedJet_idx", "SelectedJet_pt", "SelectedJet_eta",
            "SelectedJet_mass",
            f"Jet_idx{source_suffix}",
            f"Jet_pt{source_suffix}",
            f"Jet_mass{source_suffix}",
        }
        missing = sorted(required - columns)
        if missing:
            raise RuntimeError(
                f"Cannot build split JER variation '{target_suffix}'; "
                "missing columns: " + ", ".join(missing)
            )
        raw_pt = f"__splitJER_pt{target_suffix}"
        raw_mass = f"__splitJER_mass{target_suffix}"
        order = f"__splitJER_order{target_suffix}"
        for observable, output in (("pt", raw_pt), ("mass", raw_mass)):
            shifted = f"Jet_{observable}{source_suffix}"
            expression = f"""
                auto out = SelectedJet_{observable};
                for (size_t i = 0; i < SelectedJet_idx.size(); ++i) {{
                    const float pt = SelectedJet_pt[i];
                    const float aeta = std::abs(SelectedJet_eta[i]);
                    if (!({condition})) continue;
                    for (size_t j = 0; j < {shifted}.size(); ++j) {{
                        if (SelectedJet_idx[i] == Jet_idx{source_suffix}[j]) {{
                            out[i] = {shifted}[j];
                            break;
                        }}
                    }}
                }}
                return out;
            """
            rdf = rdf.Define(output, expression)
        rdf = rdf.Define(order, f"Reverse(Argsort({raw_pt}))")
        for source in nominal_properties:
            stem = source[len("SelectedJet_"):]
            target = f"SelectedJet_{stem}{target_suffix}"
            shifted_values = {"pt": raw_pt, "mass": raw_mass}
            value = shifted_values.get(stem, source)
            rdf = rdf.Define(target, f"Take({value}, {order})")
        rdf = rdf.Define(
            f"SelectedJet_p4{target_suffix}",
            f"""
                ROOT::VecOps::RVec<ROOT::Math::LorentzVector<
                    ROOT::Math::PtEtaPhiM4D<double>>> out;
                out.reserve(SelectedJet_pt{target_suffix}.size());
                for (size_t i = 0; i < SelectedJet_pt{target_suffix}.size(); ++i)
                    out.emplace_back(
                        SelectedJet_pt{target_suffix}[i],
                        SelectedJet_eta{target_suffix}[i],
                        SelectedJet_phi{target_suffix}[i],
                        SelectedJet_mass{target_suffix}[i]);
                return out;
            """,
        )
        rdf = rdf.Define(
            f"N_SelectedJets{target_suffix}",
            f"static_cast<int>(SelectedJet_idx{target_suffix}.size())",
        )
        rdf = rdf.Define(
            f"SelectedJetTagSel{target_suffix}",
            f"SelectedJet_btag_medium{target_suffix}[SelectedJet_btag_medium{target_suffix}].size() < 1 && "
            f"SelectedJet_btag_loose{target_suffix}[SelectedJet_btag_loose{target_suffix}].size() < 2",
        )
        rdf = rdf.Define(
            f"VBFJetCand{target_suffix}",
            f"FindVBFJets(SelectedJet_p4{target_suffix}, SelectedJet_IsOutsideHorn{target_suffix})",
        )
        rdf = rdf.Define(
            f"HasVBF{target_suffix}",
            f"static_cast<bool>(VBFJetCand{target_suffix}.isVBF)",
        )
        for leg in (1, 2):
            rdf = rdf.Define(
                f"VBFJetIdx_{leg}{target_suffix}",
                f"HasVBF{target_suffix} ? int(VBFJetCand{target_suffix}.leg_index[{leg - 1}]) : -1000",
            )
        columns = {str(column) for column in rdf.GetColumnNames()}
    return rdf
