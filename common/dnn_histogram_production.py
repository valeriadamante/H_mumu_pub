"""DNN-specific dataframe transformations used while producing histograms.

Keep mass remapping and DNN reevaluation in common: the histogram maker
should only decide which dataframe and output column to book.
"""

from __future__ import annotations


SIDEBAND_PAYLOADS = {
    "Z_sideband": "DNNZSidebandMassShift",
    "H_sideband": "DNNHSidebandMassShift",
}

DNN_OUTPUT_VARIABLE = "DNN_NNOutput"


def needs_sideband_mass_shift(mass_region, variable):
    """Return whether *variable* needs a sideband-specific DNN evaluation."""
    return mass_region in SIDEBAND_PAYLOADS and variable == DNN_OUTPUT_VARIABLE


def sideband_mass_expression(mass_region): 
    if mass_region == "Z_sideband":
        return "static_cast<float>(115.0 + 0.5 * (m_mumu - 70.0))"
    if mass_region == "H_sideband":
        return (
            "static_cast<float>(m_mumu < 115.0 ? "
            "115.0 + (m_mumu - 110.0) : 120.0 + (m_mumu - 135.0))"
        )
    raise ValueError(f"Unsupported DNN mass-shift region: {mass_region}")


def shifted_output_column(mass_region):
    """Name of the DNN output produced for a shifted sideband."""
    try:
        return f"{SIDEBAND_PAYLOADS[mass_region]}_NNOutput"
    except KeyError as error:
        raise ValueError(
            f"Unsupported DNN mass-shift region: {mass_region}"
        ) from error


def apply_sideband_mass_shifted_dnn(
    rdf, mass_region, *, btag_algo, era, model_set, jet_suffix="", muon_suffix=""
):
    """Remap ``m_mumu`` and evaluate only the sideband-specific DNN payload.

    With a jet or muon suffix the payload runs on the shifted inputs (the node
    is already filtered on the shifted selection); a muon variation also
    remaps its own shifted mass.
    """
    from common.dnn_application import ApplyDNN, DNNApplication

    shifted_rdf = rdf.Redefine("m_mumu", sideband_mass_expression(mass_region))
    payload = SIDEBAND_PAYLOADS[mass_region]
    suffix = jet_suffix or muon_suffix
    if not suffix:
        return ApplyDNN(shifted_rdf, [payload], btag_algo=btag_algo, era=era, model_set=model_set)
    column_map = {}
    shifted_mass = f"m_mumu{muon_suffix}" if muon_suffix else ""
    if shifted_mass and shifted_mass in {str(c) for c in rdf.GetColumnNames()}:
        remapped = f"m_mumu_sideband{muon_suffix}"
        shifted_rdf = shifted_rdf.Define(
            remapped, sideband_mass_expression(mass_region).replace("m_mumu", shifted_mass))
        column_map["m_mumu"] = remapped
    else:
        # A jet variation keeps the (remapped) nominal mass.
        column_map["m_mumu"] = "m_mumu"
    application = DNNApplication(payload_name=payload, btag_algo=btag_algo, era=era, model_set=model_set)
    return application.apply_shifted(shifted_rdf, suffix, column_map)
