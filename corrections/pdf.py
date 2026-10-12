"""PDF replica uncertainty from NanoAOD LHE PDF weights."""

import ROOT


_PDF_HELPER = r"""
#include <algorithm>
#include <cmath>
#include <ROOT/RVec.hxx>

namespace pdf_uncertainty {
inline float replicaRMS(
    const ROOT::VecOps::RVec<float>& weights,
    const unsigned int nominal_index,
    const unsigned int replica_start,
    const unsigned int replica_count
) {
    if (weights.size() <= nominal_index || replica_count < 2 ||
        weights.size() < replica_start + replica_count) {
        return 0.f;
    }

    const float nominal = static_cast<float>(weights[nominal_index]);
    if (!std::isfinite(nominal) || nominal == 0.f) {
        return 0.f;
    }

    double sum_squared = 0.;
    unsigned int valid_replicas = 0;
    for (unsigned int index = replica_start;
         index < replica_start + replica_count;
         ++index) {
        const float value = static_cast<float>(weights[index]);
        if (!std::isfinite(value)) {
            continue;
        }
        const double difference = static_cast<double>(value / nominal) - 1.;
        sum_squared += difference * difference;
        ++valid_replicas;
    }

    if (valid_replicas < 2) {
        return 0.f;
    }
    return static_cast<float>(
        std::sqrt(sum_squared / static_cast<double>(valid_replicas - 1))
    );
}

inline float replicaVariation(
    const ROOT::VecOps::RVec<float>& weights,
    const unsigned int nominal_index,
    const unsigned int replica_start,
    const unsigned int replica_count,
    const bool take_up
) {
    const float uncertainty = replicaRMS(
        weights, nominal_index, replica_start, replica_count
    );
    return take_up ? 1.f + uncertainty
                   : std::max(0.f, 1.f - uncertainty);
}
}  // namespace pdf_uncertainty
"""

if not ROOT.gInterpreter.Declare(_PDF_HELPER):
    raise RuntimeError("Failed to declare the PDF uncertainty helper")


def define_pdf_weights(df, config):
    """Define PDF up/down factors from the RMS of PDF replicas."""
    config = config or {}
    if not config.get("enabled", True):
        return (
            df.Define("weight_pdf_up", "1.f")
              .Define("weight_pdf_down", "1.f")
        )

    branch = config.get("branch", "LHEPdfWeight")
    prescription = config.get("prescription", "replicas_rms")
    missing_policy = config.get("missing_branch", "unity")
    nominal_index = int(config.get("nominal_index", 0))
    replica_start = int(config.get("replica_start", 1))
    replica_count = int(config.get("replica_count", 100))

    if prescription != "replicas_rms":
        raise ValueError(
            "Only pdf.prescription='replicas_rms' is supported"
        )
    if nominal_index < 0 or replica_start < 0 or replica_count < 2:
        raise ValueError(
            "PDF indices must be non-negative and replica_count must be >= 2"
        )

    available_columns = {str(column) for column in df.GetColumnNames()}
    if branch not in available_columns:
        message = (
            f"PDF branch '{branch}' is missing; "
            "PDFUp/Down will be set to unity."
        )
        if missing_policy == "error":
            raise RuntimeError(message)
        if missing_policy != "unity":
            raise ValueError(
                "pdf.missing_branch must be either 'unity' or 'error'"
            )
        print(f"[WARNING] {message}")
        return (
            df.Define("weight_pdf_up", "1.f")
              .Define("weight_pdf_down", "1.f")
        )

    arguments = f"{branch}, {nominal_index}u, {replica_start}u, {replica_count}u"
    return (
        df.Define(
            "weight_pdf_up",
            f"pdf_uncertainty::replicaVariation({arguments}, true)",
        )
        .Define(
            "weight_pdf_down",
            f"pdf_uncertainty::replicaVariation({arguments}, false)",
        )
    )


__all__ = ["define_pdf_weights"]
