#!/usr/bin/env python3
"""Write the agreed colon-separated synchronization event list.

Reads one `<region>_events.jsonl` produced by the sync skim and emits the
column order agreed with the other analysis:

    event:run:lumi:mu1_pt:...:jet2_mass

Muon pT are the corrected ones (BSC + KIT) and the dimuon quantities come from
the corrected muons with recovered FSR photons, because those are the values
the skim stores in `mu{1,2}_pt` and `pt_mumu`/`eta_mumu`/`phi_mumu`/`m_mumu`.

Events with fewer than two selected jets get the agreed placeholder for the
missing jet columns; the default is 999.000.
"""
import argparse
import json
from pathlib import Path
import sys

# (output column, source field in the JSONL). Jet columns are resolved
# separately because they index the per-event selected-jet vectors.
SCALAR_COLUMNS = [
    ("event", "event"),
    ("run", "run"),
    ("lumi", "luminosityBlock"),
    ("mu1_pt", "mu1_pt"),
    ("mu1_eta", "mu1_eta"),
    ("mu1_phi", "mu1_phi"),
    ("mu2_pt", "mu2_pt"),
    ("mu2_eta", "mu2_eta"),
    ("mu2_phi", "mu2_phi"),
    ("dimuon_pt", "pt_mumu"),
    ("dimuon_eta", "eta_mumu"),
    ("dimuon_phi", "phi_mumu"),
    ("dimuon_mass", "m_mumu"),
    ("njets", "N_SelectedJets"),
]

JET_FIELDS = ["pt", "eta", "phi", "mass"]

# High-level inputs of the ggH/VBF training, appended with --extra in the order
# of common/updated_DNN_configs/config.toml.
EXTRA_COLUMNS = [
    "dR_mumu", "y_mumu", "cosTheta_CS", "phi_CS", "R_pt",
    "minDeltaEtaSigned", "minDeltaPhi", "Zeppenfeld_Var", "pt_centrality",
    "vbfjet1_pt", "vbfjet1_eta", "vbfjet1_y", "vbfjet1_btagPNetQvG",
    "vbfjet2_pt", "vbfjet2_eta", "vbfjet2_y", "vbfjet2_btagPNetQvG",
    "m_jj", "delta_eta_jj", "pt_vbfj1j2",
    "SoftJetCleanedActivity_N", "SoftJetCleanedActivity_ptSum",
]

# Appended with --mc-weights; absent in data, so the exporter reports them as
# missing rather than silently dropping the column.
MC_WEIGHT_COLUMNS = ["genWeight", "puWeight", "weight_pu_Central"]

INTEGER_COLUMNS = {"event", "run", "lumi", "njets"}


def jet_value(row, index, field, missing):
    """Return selected jet `index` `field`, or the placeholder when absent."""
    values = row.get(f"SelectedJet_{field}")
    if values is None or index >= len(values):
        return missing
    return values[index]


def format_value(name, value, precision, missing_text):
    if value is None:
        return missing_text
    if name in INTEGER_COLUMNS:
        return str(int(value))
    return f"{float(value):.{precision}f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("events", type=Path,
                        help="<region>_events.jsonl from the sync skim")
    parser.add_argument("--output", type=Path,
                        help="Output file; default: stdout")
    parser.add_argument("--precision", type=int, default=3,
                        help="Decimal places for floating-point values (default: 3)")
    parser.add_argument("--missing", default=None,
                        help="Placeholder for absent jets; default: 999 at the chosen precision")
    parser.add_argument("--extra", action="store_true",
                        help="Append the high-level ggH/VBF training inputs")
    parser.add_argument("--mc-weights", action="store_true",
                        help="Append genWeight, pileup weight and muon ID/ISO scale factors")
    parser.add_argument("--no-header", action="store_true",
                        help="Omit the leading column-name line")
    args = parser.parse_args()

    if args.precision < 0:
        parser.error("--precision must be >= 0")
    missing_text = args.missing if args.missing is not None else f"{999.0:.{args.precision}f}"

    rows = []
    with args.events.open() as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    if not rows:
        sys.exit(f"No events in {args.events}")

    columns = list(SCALAR_COLUMNS)
    for jet in (1, 2):
        columns += [(f"jet{jet}_{field}", None) for field in JET_FIELDS]
    if args.extra:
        columns += [(name, name) for name in EXTRA_COLUMNS]
    if args.mc_weights:
        muon_scale_factors = sorted(
            {key for row in rows for key in row
             if key.startswith(("weight_mu1_", "weight_mu2_"))}
        )
        columns += [(name, name) for name in MC_WEIGHT_COLUMNS + muon_scale_factors]

    absent = sorted({name for name, source in columns
                     if source is not None and source not in rows[0]})
    if absent:
        print(f"[WARN] columns absent from {args.events.name}, filled with "
              f"{missing_text}: {', '.join(absent)}", file=sys.stderr)

    lines = []
    if not args.no_header:
        lines.append(":".join(name for name, _ in columns))
    for row in rows:
        values = []
        for name, source in columns:
            if source is None:
                jet = int(name[3]) - 1
                field = name.split("_", 1)[1]
                value = jet_value(row, jet, field, None)
            else:
                value = row.get(source)
            values.append(format_value(name, value, args.precision, missing_text))
        lines.append(":".join(values))

    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
        print(f"{len(rows)} events -> {args.output}")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
