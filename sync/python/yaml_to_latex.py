#!/usr/bin/env python3
"""Create a clean standalone LaTeX table from a samples YAML file.

The YAML is expected to have one top-level entry per sample and a ``nanoAOD``
field containing one or more DAS dataset paths.

Example:
    python3 samples_yaml_to_latex.py \
        --input samples.yaml \
        --output samples_table.tex \
        --title "NanoAOD 2024 datasets used in the analysis"

Optional PDF compilation:
    python3 samples_yaml_to_latex.py \
        --input samples.yaml \
        --output samples_table.tex \
        --compile-pdf
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

import yaml


LATEX_ESCAPES = {
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
    "\\": r"\textbackslash{}",
}


def escape_latex_text(value: Any) -> str:
    """Escape ordinary LaTeX text."""
    return "".join(LATEX_ESCAPES.get(ch, ch) for ch in str(value))


def breakable_monospace(value: Any) -> str:
    """Return monospace LaTeX text with explicit safe break opportunities.

    Long sample names and DAS paths otherwise tend to overflow p-columns.
    Breaks are inserted after separators while keeping their visual form.
    """
    pieces: list[str] = []
    break_after = {"/", "_", "-", ".", ":"}

    for ch in str(value):
        if ch == "_":
            pieces.append(r"\_")
        else:
            pieces.append(LATEX_ESCAPES.get(ch, ch))

        if ch in break_after:
            pieces.append(r"\allowbreak{}")

    return r"\texttt{" + "".join(pieces) + "}"


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Input YAML not found: {path}")

    with path.open("r", encoding="utf-8") as stream:
        content = yaml.safe_load(stream)

    if content is None:
        return {}
    if not isinstance(content, dict):
        raise ValueError("The top level of the YAML file must be a mapping.")
    return content


def normalise_datasets(sample_name: str, sample_info: Any) -> list[str]:
    if not isinstance(sample_info, dict):
        raise ValueError(f"Sample '{sample_name}' must contain a mapping.")

    datasets = sample_info.get("nanoAOD")
    if datasets is None:
        return []
    if isinstance(datasets, str):
        return [datasets]
    if isinstance(datasets, list):
        if not all(isinstance(item, str) for item in datasets):
            raise ValueError(
                f"Sample '{sample_name}' has a nanoAOD list containing non-string values."
            )
        return datasets

    raise ValueError(
        f"Sample '{sample_name}' has nanoAOD of type {type(datasets).__name__}; "
        "expected a string or list of strings."
    )


def split_samples(data: dict[str, Any]) -> tuple[list[tuple[str, Any]], list[tuple[str, Any]]]:
    """Split entries into data and simulation while preserving YAML order."""
    data_samples: list[tuple[str, Any]] = []
    mc_samples: list[tuple[str, Any]] = []

    for name, info in data.items():
        if isinstance(info, dict) and info.get("is_data") is True:
            data_samples.append((name, info))
        else:
            mc_samples.append((name, info))

    return data_samples, mc_samples


def emit_section(
    lines: list[str],
    heading: str,
    entries: Iterable[tuple[str, Any]],
) -> None:
    entries = list(entries)
    if not entries:
        return

    lines.extend(
        [
            r"\addlinespace[4pt]",
            rf"\multicolumn{{2}}{{@{{}}l}}{{\bfseries {escape_latex_text(heading)}}} \\",
            r"\addlinespace[2pt]",
        ]
    )

    for sample_name, sample_info in entries:
        datasets = normalise_datasets(sample_name, sample_info)
        sample_cell = breakable_monospace(sample_name)

        if not datasets:
            lines.append(
                sample_cell + r" & \textit{No NanoAOD dataset specified} \\"
            )
            lines.append(r"\addlinespace[2pt]")
            continue

        # Use one row for each NanoAOD path. The sample name is printed only once.
        for index, dataset in enumerate(datasets):
            left_cell = sample_cell if index == 0 else ""
            lines.append(
                left_cell + " & " + breakable_monospace(dataset) + r" \\"
            )

        lines.append(r"\addlinespace[2pt]")


def build_document(
    data: dict[str, Any],
    title: str,
    margin: str,
    sample_width: float,
    font_size: str,
    group_data_mc: bool,
) -> str:
    if not 0.15 <= sample_width <= 0.55:
        raise ValueError("--sample-width must be between 0.15 and 0.55.")

    dataset_width = 0.97 - sample_width
    data_samples, mc_samples = split_samples(data)

    lines = [
        # r"\documentclass[10pt,a4paper,landscape]{article}",
        # rf"\usepackage[margin={margin}]{{geometry}}",
        # r"\usepackage[T1]{fontenc}",
        # r"\usepackage{lmodern}",
        # r"\usepackage{microtype}",
        # r"\usepackage{booktabs}",
        # r"\usepackage{longtable}",
        # r"\usepackage{array}",
        # r"\usepackage{ragged2e}",
        # r"\pagestyle{plain}",
        # r"\setlength{\parindent}{0pt}",
        # r"\setlength{\tabcolsep}{4pt}",
        # r"\renewcommand{\arraystretch}{1.12}",
        # r"\setlength{\emergencystretch}{3em}",
        # r"\sloppy",
        # r"\begin{document}",
        r"\centering",
        rf"\{font_size}",
        rf"\begin{{longtable}}{{@{{}}>{{\RaggedRight\arraybackslash}}p{{{sample_width:.3f}\textwidth}} >{{\RaggedRight\arraybackslash}}p{{{dataset_width:.3f}\textwidth}}@{{}}}}",
        rf"\caption{{{escape_latex_text(title)}}} \\",
        r"\toprule",
        r"\textbf{Sample} & \textbf{NanoAOD dataset} \\",
        r"\midrule",
        r"\endfirsthead",
        rf"\caption[]{{{escape_latex_text(title)} (continued)}} \\",
        r"\toprule",
        r"\textbf{Sample} & \textbf{NanoAOD dataset} \\",
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{2}{r}{\footnotesize Continued on next page} \\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
    ]

    if group_data_mc:
        emit_section(lines, f"Data ({len(data_samples)} samples)", data_samples)
        emit_section(lines, f"Simulation ({len(mc_samples)} samples)", mc_samples)
    else:
        emit_section(lines, f"Samples ({len(data)} entries)", data.items())

    lines.extend(
        [
            r"\end{longtable}",
            # r"\end{document}",
            "",
        ]
    )
    return "\n".join(lines)


def compile_pdf(tex_path: Path) -> Path:
    command = [
        "latexmk",
        "-pdf",
        "-interaction=nonstopmode",
        "-halt-on-error",
        tex_path.name,
    ]
    result = subprocess.run(
        command,
        cwd=tex_path.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        log_tail = "\n".join((result.stdout + "\n" + result.stderr).splitlines()[-40:])
        raise RuntimeError(f"LaTeX compilation failed:\n{log_tail}")
    return tex_path.with_suffix(".pdf")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a samples YAML file into a clean standalone LaTeX table."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input YAML file")
    parser.add_argument("--output", required=True, type=Path, help="Output .tex file")
    parser.add_argument(
        "--title",
        default="NanoAOD datasets used in the analysis",
        help="Table caption",
    )
    parser.add_argument(
        "--margin",
        default="0.35cm",
        help="Page margin passed to geometry (default: 0.35cm)",
    )
    parser.add_argument(
        "--sample-width",
        type=float,
        default=0.28,
        help="Fraction of text width assigned to the sample column (default: 0.28)",
    )
    parser.add_argument(
        "--font-size",
        choices=("small", "footnotesize", "scriptsize"),
        default="footnotesize",
        help="Font size used inside the table (default: footnotesize)",
    )
    parser.add_argument(
        "--no-groups",
        action="store_true",
        help="Do not separate data and simulation samples",
    )
    parser.add_argument(
        "--compile-pdf",
        action="store_true",
        help="Compile the generated .tex file with latexmk",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        content = load_yaml(args.input)
        document = build_document(
            data=content,
            title=args.title,
            margin=args.margin,
            sample_width=args.sample_width,
            font_size=args.font_size,
            group_data_mc=not args.no_groups,
        )

        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(document, encoding="utf-8")
        print(f"Written: {args.output}")

        if args.compile_pdf:
            pdf = compile_pdf(args.output)
            print(f"Written: {pdf}")

    except (OSError, ValueError, yaml.YAMLError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())