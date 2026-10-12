#!/usr/bin/env python3
"""Build the separate PDFs published in the MkDocs Sync section."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ERAS = ("2022", "2022EE", "2023", "2023BPix", "2024", "2025", "2026")


def compile_tex(repo: Path, source: Path, output: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="hmm-sync-") as tmp_name:
        tmp = Path(tmp_name)
        command = [
            "pdflatex", "-interaction=nonstopmode", "-halt-on-error",
            f"-jobname={output.stem}", f"-output-directory={tmp}", str(source),
        ]
        for _ in range(2):
            result = subprocess.run(command, cwd=repo, text=True, capture_output=True)
            if result.returncode:
                tail = "\n".join((result.stdout + result.stderr).splitlines()[-40:])
                raise RuntimeError(f"LaTeX failed for {source}:\n{tail}")
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tmp / f"{output.stem}.pdf", output)


def dataset_wrapper(repo: Path, era: str, directory: Path) -> Path:
    wrapper = directory / f"datasets_{era}_main.tex"
    wrapper.write_text(
        "\\documentclass[8pt,a4paper,landscape]{extarticle}\n"
        "\\usepackage[margin=0.35cm]{geometry}\n"
        "\\usepackage[T1]{fontenc}\n\\usepackage{lmodern}\n"
        "\\usepackage{booktabs}\n\\usepackage{array}\n\\usepackage{url}\n"
        "\\usepackage{longtable}\n\\usepackage{caption}\n"
        "\\newcommand{\\RaggedRight}{\\raggedright}\n"
        "\\urlstyle{tt}\n\\setlength{\\tabcolsep}{2.5pt}\n"
        "\\renewcommand{\\arraystretch}{1.12}\n\\sloppy\n"
        "\\begin{document}\n\\footnotesize\n"
        f"\\input{{sync/latex/datasets_{era}.tex}}\n"
        "\\end{document}\n",
        encoding="utf-8",
    )
    return wrapper


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--datasets-only",
        action="store_true",
        help="regenerate and compile only the per-era dataset inventories",
    )
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    output_dir = args.output_dir or repo / "docs/sync"

    # Regenerate source tables before compiling them.
    if not args.datasets_only:
        subprocess.run(
            [str(repo / "sync/python/collect_corrections.py"), "--format", "latex", "--output", str(repo / "sync/latex/corrections_table.tex")],
            cwd=repo,
            check=True,
        )
        for generator in (
            "collect_skim_selections.py",
            "collect_histogram_workflow.py",
            "collect_analysis_selections.py",
        ):
            subprocess.run(
                [sys.executable, str(repo / "sync/python" / generator)],
                cwd=repo,
                check=True,
            )
        subprocess.run(
            [
                sys.executable,
                str(repo / "sync/python/cross_sections_yaml_to_latex.py"),
                "--input", str(repo / "config/crossSections13p6TeV.yaml"),
                "--output", str(repo / "sync/latex/cross_sections_table.tex"),
            ],
            cwd=repo,
            check=True,
        )
    for era in ERAS:
        subprocess.run(
            [
                sys.executable,
                str(repo / "sync/python/yaml_to_latex.py"),
                "--input", str(repo / f"config/Run3_{era}/samples.yaml"),
                "--output", str(repo / f"sync/latex/datasets_{era}.tex"),
                "--title", f"NanoAOD {era} datasets used in the analysis",
            ],
            cwd=repo,
            check=True,
        )

    if not args.datasets_only:
        compile_tex(repo, repo / "sync/latex/corrections_main.tex", output_dir / "corrections.pdf")
        compile_tex(repo, repo / "sync/latex/skim_selections_main.tex", output_dir / "skim_selections.pdf")
        compile_tex(repo, repo / "sync/latex/histogram_workflow_main.tex", output_dir / "histogram_workflow.pdf")
        # Analyst-facing compact selection flow; detailed skim/histogram PDFs are
        # retained separately as technical references.
        compile_tex(repo, repo / "sync/latex/selections_main.tex", output_dir / "selections.pdf")
        compile_tex(repo, repo / "sync/latex/cross_sections_main.tex", output_dir / "cross_sections.pdf")

    with tempfile.TemporaryDirectory(prefix="hmm-datasets-") as tmp_name:
        tmp = Path(tmp_name)
        for era in ERAS:
            compile_tex(repo, dataset_wrapper(repo, era, tmp), output_dir / f"datasets_{era}.pdf")

    label = "dataset" if args.datasets_only else "Sync"
    print(f"Wrote separate {label} PDFs to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
