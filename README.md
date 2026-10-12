# Minimal public H→μμ skim and histogram production

This directory is a standalone, publishable subset of the Run 3 H→μμ
analysis. It contains only the code, configurations, correction payloads and
ML models needed to turn NanoAOD files into skim files and skim files into ROOT
histograms.

Full public analysis documentation:
<https://valeriadamante.github.io/Hmumu_analysis/>.

Standard physics corrections and scale factors remain
enabled: luminosity and cross-section normalization, pileup, muon, trigger,
jet/JEC/JER, jet-veto-map and b-tag corrections.

## Contents

- `analysis/`: NanoAOD-to-skim code;
- `histograms/`: skim-to-histogram code;
- `common/`: shared selections, variables, DNN inference and ROOT helpers;
- `corrections/`: correction code and public payloads;
- `config/`: Run 3 era, sample, selection, systematic and histogram configs;
- `env.sh`: CMSSW/ROOT environment setup on a CMS machine.

Generated file lists, private output locations, campaign orchestration,
datacards, plots, studies and derived reweight payloads are outside the scope
of this minimal release.

## Setup

On lxplus or another machine with CVMFS and CMS software available:

```bash
source env.sh
```

The first invocation creates `soft/CMSSW_15_0_2`. A valid grid proxy may be
needed for remote NanoAOD input.

## Produce a skim

Dataset names must exist in `config/<ERA>/samples.yaml`.

```bash
python3 analysis/skim.py \
  --era Run3_2024 \
  --dataset-name DYto2Mu_MLL_105to160_amcatnloFXFX \
  --input-file root://cms-xrd-global.cern.ch//store/.../nano.root \
  --output-file output/skim.root \
  --report-file output/skim_report.json
```

Repeat `--input-file` or provide a comma-separated list to process more files.
Use `--n-events 10000` for a short test.

## Produce histograms

The report JSON is required for the full-dataset MC normalization. For one
skim/report pair:

```bash
python3 histograms/hist_maker.py \
  --era Run3_2024 \
  --dataset-name DYto2Mu_MLL_105to160_amcatnloFXFX \
  --root-input output/skim.root \
  --json-input output/skim_report.json \
  --output-file output/histograms.root \
  --systematics Central
```

Use `--systematics all` for all configured variations. Run
`python3 histograms/hist_maker.py --help` for batching, category and variable
options.

## Site-specific paths

`config/paths.yaml` contains placeholders only. Direct invocations above do not
depend on those output paths. If integrating the code into a batch workflow,
override them using `HMUMU_<KEY>` environment variables or an
`HMUMU_PATHS=/path/to/paths.yaml` file.

## Public-release checklist

The directory contains no derived DY reweight payloads, generated dataset file
lists, credentials, grid proxies, output ROOT files or user-specific storage
paths. Before publishing, choose and add the repository license required by
your collaboration or institution.
