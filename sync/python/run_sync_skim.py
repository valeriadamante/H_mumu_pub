#!/usr/bin/env python3
"""Run one exact NanoAOD file through the skim_v4 synchronization selection."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from common.paths import path  # noqa: E402
DEFAULT_OUTPUT = str(path('skim')) + '_sync'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--era', default='Run3_2024')
    parser.add_argument('--dataset', default='Muon0_Run2024C')
    parser.add_argument('--input-file', help='Exactly one agreed NanoAOD file; default: first configured file.')
    parser.add_argument('--output-dir', type=Path, default=Path(DEFAULT_OUTPUT))
    parser.add_argument('--category', choices=('baseline', 'VBF', 'ggF'), default='baseline')
    parser.add_argument('--n-events', type=int, default=-1, help='-1: entire file; positive: separate test subdirectory.')
    parser.add_argument('--run', action='store_true', help='Execute; otherwise print the exact plan.')
    args = parser.parse_args()
    if args.n_events != -1 and args.n_events <= 0:
        parser.error('--n-events must be -1 or positive')
    config_dir = REPO / 'config' / args.era
    samples = yaml.safe_load((config_dir / 'samples.yaml').read_text())
    sample = samples[args.dataset]
    input_file = args.input_file
    if not input_file:
        resolved = yaml.safe_load((config_dir / 'samples_withfiles.yaml').read_text())
        input_file = resolved[args.dataset]['filelist'][0]
    if ',' in input_file:
        parser.error('Synchronization requires exactly one input file')
    file_tag = Path(input_file).stem
    output = args.output_dir.resolve() / args.era / args.dataset / file_tag / args.category
    if args.n_events > 0:
        output = output / f'test_{args.n_events}'
    command = [sys.executable, str(REPO / 'analysis/skim.py'), '--era', args.era,
               '--dataset-name', args.dataset, '--input-file', input_file,
               '--output-file', str(output / 'events.root'), '--report-file', str(output / 'cutflow.json'),
               '--n-events', str(args.n_events), '--sync', '--sync-category', args.category]
    print('Dataset:', sample['nanoAOD'], flush=True)
    print('Input:', input_file, flush=True)
    print('Output:', output, flush=True)
    print(shlex.join(command), flush=True)
    if not args.run:
        return
    output.mkdir(parents=True, exist_ok=False)
    config_out = output / 'configuration'
    config_out.mkdir()
    hashes = {}
    paths = list(config_dir.glob('*.yaml'))
    paths += list((REPO / 'analysis').glob('*.py')) + list((REPO / 'corrections').glob('*.py'))
    paths += list((REPO / 'corrections').glob('*.h'))
    paths += [REPO / 'common/sync_skim.py', REPO / 'common/add_vars.py', REPO / 'common/jet_horn_policy.py']
    for path in paths:
        hashes[str(path.relative_to(REPO))] = hashlib.sha256(path.read_bytes()).hexdigest()
    for filename in ('maincfg.yaml', 'selections.yaml', 'triggers.yaml', 'systematics.yaml'):
        shutil.copy2(config_dir / filename, config_out / filename)
    manifest = {'command': command, 'input_file': input_file, 'dataset': args.dataset,
                'dataset_configuration': sample, 'sha256': hashes,
                'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
                'nominal_only': True, 'status': 'running'}
    manifest_path = output / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    env = dict(os.environ, ANALYSIS_PATH=str(REPO))
    with (output / 'production.log').open('w') as log:
        result = subprocess.run(command, cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT)
    manifest['status'] = 'complete' if result.returncode == 0 else 'failed'
    manifest['returncode'] = result.returncode
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    if result.returncode:
        raise SystemExit(f"Skim failed: inspect {output / 'production.log'}")
    print((output / 'summary.json').read_text())


if __name__ == '__main__':
    main()
