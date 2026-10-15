#!/bin/bash
# Hands-on reproducibility test (about 3 minutes, needs read access to /eos/cms/store/group/phys_higgs/cmshmm/run3_sync).
# Runs the TAMU-Pisa ggH comparison (Z sideband, inclusive channel) from scratch in a new directory and checks the numbers.
# usage: bash sync/reproducibility_test.sh [work_dir]     (default /tmp/$USER/sync_repro_test)
set -e
cd "$(dirname "$0")/.."
W=${1:-/tmp/$USER/sync_repro_test}; rm -rf "$W"; mkdir -p "$W"
python3 sync/python/compare_run3_sync_groups.py --group1 TAMU --group2 Pisa --sample ggH --region Z_sideband --channel inclusive \
    --output "$W/out" --cache-dir "$W/cache" --no-overlay
python3 sync/python/check_reproducibility.py "$W/out"
echo "plots: $W/out/TAMU-Pisa/ggH/Z_sideband/inclusive/all_common_weights/*.png"
