#!/bin/bash
# Condor wrapper for special_sync_skim.py: runs in the job sandbox and copies
# the tuple and its report to EOS.
proxy="$1"; analysis_path="$2"; dataset="$3"; input_file="$4"; output_dir="$5"; job_name="$6"; variant="${7:-special_20260923}"; extra="${8:-none}"; [ "$extra" = none ] && extra=""
export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY="$proxy"
sandbox="$PWD"
cd "$analysis_path"
source ./env.sh --cmssw-version CMSSW_15_0_2 || true
set -euo pipefail
python3 sync/python/special_sync_skim.py --dataset-name "$dataset" --input-file "$input_file" --variant "$variant" $extra \
  --output-file "$sandbox/$job_name.root" --report-file "$sandbox/${job_name}_report.json"
for f in "$job_name.root" "${job_name}_report.json"; do
  xrdcp -f "$sandbox/$f" "$output_dir/$f"
done
rm -f "$sandbox/$job_name.root" "$sandbox/${job_name}_report.json"
