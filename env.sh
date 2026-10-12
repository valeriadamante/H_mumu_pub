#!/usr/bin/env bash
# ==========================================================
# H_mumu env: CMSSW-only, no pip install
#
# Usage:
#   source env.sh
#   source env.sh --cmssw-version CMSSW_15_0_2
#   source env.sh --cmssw-version CMSSW_15_3_0
# ==========================================================

if [ -n "${ZSH_VERSION:-}" ]; then
    # shellcheck disable=SC2296  # zsh-only expansion; this file is sourceable in zsh.
    this_file="${(%):-%x}"
else
    this_file="${BASH_SOURCE[0]}"
fi

ANALYSIS_PATH="$(cd "$(dirname "$this_file")" && pwd)"
export ANALYSIS_PATH
export ANALYSIS_DATA_PATH="${ANALYSIS_PATH}/data"
export ANALYSIS_SOFT_PATH="${ANALYSIS_PATH}/soft"

mkdir -p "$ANALYSIS_DATA_PATH"
mkdir -p "$ANALYSIS_SOFT_PATH"

export X509_USER_PROXY="${ANALYSIS_DATA_PATH}/voms.proxy"

DEFAULT_CMSSW_VERSION="CMSSW_15_0_2"
REQUESTED_CMSSW_VERSION="$DEFAULT_CMSSW_VERSION"
FORCE_REINSTALL=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cmssw-version)
            REQUESTED_CMSSW_VERSION="$2"
            shift 2
            ;;
        --force-reinstall)
            FORCE_REINSTALL=1
            shift
            ;;
        *)
            # "source env.sh" inside a script sees the script's arguments
            # (an era, an action, --eras ...). They are not for us: failing here
            # left the environment empty, and the scripts that do ">/dev/null"
            # did not notice.
            echo "[WARNING] env.sh: ignoring argument '$1'"
            shift
            ;;
    esac
done

get_scram_arch() {
    local cmssw_version="$1"

    case "$cmssw_version" in
        CMSSW_15_0_*)
            echo "el9_amd64_gcc12"
            ;;
        CMSSW_15_1_*|CMSSW_15_2_*|CMSSW_15_3_*)
            echo "el9_amd64_gcc13"
            ;;
        CMSSW_16_*)
            echo "el8_amd64_gcc13"
            ;;
        *)
            echo "el9_amd64_gcc13"
            ;;
    esac
}

remove_path_matching() {
    local var_name="$1"
    local pattern="$2"

    eval "local old_value=\"\${${var_name}:-}\""

    if [[ -z "$old_value" ]]; then
        return 0
    fi

    local new_value
    new_value="$(echo "$old_value" | tr ':' '\n' | grep -v -E "$pattern" | paste -sd: -)"

    eval "export ${var_name}=\"${new_value}\""
}

prepend_path() {
    local var_name="$1"
    local new_path="$2"

    if [[ -z "$new_path" ]]; then
        return 0
    fi

    eval "local old_value=\"\${${var_name}:-}\""

    if [[ -z "$old_value" ]]; then
        eval "export ${var_name}=\"${new_path}\""
    else
        case ":$old_value:" in
            *":$new_path:"*) ;;
            *) eval "export ${var_name}=\"${new_path}:${old_value}\"" ;;
        esac
    fi
}

clean_env() {
    unset PYTHONHOME
    unset PYTHONSTARTUP
    unset PYTHONUSERBASE
    unset ROOTSYS

    # remove LCG
    remove_path_matching PYTHONPATH "/cvmfs/sft.cern.ch/lcg"
    remove_path_matching PATH "/cvmfs/sft.cern.ch/lcg"
    remove_path_matching LD_LIBRARY_PATH "/cvmfs/sft.cern.ch/lcg"
    remove_path_matching ROOT_INCLUDE_PATH "/cvmfs/sft.cern.ch/lcg"
    remove_path_matching CMAKE_PREFIX_PATH "/cvmfs/sft.cern.ch/lcg"

    # remove CMS LCG/ROOT fragments
    remove_path_matching PYTHONPATH "/cvmfs/cms.cern.ch/.*/lcg"
    remove_path_matching PATH "/cvmfs/cms.cern.ch/.*/lcg"
    remove_path_matching LD_LIBRARY_PATH "/cvmfs/cms.cern.ch/.*/lcg"
    remove_path_matching ROOT_INCLUDE_PATH "/cvmfs/cms.cern.ch/.*/lcg"
    remove_path_matching CMAKE_PREFIX_PATH "/cvmfs/cms.cern.ch/.*/lcg"

    remove_path_matching PYTHONPATH "/cvmfs/cms.cern.ch/.*/root"
    remove_path_matching PATH "/cvmfs/cms.cern.ch/.*/root"
    remove_path_matching LD_LIBRARY_PATH "/cvmfs/cms.cern.ch/.*/root"
    remove_path_matching ROOT_INCLUDE_PATH "/cvmfs/cms.cern.ch/.*/root"
    remove_path_matching CMAKE_PREFIX_PATH "/cvmfs/cms.cern.ch/.*/root"

    # remove local conda
    remove_path_matching PATH "/(mini|ana)conda3?/"
    remove_path_matching PYTHONPATH "/(mini|ana)conda3?/"
    remove_path_matching LD_LIBRARY_PATH "/(mini|ana)conda3?/"

    # IMPORTANT:
    # do NOT add /usr/lib64/python3.9/site-packages to PYTHONPATH.
    # It contains system cppyy/ROOT and crashes against CMSSW ROOT.
}

run_cmd() {
    echo "> $*"
    "$@"
    local rc=$?

    if [[ "$rc" != "0" ]]; then
        echo "[ERROR] Command failed:"
        echo "        $*"
        return "$rc"
    fi

    return 0
}

clean_env

source /cvmfs/cms.cern.ch/cmsset_default.sh

SCRAM_ARCH="$(get_scram_arch "$REQUESTED_CMSSW_VERSION")"
export SCRAM_ARCH
export CMSSW_AREA="${ANALYSIS_SOFT_PATH}/${REQUESTED_CMSSW_VERSION}"

if [[ "$FORCE_REINSTALL" == "1" && -d "$CMSSW_AREA" ]]; then
    echo ">>> Removing existing CMSSW area:"
    echo ">>>   $CMSSW_AREA"
    rm -rf "$CMSSW_AREA"
fi

if [[ ! -d "$CMSSW_AREA/src" ]]; then
    echo ">>> Creating CMSSW area"
    echo ">>>   CMSSW_VERSION = $REQUESTED_CMSSW_VERSION"
    echo ">>>   SCRAM_ARCH    = $SCRAM_ARCH"
    echo ">>>   path          = $CMSSW_AREA"

    mkdir -p "$ANALYSIS_SOFT_PATH"
    cd "$ANALYSIS_SOFT_PATH" || return 1

    run_cmd scramv1 project CMSSW "$REQUESTED_CMSSW_VERSION" || return 1
fi

cd "$CMSSW_AREA/src" || return 1
eval "$(scramv1 runtime -sh)"

prepend_path PYTHONPATH "$ANALYSIS_PATH"

unalias python 2>/dev/null || true
alias python=python3 2>/dev/null || true

# Combine is provided by the official pre-built CMS container on CVMFS.  Keep
# it as a shell function so `source env.sh` works in both bash and zsh and does
# not require a separate local CombinedLimit build.
export COMBINE_IMAGE="${COMBINE_IMAGE:-/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/combine-container:latest}"
export COMBINE_EXEC_SCRIPT="/cvmfs/cms.cern.ch/cat/combine_alias/exec_script.sh"

# Remove aliases left by /cvmfs/cms.cern.ch/cat/combine_env.sh: they do not
# necessarily bind the active Kerberos cache and can make EOS look absent.
unalias combine 2>/dev/null || true
run_combine_tool() {
    local combine_tool="$1"
    shift
    local -a combine_apptainer_options
    combine_apptainer_options=(-s exec --bind /cvmfs,/eos,/afs)
    if [[ "${KRB5CCNAME:-}" == FILE:* ]]; then
        local combine_ticket_path="${KRB5CCNAME#FILE:}"
        local combine_ticket_dir
        combine_ticket_dir="$(dirname "$combine_ticket_path")"
        combine_apptainer_options+=(
            --bind "${combine_ticket_dir}:${combine_ticket_dir}"
            --env "KRB5CCNAME=${KRB5CCNAME}"
        )
    fi
    apptainer "${combine_apptainer_options[@]}" \
        "$COMBINE_IMAGE" \
        "$COMBINE_EXEC_SCRIPT" \
        "$combine_tool" "$@"
}

combine() {
    run_combine_tool combine "$@"
}

combine_cards() {
    run_combine_tool combineCards.py "$@"
}

cd "$ANALYSIS_PATH" || return 1

# Storage locations (config/paths.yaml) as HMUMU_<KEY>.  Drop the ones a
# previous `source env.sh` exported, so switching checkout refreshes them;
# values set by hand are kept.
eval "unset ${HMUMU_EXPORTED_KEYS:-}"
eval "$(python3 "$ANALYSIS_PATH/common/paths.py" --shell)"

echo "=========================================================="
echo "[H_mumu env]"
echo "ANALYSIS_PATH   = $ANALYSIS_PATH"
echo "CMSSW_BASE      = $CMSSW_BASE"
echo "SCRAM_ARCH      = $SCRAM_ARCH"
echo "python3         = $(which python3)"
echo "root            = $(which root 2>/dev/null || echo not_found)"
echo "root-config     = $(which root-config 2>/dev/null || echo not_found)"
echo "combine         = official CVMFS container (${COMBINE_IMAGE})"
echo "PYTHONPATH      = $PYTHONPATH"
echo "=========================================================="

python3 - <<'PY'
import sys
print("Python =", sys.version.split()[0])

try:
    import ROOT
    print("ROOT =", ROOT.gROOT.GetVersion())
except Exception as e:
    print("[ERROR] ROOT import failed:", repr(e))
    raise SystemExit(1)

try:
    import correctionlib
    print("correctionlib =", correctionlib.__version__)
except Exception as e:
    print("[WARNING] correctionlib import failed:", repr(e))
    print("          If you need correctionlib without pip, use a CMSSW/LCG environment that already provides it.")

try:
    import htcondor
    print("htcondor =", htcondor.version())
except Exception as e:
    print("[WARNING] htcondor import failed:", repr(e))
    print("          Do not add /usr/lib64/python3.9/site-packages globally, it breaks ROOT.")
PY

echo "=========================================================="
echo ">>> Environment ready."
echo "=========================================================="
