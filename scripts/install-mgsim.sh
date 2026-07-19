#!/usr/bin/env bash
# Install the bundled MGSIM fork (third_party/MGSIM) into a conda env.
#
# MGSIM's read simulators (art, simlord, nanosim-h) are bioconda packages, not
# pip-installable, so the streamlined path is: init the submodule, create its
# pinned conda environment, then pip-install MGSIM into that env.
#
# Usage: scripts/install-mgsim.sh [env_name]   (default env: mgsim)
set -euo pipefail

env_name="${1:-mgsim}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mgsim_dir="$repo_root/third_party/MGSIM"

# 1. Make sure the submodule is checked out.
if [ ! -f "$mgsim_dir/environment.yml" ]; then
    echo "==> Fetching the MGSIM submodule ..."
    git -C "$repo_root" submodule update --init --recursive third_party/MGSIM
fi

# 2. Pick a solver: mamba is much faster than conda when available.
solver=conda
if command -v mamba >/dev/null 2>&1; then
    solver=mamba
elif ! command -v conda >/dev/null 2>&1; then
    echo "error: neither mamba nor conda found on PATH." >&2
    exit 1
fi

# 2b. On Apple Silicon, MGSIM's pinned bioconda tools (art, simlord, nanosim-h,
#     pyfastx, fqtools) have no osx-arm64 build, so install the osx-64 builds and
#     run them under Rosetta 2. CONDA_SUBDIR drives the solve to osx-64.
if [ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ]; then
    echo "==> Apple Silicon detected: installing osx-64 builds (run under Rosetta 2)."
    if ! /usr/bin/pgrep -q oahd 2>/dev/null; then
        echo "    note: Rosetta 2 may be required. Install with: softwareupdate --install-rosetta --agree-to-license" >&2
    fi
    export CONDA_SUBDIR=osx-64
fi

# 3. Create (or update) the conda environment from the pinned environment.yml,
#    patched two ways:
#    - add 'nodefaults' to the channels so conda never consults Anaconda's default
#      channel (repo.anaconda.com/pkgs/main), which is gated behind a
#      Terms-of-Service prompt (everything MGSIM needs is on conda-forge/bioconda);
#    - add 'pip' to the dependencies, since the upstream environment.yml omits it
#      and we need it to install MGSIM itself into the env.
env_yml="$(mktemp -t mgsim-env.XXXXXX)"
trap 'rm -f "$env_yml"' EXIT
have_nd=0; grep -qE '^[[:space:]]*-[[:space:]]*nodefaults[[:space:]]*$' "$mgsim_dir/environment.yml" && have_nd=1
have_pip=0; grep -qE '^[[:space:]]*-[[:space:]]*pip( |=|>|<|$)' "$mgsim_dir/environment.yml" && have_pip=1
awk -v nd="$have_nd" -v pip="$have_pip" '
    { print }
    /^channels:/     && !cd && nd  == 0 { print "  - nodefaults"; cd = 1 }
    /^dependencies:/ && !dd && pip == 0 { print "  - pip";        dd = 1 }
' "$mgsim_dir/environment.yml" > "$env_yml"

if conda env list | awk '{print $1}' | grep -qx "$env_name"; then
    echo "==> Updating existing conda env '$env_name' ..."
    "$solver" env update -n "$env_name" -f "$env_yml"
else
    echo "==> Creating conda env '$env_name' ..."
    "$solver" env create -n "$env_name" -f "$env_yml"
fi

# Pin the env to the same subdir so later `conda install`s into it stay consistent.
if [ -n "${CONDA_SUBDIR:-}" ]; then
    conda run -n "$env_name" conda config --env --set subdir "$CONDA_SUBDIR" || true
fi

# 4. Install MGSIM itself into that env (python -m pip so it works even if the
#    'pip' launcher isn't on the env's PATH).
echo "==> Installing MGSIM into '$env_name' ..."
conda run -n "$env_name" python -m pip install "$mgsim_dir"

echo
echo "Done. Activate the environment and check MGSIM with:"
echo "    conda activate $env_name"
echo "    MGSIM --list"
