#!/usr/bin/env bash
# setup_repos.sh — Clone and build external tool repositories required by emseq analysis module.
# Hand-crafted; not tangled from org.
set -euo pipefail

REPOS_DIR="${REPOS_DIR:-$HOME/repos}"

# --- Repository definitions ---
declare -A REPO_URLS=(
    [mHapTools]="https://github.com/butyuhao/mHapTools"
    [wgbs_tools]="https://github.com/nloyfer/wgbs_tools"
    [UXM_deconv]="https://github.com/nloyfer/UXM_deconv"
)

# --- Helpers ---
info()  { echo "[setup] $*"; }
err()   { echo "[ERROR] $*" >&2; }

clone_or_pull() {
    local name="$1" url="$2" dest="${REPOS_DIR}/${name}"
    if [[ -d "$dest/.git" ]]; then
        info "${name}: already cloned at ${dest}, pulling latest"
        git -C "$dest" pull --ff-only || err "${name}: pull failed (resolve manually)"
    else
        info "${name}: cloning ${url} → ${dest}"
        git clone "$url" "$dest"
    fi
}

# --- Main ---
mkdir -p "$REPOS_DIR"

for name in mHapTools wgbs_tools UXM_deconv; do
    clone_or_pull "$name" "${REPO_URLS[$name]}"
done

# --- wgbs_tools: build C extensions (skip with SKIP_BUILD=1) ---
if [[ "${SKIP_BUILD:-0}" != "1" ]] && [[ -f "${REPOS_DIR}/wgbs_tools/setup.py" ]]; then
    info "wgbs_tools: running python setup.py"
    (cd "${REPOS_DIR}/wgbs_tools" && python setup.py)
elif [[ "${SKIP_BUILD:-0}" == "1" ]]; then
    info "wgbs_tools: skipping build (SKIP_BUILD=1, build manually in conda env)"
else
    err "wgbs_tools: setup.py not found, skipping build"
fi

info "Done. Repos installed to ${REPOS_DIR}/"
