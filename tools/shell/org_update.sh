#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/emseq/emseq.org
# Author:  Jeff Szymanski
# Tangled: 2026-03-16 12:05:56
# ============================================================

set -euo pipefail

source ~/repos/basecamp/lib/basecamp_functions.sh

REPOS=(emseq)

# ---------- config ----------
REPO_DIR=~/repos/emseq
README_ORG="$REPO_DIR/emseq.org"
README_NODE=bda70cff-0713-4e32-8da1-ee83924b8f00
README_EXPORT=~/repos/emacs/scripts/emacs_export_header_to_markdown.py
# ----------------------------

save_in_emacs() {
  if command -v emacsclient >/dev/null 2>&1 && emacsclient -e "(progn t)" >/dev/null 2>&1; then
    emacsclient -e "(save-some-buffers t)" >/dev/null
    echo "💾 saved buffers via emacsclient"
  else
    echo "ℹ️ no Emacs server; skipping save"
  fi
}

update_readme() {
  echo "📝 exporting README.md from $README_ORG"
  pushd "$REPO_DIR" >/dev/null
  python3 "$README_EXPORT" --org_file "$README_ORG" --node_id "$README_NODE"
  popd >/dev/null
}

tangle_repo_org() {
  local repo=$1
  local org_file=~/repos/"$repo"/"$repo".org

  echo "=== $repo ==="
  if [[ -f "$org_file" ]]; then
    echo "⏳ tangling $repo..."
    tangle "$org_file" || { echo "❌ tangle failed for $repo"; return 1; }
  else
    echo "⚠️  no $repo.org, skipping"
  fi
}

git_update_repo() {
  local repo=$1
  echo "🔄 updating $repo..."
  pushd ~/repos/"$repo" >/dev/null || { echo "❌ cd failed for $repo"; return 1; }

  if ! output=$(git_wkflow_up 2>&1); then
    echo "❌ git_wkflow_up failed in $repo"
    echo "$output"
  else
    echo -e "$(date)\n$output"
  fi

  popd >/dev/null
}

main() {
  save_in_emacs

  for repo in "${REPOS[@]}"; do
    tangle_repo_org "$repo" || true
  done

  update_readme
  save_in_emacs

  for repo in "${REPOS[@]}"; do
    git_update_repo "$repo" || true
  done
}

main "$@"
