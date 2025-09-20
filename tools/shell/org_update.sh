#!/usr/bin/env bash
set -euo pipefail

source ~/repos/basecamp/lib/basecamp_functions.sh

REPOS=(emseq)

# ---------- config ----------
README_ORG=~/repos/emseq/emseq.org
README_NODE=bda70cff-0713-4e32-8da1-ee83924b8f00
README_EXPORT=~/repos/emacs/scripts/emacs_export_header_to_markdown.py
# ----------------------------

save_in_emacs() {
  # Save all modified buffers in a running Emacs server (non-interactive).
  if command -v emacsclient >/dev/null 2>&1 && emacsclient -e "(progn t)" >/dev/null 2>&1; then
    emacsclient -e "(save-some-buffers t)" >/dev/null
    echo "💾 saved buffers via emacsclient"
  else
    echo "ℹ️ no Emacs server; skipping save"
  fi
}

update_readme() {
    python3 ~/repos/emacs/scripts/emacs_export_header_to_markdown.py --org_file ~/repos/emseq/emseq.org --node_id bda70cff-0713-4e32-8da1-ee83924b8f00
}



tangle_repo_org() {
  local repo=$1
  local org_file=~/repos/"$repo"/"$repo".org

  echo "=== $repo ==="
  if [[ -f "$org_file" ]]; then
    echo "⏳ tangling $repo..."
    tangle "$org_file" || echo "❌ tangle failed for $repo"
  else
    echo "⚠️  no $repo.org, skipping"
  fi
}

git_update_repo() {
  local repo=$1
  echo "🔄 updating $repo..."
  pushd ~/repos/"$repo" >/dev/null || { echo "❌ cd failed for $repo"; return; }

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
    tangle_repo_org "$repo"
  done

  update_readme
  # Save again in case tangling opened/modified buffers
  save_in_emacs

  for repo in "${REPOS[@]}"; do
    git_update_repo "$repo"
  done
}

main "$@"
