#!/usr/bin/env bash
#set -euo pipefail

# assume tangle() is already defined
source ~/repos/basecamp/lib/basecamp_functions.sh

# Remake readme
python3 ~/repos/emacs/scripts/emacs_export_header_to_markdown.py --org_file ~/repos/emseq/emseq.org --node_id bda70cff-0713-4e32-8da1-ee83924b8f00

#repos=(mpnst emseq cfdna chaudhuri-lab)
repos=(emseq)

for repo in "${repos[@]}"; do
  echo "=== $repo ==="
  cd ~/repos/"$repo" || { echo "❌ cd fail"; continue; }

  org_file="${PWD}/${repo}.org"
  if [[ -f "$org_file" ]]; then
    echo "⏳ tangling ${repo}.org..."
    tangle "$org_file" || { echo "❌ tangle failed"; continue; }
  else
    echo "⚠️  no ${repo}.org, skipping tangle"
  fi
done

for repo in "${repos[@]}"; do
    echo "Updating $repo..."
    cd ~/repos/"$repo" || { echo "Failed to cd into $repo"; continue; }

    set +e
    output=$(git_wkflow_up 2>&1)
    status=$?
    set -e

    if [ $status -ne 0 ]; then
        echo "Error in $repo:"
        echo "$output"
        continue
    fi

    echo -e "$(date)\n$output"
done
