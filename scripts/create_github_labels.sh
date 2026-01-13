#!/usr/bin/env bash
set -euo pipefail

# Script to create all required GitHub labels for MSE_LAB
# Run this BEFORE creating issues

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if ! command -v gh &> /dev/null; then
    echo "ERROR: gh CLI not found. Install from https://cli.github.com/"
    exit 1
fi

echo "🏷️  Creating GitHub labels..."
echo ""

# Label definitions: name, color, description
declare -a LABELS=(
    "capability|0e8a16|Adds new lab capability (always commit-worthy)"
    "validation|1d76db|Proves existing capabilities work"
    "docs|0075ca|Documentation only"
    "infra|fbca04|GitHub/CI/environment setup"
    "stageA|c2e0c6|Builder-related"
    "stageB|bfdadc|Simulator-related"
    "stageC|f9d0c4|Archive-related"
    "stageD|fef2c0|Reducer-related"
    "composer|d4c5f9|Ledger/CLI-related"
    "day1|d93f0b|Must complete by end of Day 1"
    "day2|fbca04|Must complete by end of Day 2"
    "day3|0e8a16|Must complete by end of Day 3"
    "critical|b60205|Blocking other work"
)

success_count=0
skip_count=0
total=${#LABELS[@]}

for label_spec in "${LABELS[@]}"; do
    IFS='|' read -r name color description <<< "$label_spec"
    
    # Check if label already exists
    if gh label list | grep -q "^$name"; then
        echo "⏭️  Skipping '$name' (already exists)"
        skip_count=$((skip_count + 1))
    else
        if gh label create "$name" --color "$color" --description "$description" > /dev/null 2>&1; then
            echo "✅ Created '$name'"
            success_count=$((success_count + 1))
        else
            echo "❌ Failed to create '$name'"
        fi
    fi
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Created $success_count labels"
echo "⏭️  Skipped $skip_count labels (already existed)"
echo ""
echo "View labels: https://github.com/joewa6/MSE_LAB/labels"
echo ""
