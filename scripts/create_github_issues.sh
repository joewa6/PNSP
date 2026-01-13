#!/usr/bin/env bash
set -euo pipefail

# Script to create all GitHub issues from docs/PROJECT_ISSUES.md
# Uses gh CLI to automate issue creation with proper labels and bodies

BLUEPRINT="docs/PROJECT_ISSUES.md"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$REPO_ROOT"

if ! command -v gh &> /dev/null; then
    echo "ERROR: gh CLI not found. Install from https://cli.github.com/"
    exit 1
fi

if [ ! -f "$BLUEPRINT" ]; then
    echo "ERROR: $BLUEPRINT not found"
    exit 1
fi

echo "🚀 Creating GitHub issues from $BLUEPRINT..."
echo ""

# Issue definitions: title, labels, start_marker, end_marker
declare -a ISSUES=(
    # DAY 1
    "Port Stage A builder into Orchestra|capability,stageA,day1,critical|### Issue #1:|### Issue #2:"
    "Port Stage B local simulator into Orchestra|capability,stageB,day1,critical|### Issue #2:|### Issue #3:"
    "Create Stage C pointer-based archival|capability,stageC,day1,critical|### Issue #3:|### Issue #4:"
    "End-to-end smoke test (A → B → C)|validation,day1,critical|### Issue #4:|### Issue #5:"
    "Create run_manifest.json schema definition|docs,day1|### Issue #5:|## DAY 2 ISSUES"
    
    # DAY 2
    "Stage A batch mode (read test.csv)|capability,stageA,day2|### Issue #6:|### Issue #7:"
    "Stage C legacy WD8 indexer|capability,stageC,day2,critical|### Issue #7:|### Issue #8:"
    "Stage B rsync archival after simulation|capability,stageB,day2|### Issue #8:|### Issue #9:"
    "Composer ledger scan (multi-archive support)|capability,composer,day2,critical|### Issue #9:|### Issue #10:"
    "Composer origin column logic|capability,composer,day2|### Issue #10:|## DAY 3 ISSUES"
    
    # DAY 3
    "Stage D metric extraction (first reducer)|capability,stageD,day3,critical|### Issue #11:|### Issue #12:"
    "Composer rank table from Stage D metrics|capability,composer,day3|### Issue #12:|### Issue #13:"
    "Stage B HPC skeleton|capability,stageB,day3|### Issue #13:|### Issue #14:"
    "Repository README with quickstart|docs,day3|### Issue #14:|### Issue #15:"
    "Pin Python environment (environment.yml)|infra,day3|### Issue #15:|### Issue #16:"
    "GitHub branch protection rules|infra,day3|### Issue #16:|## 🏷️ LABEL TAXONOMY"
)

issue_count=0
success_count=0

for issue_spec in "${ISSUES[@]}"; do
    IFS='|' read -r title labels start_marker end_marker <<< "$issue_spec"
    issue_count=$((issue_count + 1))
    
    echo "Creating issue #$issue_count: $title"
    echo "  Labels: $labels"
    
    # Extract issue body from blueprint
    # Find content between start_marker and end_marker, excluding the markers themselves
    body=$(awk "/$start_marker/,/$end_marker/" "$BLUEPRINT" | sed '1d;$d')
    
    if [ -z "$body" ]; then
        echo "  ⚠️  WARNING: Empty body for issue #$issue_count"
        continue
    fi
    
    # Create the issue
    if gh issue create \
        --title "$title" \
        --label "$labels" \
        --body "$body" > /dev/null 2>&1; then
        echo "  ✅ Created"
        success_count=$((success_count + 1))
    else
        echo "  ❌ FAILED"
    fi
    
    echo ""
    
    # Rate limit protection (GitHub API has rate limits)
    sleep 1
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Created $success_count/$issue_count issues"
echo ""
echo "View issues: https://github.com/joewa6/MSE_LAB/issues"
echo ""
echo "Next steps:"
echo "  1. Create GitHub Project Board: https://github.com/joewa6/MSE_LAB/projects"
echo "  2. Add columns: Backlog, Day 1, Day 2, Day 3, Review, Done"
echo "  3. Assign issues to columns per docs/PROJECT_ISSUES.md"
echo ""
