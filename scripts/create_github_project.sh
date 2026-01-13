#!/usr/bin/env bash
set -euo pipefail

# Script to create GitHub Project Board (v2) and populate with issues
# Requires: gh CLI with 'project' scope

OWNER="joewa6"
REPO="MSE_LAB"
PROJECT_TITLE="MSE_LAB Bootstrap"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "🎯 Creating GitHub Project Board..."
echo ""

# Check gh CLI
if ! command -v gh &> /dev/null; then
    echo "❌ ERROR: gh CLI not found"
    exit 1
fi

# Check auth and scopes
echo "Checking authentication..."
if ! gh auth status 2>&1 | grep -q "Logged in"; then
    echo "❌ ERROR: Not logged in to GitHub"
    echo "Run: gh auth login"
    exit 1
fi

# Check for project scope
if ! gh auth status 2>&1 | grep -q "'project'"; then
    echo "⚠️  Missing 'project' scope. Refreshing token..."
    if gh auth refresh -s project; then
        echo "✅ Token refreshed with project scope"
    else
        echo "❌ Failed to refresh token"
        echo "Please run: gh auth refresh -s project"
        exit 1
    fi
fi

echo ""

# Create project
echo "Creating project: $PROJECT_TITLE"
if PROJECT_URL=$(gh project create --owner "$OWNER" --title "$PROJECT_TITLE" --format json | jq -r '.url' 2>/dev/null); then
    echo "✅ Project created: $PROJECT_URL"
    PROJECT_NUMBER=$(echo "$PROJECT_URL" | grep -oP '\d+$')
    echo "   Project number: $PROJECT_NUMBER"
else
    echo "❌ Failed to create project"
    echo ""
    echo "Manual fallback:"
    echo "1. Go to https://github.com/$OWNER?tab=projects"
    echo "2. Click 'New project'"
    echo "3. Choose 'Board' template"
    echo "4. Name it: $PROJECT_TITLE"
    exit 1
fi

echo ""

# Get list of all issues
echo "Fetching issues..."
ISSUE_COUNT=$(gh issue list --limit 200 --json number -q 'length')
echo "Found $ISSUE_COUNT issues"

if [ "$ISSUE_COUNT" -eq 0 ]; then
    echo "⚠️  No issues found. Run ./scripts/create_github_issues.sh first"
    exit 1
fi

echo ""

# Add issues to project
echo "Adding issues to project..."
added=0
failed=0

for issue_num in $(gh issue list --limit 200 --json number -q '.[].number' | sort -n); do
    issue_url="https://github.com/$OWNER/$REPO/issues/$issue_num"
    
    if gh project item-add --owner "$OWNER" --number "$PROJECT_NUMBER" --url "$issue_url" > /dev/null 2>&1; then
        echo "  ✅ Added issue #$issue_num"
        added=$((added + 1))
    else
        echo "  ⚠️  Failed to add issue #$issue_num"
        failed=$((failed + 1))
    fi
    
    sleep 0.5  # Rate limit protection
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Added $added/$ISSUE_COUNT issues to project"

if [ $failed -gt 0 ]; then
    echo "⚠️  Failed to add $failed issues (may already exist)"
fi

echo ""
echo "📊 Project Board: $PROJECT_URL"
echo ""
echo "Next steps:"
echo "  1. Open project board in browser"
echo "  2. Configure 'Status' field with columns:"
echo "     - Backlog"
echo "     - Day 1"
echo "     - Day 2"
echo "     - Day 3"
echo "     - Review"
echo "     - Done"
echo "  3. Assign issues to columns:"
echo "     Day 1: issues #1-#5"
echo "     Day 2: issues #6-#10"
echo "     Day 3: issues #11-#16"
echo ""
echo "Or use the web UI to bulk-assign statuses."
echo ""
