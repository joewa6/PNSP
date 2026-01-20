# GitHub Project Board - Quick Start

## ⚡ You have two options:

### Option A: Manual (5 minutes, guaranteed to work)

Follow the step-by-step guide in:
**`BuildingPlan/GITHUB_PROJECT.plan`**

Quick summary:
1. Go to https://github.com/joewa6?tab=projects
2. Create new Board project: "MSE_LAB Bootstrap"
3. Configure Status columns: Backlog, Day 1, Day 2, Day 3, Review, Done
4. Add all issues from repo
5. Drag issues to columns per `docs/PROJECT_ISSUES.md`

**Do this.** It's simple and works.

---

### Option B: Automated (requires token refresh)

```bash
# 1. Refresh GitHub token (requires browser auth)
gh auth refresh -s project

# 2. Run automation script
./scripts/create_github_project.sh
```

**Note**: The token refresh opens a browser. If you're on a headless server or the browser auth fails, use Option A.

---

## Current state

✅ **16 issues created** with full acceptance criteria
✅ **13 labels created** (capability, stageA-D, day1-3, critical, etc.)
✅ **Automation scripts committed** (can recreate everything)

**Missing**: GitHub Project Board (do Option A or B above)

---

## What the board is for

**Execution view** of the 3-day bootstrap:

- Day 1 column: Issues #1-5 (end-to-end smoke test)
- Day 2 column: Issues #6-10 (spine integration)
- Day 3 column: Issues #11-16 (reduce + hardening)

Move issues to "Done" only when acceptance criteria are met.

---

**After board is created, start Issue #1: Port Stage A builder**
