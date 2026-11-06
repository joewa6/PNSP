"""Provenance for new pipeline executions (not historical reconstruction)."""

from datetime import datetime, timezone
from pathlib import Path
import subprocess


def execution_provenance(stage, contaminant, replica, inputs, output):
    """Record the executing code and explicit paths; unavailable Git data is null."""
    root = Path(__file__).resolve().parents[1]
    commit = dirty = None
    # Do not accidentally identify a parent repository when running a ZIP export.
    if (root / ".git").exists():
        try:
            commit = subprocess.check_output(
                ["git", "-C", str(root), "rev-parse", "HEAD"], text=True,
                stderr=subprocess.DEVNULL, timeout=5,
            ).strip()
            dirty = bool(subprocess.check_output(
                ["git", "-C", str(root), "status", "--porcelain"], text=True,
                stderr=subprocess.DEVNULL, timeout=5,
            ).strip())
        except (OSError, subprocess.SubprocessError):
            commit = dirty = None
    return {
        "schema_version": 1,
        "record_type": "pipeline_execution",
        "stage": stage,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "contaminant_id": contaminant,
        "replica_id": replica,
        "pipeline_git_commit": commit,
        "pipeline_worktree_dirty": dirty,
        "input_paths": [str(Path(path).resolve()) for path in inputs],
        "output_path": str(Path(output).resolve()),
    }
