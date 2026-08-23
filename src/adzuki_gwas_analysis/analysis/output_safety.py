"""Shared "is this output directory safe to write into" guard.

Extracted from :mod:`adzuki_gwas_analysis.analysis.batch` (Issue #9/PR #16) so
:mod:`adzuki_gwas_analysis.analysis.report` (Issue #11) can apply the exact same
never-merge-into-an-existing-directory contract to its own ``--output-dir`` without a second,
drifting copy of the same logic.
"""

from __future__ import annotations

from pathlib import Path

from adzuki_gwas_analysis.errors import BatchOutputDirectoryUnsafeError


def check_output_dir_is_safe(output_dir: Path) -> None:
    """Fail fast unless ``output_dir`` is a plain, empty (or not-yet-existing) directory.

    Never silently reused, merged into, or replaced: every caller of this function
    publishes its full output as one all-or-nothing unit, so anything already there --
    or anything this check cannot positively confirm is an ordinary real directory -- is
    rejected before any work is done.
    """
    if output_dir.is_symlink():
        raise BatchOutputDirectoryUnsafeError(
            output_dir=str(output_dir), reason="path is a symlink, not a plain directory"
        )
    if output_dir.exists():
        if not output_dir.is_dir():
            raise BatchOutputDirectoryUnsafeError(
                output_dir=str(output_dir), reason="path exists and is not a directory"
            )
        if any(output_dir.iterdir()):
            raise BatchOutputDirectoryUnsafeError(
                output_dir=str(output_dir),
                reason=(
                    "directory already exists and is not empty -- output is published as "
                    "one all-or-nothing unit and never merges into an existing directory"
                ),
            )
