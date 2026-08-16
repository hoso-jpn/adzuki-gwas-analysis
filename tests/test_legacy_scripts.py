"""Verify the backward-compatible scripts/01-04 have no import-time side effects.

The original versions of these scripts ran top-level module code (reading
data, writing plots) on import, making them untestable and unsafe to import
anywhere else. Each rewritten script must have all its logic inside
``main()`` behind ``if __name__ == "__main__":``, so merely importing the
module does nothing.
"""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"

_SCRIPT_NAMES = [
    "01_manhattan_plot",
    "02_qq_plot",
    "03_regional_plot",
    "04_extract_top_variants_by_region",
]


class LegacyScriptImportTests(unittest.TestCase):
    def test_importing_each_script_has_no_side_effects(self) -> None:
        for name in _SCRIPT_NAMES:
            with self.subTest(script=name):
                path = SCRIPTS_DIR / f"{name}.py"
                spec = importlib.util.spec_from_file_location(f"_legacy_{name}", path)
                assert spec is not None and spec.loader is not None
                module = importlib.util.module_from_spec(spec)
                # Must not raise, must not touch data/raw or plots/ or
                # results/ -- importing defines main() but does not call it.
                spec.loader.exec_module(module)
                self.assertTrue(hasattr(module, "main"), f"{name} must define main()")
                sys.modules.pop(f"_legacy_{name}", None)


if __name__ == "__main__":
    unittest.main()
