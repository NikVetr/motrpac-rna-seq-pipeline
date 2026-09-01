from __future__ import annotations

import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import validate_jsons  # noqa: E402


class ValidateJsonsTests(unittest.TestCase):
    def test_comparison_reads_both_files_and_returns_meaningful_status(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            left = root / "left.json"
            right = root / "right.json"
            left.write_text(json.dumps({"value": 1}), encoding="utf-8")
            right.write_text(json.dumps({"value": 1}), encoding="utf-8")
            with redirect_stdout(StringIO()):
                self.assertEqual(0, validate_jsons.main([str(left), str(right)]))

            right.write_text(json.dumps({"value": 2}), encoding="utf-8")
            with redirect_stdout(StringIO()):
                self.assertEqual(1, validate_jsons.main([str(left), str(right)]))

            right.write_text("not-json", encoding="utf-8")
            with redirect_stderr(StringIO()):
                self.assertEqual(2, validate_jsons.main([str(left), str(right)]))


if __name__ == "__main__":
    unittest.main()
