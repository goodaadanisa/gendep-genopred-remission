#!/usr/bin/env python3
"""Compatibility wrapper for the installed ``gendep`` command-line interface."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
# Executing this file directly places ``scripts/`` first on ``sys.path``.
# Always place ``src/`` ahead of it so this wrapper cannot shadow the
# installed ``gendep`` package with its own filename.
sys.path.insert(0, str(SRC))

from gendep.cli import main  # noqa: E402


if __name__ == "__main__":
    main()
