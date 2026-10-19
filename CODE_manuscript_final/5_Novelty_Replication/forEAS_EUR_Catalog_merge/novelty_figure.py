#!/usr/bin/env python3
"""Compatibility wrapper: prefer trait-specific scripts.

  python3 Figure4A_SuppleTable6.py
  python3 SuppleFigure7_SuppleTable7.py

This entry point still accepts ``--trait-group Cont|Cate`` and delegates.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from Figure4A_SuppleTable6 import CFG as CONT_CFG
from SuppleFigure7_SuppleTable7 import CFG as CATE_CFG
from novelty_figure_common import run


def resolve_trait_group(cli_value: str | None) -> str:
    raw = cli_value if cli_value is not None else os.environ.get("TRAIT_GROUP")
    if raw is None or str(raw).strip() == "":
        raise SystemExit(
            "Required: --trait-group Cont|Cate (or set env TRAIT_GROUP).\n"
            "Prefer: python3 Figure4A_SuppleTable6.py  or  "
            "python3 SuppleFigure7_SuppleTable7.py"
        )
    value = str(raw).strip()
    if value not in ("Cont", "Cate"):
        raise SystemExit(f"--trait-group must be Cont or Cate; got {value!r}.")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--trait-group",
        choices=["Cont", "Cate"],
        default=None,
        help="Trait group Cont|Cate (required unless TRAIT_GROUP is set).",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Merged novelty .xlsx (default: trait-specific ExtendedData2_*.xlsx).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Figure output directory "
            "(default: Cont -> Figure4; Cate -> SuppleFigure7)."
        ),
    )
    args = parser.parse_args()
    trait_group = resolve_trait_group(args.trait_group)
    cfg = CONT_CFG if trait_group == "Cont" else CATE_CFG
    run(cfg, args.input, args.output_dir)


if __name__ == "__main__":
    main()
