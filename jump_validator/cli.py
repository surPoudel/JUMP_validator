"""Command-line interface for JUMP Validator."""

from __future__ import annotations

import argparse
import sys

from .validator import run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="jump-validator",
        description="Validate spectra and generate a manual validation report (Excel + plots).",
    )
    parser.add_argument(
        "params_file",
        help="Path to the INI-style .params file (see jump_validator.params for an example).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    run(args.params_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
