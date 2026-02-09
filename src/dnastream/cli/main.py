from __future__ import annotations

import argparse

from . import create


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dnastream",
        description="DNAStream command line interface",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)
    create.add_parser(sub)
    return p


def main(argv=None) -> int:
    p = build_parser()
    args = p.parse_args(argv)
    return args.func(args)
