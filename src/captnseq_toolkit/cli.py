"""Top-level `captnseq-toolkit` command: one argparse subcommand per tool."""

import argparse
import sys

from captnseq_toolkit import __version__, essentiality, sam2sites, sites2genes

# name -> (module, one-line help). Each module provides add_arguments(parser) and run_args(args).
SUBCOMMANDS = {
    "sam2sites": (sam2sites, "SAM -> insertion sites (BED + bedGraph tracks)"),
    "sites2genes": (sites2genes, "insertion sites -> per-gene insertion counts / index"),
    "essentiality": (essentiality, "per-gene counts -> essential / ambiguous / non-essential calls"),
}


def build_parser():
    parser = argparse.ArgumentParser(
        prog="captnseq-toolkit",
        description="Analysis tools for CAPTn-seq transposon-insertion sequencing data.",
    )
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="<subcommand>", required=True)
    for name, (module, help_text) in SUBCOMMANDS.items():
        sub = subparsers.add_parser(
            name,
            help=help_text,
            description=module.__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        module.add_arguments(sub)
        sub.set_defaults(func=module.run_args)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
