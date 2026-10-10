"""Top-level `captnseq-toolkit` command: dispatches to one subcommand per tool."""

import sys

from captnseq_toolkit import __version__, essentiality, sam2sites, sites2genes

SUBCOMMANDS = {
    "sam2sites": (sam2sites.main, "SAM -> insertion sites (BED + bedGraph tracks)"),
    "sites2genes": (sites2genes.main, "insertion sites -> per-gene insertion counts / index"),
    "essentiality": (essentiality.main, "per-gene counts -> essential / ambiguous / non-essential calls"),
}


def _usage():
    lines = [f"captnseq-toolkit {__version__}", "", "usage: captnseq-toolkit <subcommand> [options]", "", "subcommands:"]
    lines += [f"  {name:<14}{desc}" for name, (_, desc) in SUBCOMMANDS.items()]
    lines += ["", "Run `captnseq-toolkit <subcommand> --help` for subcommand options."]
    return "\n".join(lines)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(_usage())
        return 0
    if argv[0] in ("-V", "--version"):
        print(f"captnseq-toolkit {__version__}")
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd not in SUBCOMMANDS:
        print(f"captnseq-toolkit: unknown subcommand '{cmd}'\n\n{_usage()}", file=sys.stderr)
        return 2
    SUBCOMMANDS[cmd][0](rest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
