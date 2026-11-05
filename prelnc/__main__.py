from __future__ import absolute_import

"""``python -m prelnc`` entry point.

Subcommands mirror the installed console scripts so the package can be used from
a source checkout without installing it::

    python -m prelnc -i in.fa --species human      # same as `prelnc`
    python -m prelnc predict -i in.fa ...
    python -m prelnc train -l lnc.fa -p coding.fa ...
    python -m prelnc hexamer -c cds.fa -n nc.fa ...
"""

import sys

from prelnc.cli import hexamer_main, main, train_main

SUBCOMMANDS = {
    "predict": main,
    "train": train_main,
    "hexamer": hexamer_main,
}


def run(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in SUBCOMMANDS:
        return SUBCOMMANDS[argv[0]](argv[1:])
    return main(argv)


if __name__ == "__main__":
    sys.exit(run())
