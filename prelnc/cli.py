from __future__ import absolute_import

import argparse
import os
import sys
import textwrap

from prelnc import __version__
from prelnc.config import available_species, resolve_profile
from prelnc.hexamer import build_hexamer_table, write_hexamer_file, write_hexamer_table
from prelnc.outputs import OutputError, delimiter_for, describe_table, resolve_outputs
from prelnc.pipelines import run_prediction_pipeline, run_training_pipeline
from prelnc.progress import PREFIX, Progress, info, print_summary

VERSION_STRING = "PreLnc %s" % __version__

MACOS_PRIVACY_HINT = (
    "macOS blocks Documents, Desktop and Downloads until the app running this command "
    "is allowed to read them. Grant access under System Settings > Privacy & Security > "
    "Files and Folders (or Full Disk Access), then run the command again. To skip that "
    "step, copy the input to a folder outside those three, for example /tmp or your "
    "home directory."
)


def _permission_message(error):
    """Explain a denied read/write, which on macOS is usually a privacy setting."""
    target = getattr(error, "filename", None)
    if target:
        headline = "Permission denied: %s" % os.path.abspath(target)
    else:
        headline = "Permission denied."
    lines = ["%s ERROR: %s" % (PREFIX, headline)]
    for line in textwrap.wrap(MACOS_PRIVACY_HINT, width=70):
        lines.append("%s %s" % (PREFIX, line))
    return "\n".join(lines)


def _guard(func):
    """Turn Ctrl-C and expected errors into a short message and a clean exit code."""

    def wrapper(argv=None):
        try:
            return func(argv)
        except KeyboardInterrupt:
            sys.stderr.write("\n%s Interrupted. No files were written.\n" % PREFIX)
            return 130
        except PermissionError as error:
            # A denied read is almost always a macOS privacy setting, not a bug.
            sys.stderr.write(_permission_message(error) + "\n")
            return 1
        except (OutputError, ValueError, RuntimeError, OSError) as error:
            sys.stderr.write("%s ERROR: %s\n" % (PREFIX, error))
            return 1

    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def _add_common_output_arguments(parser, default_hint):
    parser.add_argument(
        "-o",
        "--output",
        help="Output file (.tsv or .csv). Default: %s in the current directory" % default_hint,
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists (default: write a numbered copy)",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Hide progress bars (the final summary is still printed)",
    )


def _validate_output(path):
    """Fail on an unsupported extension before spending time on the computation."""
    if path:
        delimiter_for(path)
    return path


def build_predict_parser():
    parser = argparse.ArgumentParser(
        prog="prelnc",
        description="Predict lncRNA probability from transcript FASTA.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    _add_common_output_arguments(parser, "./prelnc_<species>_prediction_<timestamp>.tsv")
    parser.add_argument("--species", choices=available_species(), help="Built-in species profile")
    parser.add_argument("--type", choices=["animal", "plant"], dest="feature_type", help="Feature set for custom runs")
    parser.add_argument("-m", "--model", help="Custom model path")
    parser.add_argument("-r", "--hexamer", help="Custom hexamer table path")
    parser.add_argument("--version", action="version", version=VERSION_STRING)
    return parser


def build_train_parser():
    parser = argparse.ArgumentParser(
        prog="prelnc-train",
        description="Train a PreLnc model from lncRNA and protein-coding FASTA files.",
    )
    parser.add_argument("-l", "--lnc", required=True, help="lncRNA FASTA file")
    parser.add_argument("-p", "--pc", required=True, help="protein-coding FASTA file")
    _add_common_output_arguments(parser, "./prelnc_<species>_model_<timestamp>.pkl")
    parser.add_argument("--species", choices=available_species(), help="Built-in species profile")
    parser.add_argument("--type", choices=["animal", "plant"], dest="feature_type", help="Feature set for custom runs")
    parser.add_argument("-r", "--hexamer", help="Hexamer table path")
    parser.add_argument("--version", action="version", version=VERSION_STRING)
    return parser


def build_hexamer_parser():
    parser = argparse.ArgumentParser(
        prog="prelnc-hexamer",
        description="Build a hexamer frequency table for a new species.",
    )
    parser.add_argument(
        "-c", "--coding", required=True,
        help="Coding sequences in FASTA format (CDS only, i.e. start codon to stop codon, no UTR)",
    )
    parser.add_argument("-n", "--noncoding", required=True, help="Non-coding sequences in FASTA format")
    parser.add_argument(
        "-o",
        "--output",
        help="Output table (.tsv or .csv). Omit to write to stdout",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists (default: write a numbered copy)",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Hide progress bars (the final summary is still printed)",
    )
    parser.add_argument("--version", action="version", version=VERSION_STRING)
    return parser


@_guard
def predict_main(argv=None):
    parser = build_predict_parser()
    args = parser.parse_args(argv)
    _validate_output(args.output)
    profile = resolve_profile(
        species=args.species,
        feature_type=args.feature_type,
        model_path=args.model,
        hexamer_path=args.hexamer,
    )
    if not profile["model"]:
        parser.error("A model path is required when no built-in species profile is selected.")
    if not profile["hexamer"]:
        parser.error("A hexamer table is required when no built-in species profile is selected.")

    stem_key = args.species or args.feature_type
    run_prediction_pipeline(
        input_path=args.input,
        model_path=profile["model"],
        hexamer_path=profile["hexamer"],
        feature_names=profile["features"],
        output_path=args.output,
        default_stem="prelnc_%s_prediction" % stem_key,
        force=args.force,
        quiet=args.quiet,
    )
    return 0


@_guard
def train_main(argv=None):
    parser = build_train_parser()
    args = parser.parse_args(argv)
    profile = resolve_profile(
        species=args.species,
        feature_type=args.feature_type,
        hexamer_path=args.hexamer,
    )
    if not profile["hexamer"]:
        parser.error("A hexamer table is required when no built-in species profile is selected.")

    stem_key = args.species or args.feature_type
    run_training_pipeline(
        lnc_path=args.lnc,
        pc_path=args.pc,
        hexamer_path=profile["hexamer"],
        feature_names=profile["features"],
        output_path=args.output,
        default_stem="prelnc_%s_model" % stem_key,
        force=args.force,
        quiet=args.quiet,
    )
    return 0


@_guard
def hexamer_main(argv=None):
    parser = build_hexamer_parser()
    args = parser.parse_args(argv)
    _validate_output(args.output)

    def progress_factory(label, total):
        return Progress(total, label, enabled=not args.quiet)

    rows = build_hexamer_table(
        args.coding, args.noncoding, progress_factory=progress_factory
    )

    if not args.output:
        write_hexamer_table(sys.stdout, rows)
        print_summary(
            [
                ("Output", "stdout (%s)" % describe_table(len(rows), 3)),
                ("Hint", "pass -o table.tsv to write the table to a file"),
            ],
            title="Hexamer table written",
        )
        return 0

    plan = resolve_outputs(args.output, companions=(), force=args.force)
    if plan["renamed"]:
        info("%s already exists, writing to %s instead" % (plan["requested"], plan["main"]))
    write_hexamer_file(plan["main"], rows)
    print_summary(
        [
            ("Coding input", "%s" % args.coding),
            ("Noncoding input", "%s" % args.noncoding),
            ("Output", "%s (%s)" % (plan["main"], describe_table(len(rows), 3))),
        ],
        title="Hexamer table written",
    )
    return 0


def main(argv=None):
    return predict_main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
