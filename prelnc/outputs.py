from __future__ import absolute_import

"""Output path resolution and tabular writers.

Every table PreLnc writes goes through this module so that:

* the delimiter always matches the file extension (``.csv`` -> comma,
  ``.tsv``/``.txt`` -> tab);
* an existing file is never overwritten silently (a numbered suffix is added
  instead, unless the caller explicitly asks to replace it);
* files are written to a temporary sibling and moved into place, so an
  interrupted run cannot leave a half-written table behind.
"""

import contextlib
import csv
import os
import tempfile
import time

TSV_EXTENSIONS = (".tsv", ".tab", ".txt")
CSV_EXTENSIONS = (".csv",)
SUPPORTED_EXTENSIONS = CSV_EXTENSIONS + TSV_EXTENSIONS

DEFAULT_EXTENSION = ".tsv"


class OutputError(RuntimeError):
    """Raised when an output path cannot be used."""


def delimiter_for(path):
    """Return the delimiter implied by ``path``'s extension."""
    extension = os.path.splitext(path)[1].lower()
    if extension in CSV_EXTENSIONS:
        return ","
    if extension in TSV_EXTENSIONS:
        return "\t"
    raise OutputError(
        "Unsupported output extension %r for %s. Use .csv for comma separated "
        "or .tsv/.txt for tab separated output." % (extension or "(none)", path)
    )


def plural(count, word, plural_form=None):
    """Return ``"3 sequences"`` / ``"1 sequence"``."""
    if count == 1:
        return "%d %s" % (count, word)
    return "%d %s" % (count, plural_form or word + "s")


def describe_table(rows, columns):
    """Return ``"100 rows x 3 columns"`` for use in progress summaries."""
    return "%s x %s" % (plural(rows, "row"), plural(columns, "column"))


def timestamped_stem(prefix):
    """Build a default file stem such as ``prelnc_human_prediction_20261105_132500``."""
    return "%s_%s" % (prefix, time.strftime("%Y%m%d_%H%M%S"))


def _ensure_directory(path):
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        raise OutputError("Output directory does not exist: %s" % directory)
    if not os.access(directory, os.W_OK):
        raise OutputError("Output directory is not writable: %s" % directory)
    return directory


def _numbered(path, index):
    """Return ``path`` with a ``_<index>`` suffix inserted before the extension."""
    if index <= 0:
        return path
    stem, extension = os.path.splitext(path)
    return "%s_%d%s" % (stem, index, extension)


def _candidates(output, companions, index):
    main = _numbered(output, index)
    stem = os.path.splitext(main)[0]
    return [main] + [stem + suffix for suffix in companions]


def resolve_outputs(
    output=None,
    default_stem=None,
    companions=(),
    force=False,
    cwd=None,
    table=True,
    default_extension=DEFAULT_EXTENSION,
):
    """Decide the final paths for one run.

    ``output`` may be ``None``, in which case a timestamped name is generated in
    the current working directory (or ``cwd``). ``companions`` are suffixes such
    as ``".features.tsv"`` that share the output stem.

    ``table=False`` skips the tabular-extension check, for outputs such as a
    pickled model that are not tables themselves.

    Returns a dict with ``main``, ``companions`` (list of paths), ``requested``
    and ``renamed`` (True when a suffix had to be appended).
    """
    if output is None:
        if not default_stem:
            raise OutputError("Either an output path or a default stem is required.")
        cwd = cwd or os.getcwd()
        output = os.path.join(cwd, timestamped_stem(default_stem) + default_extension)

    output = os.path.abspath(os.path.expanduser(output))
    if table:
        delimiter_for(output)
    _ensure_directory(output)

    if force:
        index = 0
    else:
        index = 0
        while any(os.path.exists(path) for path in _candidates(output, companions, index)):
            index += 1
            if index > 9999:
                raise OutputError(
                    "Could not find a free file name near %s after 10000 attempts." % output
                )

    paths = _candidates(output, companions, index)
    return {
        "main": paths[0],
        "companions": paths[1:],
        "requested": output,
        "renamed": index > 0,
    }


@contextlib.contextmanager
def atomic_open(path, mode="w", encoding="utf-8", newline=""):
    """Write ``path`` through a temporary file and move it into place on success."""
    directory = _ensure_directory(path)
    handle, temp_path = tempfile.mkstemp(prefix=".prelnc-", suffix=".tmp", dir=directory)
    try:
        if "b" in mode:
            stream = os.fdopen(handle, mode)
        else:
            stream = os.fdopen(handle, mode, encoding=encoding, newline=newline)
        with stream:
            yield stream
        os.chmod(temp_path, _default_file_mode())
        os.replace(temp_path, path)
        temp_path = None
    finally:
        if temp_path is not None and os.path.exists(temp_path):
            os.unlink(temp_path)


def _default_file_mode():
    """Respect the process umask so outputs are not always 0600."""
    current = os.umask(0)
    os.umask(current)
    return 0o666 & ~current


def write_table(path, header, rows, delimiter=None):
    """Write a header plus rows to ``path`` in the delimiter its extension implies."""
    delimiter = delimiter or delimiter_for(path)
    with atomic_open(path, "w", newline="") as stream:
        writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
        writer.writerow(list(header))
        for row in rows:
            writer.writerow(list(row))
    return path
