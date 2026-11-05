from __future__ import absolute_import

from prelnc.io_utils import count_fasta_records
from prelnc.legacy.frame_kmer import kmer_freq_file
from prelnc.outputs import write_table


WORD_SIZE = 6

CODING_STEP = 3
NONCODING_STEP = 1

HEADER = ["hexamer", "coding", "noncoding"]


class HexamerTableError(RuntimeError):
    pass


def build_hexamer_table(coding_path, noncoding_path, word_size=WORD_SIZE, progress_factory=None):
    """Return [(hexamer, coding_freq, noncoding_freq), ...].

    Coding kmers are counted every 3 bases (in-frame codons), non-coding kmers
    every base, which is what the original CPAT make_hexamer_tab.py did.

    ``progress_factory`` is an optional callable ``(label, total) -> Progress``
    used to report progress for the two counting passes.
    """
    coding_progress = _start_progress(progress_factory, "coding k-mers", coding_path)
    coding_counts = kmer_freq_file(
        coding_path, word_size=word_size, step_size=CODING_STEP, frame=0, progress=coding_progress
    )
    _finish_progress(coding_progress)

    noncoding_progress = _start_progress(progress_factory, "noncoding k-mers", noncoding_path)
    noncoding_counts = kmer_freq_file(
        noncoding_path, word_size=word_size, step_size=NONCODING_STEP, frame=0, progress=noncoding_progress
    )
    _finish_progress(noncoding_progress)

    coding_total = float(sum(coding_counts.values()))
    noncoding_total = float(sum(noncoding_counts.values()))
    if coding_total <= 0:
        raise HexamerTableError("No %d-mer found in coding FASTA: %s" % (word_size, coding_path))
    if noncoding_total <= 0:
        raise HexamerTableError("No %d-mer found in noncoding FASTA: %s" % (word_size, noncoding_path))

    rows = []
    for kmer in sorted(coding_counts):
        if "N" in kmer:
            continue
        rows.append((
            kmer,
            coding_counts[kmer] / coding_total,
            noncoding_counts.get(kmer, 0) / noncoding_total,
        ))
    return rows


def _start_progress(progress_factory, label, path):
    if progress_factory is None:
        return None
    return progress_factory(label, count_fasta_records(path))


def _finish_progress(progress):
    if progress is not None:
        progress.finish()


def write_hexamer_table(handle, rows):
    handle.write("\t".join(HEADER) + "\n")
    for kmer, coding_freq, noncoding_freq in rows:
        handle.write("%s\t%.10f\t%.10f\n" % (kmer, coding_freq, noncoding_freq))


def write_hexamer_file(path, rows):
    """Write the table atomically, honouring the delimiter implied by ``path``."""
    table_rows = [[kmer, "%.10f" % coding_freq, "%.10f" % noncoding_freq] for kmer, coding_freq, noncoding_freq in rows]
    return write_table(path, HEADER, table_rows)
