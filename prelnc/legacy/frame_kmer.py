from __future__ import absolute_import

import gzip
import itertools
import math
import re
from collections import Counter


DNA_PATTERN = re.compile(r"^[ACGTN]+$")


def word_generator(seq, word_size, step_size, frame=0):
    for idx in range(frame, len(seq), step_size):
        word = seq[idx:idx + word_size]
        if len(word) == word_size:
            yield word


def all_possible_kmer(length):
    for item in itertools.product(["A", "C", "G", "T", "N"], repeat=length):
        yield "".join(item)


def seq_generator(lines):
    tmpseq = ""
    name = ""
    for raw_line in lines:
        line = raw_line.strip().upper()
        if line.startswith(("#", " ", "\n")):
            continue
        if line.startswith((">", "@")):
            if tmpseq:
                yield [name, tmpseq]
                tmpseq = ""
            name = line.split()[0][1:]
        elif DNA_PATTERN.match(line):
            tmpseq += line
    yield [name, tmpseq]


def kmer_freq_lines(lines, word_size, step_size=1, frame=0, min_count=0, progress=None):
    ret_dict = {}
    count_table = Counter()
    for _, seq in seq_generator(lines):
        if seq:
            count_table.update(word_generator(seq, word_size=word_size, step_size=step_size, frame=frame))
            if progress is not None:
                progress.advance()

    for kmer in all_possible_kmer(word_size):
        if kmer not in count_table:
            count_table[kmer] = 0
        if count_table[kmer] >= min_count and "N" not in kmer:
            ret_dict[kmer] = count_table[kmer]
    return ret_dict


def kmer_freq_file(path, word_size, step_size=1, frame=0, min_count=0, progress=None):
    """Count kmers in a FASTA file. Gzipped input is detected by suffix."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        return kmer_freq_lines(
            handle,
            word_size=word_size,
            step_size=step_size,
            frame=frame,
            min_count=min_count,
            progress=progress,
        )


def kmer_ratio(seq, word_size, step_size, coding, noncoding):
    if len(seq) < word_size:
        return 0.0

    sum_of_log_ratio = 0.0
    frame_count = 0.0
    for kmer in word_generator(seq=seq, word_size=word_size, step_size=step_size, frame=0):
        if (kmer not in coding) or (kmer not in noncoding):
            continue
        if coding[kmer] > 0 and noncoding[kmer] > 0:
            sum_of_log_ratio += math.log(coding[kmer] / noncoding[kmer])
        elif coding[kmer] > 0 and noncoding[kmer] == 0:
            sum_of_log_ratio += 1
        elif coding[kmer] == 0 and noncoding[kmer] > 0:
            sum_of_log_ratio -= 1
        else:
            continue
        frame_count += 1

    if not frame_count:
        return -1.0
    return sum_of_log_ratio / frame_count

