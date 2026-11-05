from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import platform

STOP_CODONS = {"tag", "taa", "tga"}


def find_txcdspredict():
    """Locate the optional external txCdsPredict binary.

    The package no longer ships the Linux-only ELF binary: the default CDS
    prediction is a pure Python implementation, and the external binary is an
    opt-in escape hatch (see should_use_external_txcds).
    """
    binary = shutil.which("txCdsPredict")
    if binary:
        return binary

    raise RuntimeError(
        "txCdsPredict was not found in PATH. Either unset PRELNC_USE_EXTERNAL_TXCDS "
        "to use the built-in Python implementation, or put the txCdsPredict binary "
        "(Linux x86-64 only) on your PATH."
    )


def parse_txcds_output(path):
    cds_len = {}
    cds_score = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            transcript_id = fields[0].split("|")[0]
            pred_start = int(fields[1])
            pred_end = int(fields[2])
            cds_len[transcript_id] = pred_end - pred_start
            cds_score[transcript_id] = float(fields[5])
    return cds_len, cds_score


def _iter_fasta_records(path):
    name = None
    parts = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def _find_orf_end(sequence, start_pos):
    seq_len = len(sequence)
    pos = start_pos + 3
    while pos + 3 <= seq_len:
        codon = sequence[pos:pos + 3]
        if codon in STOP_CODONS:
            return pos + 3
        pos += 3
    return seq_len


def _is_kozak(sequence, start_pos):
    if sequence[start_pos:start_pos + 3] != "atg":
        return False
    if start_pos + 3 < len(sequence) and sequence[start_pos + 3] == "g":
        return True
    if start_pos >= 3 and sequence[start_pos - 3] in ("a", "g"):
        return True
    return False


def _covered_bases_before(length, intervals):
    diff = [0] * (length + 1)
    for start, end in intervals:
        if start >= length:
            continue
        end = min(end, length)
        diff[start] += 1
        diff[end] -= 1

    covered_before = [0] * length
    active = 0
    covered = 0
    for idx in range(length):
        active += diff[idx]
        covered_before[idx] = covered
        if active > 0:
            covered += 1
    return covered_before


def _candidate_starts(sequence, any_start=True):
    starts = set()
    if any_start:
        for idx in range(min(3, len(sequence))):
            if idx + 3 <= len(sequence):
                starts.add(idx)
    start = sequence.find("atg")
    while start >= 0:
        starts.add(start)
        start = sequence.find("atg", start + 1)
    return sorted(starts)


def _build_upstream_penalties(sequence):
    atg_intervals = []
    kozak_intervals = []
    for start_pos in _candidate_starts(sequence, any_start=False):
        end_pos = _find_orf_end(sequence, start_pos)
        atg_intervals.append((start_pos, end_pos))
        if _is_kozak(sequence, start_pos):
            kozak_intervals.append((start_pos, end_pos))
    return _covered_bases_before(len(sequence), atg_intervals), _covered_bases_before(len(sequence), kozak_intervals)


def _score_candidate(sequence, start_pos, end_pos, up_atg_count, up_kozak_count):
    size = end_pos - start_pos
    size -= size % 3
    end_pos = start_pos + size
    score = float(size)
    if sequence[start_pos:start_pos + 3] == "atg":
        score += 50.0
        if _is_kozak(sequence, start_pos):
            score += 100.0
    if end_pos <= len(sequence) and end_pos - 3 >= 0 and sequence[end_pos - 3:end_pos] in STOP_CODONS:
        score += 50.0
    score -= 0.5 * up_atg_count[start_pos]
    score -= 0.5 * up_kozak_count[start_pos]
    return score, size


def _predict_single_record(sequence):
    sequence = sequence.lower()
    up_atg_count, up_kozak_count = _build_upstream_penalties(sequence)
    best = None
    for start_pos in _candidate_starts(sequence, any_start=True):
        raw_end = _find_orf_end(sequence, start_pos)
        score, size = _score_candidate(sequence, start_pos, raw_end, up_atg_count, up_kozak_count)
        candidate = {
            "start": start_pos,
            "end": start_pos + size,
            "length": size,
            "score": score,
            "is_atg": sequence[start_pos:start_pos + 3] == "atg",
        }
        key = candidate["score"]
        if best is None or key > best[0]:
            best = (key, candidate)
    return best[1]


def run_python_txcds_predict(fasta_path, progress=None):
    cds_len = {}
    cds_score = {}
    for transcript_id, sequence in _iter_fasta_records(fasta_path):
        prediction = _predict_single_record(sequence)
        cds_len[transcript_id] = prediction["length"]
        cds_score[transcript_id] = int(prediction["score"])
        if progress is not None:
            progress.advance()
    return cds_len, cds_score


def should_use_external_txcds():
    return os.environ.get("PRELNC_USE_EXTERNAL_TXCDS", "").lower() in ("1", "true", "yes")


def run_txcds_predict(fasta_path, workdir=None, progress=None):
    if not should_use_external_txcds():
        return run_python_txcds_predict(fasta_path, progress=progress)

    system_name = platform.system()
    if system_name != "Linux":
        raise RuntimeError(
            "Bundled txCdsPredict only supports Linux, but the current system is %s. "
            "Use a Linux environment, container, or replace txCdsPredict with a compatible build."
            % system_name
        )

    binary = find_txcdspredict()
    fd, output_path = tempfile.mkstemp(prefix="prelnc_", suffix=".cds", dir=workdir)
    os.close(fd)
    command = [binary, fasta_path, "-anyStart", output_path]
    try:
        subprocess.check_call(command)
        return parse_txcds_output(output_path)
    except OSError as error:
        raise RuntimeError("Failed to execute txCdsPredict: %s" % error)
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)
