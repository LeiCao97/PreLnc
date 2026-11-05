from __future__ import absolute_import

import gzip

from prelnc.outputs import write_table

PREDICTION_THRESHOLD = 0.5


def load_hexamer_table(path):
    coding = {}
    noncoding = {}
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if fields[0] == "hexamer":
                continue
            coding[fields[0]] = float(fields[1])
            noncoding[fields[0]] = float(fields[2])
    return coding, noncoding


def count_fasta_records(path):
    """Count ``>`` headers without parsing sequences (fast, used for progress)."""
    opener = gzip.open if path.endswith(".gz") else open
    total = 0
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                total += 1
    return total


def write_feature_file(path, transcript_ids, rows, feature_names, include_label=False):
    """Write the feature matrix with ``transcript_id`` as the first column."""
    header = ["transcript_id"] + list(feature_names)
    if include_label:
        header.append("label")
    table_rows = []
    for transcript, row in zip(transcript_ids, rows):
        table_rows.append([transcript] + list(row))
    return write_table(path, header, table_rows)


def prediction_label(score, threshold=PREDICTION_THRESHOLD):
    return "lncRNA" if score >= threshold else "coding"


def write_prediction_output(path, rows, probabilities):
    table_rows = []
    for row, score in zip(rows, probabilities):
        table_rows.append([row["transcript_id"], "%.8f" % score, prediction_label(score)])
    return write_table(path, ["transcript_id", "prediction_score", "prediction"], table_rows)


def rows_to_matrix(rows, feature_names):
    matrix = []
    for row in rows:
        matrix.append([row[name] for name in feature_names])
    return matrix


def rows_to_feature_values(row, feature_names):
    return [row[name] for name in feature_names]
