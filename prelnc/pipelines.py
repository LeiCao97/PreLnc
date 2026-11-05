from __future__ import absolute_import

import os
import platform
import time

from prelnc.external import run_txcds_predict
from prelnc.feature_extractors import extract_feature_rows
from prelnc.io_utils import (
    PREDICTION_THRESHOLD,
    count_fasta_records,
    load_hexamer_table,
    rows_to_feature_values,
    rows_to_matrix,
    write_feature_file,
    write_prediction_output,
)
from prelnc.model_ops import (
    load_model,
    predict_probabilities,
    save_model,
    train_model,
    training_accuracy,
)
from prelnc.outputs import (
    OutputError,
    describe_table,
    plural,
    resolve_outputs,
    write_table,
)
from prelnc.progress import Progress, format_duration, info, print_summary, warn


def _extract_rows(fasta_path, hexamer_tables, feature_names, quiet, stage):
    """Run CDS prediction and feature extraction with one progress bar each."""
    total = count_fasta_records(fasta_path)
    if total == 0:
        raise ValueError("No sequences found in FASTA file: %s" % fasta_path)

    info("%s: %s (%s)" % (stage, fasta_path, plural(total, "sequence")))

    cds_progress = Progress(total, "%s CDS" % stage, enabled=not quiet)
    cds_len, cds_score = run_txcds_predict(fasta_path, progress=cds_progress)
    cds_progress.finish()

    feature_progress = Progress(total, "%s features" % stage, enabled=not quiet)
    rows = extract_feature_rows(
        fasta_path, cds_len, cds_score, hexamer_tables, feature_names, progress=feature_progress
    )
    feature_progress.finish()
    return rows, total


def _prepare_outputs(plan):
    """Warn about renames or forced overwrites before anything is written."""
    existed = os.path.exists(plan["main"])
    if plan["renamed"]:
        info("%s already exists, writing to %s instead" % (plan["requested"], plan["main"]))
    elif existed:
        warn("%s already exists and will be overwritten (--force)" % plan["main"])
    return existed


def run_prediction_pipeline(
    input_path,
    model_path,
    hexamer_path,
    feature_names,
    output_path=None,
    default_stem="prelnc_prediction",
    force=False,
    quiet=False,
):
    started = time.time()
    hexamer_tables = load_hexamer_table(hexamer_path)
    rows, total = _extract_rows(input_path, hexamer_tables, feature_names, quiet, "Predict")

    if len(rows) != total:
        warn(
            "Header scan counted %d sequences but %d were parsed; progress is an estimate"
            % (total, len(rows))
        )
        total = len(rows)

    info("Scoring %s x %s ..." % (plural(total, "sequence"), plural(len(feature_names), "feature")))
    model = load_model(model_path)
    probabilities = predict_probabilities(model, rows_to_matrix(rows, feature_names))

    plan = resolve_outputs(
        output_path, default_stem=default_stem, companions=(".features.tsv",), force=force
    )
    _prepare_outputs(plan)

    write_prediction_output(plan["main"], rows, probabilities)
    write_feature_file(
        plan["companions"][0],
        [row["transcript_id"] for row in rows],
        [rows_to_feature_values(row, feature_names) for row in rows],
        feature_names,
        include_label=False,
    )

    scores = list(probabilities)
    lnc_count = sum(1 for score in scores if score >= PREDICTION_THRESHOLD)
    coding_count = total - lnc_count
    elapsed = time.time() - started

    print_summary(
        [
            ("Input", "%s (%s)" % (input_path, plural(total, "sequence"))),
            ("lncRNA", "%d (%.1f%%)" % (lnc_count, _percent(lnc_count, total))),
            ("coding", "%d (%.1f%%)" % (coding_count, _percent(coding_count, total))),
            ("Elapsed", format_duration(elapsed)),
            ("Predictions", "%s (%s)" % (plan["main"], describe_table(total, 3))),
            (
                "Features",
                "%s (%s)" % (plan["companions"][0], describe_table(total, len(feature_names) + 1)),
            ),
        ],
        title="Prediction finished",
    )

    return {
        "output": plan["main"],
        "features": plan["companions"][0],
        "rows": total,
        "lncrna": lnc_count,
        "coding": coding_count,
        "elapsed": elapsed,
    }


def run_training_pipeline(
    lnc_path,
    pc_path,
    hexamer_path,
    feature_names,
    output_path=None,
    default_stem="prelnc_model",
    force=False,
    quiet=False,
):
    started = time.time()
    hexamer_tables = load_hexamer_table(hexamer_path)

    lnc_rows, lnc_total = _extract_rows(lnc_path, hexamer_tables, feature_names, quiet, "lncRNA")
    pc_rows, pc_total = _extract_rows(pc_path, hexamer_tables, feature_names, quiet, "coding")

    if not lnc_rows or not pc_rows:
        raise ValueError("Both the lncRNA and the coding FASTA need at least one sequence.")

    all_rows = lnc_rows + pc_rows
    transcript_ids = [row["transcript_id"] for row in all_rows]
    feature_values = [rows_to_feature_values(row, feature_names) for row in all_rows]
    labels = [1] * len(lnc_rows) + [0] * len(pc_rows)
    matrix = rows_to_matrix(all_rows, feature_names)

    info(
        "Training random forest on %s x %s (10 trees) ..."
        % (plural(len(labels), "sample"), plural(len(feature_names), "feature"))
    )
    model = train_model(matrix, labels)
    accuracy = training_accuracy(model, matrix, labels)

    plan = resolve_outputs(
        output_path,
        default_stem=default_stem,
        companions=(".feature.tsv", ".report.tsv"),
        force=force,
        table=False,
        default_extension=".pkl",
    )
    _prepare_outputs(plan)

    save_model(model, plan["main"])
    write_feature_file(
        plan["companions"][0], transcript_ids, feature_values, feature_names, include_label=True
    )

    elapsed = time.time() - started
    report_rows = [
        ("created_at", time.strftime("%Y-%m-%d %H:%M:%S")),
        ("input_lncrna", os.path.abspath(lnc_path)),
        ("input_coding", os.path.abspath(pc_path)),
        ("hexamer_table", os.path.abspath(hexamer_path)),
        ("samples_lncrna", lnc_total),
        ("samples_coding", pc_total),
        ("samples_total", len(labels)),
        ("features", len(feature_names)),
        ("n_estimators", 10),
        ("train_accuracy", "%.6f" % accuracy),
        ("sklearn_version", _sklearn_version()),
        ("python_version", platform.python_version()),
        ("elapsed_sec", "%.2f" % elapsed),
    ]
    write_table(plan["companions"][1], ["metric", "value"], report_rows)

    print_summary(
        [
            (
                "Samples",
                "lncRNA %d, coding %d, total %d" % (lnc_total, pc_total, len(labels)),
            ),
            ("Features", "%d" % len(feature_names)),
            ("Train accuracy", "%.3f" % accuracy),
            ("Elapsed", format_duration(elapsed)),
            ("Model", "%s (%s)" % (plan["main"], _human_size(plan["main"]))),
            (
                "Feature table",
                "%s (%s)" % (plan["companions"][0], describe_table(len(labels), len(feature_names) + 2)),
            ),
            ("Report", "%s (%s)" % (plan["companions"][1], plural(len(report_rows), "metric"))),
        ],
        title="Training finished",
    )

    return {
        "model": plan["main"],
        "features": plan["companions"][0],
        "report": plan["companions"][1],
        "rows": len(labels),
        "accuracy": accuracy,
        "elapsed": elapsed,
    }


def _percent(part, total):
    if not total:
        return 0.0
    return 100.0 * part / float(total)


def _human_size(path):
    try:
        size = float(os.path.getsize(path))
    except OSError:
        return "unknown size"
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return "%.1f %s" % (size, unit) if unit != "B" else "%d B" % size
        size /= 1024.0
    return "%d B" % size


def _sklearn_version():
    try:
        import sklearn

        return sklearn.__version__
    except Exception:
        return "unknown"


__all__ = [
    "OutputError",
    "run_prediction_pipeline",
    "run_training_pipeline",
]
