# PreLnc 2.0

PreLnc distinguishes long non-coding transcripts from protein-coding transcripts for several plant and animal species, using a random-forest model over sequence-derived features (Fickett score, hexamer ratio, ORF integrity, CDS prediction score, GC content and trinucleotide frequencies).

Original method: [Genes 2020, 11(9), 981](https://www.mdpi.com/2073-4425/11/9/981/htm).

Version 2.0 is a rewrite of the published 1.x code (`1.x` branch, tag `v1.0`): it is a real installable package with a single command-line interface, it no longer needs the Linux-only `txCdsPredict` binary, it writes proper CSV/TSV tables without overwriting your files, and it reports progress while it runs. Feature computation and the bundled models are unchanged, so `prediction_score` values are identical to the published version.

---

## Which version do I want?

| Branch / tag | Contents |
|---|---|
| `master` (default) | PreLnc 2.0 — Python 3, installable, maintained. This README. |
| `1.x` | The 2020 Python 2 release, frozen and no longer maintained. Kept unchanged so published results can be reproduced; its README points back here. |
| `v2.0.0` / `v1.0` | Immutable tags of the two releases, for a fixed snapshot. |

Stay on `master` unless you specifically need the original scripts: `1.x` requires Python 2.7, a Linux x86-64 `txCdsPredict` binary, and manual `-m` / `-r` model paths.

---

## Changes since 1.x

`1.x` was a set of Python 2 scripts; 2.0 is one Python 3 package. What that means in practice:

| Area | 1.x | PreLnc 2.0 |
|---|---|---|
| Python | 2.7 only (`from __future__ import print_function`, `import Bio.Alphabet`) | 3.7+ (tests run on 3.7 and 3.8) |
| Install | `git clone`, `pip install biopython sklearn`, then edit `~/.bashrc` so `txCdsPredict` is on `PATH` | `python -m pip install .`; no manual `PATH` editing |
| Commands | 4 scripts inside `bin/`, pick one by species: `prelnc.py` / `prelnc_plant.py` / `createmodel.py` / `createmodel_plant.py` | 3 console commands: `prelnc`, `prelnc-train`, `prelnc-hexamer` |
| Species | encoded in the script name | `--species human\|mouse\|cow\|ath\|osa\|zma` (or `--type animal\|plant` with your own model) |
| Model / hexamer paths | passed by hand with `-m` / `-r`, relative to the checkout | bundled in the package and resolved from `--species`; `-m` / `-r` only to override |
| CDS prediction | required the compiled `txCdsPredict` (Linux x86-64 ELF) on `PATH` | pure Python implementation by default, no external binary; the binary is still usable via `PRELNC_USE_EXTERNAL_TXCDS=1` |
| Output | one tab-separated table with no header and no transcript id, assembled by `paste`-ing a features file with the score column | two header-ed CSV/TSV tables: predictions (`transcript_id`, `prediction_score`, `prediction`) and features (`transcript_id` + one column per feature) |
| Overwriting | silently replaced the target file | writes `out_1.tsv`, `out_2.tsv`, … instead; `--force` to overwrite; companion files renamed with the same stem |
| Failed runs | scratch files (`<output>.feature`, a fixed-name `TMP.txt` in the working directory) could be left behind; assembly shelled out to `paste`/`rm`, so it was Unix-only | written to a temporary sibling and moved into place only on success, no shell calls |
| Feedback | two "start to ..." lines on stderr, then total cost | per-stage progress (percent, done/total, ETA) on stderr and an end-of-run summary listing counts and every file written |
| Training | `createmodel.py` wrote only a `.pkl`; feature order was implicit | `prelnc-train` writes `.pkl` + `.feature.tsv` + `.report.tsv`, and reports training accuracy |
| Hexamer table | `bin/make_hexamer_tab.py`, an unmodified copy of the CPAT script that needs CPAT's `cpmodule` on `sys.path` and dies with `NameError` on empty input | `prelnc-hexamer`, same algorithm, explicit error instead of a traceback |
| Code duplication | `prelnc.py` / `prelnc_plant.py` and `createmodel.py` / `createmodel_plant.py` were near-identical; the CPAT/CPC2 primitives lived in `bin/cpmodule/` | one code path in `prelnc/`; the same primitives kept verbatim under `prelnc/legacy/` |
| Repository size | 18 MB, half of it an unused `models/modelFor8CPUs/` bundle | 9.6 MB, one copy of each bundled model |
| Tests | none | 27 unittest cases: hexamer parity with the original script, a full prediction run on `example/example.fa`, overwrite/delimiter handling, progress formatting |
| Version | no version string in the code (tag `v1.0` added retroactively) | `2.0.0`, `prelnc --version` |

What did **not** change: the feature definitions, the hexamer ratio computation, and the six bundled random-forest models. On `example/example.fa` and the first 100 *O. sativa* transcripts, `prediction_score` and all feature values are bit-for-bit identical to the published output.

Output-format changes are the one thing that can break a downstream script: 1.x wrote a single headerless table (all features plus the score) to `<output>`. 2.0 writes two header-ed tables — `<output>` with `transcript_id`, `prediction_score` and `prediction`, and `<stem>.features.tsv` with `transcript_id` plus one column per feature. Existing files are also renamed instead of overwritten. Everything else is additive.

---

## Install

```bash
python -m pip install .
```

Use `--user` if the environment is not writable:

```bash
python -m pip install --user .
```

This installs three commands:

| Command | Purpose |
|---|---|
| `prelnc` | predict lncRNA probability |
| `prelnc-train` | train a model for a new species |
| `prelnc-hexamer` | build a hexamer frequency table |

Dependencies (`biopython`, `numpy`, `scikit-learn`) are installed automatically.

## Quick start

```bash
prelnc -i example/example.fa --species human -o example.tsv
```

`-o` is optional. Without it the result goes to the current directory as
`./prelnc_human_prediction_<timestamp>.tsv`, so a run never clobbers a file you
already have.

Output:

- `example.tsv` — `transcript_id`, `prediction_score` (probability of being non-coding) and `prediction` (`lncRNA` above 0.5, otherwise `coding`)
- `example.features.tsv` — `transcript_id` plus one column per extracted feature, for inspection

From a source checkout without installing, use `python -m prelnc`; the same
entry point accepts `predict`, `train` and `hexamer` subcommands:

```bash
python -m prelnc train -l lnc.fa -p coding.fa --species human -o human_model.pkl
```

## Output files

### Format

The extension decides the delimiter, checked before the run starts:

| Extension | Delimiter |
|---|---|
| `.csv` | comma |
| `.tsv`, `.tab`, `.txt` | tab |

Any other extension is rejected with an error. All tables are UTF-8 with a
header row.

### Existing files are never overwritten

If the target exists, a numbered copy is written instead (`out.tsv` →
`out_1.tsv`, `out_2.tsv`, …) and the actual path is printed. Pass `--force` to
overwrite on purpose. A companion file (feature table, report) is renamed
together with its result file so the pair always keeps the same stem.

Files are written to a temporary sibling and moved into place only after the run
succeeds, so an interrupted run leaves no half-written table behind.

## Progress and summary

Every command reports what it is doing on stderr. On a terminal the line is
refreshed in place; when stderr is redirected it prints one line per 5%, which
keeps log files readable:

```
[PreLnc] Predict: transcripts.fa (5000 sequences)
[PreLnc] Predict CDS [ 45%] 2250/5000 items · elapsed 1m38s · eta 2m00s
```

At the end a summary lists how many sequences were processed, how many were
called lncRNA, the elapsed time, and exactly which files were written:

```
[PreLnc] Prediction finished
[PreLnc] Input       : transcripts.fa (5000 sequences)
[PreLnc] lncRNA      : 2140 (42.8%)
[PreLnc] coding      : 2860 (57.2%)
[PreLnc] Elapsed     : 3m12s
[PreLnc] Predictions : ./prelnc_human_prediction_20261105_132500.tsv (5000 rows x 3 columns)
[PreLnc] Features    : ./prelnc_human_prediction_20261105_132500.features.tsv (5000 rows x 21 columns)
```

`-q/--quiet` hides the progress lines; the summary is still printed.

## Built-in species

| `--species` | Feature set |
|---|---|
| `human`, `mouse`, `cow` | animal (20 features) |
| `ath`, `osa`, `zma` | plant (33 features) |

`--species` resolves the feature set, the bundled model and the bundled hexamer table automatically. Pass `-m/--model` or `-r/--hexamer` only to override them.

### Custom model or species

```bash
prelnc -i transcripts.fa --type animal -m my_model.pkl -r my_hexamer.tsv -o out.txt
```

`--type` selects the feature set (`animal` or `plant`) when you supply your own model.

## Training a model

```bash
prelnc-train -l lnc.fa -p coding.fa --species human -o human_model.pkl
```

Writes three files, all sharing the same stem (rename on conflict applies to the
whole set):

| File | Content |
|---|---|
| `human_model.pkl` | the trained model (binary, not a table) |
| `human_model.feature.tsv` | labelled feature matrix: `transcript_id`, features, `label` (1 = lncRNA, 0 = coding) |
| `human_model.report.tsv` | key/value training report: sample counts, feature count, tree count, training accuracy, versions |

The summary printed at the end includes the training accuracy; note that it is
measured on the training set, so it is a sanity check rather than a
generalisation estimate.

## Building a hexamer table for a new species

The hexamer table is the coding vs non-coding 6-mer frequency ratio table, and it is species specific.

```bash
prelnc-hexamer -c coding_cds.fa -n noncoding.fa -o my_hexamer.tsv
```

- `-c` must be **CDS only** (start codon to stop codon, no UTR)
- `-n` is the non-coding transcript set
- omit `-o` to write the table to stdout, as the original `make_hexamer_tab.py` did
- `-o` accepts `.tsv` or `.csv`, and refuses to overwrite an existing table unless `--force` is given

Then train:

```bash
prelnc-train -l lnc.fa -p coding_cds.fa --type animal -r my_hexamer.tsv -o my_model.pkl
```

## Column reference

Prediction table (`-o`):

```
transcript_id	prediction_score	prediction
ENST00000658521.1	1.00000000	lncRNA
```

- `prediction_score` — probability that the transcript is non-coding
- `prediction` — `lncRNA` when the score is >= 0.5, otherwise `coding`

Feature table (`<stem>.features.tsv`) — `transcript_id` plus one column per
feature, in the order defined by the selected feature set (`prelnc/config.py`).
The id column is included so rows can be joined back to the result table.

## CDS prediction backend

By default `prelnc` uses a pure Python implementation of `txCdsPredict -anyStart`, so it runs on any platform with no external binary.

The bundled C++ binary was Linux x86-64 only and is no longer shipped. To use it instead:

```bash
export PRELNC_USE_EXTERNAL_TXCDS=1
# txCdsPredict must be on PATH
```

Note that the Python implementation is a re-implementation: scores are close but not guaranteed identical to the C++ original. Use the external binary if you need to reproduce published numbers exactly.

## Project layout

```text
prelnc/
  cli.py              command-line entry points
  config.py           species profiles and feature sets (single source of truth)
  external.py         CDS prediction (Python implementation + optional binary)
  feature_extractors.py
  hexamer.py          hexamer table construction
  io_utils.py         FASTA/table read-write helpers
  model_ops.py        model load / predict / train
  outputs.py          output path resolution, overwrite protection, atomic writers
  pipelines.py        prediction and training workflows
  progress.py         progress bars and end-of-run summaries
  legacy/             Fickett, ORF and k-mer primitives inherited from CPAT/CPC2
  data/models/        bundled models and hexamer tables
example/              example input
tests/                unit tests and fixtures
```

## Tests

```bash
python -m unittest discover -s tests
```

Covers hexamer table construction against a reference output from the original
`make_hexamer_tab.py`, full prediction runs on `example/example.fa`, overwrite
protection and delimiter handling, progress/summary formatting, and version
consistency.

## Notes and limitations

- The bundled `.pkl` models were trained with an old scikit-learn. `pyproject.toml` pins `numpy<2` and `scikit-learn<1.4` so they remain loadable; retraining is recommended for long-term use.
- `pip install .` needs setuptools 61+ (or network access, so pip can build in an isolated environment). On an offline machine with an old toolchain, upgrade pip/setuptools first.
- Requires Python 3.7+.
- `prediction_score` values are unchanged from 1.x for the same input and models (verified on `example/example.fa` and 100 *O. sativa* transcripts). What is new in 2.0 is the split into a prediction table (`transcript_id`, `prediction_score`, `prediction`) and a feature table (`<stem>.features.tsv`, with `transcript_id` and a header), plus the refusal to overwrite existing files.
- The `hexamer` column of a historical development snapshot (`test_osa5000-output2.txt`, not part of the 1.x release) cannot be reproduced by any hexamer table bundled here: all 48 combinations (6 tables x 4 sequence segments x 2 steps) were tried, none matched. That file most likely came from a table that was never committed, so regenerate baselines with the current code rather than using it for regression.
- Progress ETAs are extrapolated from elapsed time, so they are rough for inputs whose transcripts vary a lot in length.

## License

GPL-3.0. See [LICENSE](LICENSE).
