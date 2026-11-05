from __future__ import absolute_import

import importlib
import pickle
import sys
import warnings

import numpy


def _legacy_module_aliases():
    aliases = {}
    module_map = {
        "sklearn.ensemble.forest": "sklearn.ensemble._forest",
        "sklearn.tree.tree": "sklearn.tree._classes",
    }
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"A NumPy version >=1\.16\.5 and <1\.23\.0 is required for this version of SciPy",
            category=UserWarning,
        )
        for legacy_name, current_name in module_map.items():
            if legacy_name in sys.modules:
                continue
            try:
                aliases[legacy_name] = importlib.import_module(current_name)
            except ImportError:
                continue
    return aliases


def load_model(path):
    aliases = _legacy_module_aliases()
    with open(path, "rb") as handle:
        try:
            for legacy_name, module in aliases.items():
                sys.modules[legacy_name] = module
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r"A NumPy version >=1\.16\.5 and <1\.23\.0 is required for this version of SciPy",
                    category=UserWarning,
                )
                warnings.filterwarnings(
                    "ignore",
                    message=r"Trying to unpickle estimator .* from version .* when using version .*",
                    category=UserWarning,
                )
                return pickle.load(handle, encoding="latin1")
        finally:
            for legacy_name in aliases:
                sys.modules.pop(legacy_name, None)


def save_model(model, path):
    from prelnc.outputs import atomic_open

    with atomic_open(path, "wb") as handle:
        pickle.dump(model, handle)
    return path


def predict_probabilities(model, matrix):
    array = numpy.array(matrix, dtype=numpy.float32)
    return model.predict_proba(array)[:, 1]


def training_accuracy(model, matrix, labels):
    """Accuracy on the training set, as a quick sanity check of the fit."""
    array = numpy.array(matrix, dtype=numpy.float32)
    return float(model.score(array, numpy.array(labels)))


def train_model(matrix, labels):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"A NumPy version >=1\.16\.5 and <1\.23\.0 is required for this version of SciPy",
            category=UserWarning,
        )
        from sklearn.ensemble import RandomForestClassifier

    array = numpy.array(matrix, dtype=numpy.float32)
    y = numpy.array(labels)
    classifier = RandomForestClassifier(n_estimators=10, class_weight="balanced")
    return classifier.fit(array, y)
