from __future__ import absolute_import

import os


PACKAGE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(PACKAGE_DIR, "data")
MODELS_DIR = os.path.join(DATA_DIR, "models")

START_CODONS = "ATG"
STOP_CODONS = "TAG,TAA,TGA"

ANIMAL_FEATURES = [
    "cds_percentage",
    "fickett_score",
    "hexamer_score",
    "cds_score",
    "CGG",
    "TAG",
    "gc_content",
    "std_stop",
    "isoelectric_point",
    "ACG",
    "GGC",
    "sequence_length",
    "CGT",
    "AGC",
    "GAC",
    "GGG",
    "TCA",
    "CAT",
    "orf_integrity",
    "CAG",
]

PLANT_FEATURES = [
    "sequence_length",
    "hexamer_score",
    "cds_score",
    "std_stop",
    "fickett_score",
    "isoelectric_point",
    "GGG",
    "TAG",
    "CGA",
    "CAA",
    "cds_percentage",
    "GAT",
    "ACT",
    "gc_content",
    "ATC",
    "orf_integrity",
    "AAA",
    "AAC",
    "CTA",
    "AGG",
    "GAG",
    "AAT",
    "GAC",
    "CAG",
    "GGA",
    "CAT",
    "AGC",
    "TAT",
    "TAA",
    "ACC",
    "GTA",
    "ATG",
    "TTA",
]

SPECIES_CONFIG = {
    "human": {
        "type": "animal",
        "model": "human_model.pkl",
        "hexamer": "human_hexamer.tsv",
    },
    "mouse": {
        "type": "animal",
        "model": "mouse_model.pkl",
        "hexamer": "mouse_hexamer.tsv",
    },
    "cow": {
        "type": "animal",
        "model": "cow_model.pkl",
        "hexamer": "cow_hexamer.tsv",
    },
    "ath": {
        "type": "plant",
        "model": "ath_model.pkl",
        "hexamer": "ath_hexamer.tsv",
    },
    "osa": {
        "type": "plant",
        "model": "osa_model.pkl",
        "hexamer": "osa_hexamer.tsv",
    },
    "zma": {
        "type": "plant",
        "model": "zma_model.pkl",
        "hexamer": "zma_hexamer.tsv",
    },
}

FEATURE_SETS = {
    "animal": ANIMAL_FEATURES,
    "plant": PLANT_FEATURES,
}


def available_species():
    return sorted(SPECIES_CONFIG.keys())


def feature_names(feature_type):
    try:
        return FEATURE_SETS[feature_type]
    except KeyError:
        raise ValueError("Unknown feature type: %s" % feature_type)


def default_model_path(species):
    return os.path.join(MODELS_DIR, SPECIES_CONFIG[species]["model"])


def default_hexamer_path(species):
    return os.path.join(MODELS_DIR, SPECIES_CONFIG[species]["hexamer"])


def resolve_profile(species=None, feature_type=None, model_path=None, hexamer_path=None):
    if species:
        if species not in SPECIES_CONFIG:
            raise ValueError("Unsupported species: %s" % species)
        config = SPECIES_CONFIG[species]
        resolved_type = config["type"]
        resolved_model = model_path or default_model_path(species)
        resolved_hexamer = hexamer_path or default_hexamer_path(species)
    else:
        if not feature_type:
            raise ValueError("Either species or feature type must be provided")
        resolved_type = feature_type
        resolved_model = model_path
        resolved_hexamer = hexamer_path

    if resolved_type not in FEATURE_SETS:
        raise ValueError("Unsupported feature type: %s" % resolved_type)

    return {
        "type": resolved_type,
        "model": resolved_model,
        "hexamer": resolved_hexamer,
        "features": feature_names(resolved_type),
    }

