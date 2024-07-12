"""Configuration variables."""

import importlib.resources as resourcer

PACKAGE_DIR = resourcer.files("msp_genomes")
DATABASES_DIR = PACKAGE_DIR / "databases"
MSP_COLLECTIONS = DATABASES_DIR / "filtered_collections.xlsx"
PHENOTYPES_TABLE = DATABASES_DIR / "phenotypes.txt"
DISINFINDER_DB = DATABASES_DIR / "disinfinder_db"
POINTFINDER_DB = DATABASES_DIR / "pointfinder_db"
PLASMIDFINDER_DB = DATABASES_DIR / "plasmidfinder_db"
RESFINDER_DB = DATABASES_DIR / "resfinder_db"
VIRULENCEFINDER_DB = DATABASES_DIR / "virulencefinder_db"
MLST_DB = DATABASES_DIR / "mlst_db"
MLST_CONFIG = MLST_DB / "config"
TMP_DIR = PACKAGE_DIR / "tmp"
OUTPUT_NAMES_COMPILATIONS = {
    "mlst": "mlst_compilation.xlsx",
    "plasmidfinder": "plasmidfinder_compilation.xlsx",
    "resfinder": "resfinder_compilation.xlsx",
    "virulencefinder": "virulencefinder_compilation.xlsx",
}
