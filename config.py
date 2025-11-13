#!/usr/bin/env python3
"""
Global configuration for MD-FoldX pipeline
"""

import os
from pathlib import Path

# FoldX paths
FOLDX_PATH = os.path.expanduser("~/foldx_4.1/foldx_20251231")
ROTABASE_PATH = os.path.expanduser("~/foldx_4.1/rotabase.txt")

# Default parallel workers
DEFAULT_WORKERS = None  # None = auto-detect (CPU count - 2)

# GROMACS settings
GMX_SELECTION = "Protein"  # Default group for trjconv

# Residue name mappings (CHARMM -> PDB)
HIS_RENAME_MAP = {
    'HSE': 'HIS',
    'HSD': 'HIS',
    'HSP': 'HIS',
    'HIE': 'HIS',
    'HID': 'HIS',
    'HIP': 'HIS',
    'HIS': 'HIS'
}

FIX_ILE_CD = True  # Fix ILE CD -> CD1
