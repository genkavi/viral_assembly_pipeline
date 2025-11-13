# Viral Assembly FoldX Workflow

A computational pipeline for analyzing complete viral assemblies with icosahedral symmetries using molecular dynamics and FoldX protein stability calculations.

## Overview

This pipeline integrates molecular dynamics (MD) simulations with FoldX stability analysis to study viral envelope proteins and their mutations. It handles complex icosahedral viral structures with multiple assembly types (3-fold, 5-fold, dimer configurations) and performs comprehensive stability analysis through alanine scanning mutagenesis.

The pipeline offers two optimization strategies: FoldX RepairPDB for fast side chain optimization, and MODELLER for more thorough refinement with backbone flexibility to resolve assembly clashes.

- Extract frames from MD trajectories
- Assemble complete viral envelopes from asymmetric units using BIOMT transformations
- Optimize structures with FoldX RepairPDB (side chains only) or MODELLER (backbone + side chains)
- Calculate stability metrics (ΔΔG) for viral assemblies
- Perform alanine scanning mutagenesis for all residues
- Statistical analysis with proper error propagation for hierarchical data

## Installation

### Prerequisites

- Python 3.10+
- FoldX (tested with version 4.1)
- MODELLER (for structure refinement with backbone flexibility)
- GROMACS (for trajectory processing)
- Conda or Mamba (recommended)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/viral-assembly-foldx.git
cd viral-assembly-foldx
```

2. Create conda environment:
```bash
conda env create -f environment_minimal.yml
conda activate md_foldx

# Install MODELLER (requires free academic license from https://salilab.org/modeller/)
conda install -c salilab modeller
```

3. Configure FoldX paths in `config.py`:
```python
FOLDX_PATH = "/path/to/foldx"
ROTABASE_PATH = "/path/to/rotabase.txt"
```

## Project Structure

```
.
├── config.py                      # Global configuration (FoldX paths, parallel workers)
├── extract_frames.py              # Extract frames from MD trajectories
├── assemble_structure.py          # Build viral assemblies from templates
├── optimize_structures.py         # Run FoldX RepairPDB optimization
├── refine_assembly_modeller.py    # MODELLER-based refinement with backbone flexibility
├── foldx_wrapper.py              # FoldX command wrapper (Stability, AlaScan)
├── analyze_foldx.py              # Statistical analysis of FoldX results
├── run_pipeline_stability.sh     # Complete workflow automation
├── alascan.ipynb                 # Alanine scan analysis and visualization
├── mutations.ipynb               # Mutation impact analysis
└── environment_minimal.yml       # Conda environment specification
```

## Usage

### Quick Start

Run the complete pipeline:
```bash
bash run_pipeline_stability.sh
```

### Individual Components

#### 1. Extract MD Frames
```bash
python extract_frames.py \
    -t trajectory.xtc \
    -s topology.pdb \
    -o output_dir/ \
    --start 500000 --end 1000000 --dt 100000
```

#### 2. Assemble Viral Structures
```bash
python assemble_structure.py \
    --input frames/*.pdb \
    --template template_assembly.pdb \
    --output-dir assemblies/ \
    --workers 96
```

#### 3. Optimize Structures (FoldX)
```bash
python optimize_structures.py \
    --input structures/*.pdb \
    --output-dir optimized/ \
    --repair \
    --workers 96
```

#### 3b. Refine Assemblies (MODELLER - Alternative/Additional)
```bash
python refine_assembly_modeller.py \
    --input assemblies/*.pdb \
    --output-dir assemblies-refined/ \
    --cg-steps 50 --md-steps 100 --md-temp 300 \
    --workers 96
```

#### 4. Stability Analysis
```bash
python foldx_wrapper.py \
    --command Stability \
    --input optimized/*.pdb \
    --output-dir stability_results/ \
    --workers 96
```

#### 5. Alanine Scanning
```bash
python foldx_wrapper.py \
    --command AlaScan \
    --input optimized/*.pdb \
    --chains A B C D E F \
    --output-dir alascan_results/ \
    --workers 96
```

#### 6. Analyze Results
```bash
python analyze_foldx.py \
    --input stability_results/ \
    --output all_stability.csv \
    --assembly-totals assembly_totals.csv \
    --viral-envelope viral_envelope_stability.csv \
    --energy-column total_energy
```

### Output Files

The pipeline generates several output directories:

```
viral_assembly_pipeline_results/
├── frames/                     # Extracted MD frames
├── frames-optimized/           # Optimized individual frames
├── assembly/                   # Assembled viral envelopes (by template)
├── assembly-optimized/         # Optimized assemblies (FoldX RepairPDB)
├── assembly-refined/           # MODELLER-refined assemblies (optional)
│   └── *.out                  # Energy reports for each structure
├── stability/                  # Stability analysis results
│   ├── all_stability.csv      # Per-chain energy values
│   ├── assembly_totals.csv    # Total assembly energies
│   └── viral_envelope_stability.csv  # Averaged results
└── alascan/                    # Alanine scan mutation data
```

### Statistical Analysis

The analysis pipeline properly handles hierarchical data structure:
- **Independent samples**: MD trajectory frames (averaged with standard error)
- **Correlated measurements**: Symmetry-related chains (averaged without inflating variance)
- **Error propagation**: Calculates uncertainties at each hierarchical level

## Configuration

Edit `config.py` to customize:

```python
# FoldX executable paths
FOLDX_PATH = "/path/to/foldx"
ROTABASE_PATH = "/path/to/rotabase.txt"

# Parallel processing
DEFAULT_WORKERS = None  # Auto-detect CPU count

# Residue naming (for CHARMM force fields)
HIS_RENAME_MAP = {'HSE': 'HIS', 'HSD': 'HIS', ...}
FIX_ILE_CD = True  # Fix ILE CD -> CD1 naming
```

### MODELLER Refinement Parameters

The MODELLER refinement uses a three-stage protocol: **CG → MD → CG**

```bash
--cg-steps 50      # Conjugate gradient steps per round (20-100)
--md-steps 100     # Molecular dynamics steps (50-200)
--md-temp 300      # MD temperature in Kelvin (300-400)
```

**Parameter guidelines:**
- **Quick refinement**: `--cg-steps 20 --md-steps 50 --md-temp 300` (fast, mild relaxation)
- **Standard refinement**: `--cg-steps 50 --md-steps 100 --md-temp 300` (default, balanced)
- **Thorough refinement**: `--cg-steps 100 --md-steps 200 --md-temp 400` (slow, aggressive relaxation)

Higher temperatures and more MD steps provide stronger relaxation but increase computational cost.

## Analysis Notebooks

- **alascan.ipynb**: Interactive visualization and analysis of alanine scanning results
- **mutations.ipynb**: Analysis of specific mutation impacts (e.g., Asibi→17D mutations)

Both notebooks include:
- ΔΔG distribution analysis
- Identification of stabilizing/destabilizing residues
- Structure coloring by energy contribution
- Statistical comparisons across frames and assemblies

