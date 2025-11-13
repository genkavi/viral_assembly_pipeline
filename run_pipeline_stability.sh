#!/bin/bash

set -e  # Exit if any command fails

################################################################
# Configuration
################################################################

OUTPUT_DIR="viral_assembly_pipeline_results" # Output directory
WORKERS=96 # Number of parallel workers


###############################################################
# Path to the pipeline directory
############################################################### 
PIPELINE_DIR="$(dirname "$(readlink -f "$0")")" # Path to the pipeline directory

echo "=========================================="
echo "Extract frames from trajectories"
echo "=========================================="

# Extract frames from each trajectory

for traj in 0 1 2; do
python ${PIPELINE_DIR}/extract_frames.py \
    -t $traj/nW.xtc -s $traj/nW.pdb \
    -o ${OUTPUT_DIR}/frames/ -p ${traj}_ \
    --start 500000 --end 1000000 --dt 100000
done

echo ""
echo "=========================================="
echo "Optimize extracted structures"
echo "=========================================="

python ${PIPELINE_DIR}/optimize_structures.py \
    --input ${OUTPUT_DIR}/frames/*.pdb \
    --output-dir ${OUTPUT_DIR}/frames-foldx_optimized \
    --repair \
    --workers ${WORKERS}


echo ""
echo "=========================================="
echo "Assemble structures with different templates"
echo "=========================================="

for template in ${PIPELINE_DIR}/templates/*.pdb; do 
    template_name=$(basename ${template} .pdb)
    python ${PIPELINE_DIR}/assemble_structure.py  \
        --input ${OUTPUT_DIR}/frames-foldx_optimized/*.pdb  \
        --template ${template}  \
        --output-dir ${OUTPUT_DIR}/assembly/${template_name}\
        --workers 96
done

echo ""
echo "=========================================="
echo "Optimize all assembled structures with Modeller: "
echo "  Resolves Major Clashes and Optimizes Structure"
echo "=========================================="

python ${PIPELINE_DIR}/refine_assembly_modeller.py \
    --input ${OUTPUT_DIR}/assembly/*/*.pdb \
    --output-dir ${OUTPUT_DIR}/assembly-modeller_optimized \
    --cg-steps 100 --md-steps 200 --md-temp 400 \
    --workers ${WORKERS}   

echo ""
echo "=========================================="
echo "Optimize all Modeller-optimized assembled structures with Foldx"
echo "=========================================="

python ${PIPELINE_DIR}/optimize_structures.py\
    --input ${OUTPUT_DIR}/assembly-modeller_optimized/*/*.pdb\
    --output-dir ${OUTPUT_DIR}/assembly-foldx_optimized\
    --repair\
    --workers ${WORKERS}

echo ""
echo "=========================================="
echo "Stability analysis of all optimized assemblies"
echo "=========================================="

python ${PIPELINE_DIR}/foldx_wrapper.py \
    --command Stability \
    --input ${OUTPUT_DIR}/assembly-foldx_optimized/*/*.pdb \
    --output-dir ${OUTPUT_DIR}/stability \
    --workers 96


echo ""
echo "=========================================="
echo "Analyze stability results"
echo "=========================================="

python ${PIPELINE_DIR}/analyze_foldx.py \
    --input ${OUTPUT_DIR}/stability \
    --output ${OUTPUT_DIR}/stability/all_stability.csv \
    --assembly-totals ${OUTPUT_DIR}/stability/assembly_totals.csv \
    --viral-envelope ${OUTPUT_DIR}/stability/viral_envelope_stability.csv \
    --energy-column total_energy

python ${PIPELINE_DIR}/foldx_wrapper.py \
    --command AlaScan \
    --input ${OUTPUT_DIR}/assembly-foldx_optimized/*/*.pdb \
    --output-dir ${OUTPUT_DIR}/alascan \
    --workers ${WORKERS}



python ${PIPELINE_DIR}/foldx_wrapper.py \
    --command BuildModel \
    --input ${OUTPUT_DIR}/assembly-foldx_optimized/*/*.pdb \
    --mutations G52R A56V A170V T173I K200T M299I S305F P325S K331R T380R A407V A416T \
    --chains A B E F I J M N Q R   \
    --output-dir ${OUTPUT_DIR}/asibi_to_17d_mutations \
    --workers 96

PIPELINE_DIR=../viral_assembly_pipeline
OUTPUT_DIR=./viral_assembly_pipeline_results
python ${PIPELINE_DIR}/foldx_wrapper.py \
    --command BuildModel \
    --input ${OUTPUT_DIR}/assembly-foldx_optimized/*/*.pdb \
    --mutations R52G V56A V170A I173T T200K I299M F305S S325P R331K R380T V407A T416A \
    --chains A B E F I J M N Q R \
    --output-dir ${OUTPUT_DIR}/17d_to_asibi_mutations \
    --workers 96
