#!/usr/bin/env bash
# ===========================================================================
# run_batch_slurm.sh -- Submit chimera_scan SLURM jobs for multiple species
#
# Loops through species results directories and submits one sbatch job per
# species. Each job runs chimera_scan.py with appropriate resources.
#
# Usage:
#   bash examples/run_batch_slurm.sh                   # all species
#   bash examples/run_batch_slurm.sh Gmax Athaliana    # specific species only
# ===========================================================================
set -euo pipefail

# ---- Paths (edit these) ---------------------------------------------------

# Base directory containing per-species results subdirectories
# Each subdirectory should have a reference_proteins.fasta file.
RESULTS_DIR="/path/to/results"

# chimera_scan.py location
SCRIPT="/path/to/chimera_scan/chimera_scan.py"

# SwissProt MMseqs2 database prefix (no extensions)
# Download from: https://doi.org/10.6084/m9.figshare.28236284
SWISSPROT_DB="/path/to/mmseqs_db/uniprot_sprot_odb10_plants"

# SwissProt FASTA with full headers (for Stage 4 annotation)
SWISSPROT_FASTA="/path/to/swissprot_plants.fasta"

# Directory for SLURM log files
LOGDIR="${RESULTS_DIR}/logs"

# ---- SLURM resource settings ----------------------------------------------

CPUS=16
MEMORY="64G"
TIME_LIMIT="7-00:00:00"     # 7 days

# ---- Environment activation command ---------------------------------------
# This command is prepended to each SLURM job to activate the conda env.
CONDA_INIT='eval "$(micromamba shell hook -s bash)" && micromamba activate AnnotationSplitter'

# ---- MMseqs2 settings -----------------------------------------------------

MMSEQS_PATH="mmseqs"
MMSEQS_THREADS="${CPUS}"

# ---- Species list ----------------------------------------------------------
# Default: all 13 plant species. Override by passing species names as arguments.

if [[ $# -gt 0 ]]; then
    SPECIES=("$@")
else
    SPECIES=(
        Athaliana
        Bdistachyon
        BrapaO
        Gmax
        Lsativa
        Osativa
        Ppatens
        Ptrichocarpa
        Sbicolor
        Sitalica
        Slycopersicum
        Vvinifera
        Zmays
    )
fi

# ---- Validation ------------------------------------------------------------

if [[ "${RESULTS_DIR}" == "/path/to/"* ]]; then
    echo "ERROR: Please edit the paths in this script before running." >&2
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: chimera_scan.py not found at ${SCRIPT}" >&2
    exit 1
fi

# ---- Create log directory --------------------------------------------------

mkdir -p "${LOGDIR}"

# ---- Submit jobs -----------------------------------------------------------

echo "Submitting chimera_scan jobs for ${#SPECIES[@]} species ..."
echo "  SLURM resources: ${CPUS} CPUs, ${MEMORY} RAM, ${TIME_LIMIT}"
echo "  Logs: ${LOGDIR}/"
echo ""

SUBMITTED=0
SKIPPED=0

for sp in "${SPECIES[@]}"; do
    SPD="${RESULTS_DIR}/${sp}"

    # Check that the species directory and reference FASTA exist
    if [[ ! -d "${SPD}" ]]; then
        echo "SKIP ${sp}: directory not found (${SPD})"
        ((SKIPPED++))
        continue
    fi

    REF_FASTA="${SPD}/reference_proteins.fasta"
    if [[ ! -f "${REF_FASTA}" ]]; then
        echo "SKIP ${sp}: no reference_proteins.fasta in ${SPD}"
        ((SKIPPED++))
        continue
    fi

    # Build the chimera_scan command
    RUN_CMD="${CONDA_INIT} && python ${SCRIPT}"
    RUN_CMD+=" --ref-fasta ${REF_FASTA}"
    RUN_CMD+=" --swissprot-db ${SWISSPROT_DB}"
    RUN_CMD+=" --swissprot-fasta ${SWISSPROT_FASTA}"
    RUN_CMD+=" --mmseqs-path ${MMSEQS_PATH}"
    RUN_CMD+=" --mmseqs-threads ${MMSEQS_THREADS}"
    RUN_CMD+=" --output-dir ${SPD}"
    RUN_CMD+=" --verbose"

    # If a pre-existing MMseqs2 result exists, reuse it
    MMSEQS_OUT="${SPD}/reference_proteins.mmseqs.out"
    if [[ -f "${MMSEQS_OUT}" ]]; then
        RUN_CMD+=" --mmseqs-out ${MMSEQS_OUT}"
    fi

    # Submit via sbatch
    JOB=$(sbatch \
        --job-name="CS_${sp}" \
        --output="${LOGDIR}/CS_${sp}_%j.out" \
        --error="${LOGDIR}/CS_${sp}_%j.err" \
        --cpus-per-task="${CPUS}" \
        --mem="${MEMORY}" \
        --time="${TIME_LIMIT}" \
        --parsable \
        --wrap "${RUN_CMD}")

    echo "  ${sp}  jobid=${JOB}"
    ((SUBMITTED++))
done

echo ""
echo "Submitted: ${SUBMITTED}  Skipped: ${SKIPPED}"
echo "Monitor with: squeue -u \$(whoami) -n CS_"
