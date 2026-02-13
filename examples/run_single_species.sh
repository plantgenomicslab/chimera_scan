#!/usr/bin/env bash
# ===========================================================================
# run_single_species.sh -- Run chimera_scan on a single species
#
# This script demonstrates how to run chimera_scan.py for one species.
# Adjust the paths below to match your local setup.
#
# Usage:
#   bash examples/run_single_species.sh
# ===========================================================================
set -euo pipefail

# ---- Paths (edit these) ---------------------------------------------------

# Directory containing chimera_scan.py
CHIMERA_SCAN_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# chimera_scan.py location
SCRIPT="${CHIMERA_SCAN_DIR}/chimera_scan.py"

# Reference proteome FASTA (from AnnotationSplitter output)
# This file contains the translated longest-isoform proteins.
REF_FASTA="/path/to/results/Gmax/reference_proteins.fasta"

# MMseqs2 SwissProt database prefix (no file extensions)
# Download from: https://doi.org/10.6084/m9.figshare.28236284
SWISSPROT_DB="/path/to/mmseqs_db/uniprot_sprot_odb10_plants"

# SwissProt FASTA with full headers (optional but recommended for Stage 4)
SWISSPROT_FASTA="/path/to/swissprot_plants.fasta"

# Output directory for this species
OUTPUT_DIR="/path/to/results/Gmax/chimera_scan_output"

# (Optional) Pre-existing MMseqs2 result from a previous run.
# Uncomment to skip Stage 1 if the search was already done.
# MMSEQS_OUT="/path/to/results/Gmax/reference_proteins.mmseqs.out"

# ---- MMseqs2 settings -----------------------------------------------------

MMSEQS_PATH="mmseqs"       # Path to mmseqs binary (or just "mmseqs" if on PATH)
MMSEQS_THREADS=16           # Number of threads for MMseqs2 searches

# ---- Analysis parameters (defaults shown) ---------------------------------

EVALUE=1e-5                 # E-value cutoff for Pattern B filtering
MIN_HITS=5                  # Minimum hits per query for Pattern B analysis
MIN_COVERAGE=0.80           # Minimum best-target subject coverage per domain
OUTPUT_PREFIX="chimera_fullscan"  # Prefix for output filenames

# ---- Environment activation -----------------------------------------------
# Activate the conda/micromamba environment containing pandas, biopython, mmseqs2.
# Uncomment and adjust one of the following:

# eval "$(micromamba shell hook -s bash)" && micromamba activate AnnotationSplitter
# conda activate AnnotationSplitter

# ---- Validation ------------------------------------------------------------

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: chimera_scan.py not found at ${SCRIPT}" >&2
    exit 1
fi

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "ERROR: Reference FASTA not found at ${REF_FASTA}" >&2
    exit 1
fi

if [[ "${SWISSPROT_DB}" == "/path/to/"* ]]; then
    echo "ERROR: Please edit the paths in this script before running." >&2
    exit 1
fi

# ---- Run chimera_scan ------------------------------------------------------

echo "Running chimera_scan for a single species ..."
echo "  Reference FASTA : ${REF_FASTA}"
echo "  SwissProt DB    : ${SWISSPROT_DB}"
echo "  Output dir      : ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"

# Build the command
CMD=(
    python "${SCRIPT}"
    --ref-fasta "${REF_FASTA}"
    --swissprot-db "${SWISSPROT_DB}"
    --swissprot-fasta "${SWISSPROT_FASTA}"
    --mmseqs-path "${MMSEQS_PATH}"
    --mmseqs-threads "${MMSEQS_THREADS}"
    --evalue "${EVALUE}"
    --min-hits "${MIN_HITS}"
    --min-coverage "${MIN_COVERAGE}"
    --output-prefix "${OUTPUT_PREFIX}"
    --output-dir "${OUTPUT_DIR}"
    --verbose
)

# If a pre-existing MMseqs2 result is available, add it
if [[ -n "${MMSEQS_OUT:-}" && -f "${MMSEQS_OUT}" ]]; then
    CMD+=(--mmseqs-out "${MMSEQS_OUT}")
    echo "  Reusing MMseqs2 result: ${MMSEQS_OUT}"
fi

echo "Command:"
echo "  ${CMD[*]}"
echo ""

"${CMD[@]}"

echo ""
echo "Done. Output files are in: ${OUTPUT_DIR}/"
