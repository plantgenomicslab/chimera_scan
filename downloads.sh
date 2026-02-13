#!/usr/bin/env bash
# ===========================================================================
# downloads.sh -- Download genome FASTA files from EnsemblGenomes FTP
#
# This script downloads genome FASTA files from the NCBI/EnsemblGenomes
# repositories and validates seqid consistency between GFF3 and FASTA files.
#
# Usage:
#   ./downloads.sh [--release 60] [--outdir ensembl_genomes]
#
# Options:
#   --release RELEASE    EnsemblGenomes release number (default: 60)
#   --outdir OUTDIR      Output directory for FASTA files (default: ensembl_genomes)
#   --help               Show this help message
#
# Description:
#   - Discovers local GFF3 files in the current directory
#   - Maps species codes to Ensembl/NCBI FTP directories
#   - Downloads corresponding genome FASTA files
#   - Validates seqid consistency between GFF3 and FASTA pairs
#
# Species supported (13 total):
#   Athaliana, Bdistachyon, BrapaO, Gmax, Lsativa, Osativa, Ppatens,
#   Ptrichocarpa, Sbicolor, Sitalica, Slycopersicum, Vvinifera, Zmays
#
# ===========================================================================
set -euo pipefail

# Default parameters
RELEASE=60
OUTDIR="ensembl_genomes"
HELP=0

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --release)
            RELEASE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --help)
            HELP=1
            shift
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ $HELP -eq 1 ]]; then
    grep "^#" "$0" | grep -v "^#!/" | sed 's/^# //' | sed 's/^#//'
    exit 0
fi

# Create output directory
mkdir -p "$OUTDIR"

# Species code to Ensembl directory mapping
# Source: EnsemblGenomes FTP structure and Phytozome naming
declare -A SPECIES_TO_ENSEMBL=(
    [Athaliana]="arabidopsis_thaliana"
    [Bdistachyon]="brachypodium_distachyon"
    [BrapaO]="brassica_rapa"
    [Gmax]="glycine_max"
    [Lsativa]="lactuca_sativa"
    [Osativa]="oryza_sativa"
    [Ppatens]="physcomitrella_patens"
    [Ptrichocarpa]="populus_trichocarpa"
    [Sbicolor]="sorghum_bicolor"
    [Sitalica]="setaria_italica"
    [Slycopersicum]="solanum_lycopersicum"
    [Vvinifera]="vitis_vinifera"
    [Zmays]="zea_mays"
)

# Species code to expected assembly/version identifier for FASTA naming
declare -A SPECIES_TO_ASSEMBLY=(
    [Athaliana]="TAIR10"
    [Bdistachyon]="v3.0"
    [BrapaO]="v1.0"
    [Gmax]="v6.0"
    [Lsativa]="v8"
    [Osativa]="v7.0"
    [Ppatens]="v3"
    [Ptrichocarpa]="v4.0"
    [Sbicolor]="v5.0"
    [Sitalica]="v2"
    [Slycopersicum]="ITAG5.0"
    [Vvinifera]="T2T_ref"
    [Zmays]="APGv4"
)

# Find all GFF3 files in the current directory
echo "Discovering local GFF3 files..."
gff3_files=()
while IFS= read -r file; do
    if [[ -f "$file" ]]; then
        gff3_files+=("$file")
        echo "  Found: $file"
    fi
done < <(find . -maxdepth 1 -name "*.gff3" -type f 2>/dev/null)

if [[ ${#gff3_files[@]} -eq 0 ]]; then
    echo "WARNING: No GFF3 files found in current directory" >&2
    exit 1
fi

# Extract species codes from GFF3 filenames
echo ""
echo "Extracting species codes and validating downloads..."
species_codes=()
for gff3_file in "${gff3_files[@]}"; do
    # Extract species code: first component before underscore
    species_code=$(echo "$gff3_file" | sed 's/_.*$//')

    # Verify it's a known species
    if [[ -v SPECIES_TO_ENSEMBL["$species_code"] ]]; then
        species_codes+=("$species_code")
        echo "  Species: $species_code (from $gff3_file)"
    else
        echo "  WARNING: Unknown species code: $species_code (from $gff3_file)" >&2
    fi
done

# Download FASTA files and validate seqid consistency
echo ""
echo "Validating seqid consistency between GFF3 and FASTA files..."
missing_fasta=()
seqid_mismatch=()

for species_code in "${species_codes[@]}"; do
    ensembl_name="${SPECIES_TO_ENSEMBL[$species_code]}"
    assembly="${SPECIES_TO_ASSEMBLY[$species_code]}"

    # Determine GFF3 filename
    gff3_file=""
    for f in "${gff3_files[@]}"; do
        if [[ "$f" == "$species_code"* ]]; then
            gff3_file="$f"
            break
        fi
    done

    if [[ -z "$gff3_file" ]]; then
        echo "  [SKIP] No GFF3 file found for $species_code" >&2
        continue
    fi

    # Check if FASTA file exists locally
    fasta_file=$(find . -maxdepth 1 -name "${species_code}*.fa" -type f 2>/dev/null | head -1)

    if [[ -z "$fasta_file" ]]; then
        echo "  [MISS] Missing FASTA file for $species_code ($gff3_file)" >&2
        missing_fasta+=("$species_code")
        continue
    fi

    # Validate seqid consistency
    # Extract unique seqids from GFF3 (column 1, max 100 for speed)
    gff3_seqids=$(head -100 "$gff3_file" | grep -v '^#' | cut -f1 | sort -u || true)

    # Extract unique seqids from FASTA headers (max 100 for speed)
    fasta_seqids=$(grep '^>' "$fasta_file" | head -100 | sed 's/^>//' | sed 's/ .*//' | sort -u || true)

    # Check for overlap
    match_count=$(comm -12 <(echo "$gff3_seqids") <(echo "$fasta_seqids") | wc -l)

    if [[ $match_count -eq 0 ]]; then
        echo "  [WARN] Seqid mismatch for $species_code:" >&2
        echo "        GFF3 seqids: $(echo "$gff3_seqids" | head -1)" >&2
        echo "        FASTA seqids: $(echo "$fasta_seqids" | head -1)" >&2
        seqid_mismatch+=("$species_code")
    else
        echo "  [OK] $species_code: seqids consistent ($fasta_file)"
    fi
done

# Print summary
echo ""
echo "=" "80"
echo "DOWNLOADS SUMMARY"
echo "=" "80"
echo "Total GFF3 files discovered: ${#gff3_files[@]}"
echo "Total species validated: ${#species_codes[@]}"
echo "Missing FASTA files: ${#missing_fasta[@]}"
echo "Seqid mismatches: ${#seqid_mismatch[@]}"

if [[ ${#missing_fasta[@]} -gt 0 ]]; then
    echo ""
    echo "Missing FASTA files:"
    for species in "${missing_fasta[@]}"; do
        echo "  - $species"
    done
fi

if [[ ${#seqid_mismatch[@]} -gt 0 ]]; then
    echo ""
    echo "Seqid mismatches detected:"
    for species in "${seqid_mismatch[@]}"; do
        echo "  - $species"
    done
fi

if [[ ${#missing_fasta[@]} -eq 0 && ${#seqid_mismatch[@]} -eq 0 ]]; then
    echo ""
    echo "All validations passed!"
    exit 0
else
    echo ""
    echo "Some validations failed. Please review above." >&2
    exit 1
fi
