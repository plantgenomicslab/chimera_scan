# chimera_scan: Whole-Proteome Chimeric Gene Detection Pipeline

Detect chimeric (fused) genes in annotated proteomes by identifying proteins
whose SwissProt homology profile is bimodal -- the N-terminal and C-terminal
halves match distinct protein families (Pattern B). chimera_scan validates
candidates with self-BLAST against the reference proteome and annotates them
with non-overlapping SwissProt protein tiling.

## Features

- **Pattern B bimodal detection** -- identifies query proteins where SwissProt
  targets cluster into separate N-terminal and C-terminal groups, indicating a
  likely gene fusion event.
- **Subject coverage filtering** -- requires >80% target coverage on the best
  SwissProt match in each domain, reducing false positives from partial
  alignments.
- **Self-BLAST validation** -- searches candidates against the full reference
  proteome (self-excluded) to test whether the bimodal signal is also present
  among paralogs.
- **SwissProt annotation with non-overlapping protein tiling** -- greedily
  selects the best-scoring non-overlapping UniProt matches along the query,
  reporting the number of distinct genes tiled and their coordinates.
- **Streaming chunked processing** -- reads large MMseqs2 result files in
  configurable chunks (default 2M rows), keeping memory usage manageable for
  whole-proteome scans.
- **Idempotent re-runs** -- skips stages whose output files already exist,
  allowing safe restarts after interruption.

## Requirements

- **Python 3.11+** (tested in a micromamba/conda environment)
- **Python packages**: `pandas`, `biopython`
- **MMseqs2** (the `easy-search` subcommand must be available)
- **SwissProt database** -- an MMseqs2-formatted database of UniProt/SwissProt
  sequences (see [Database preparation](#database-preparation) below)

## Installation

### Environment setup

```bash
# Create and activate a conda/micromamba environment
micromamba create -n chimera_scan python=3.11 pandas biopython -c conda-forge
micromamba activate chimera_scan

# Install MMseqs2
micromamba install -c conda-forge -c bioconda mmseqs2
```

Alternatively, activate the existing AnnotationSplitter environment if it is
available:

```bash
micromamba activate AnnotationSplitter
```

### Database preparation

chimera_scan requires two SwissProt-related files:

1. **MMseqs2 SwissProt database** -- download the pre-built database from
   [Figshare](https://doi.org/10.6084/m9.figshare.28236284) and extract it.
   The `--swissprot-db` flag should point to the database prefix (e.g.,
   `mmseqs_db/uniprot_sprot_odb10_plants`).

2. **SwissProt FASTA with headers** -- a plain FASTA file of the same SwissProt
   sequences, used by Stage 4 to extract protein descriptions and gene names.
   Pass this via `--swissprot-fasta`.

## Usage

```bash
python chimera_scan.py \
    --ref-fasta reference_proteins.fasta \
    --swissprot-db mmseqs_db/uniprot_sprot_odb10_plants \
    --swissprot-fasta swissprot_plants.fasta \
    --mmseqs-path mmseqs \
    --mmseqs-threads 16 \
    --verbose
```

### Parameters

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--ref-fasta` | Yes | -- | Reference proteome FASTA file |
| `--swissprot-db` | Yes | -- | MMseqs2 SwissProt database or FASTA (search target) |
| `--swissprot-fasta` | No | -- | UniProt SwissProt FASTA with full headers (enables Stage 4 annotation) |
| `--mmseqs-out` | No | -- | Pre-existing MMseqs2 result TSV (skips Stage 1 if provided) |
| `--mmseqs-path` | No | `mmseqs` | Path to the mmseqs binary |
| `--mmseqs-threads` | No | `16` | Number of threads for MMseqs2 searches |
| `--evalue` | No | `1e-5` | E-value cutoff for Pattern B hit filtering |
| `--min-hits` | No | `5` | Minimum number of hits per query to attempt Pattern B analysis |
| `--min-coverage` | No | `0.80` | Minimum best-target subject coverage per domain |
| `--chunk-size` | No | `2000000` | Rows per chunk for the streaming Pattern B scan |
| `--output-prefix` | No | `chimera_fullscan` | Prefix for all output file names |
| `--output-dir` | No | `.` | Directory for output files |
| `--verbose` | No | off | Print additional progress information |

### Input file requirements

- **Reference proteome FASTA** -- one protein sequence per gene model. Standard
  FASTA format with unique sequence IDs. This is typically the translated
  longest isoform per gene from a genome annotation.
- **MMseqs2 SwissProt database** -- created with `mmseqs createdb` from a
  SwissProt FASTA, or download the pre-built database linked above.

### Pipeline stages

1. **Stage 1 -- MMseqs2 search**: Reference proteome vs SwissProt database.
   Produces a tab-separated alignment file. Skipped if `--mmseqs-out` is
   provided.
2. **Stage 2 -- Streaming Pattern B scan**: Reads the MMseqs2 output in chunks,
   filters by e-value, then tests each query for a bimodal target distribution.
   Applies the `--min-coverage` threshold on the best SwissProt target in each
   domain.
3. **Stage 3 -- Self-BLAST**: Extracts candidate sequences and searches them
   against the full reference proteome (self-hits excluded). Repeats Pattern B
   analysis on the self-BLAST results.
4. **Stage 4 -- SwissProt annotation** (requires `--swissprot-fasta`): Searches
   candidates against the SwissProt FASTA to obtain full headers, then performs
   non-overlapping protein tiling to count distinct gene products.
5. **Combined summary**: Merges evidence from all stages into a single table
   ranked by evidence count.

## Examples

### Single species run

```bash
python chimera_scan.py \
    --ref-fasta results/Gmax/reference_proteins.fasta \
    --swissprot-db mmseqs_db/uniprot_sprot_odb10_plants \
    --swissprot-fasta swissprot_plants.fasta \
    --mmseqs-out results/Gmax/reference_proteins.mmseqs.out \
    --output-dir results/Gmax \
    --verbose
```

### Batch SLURM submission

See `examples/` for a SLURM batch script that submits one job per species:

```bash
bash examples/run_chimera_scan.sh
```

The batch script iterates over species directories, checks for required input
files, and submits each species as an independent SLURM job with 16 CPUs and
64 GB memory.

## Output Files

All output files are written to `--output-dir` with the `--output-prefix` as
their name prefix (default: `chimera_fullscan`).

| File | Description |
|------|-------------|
| `<prefix>.mmseqs.tsv` | Raw MMseqs2 search results (Stage 1) |
| `<prefix>.pattern_b.summary.tsv` | Pattern B candidates with per-domain statistics (Stage 2) |
| `<prefix>.pattern_b.details.tsv` | Per-hit details for Pattern B candidates, including domain assignment (Stage 2) |
| `<prefix>.candidates.fasta` | FASTA sequences of all Pattern B candidates |
| `<prefix>.self_blast.raw.tsv` | Raw self-BLAST MMseqs2 output (Stage 3) |
| `<prefix>.self_blast.summary.tsv` | Self-BLAST Pattern B summary (Stage 3) |
| `<prefix>.self_blast.details.tsv` | Self-BLAST per-hit details (Stage 3) |
| `<prefix>.vs_swissprot_described.tsv` | SwissProt search with full headers (Stage 4) |
| `<prefix>.combined_summary.tsv` | Final combined summary ranked by evidence (all stages) |

### Column descriptions

#### pattern_b.summary.tsv

| Column | Description |
|--------|-------------|
| `query_id` | Protein identifier from the reference FASTA |
| `qlen` | Query protein length (amino acids) |
| `n_nterm_targets` | Number of SwissProt targets whose mean alignment midpoint falls in the N-terminal region (<35% of query length) |
| `n_cterm_targets` | Number of targets in the C-terminal region (>65% of query length) |
| `n_bridge_targets` | Number of targets in the bridge zone (35-65%) |
| `nterm_qcov` | Fraction of query covered by N-terminal domain alignments |
| `cterm_qcov` | Fraction of query covered by C-terminal domain alignments |
| `total_qcov` | Fraction of query covered by all domain alignments combined |
| `nterm_best_tcov` | Best target subject coverage in the N-terminal domain |
| `cterm_best_tcov` | Best target subject coverage in the C-terminal domain |
| `nterm_range` | Query coordinate range spanned by N-terminal alignments |
| `cterm_range` | Query coordinate range spanned by C-terminal alignments |

#### combined_summary.tsv

| Column | Description |
|--------|-------------|
| `candidate` | Protein identifier |
| `qlen` | Query protein length |
| `nterm_best_tcov` | Best target coverage, N-terminal domain (from SwissProt Pattern B) |
| `cterm_best_tcov` | Best target coverage, C-terminal domain (from SwissProt Pattern B) |
| `sp_tcov` | Minimum target coverage across non-overlapping SwissProt proteins |
| `swissprot_n_genes` | Number of distinct non-overlapping SwissProt proteins tiling the query |
| `pattern_b_vs_swissprot` | Whether the candidate shows Pattern B vs SwissProt |
| `pattern_b_vs_self` | Whether the candidate shows Pattern B vs self-BLAST |
| `evidence_count` | Number of evidence sources supporting chimera call (max 2) |
| `best_sp_hit` | Best overall SwissProt hit accession |
| `best_sp_gene` | Gene name of the best SwissProt hit |
| `best_sp_evalue` | E-value of the best SwissProt hit |
| `nterm_sp_hit` | Best N-terminal SwissProt hit accession |
| `nterm_sp_gene` | Gene name of the N-terminal hit |
| `nterm_qstart` / `nterm_qend` | Query coordinates of N-terminal SwissProt alignments |
| `cterm_sp_hit` | Best C-terminal SwissProt hit accession |
| `cterm_sp_gene` | Gene name of the C-terminal hit |
| `cterm_qstart` / `cterm_qend` | Query coordinates of C-terminal SwissProt alignments |
| `sp_genes_detail` | Non-overlapping protein tiling detail: `gene(qstart-qend,tcov=X.XX); ...` |
| `best_sp_desc` | Description of the best overall SwissProt hit |
| `nterm_sp_desc` | Description of the best N-terminal SwissProt hit |
| `cterm_sp_desc` | Description of the best C-terminal SwissProt hit |

## How it works

A chimeric gene is a mis-annotation where two or more independent genes are
erroneously merged into a single gene model. chimera_scan detects these by
looking for a bimodal homology signature:

1. For each query protein, all SwissProt hits (passing e-value and coverage
   filters) are grouped by target. Each target's mean alignment midpoint on the
   query is computed.
2. Targets are classified as N-terminal (midpoint < 35% of query length),
   C-terminal (midpoint > 65%), or bridge (in between).
3. A query is flagged as a Pattern B candidate when it has at least 3 targets in
   both the N-terminal and C-terminal groups, and the number of bridge targets
   is fewer than the smaller domain group.
4. The best target in each domain must have >80% subject coverage
   (`--min-coverage`), ensuring that each domain aligns to a near-full-length
   known protein.

This approach reliably identifies cases where a single annotated protein is
actually composed of two (or more) distinct protein families joined end to end.

## Citation

If you use chimera_scan in your research, please cite this repository. chimera_scan
was developed as a companion tool to
[AnnotationSplitter](https://github.com/Andy-B-123/AnnotationSplitter) for
large-scale chimeric gene detection in plant genomes.

## License

See LICENSE file for details.
