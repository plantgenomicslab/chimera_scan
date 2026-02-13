# Requirements

## Python

- **Python 3.11+** (tested with Python 3.11.5)

### Python packages

| Package    | Minimum version | Purpose                          |
|------------|-----------------|----------------------------------|
| pandas     | 2.0+            | Tabular data manipulation        |
| biopython  | 1.80+           | FASTA parsing and sequence I/O   |

Install with pip:

```bash
pip install pandas biopython
```

Or use the pre-configured conda environment:

```bash
micromamba activate AnnotationSplitter
```

## External tools

### MMseqs2

MMseqs2 is required for sequence similarity searches. It must be available on
your `PATH` or specified via `--mmseqs-path`.

Install MMseqs2:

```bash
# conda / micromamba
conda install -c bioconda mmseqs2

# or from source: https://github.com/soedinglab/MMseqs2
```

## Databases

### SwissProt MMseqs2 database (required)

The pre-built MMseqs2 SwissProt database is used for the primary search
(Stage 1) and is specified with `--swissprot-db`.

Download from: https://doi.org/10.6084/m9.figshare.28236284

After downloading, the database files should be in a single directory. Pass the
database prefix (without file extensions) to `--swissprot-db`. For example, if
the files are:

```
mmseqs_db/uniprot_sprot_odb10_plants
mmseqs_db/uniprot_sprot_odb10_plants.index
mmseqs_db/uniprot_sprot_odb10_plants.dbtype
mmseqs_db/uniprot_sprot_odb10_plants_h
mmseqs_db/uniprot_sprot_odb10_plants_h.index
mmseqs_db/uniprot_sprot_odb10_plants_h.dbtype
```

Then use `--swissprot-db mmseqs_db/uniprot_sprot_odb10_plants`.

### SwissProt FASTA file (optional, recommended)

A UniProt/SwissProt FASTA file with full headers is needed for Stage 4
(description annotation). This enables gene name extraction and protein
description reporting in the final summary.

Specify with `--swissprot-fasta`. The file should contain UniProt-format headers
with `GN=` and `OS=` fields, e.g.:

```
>sp|Q9FIB0|... Protein name OS=Arabidopsis thaliana ... GN=GENE1 ...
```

A suitable file (`swissprot_plants.fasta` or `uniprot_sprot_odb10_plants.fasta`)
can be obtained from the same Figshare repository linked above.

## Hardware recommendations

- **CPU**: 4+ cores recommended (MMseqs2 uses `--mmseqs-threads`, default 16)
- **RAM**: 8 GB minimum; 16+ GB recommended for large proteomes
- **Disk**: ~5 GB for databases + working space for intermediate MMseqs2 files
