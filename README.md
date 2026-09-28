# Sma3s v3

[![Conda](https://img.shields.io/conda/vn/PMC_fps/sma3s.svg)](https://anaconda.org/PMC_fps/sma3s)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

**Sequence Massive Annotator using 3 modules**

Sma3s is a command-line suite for large-scale functional annotation, inspection
of high-confidence sequence matches, taxonomic reporting and functional
enrichment.

The Conda package installs four coordinated commands:

| Command | Purpose |
|---|---|
| `sma3s` | Run the three-module functional annotation workflow. |
| `download_sma3s_db` | Download the ready-to-use Sma3s reference database and SQLite caches. |
| `sma3s-hit-filter` | Select hits by identity/coverage thresholds and optionally attach UniProt species and TaxID information. |
| `sma3s-enrichment` | Perform GO and UniProt Keyword enrichment and optional genus-abundance summaries. |

> [!NOTE]
> This README assumes that the Conda recipe exposes the scripts with the command
> names shown above. The recommended source filename for the current
> `filter_perfect_hits_uniprot_sqlite.py` script is `sma3s_hit_filter.py`, with
> the public command `sma3s-hit-filter`.

---

## Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [End-to-end workflow](#end-to-end-workflow)
- [Sma3s annotation](#sma3s-annotation)
- [Database download](#database-download)
- [Hit filtering and taxonomy](#hit-filtering-and-taxonomy)
- [Functional enrichment](#functional-enrichment)
- [Input files](#input-files)
- [Output files](#output-files)
- [Reference databases](#reference-databases)
- [Performance and parallelization](#performance-and-parallelization)
- [Caches and restart behavior](#caches-and-restart-behavior)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Issues](#issues)

---

## Overview

Sma3s combines three complementary annotation strategies:

1. **Annotator 1 (A1)** transfers annotation from a direct,
   high-identity and high-coverage reference hit.
2. **Annotator 2 (A2)** identifies probable orthologues through reciprocal
   best-hit analysis.
3. **Annotator 3 (A3)** detects functional terms enriched among multiple
   homologous reference sequences.

The main annotation workflow can recover:

- gene names;
- protein descriptions;
- Enzyme Commission numbers;
- Gene Ontology terms;
- UniProt Keywords;
- UniProt pathways;
- GO Slim terms;
- annotation provenance;
- supporting hit identities, coverages, alignment lengths and E-values.

The auxiliary commands extend the core workflow:

<img width="1672" height="941" alt="sma3s_workflow" src="https://github.com/user-attachments/assets/711d9448-36ed-43c1-bc83-268c20b0dd07" />


## Features

### Annotation

- Protein and translated nucleotide searches with MMseqs2.
- Independent or combined execution of A1, A2 and A3.
- Direct support for UniProt `.dat`, `.dat.gz`, and prepared
  FASTA/`.annot` pairs.
- Taxonomic inclusion and exclusion by genus, family or order.
- Optional GO namespace columns, GO Slim terms and evidence filtering.
- Optional reporting of the exact reference hits supporting each annotation.

### Scalability

- Parallel MMseqs2 searches.
- Parallel UniProt parsing.
- Parallel per-query annotation.
- Parallel `.gz` decompression through `rapidgzip`.
- SQLite-backed FASTA indexes, annotations, hits and reciprocal-best-hit maps.
- Sequential or indexed FASTA extraction for reciprocal searches.
- Reuse of completed searches and current caches.

### Downstream analysis

- Resumable parallel download of the prepared Sma3s database.
- Configurable filtering by identity, query coverage and target coverage.
- UniProt species and NCBI TaxID recovery through a local SQLite index.
- GO and Keyword enrichment using the complete annotation table as universe.
- Benjamini-Hochberg false-discovery-rate correction.
- GO Biological Process, Molecular Function and Cellular Component summaries.
- Dot plots, bar plots and optional Keyword treemaps.
- Optional genus-abundance and per-gene genus-richness summaries.

---

## Installation

Sma3s is distributed through the `pmc_fps` Conda channel.

### Mamba

```bash
mamba create -n sma3s \
    -c pmc_fps \
    -c conda-forge \
    -c bioconda \
    sma3s=3
```

### Conda

```bash
conda create -n sma3s \
    -c pmc_fps \
    -c conda-forge \
    -c bioconda \
    sma3s=3
```

Activate the environment:

```bash
conda activate sma3s
```

### Verify the installation

```bash
sma3s --check-install
```

The diagnostic checks:

- Python;
- MMseqs2 availability and execution;
- SQLite;
- write access to the temporary directory;
- standard gzip support;
- optional parallel gzip support through `rapidgzip`.

A successful core installation ends with:

```text
Installation check: OK
```

A custom MMseqs2 executable or temporary directory can be tested with:

```bash
sma3s --check-install \
    --mmseqs-bin /path/to/mmseqs \
    --tmpdir /scratch/$USER/sma3s
```

### Package dependencies

The Conda recipe should provide the dependencies required by the complete suite:

```text
python
mmseqs2
rapidgzip
pandas
scipy
matplotlib
squarify
```

`rapidgzip` is optional in the source script because standard gzip is available
as a fallback, but including it in the Conda package enables parallel
decompression by default.

`squarify` is only required for Keyword treemaps, but including it ensures that
all documented enrichment outputs are available.

---

## Quick start

### 1. Download a prepared database

For example, download the bacterial reference:

```bash
download_sma3s_db --db bacteria
```

By default, taxon-specific databases are stored under:

```text
db_sma3s/
```

### 2. Annotate a protein FASTA

```bash
sma3s \
    -i proteins.faa \
    -d db_sma3s/bacteria/uniprot_bacteria.fasta \
    -num_threads 16
```

### 3. Record the annotation-supporting hits

Use `--report-hit-metrics` when the output will be inspected with
`sma3s-hit-filter`:

```bash
sma3s \
    -i proteins.faa \
    -d db_sma3s/bacteria/uniprot_bacteria.fasta \
    -num_threads 16 \
    -source \
    --report-hit-metrics
```

### 4. Select strict 100/100/100 matches

```bash
sma3s-hit-filter filter \
    -i proteins_uniprot_sprot_trembl_all_source_hitmetrics.tsv \
    -o strict_hits
```

### 5. Run functional enrichment

Create a text file containing one query gene ID per line:

```text
gene_001
gene_014
gene_205
```

Then run:

```bash
sma3s-enrichment \
    -a proteins_uniprot_sprot_trembl_all.tsv \
    -i genes_of_interest.txt \
    -d enrichment \
    -o genes_of_interest \
    --go-obo go-basic.obo
```

---

## End-to-end workflow

This example annotates a pangenome, identifies strict cross-reference matches,
assigns species to those hits and performs functional and taxonomic summaries.

### Step 1: download the prepared reference

```bash
download_sma3s_db \
    --db bacteria \
    --workers 2
```

### Step 2: annotate the protein catalog

```bash
sma3s \
    -i pangenome_proteins.faa \
    -d db_sma3s/bacteria/uniprot_bacteria.fasta \
    -num_threads 32 \
    --annotation-workers 8 \
    -source \
    --report-hit-metrics
```

Do not use `-go` in this run when the same table will be used directly by
`sma3s-enrichment`. The enrichment command expects one combined `GO` column and
can split GO terms into namespaces using `--go-obo`.

### Step 3: build a UniProt species index

A species index is only required when species and TaxID information should be
added to selected hits.

```bash
sma3s-hit-filter build-index \
    --uniprot-dat uniprot_sprot_trembl_all.dat.gz \
    --species-db uniprot_sprot_trembl_all.species.sqlite
```

The index stores both UniProt entry names and accessions.

### Step 4: select strict hits

```bash
sma3s-hit-filter filter \
    -i pangenome_proteins_uniprot_sprot_trembl_all_source_hitmetrics.tsv \
    -o pangenome_100_100_100 \
    --species-db uniprot_sprot_trembl_all.species.sqlite --min-identity 100 --min-qcov 100 --min-tcov 100
```

The default criteria are:

```text
identity >= 100%
query coverage >= 100%
target coverage >= 100%
```

### Step 5: perform enrichment and genus summaries

```bash
sma3s-enrichment \
    -a pangenome_proteins_uniprot_sprot_trembl_all_source_hitmetrics.tsv \
    -i candidate_genes.txt \
    -d candidate_enrichment \
    -o candidates \
    --go-obo go-basic.obo \
    --taxa-hits pangenome_100_100_100.perfect_hits.long.tsv \
    --min-term-size 2 \
    --alpha 0.05 \
    --top 25 \
    --use-namespace-fdr-for-namespace-plots \
    --make-keyword-treemap
```

The taxonomic component reports genus abundance and per-gene genus diversity.
It does not perform a genus-enrichment test.

---

# Sma3s annotation

## Basic command

```bash
sma3s -i QUERY_FASTA -d REFERENCE -num_threads CPUS
```

Protein input is assumed unless `-nucl` is supplied.

## Reference input

`-d` accepts:

| Input | Behavior |
|---|---|
| `reference.dat.gz` | Decompress or reuse `reference.dat`, then create or reuse `reference.fasta` and `reference.annot`. |
| `reference.dat` | Create or reuse `reference.fasta` and `reference.annot`. |
| `reference.fasta` | Use the FASTA directly; `reference.annot` must exist beside it. |

## Annotation modules

Choose modules with `-a`:

| Value | Modules |
|---|---|
| `1` | A1 only |
| `2` | A2 only |
| `3` | A3 only |
| `12` | A1 + A2 |
| `13` | A1 + A3 |
| `23` | A2 + A3 |
| `123` | A1 + A2 + A3; default |

## Default thresholds

### Initial MMseqs2 search

```text
E-value <= 1e-6
Sensitivity = 7.5
Maximum retained hits = 250
Low-complexity masking = disabled
```

### Annotator 1

```text
Identity >= 90%
Target/reference coverage >= 90%
P value <= 0.1
```

### Annotator 2

```text
Identity >= Rost(20, alignment length)
Target/reference coverage >= 80%
P value <= 0.1
Reciprocal best hit required
```

If any identity or coverage threshold is supplied manually, A2 changes to fixed
thresholds. Its fixed identity default becomes 75% unless `-id2` is provided.

### Annotator 3

```text
Identity > Rost(20, alignment length)
No explicit coverage threshold
Hypergeometric enrichment P value <= 0.1
```

## Common examples

### All annotation modules

```bash
sma3s \
    -i proteins.faa \
    -d reference.dat.gz \
    -a 123 \
    -num_threads 16
```

### Direct high-confidence transfer only

```bash
sma3s \
    -i proteins.faa \
    -d reference.dat.gz \
    -a 1 \
    -id1 90 \
    -cov1 90 \
    -num_threads 16
```

### Translated nucleotide search

```bash
sma3s \
    -i transcripts.fna \
    -d reference.dat.gz \
    -nucl \
    -num_threads 16
```

When no manual thresholds are supplied, nucleotide-mode identity and coverage
thresholds are reduced to 51%.

### Include detailed support information

```bash
sma3s \
    -i proteins.faa \
    -d reference.dat.gz \
    -source \
    --report-hit-metrics \
    -num_threads 16
```

### Select a taxonomic subset of UniProt

```bash
sma3s \
    -i proteins.faa \
    -d uniprot_trembl.dat.gz \
    --family Enterobacteriaceae \
    -num_threads 32
```

### Exclude a genus

```bash
sma3s \
    -i proteins.faa \
    -d uniprot_trembl.dat.gz \
    --exclude-genus Vibrio \
    -num_threads 32
```

Taxonomic exclusion automatically enables hit-metric columns.

### Large search with memory control

```bash
sma3s \
    -i catalog.faa \
    -d uniprot_trembl.dat.gz \
    -num_threads 64 \
    --annotation-workers 16 \
    --split-memory-limit 180G \
    --tmpdir /scratch/$USER/sma3s
```

## Main options

### Input and workflow

| Option | Default | Description |
|---|---:|---|
| `--check-install` | disabled | Verify the installation and exit. |
| `-i FILE` | required | Query FASTA. |
| `-d FILE` | required | `.dat.gz`, `.dat`, or FASTA/`.annot` reference. |
| `-a VALUE` | `123` | Annotation-module combination. |
| `-b FILE` | generated | Existing or target MMseqs2 tabular file. |

### Output and annotation content

| Option | Default | Description |
|---|---:|---|
| `-source` | disabled | Add `ANNOTATOR` and `USED_UNIPROT_SEQUENCES`. |
| `--report-hit-metrics` | disabled | Add hit ID, identity, query coverage, target coverage, alignment length and E-value. |
| `-go` | disabled | Split GO into BP, MF and CC identifier/name columns. |
| `-goslim` | disabled | Add a GO Slim column and summary section. |
| `-noempty` | disabled | Omit unannotated queries. |
| `-quality` | disabled | Apply the Sma3s ECO/IEA evidence filter. |
| `-nopred` | disabled | Discard UniProt entries with `PE=4` while building the reference. |

### Similarity thresholds

| Option | Default | Description |
|---|---:|---|
| `-r FLOAT` | `20` | N parameter of the Rost identity curve. |
| `-p FLOAT` | `0.1` | Maximum annotation/enrichment P value. |
| `-id1 FLOAT` | `90` | Minimum A1 identity. |
| `-cov1 FLOAT` | `90` | Minimum A1 target coverage. |
| `-id2 FLOAT` | Rost curve | Minimum A2 identity; fixed default is 75 after manual threshold customization. |
| `-cov2 FLOAT` | `80` | Minimum A2 target coverage. |

### Reference interpretation and taxonomy

| Option | Default | Description |
|---|---:|---|
| `-uniprot` | disabled | Treat the reference as UniProt rather than UniRef and enable A3 clustering. |
| `--genus NAME` | disabled | Keep an exact genus lineage node. |
| `--exclude-genus NAME` | disabled | Remove an exact genus lineage node. |
| `--family NAME` | disabled | Keep an exact family lineage node. |
| `--exclude-family NAME` | disabled | Remove an exact family lineage node. |
| `--order NAME` | disabled | Keep an exact order lineage node. |
| `--exclude-order NAME` | disabled | Remove an exact order lineage node. |

Taxonomic filters require `.dat` or `.dat.gz` input because they are evaluated
from UniProt `OC` lineage records.

### MMseqs2 and parallelism

| Option | Default | Description |
|---|---:|---|
| `-num_threads INT` | `1` | MMseqs2 threads and upper limit for derived worker counts. |
| `--annotation-workers INT` | `max(1, threads // 4)` | Python processes for database parsing and final annotation. |
| `--cluster-threads INT` | derived | Threads per A3 clustering job. |
| `--decompression-threads INT` | `-num_threads` | Threads used by `rapidgzip`. |
| `--max-seqs INT` | `250` | Maximum initial hits retained per query. |
| `--sensitivity FLOAT` | `7.5` | MMseqs2 sensitivity. |
| `--split-memory-limit VALUE` | unset | MMseqs2 memory limit, for example `80G`. |
| `--db-load-mode {0,1,2}` | unset | Explicit MMseqs2 database load mode. |
| `--mmseqs-bin PATH` | `mmseqs` | Executable name or path. |
| `-filter` | disabled | Enable MMseqs2 low-complexity masking. |
| `--force-search` | disabled | Repeat the initial search instead of reusing its output. |

### Reciprocal search

| Option | Default | Description |
|---|---:|---|
| `--reciprocal-fasta-mode {auto,index,stream}` | `auto` | Construct the reciprocal FASTA through indexed access or sequential scanning. |
| `--reciprocal-stream-threshold INT` | `50000` | Candidate count at which `auto` switches to sequential scanning. |

### Temporary data and cleanup

| Option | Default | Description |
|---|---:|---|
| `--tmpdir DIR` | `./sma3s_mmseqs_tmp` | Temporary directory. |
| `--keep-tmp` | disabled | Retain temporary files. |
| `--force-clean-outputs` | disabled | Remove outputs and run-specific caches. |
| `--force-clean-db` | disabled | Remove files and caches derived from `.dat` or `.dat.gz`. |
| `--force-clean` | disabled | Apply both cleanup modes. |
| `--clean-only` | disabled | Clean and exit without annotation. |

---

# Database download

Sma3s requires a UniProt-based reference database containing protein sequences,
annotations, and two pre-built SQLite cache files. Pre-built databases can be
downloaded directly using `download_sma3s_db.py`.

## Available databases

The downloader supports taxon-specific UniProt databases. Database availability
is checked directly against the Sma3s database server.

Currently prepared databases include:

| Database | Description |
|---|---|
| `bacteria` | UniProt bacterial proteins |
| `archaea` | UniProt archaeal proteins |

Additional databases, including viral, fungal, plant, and metazoan datasets,
will be progressively added.

The current status of all supported databases can be checked with:

```bash
download_sma3s_db.py --list-dbs
```

Databases for which one or more required files are not yet available are
reported as:

```text
UNDER PREPARATION
```

## Download a database

To download the bacterial database:

```bash
download_sma3s_db.py --db bacteria
```

To download the archaeal database:

```bash
download_sma3s_db.py --db archaea
```

By default, databases are stored in:

```text
db_sma3s/
```

with a separate directory for each taxonomic group:

```text
db_sma3s/
├── bacteria/
│   ├── uniprot_bacteria.fasta
│   ├── uniprot_bacteria.annot
│   ├── uniprot_bacteria.fasta.sma3s_fasta_index.sqlite
│   └── uniprot_bacteria.annot.q0_go0_goslim0.sma3s_annot.sqlite
│
└── archaea/
    ├── uniprot_archaea.fasta
    ├── uniprot_archaea.annot
    ├── uniprot_archaea.fasta.sma3s_fasta_index.sqlite
    └── uniprot_archaea.annot.q0_go0_goslim0.sma3s_annot.sqlite
```

A custom output directory can be specified with:

```bash
download_sma3s_db.py --db bacteria --output /path/to/databases
```

## Required database files

Each Sma3s reference database consists of four files:

```text
<prefix>.fasta
<prefix>.annot
<prefix>.fasta.sma3s_fasta_index.sqlite
<prefix>.annot.q0_go0_goslim0.sma3s_annot.sqlite
```

The FASTA and annotation files contain the UniProt reference data, whereas the
SQLite files contain pre-built Sma3s caches used to accelerate sequence
retrieval and annotation processing.

A database is considered available only when all four files are present on the
server. If any required file is missing, the downloader reports the database as
being under preparation and does not start an incomplete download.

## Parallel downloads

Database files are downloaded in parallel using two simultaneous workers by
default:

```bash
download_sma3s_db.py --db bacteria --workers 2
```

The number of simultaneous downloads can be modified using `--workers`.

Interrupted downloads are stored temporarily as `.part` files and can be
resumed by running the same command again.

## Additional options

The complete command-line help can be displayed with:

```bash
download_sma3s_db.py --help
```

Main options include:

```text
--db DATABASE       Reference database to download
--list-dbs          List databases and their current availability
-o, --output DIR    Output directory [default: db_sma3s]
-w, --workers INT   Number of simultaneous downloads [default: 2]
--overwrite         Replace previously downloaded files
--retries INT       Maximum number of download attempts per file
--timeout INT       Connection timeout in seconds
```

---

# Hit filtering and taxonomy

## Recommended name

The recommended replacement for
`filter_perfect_hits_uniprot_sqlite.py` is:

```text
Source file: sma3s_hit_filter.py
Conda command: sma3s-hit-filter
```

This name is preferable to `perfect_hits` because the identity and coverage
thresholds are configurable. The default is a perfect 100/100/100 match, but
the program can also select less restrictive matches.

Other possible names are:

| Name | Comment |
|---|---|
| `sma3s-exact-hits` | Clear for the default strict use, but less accurate when thresholds are reduced. |
| `sma3s-hit-taxonomy` | Emphasizes species/TaxID reporting, although taxonomy is optional. |
| `sma3s-hit-report` | Broad and future-proof, but less explicit about filtering. |

## Subcommands

```bash
sma3s-hit-filter build-index [OPTIONS]
sma3s-hit-filter filter [OPTIONS]
```

## Build a UniProt species index

```bash
sma3s-hit-filter build-index \
    --uniprot-dat uniprot_sprot_trembl_all.dat.gz \
    --species-db uniprot_sprot_trembl_all.species.sqlite
```

The index stores:

- UniProt entry name;
- UniProt primary and secondary accessions;
- organism name from `OS`;
- NCBI TaxID from `OX`;
- identifier type (`ID` or `AC`).

The input can be `.dat` or `.dat.gz`.

### Build-index options

| Option | Default | Description |
|---|---:|---|
| `--uniprot-dat FILE` | required | UniProt `.dat` or `.dat.gz`. |
| `--species-db FILE` | required | Output SQLite database. |
| `--overwrite` | disabled | Replace an existing index. |
| `--batch-size INT` | `50000` | SQLite insertion batch size. |
| `--progress-every INT` | `100000` | Progress-report interval in UniProt records. |

## Filter annotation hits

The annotation table must contain:

```text
ANNOTATION_HIT_SEQUENCE
ANNOTATION_HIT_IDENTITY
ANNOTATION_HIT_QCOV
ANNOTATION_HIT_TCOV
```

Generate these columns with:

```bash
sma3s ... --report-hit-metrics
```

Then filter:

```bash
sma3s-hit-filter filter \
    -i annotation_hitmetrics.tsv \
    -o perfect_uniprot_hits
```

### Add species and TaxID

```bash
sma3s-hit-filter filter \
    -i annotation_hitmetrics.tsv \
    -o perfect_uniprot_hits \
    --species-db uniprot.species.sqlite
```

### Use relaxed thresholds

```bash
sma3s-hit-filter filter \
    -i annotation_hitmetrics.tsv \
    -o hits_90_90_90 \
    --min-identity 90 \
    --min-qcov 90 \
    --min-tcov 90
```

### Filter options

| Option | Default | Description |
|---|---:|---|
| `-i`, `--input FILE` | required | Sma3s annotation table; `.gz` is accepted. |
| `-o`, `--output-prefix PREFIX` | `perfect_uniprot_hits` | Prefix for the three output files. |
| `--species-db FILE` | unset | Optional index generated by `build-index`. |
| `--sep {auto,tab,comma}` | `auto` | Input separator. |
| `--min-identity FLOAT` | `100` | Minimum hit identity. |
| `--min-qcov FLOAT` | `100` | Minimum query coverage. |
| `--min-tcov FLOAT` | `100` | Minimum target coverage. |

## Hit-filter outputs

### Gene-level table

```text
<PREFIX>.perfect_genes.tsv
```

Contains the original annotation row plus:

```text
N_PERFECT_HITS
PERFECT_HIT_IDS
PERFECT_HIT_SPECIES
PERFECT_HIT_TAXIDS
```

The output filename retains `perfect` even when relaxed thresholds are used.
A future revision could rename this suffix to `.selected_genes.tsv`.

### Long hit table

```text
<PREFIX>.perfect_hits.long.tsv
```

One row is written per selected query-hit pair:

```text
ID
GENENAME
DESCRIPTION
HIT_ID
HIT_IDENTITY
HIT_QCOV
HIT_TCOV
HIT_ALNLEN
HIT_EVALUE
HIT_SPECIES
HIT_TAXID
```

### Summary

```text
<PREFIX>.summary.tsv
```

Reports:

- total input genes;
- genes with at least one selected hit;
- total selected hits;
- selected hits without a species assignment;
- rows with mismatched metric-field lengths;
- identity and coverage thresholds.

---

# Functional enrichment

## Command

```bash
sma3s-enrichment \
    -a ANNOTATION_TSV \
    -i QUERY_IDS
```

The annotation table defines the GO and Keyword universe. The query ID file
contains the subset to test.

## Important GO-column compatibility

By default, `sma3s-enrichment` expects:

```text
#ID
GENENAME
GO
KEYWORD
```

The simplest compatible annotation run is therefore:

```bash
sma3s ...              # without -go
```

If the Sma3s annotation was generated with `-go`, specify one of the separated
columns explicitly, for example:

```bash
sma3s-enrichment \
    -a annotation_go.tsv \
    -i genes.txt \
    --go-column 'GO(P)ID'
```

That run tests only the selected GO column. For a combined enrichment followed
by BP/MF/CC splitting, retain the single `GO` column and provide
`--go-obo go-basic.obo`.

## Basic enrichment

```bash
sma3s-enrichment \
    -a annotation.tsv \
    -i genes_of_interest.txt
```

Default output directory:

```text
enrichment/
```

Default prefix:

```text
functional_enrichment
```

## GO names and namespaces

```bash
sma3s-enrichment \
    -a annotation.tsv \
    -i genes_of_interest.txt \
    --go-obo go-basic.obo
```

This enables:

- GO names;
- GO namespaces;
- BP, MF and CC tables;
- BP, MF and CC plots;
- namespace-specific FDR values.

## Keyword treemap

```bash
sma3s-enrichment \
    -a annotation.tsv \
    -i genes_of_interest.txt \
    --make-keyword-treemap
```

## Taxonomic summaries

Use the long table generated by `sma3s-hit-filter`:

```bash
sma3s-enrichment \
    -a annotation.tsv \
    -i genes_of_interest.txt \
    --taxa-hits perfect_uniprot_hits.perfect_hits.long.tsv
```

This produces:

- genus abundance among query genes;
- number of hit rows per genus;
- genus diversity per query gene;
- genes with the highest number of distinct genus hits.

Genus is inferred from the first meaningful token of `HIT_SPECIES`.
`Candidatus <Genus>` names receive special handling.

These outputs are abundance and diversity summaries, not genus enrichment.

## Statistical method

For each GO or Keyword term, the program performs a one-sided hypergeometric
over-representation test:

```text
Universe = all valid IDs in the annotation table
Study set = supplied query IDs found in the universe
Term set = universe IDs annotated with the term
```

P values are adjusted using the Benjamini-Hochberg procedure.

When GO namespaces are available, the output includes:

- global FDR across all tested GO terms;
- namespace-specific FDR calculated independently within BP, MF or CC.

## Enrichment options

### Required and output options

| Option | Default | Description |
|---|---:|---|
| `-a`, `--annotations FILE` | required | Annotation universe. |
| `-i`, `--ids FILE` | required | Query IDs, one per line; first column is used. |
| `-o`, `--out-prefix PREFIX` | `functional_enrichment` | Output prefix inside the output directory. |
| `-d`, `--outdir DIR` | `enrichment` | Output directory. |

### Column options

| Option | Default | Description |
|---|---:|---|
| `--id-column NAME` | `#ID` | Annotation ID column. |
| `--genename-column NAME` | `GENENAME` | Gene label column. |
| `--go-column NAME` | `GO` | GO term column. |
| `--keyword-column NAME` | `KEYWORD` | Keyword column. |
| `--sep STRING` | tab | Annotation and taxonomy table separator. |
| `--term-delimiter-regex REGEX` | `;` | GO/Keyword term separator. |

### GO metadata

| Option | Default | Description |
|---|---:|---|
| `--go-obo FILE` | unset | `go-basic.obo` file providing names and namespaces. |
| `--go-names-tsv FILE` | unset | Alternative GO ID/name mapping. |
| `--no-split-go-namespaces` | disabled | Do not create separate namespace outputs. |
| `--use-namespace-fdr-for-namespace-plots` | disabled | Use BP/MF/CC-specific FDR in namespace plots. |

### Taxonomy summaries

| Option | Default | Description |
|---|---:|---|
| `--taxa-hits FILE` | unset | Long hit/taxonomy table. |
| `--taxa-id-column NAME` | `ID` | Query ID column in the taxonomy table. |
| `--taxa-species-column NAME` | `HIT_SPECIES` | Species column. |
| `--write-species-table` | disabled | Also calculate a species-level enrichment table. |

### Filtering and plots

| Option | Default | Description |
|---|---:|---|
| `--min-term-size INT` | `1` | Ignore terms assigned to fewer universe IDs. |
| `--max-term-size INT` | unlimited | Ignore terms assigned to more universe IDs. |
| `--alpha FLOAT` | `0.05` | FDR cutoff used to prioritize plotted terms. |
| `--top INT` | `20` | Maximum terms or taxa shown per plot. |
| `--plot-formats LIST` | `png,pdf` | Comma-separated formats. |
| `--make-keyword-treemap` | disabled | Generate a Keyword treemap. |

## Enrichment outputs

With prefix `functional_enrichment`, the command always writes:

```text
functional_enrichment.GO.enrichment.tsv
functional_enrichment.KEYWORD.enrichment.tsv
functional_enrichment.GO.dotplot.png
functional_enrichment.GO.dotplot.pdf
functional_enrichment.GO.barplot.png
functional_enrichment.GO.barplot.pdf
functional_enrichment.KEYWORD.dotplot.png
functional_enrichment.KEYWORD.dotplot.pdf
functional_enrichment.KEYWORD.barplot.png
functional_enrichment.KEYWORD.barplot.pdf
functional_enrichment.summary.txt
```

With `--go-obo`, namespace-specific tables and plots are also written:

```text
functional_enrichment.GO.biological_process.enrichment.tsv
functional_enrichment.GO.molecular_function.enrichment.tsv
functional_enrichment.GO.cellular_component.enrichment.tsv
```

With `--make-keyword-treemap`:

```text
functional_enrichment.KEYWORD.treemap.<FORMAT>
```

With `--taxa-hits`:

```text
functional_enrichment.GENUS.abundance.tsv
functional_enrichment.GENE.genus_richness.tsv
functional_enrichment.GENUS.abundance.<FORMAT>
functional_enrichment.GENE.genus_richness.<FORMAT>
```

With `--write-species-table`:

```text
functional_enrichment.SPECIES.enrichment.tsv
```

---

## Input files

### Query FASTA

Protein FASTA:

```text
>gene_001
MKKIGYSAPRQT...
>gene_002
MNNNKDLTQLAE...
```

Nucleotide FASTA requires `-nucl`.

### UniProt text database

Accepted forms:

```text
reference.dat
reference.dat.gz
```

### Prepared reference pair

```text
reference.fasta
reference.annot
```

Reference FASTA identifiers are read up to the first whitespace and must match
the first field of the `.annot` file.

### Annotation table for hit filtering

Required columns:

```text
ANNOTATION_HIT_SEQUENCE
ANNOTATION_HIT_IDENTITY
ANNOTATION_HIT_QCOV
ANNOTATION_HIT_TCOV
```

Optional but propagated to the long table:

```text
#ID
GENENAME
DESCRIPTION
ANNOTATION_HIT_ALNLEN
ANNOTATION_HIT_EVALUE
```

### Query ID list for enrichment

One identifier per line is recommended:

```text
gene_001
gene_002
gene_003
```

When a line contains several columns, only the first token is used.

---

## Output files

## Main annotation table

Default columns:

```text
#ID
GENENAME
DESCRIPTION
ENZYME
GO
KEYWORD
PATHWAY
```

`-go` replaces `GO` with:

```text
GO(P)ID
GO(P)NAME
GO(F)ID
GO(F)NAME
GO(C)ID
GO(C)NAME
```

`-source` adds:

```text
ANNOTATOR
USED_UNIPROT_SEQUENCES
```

`--report-hit-metrics` adds:

```text
ANNOTATION_HIT_SEQUENCE
ANNOTATION_HIT_IDENTITY
ANNOTATION_HIT_QCOV
ANNOTATION_HIT_TCOV
ANNOTATION_HIT_ALNLEN
ANNOTATION_HIT_EVALUE
```

Multiple supporting hits are represented as positionally matched,
semicolon-separated values.

## Main summary table

The `_summary.tsv` file reports:

- total query sequences;
- annotated query sequences;
- counts for each annotation type;
- counts for A1, A2, A3 and combined branches;
- GO Slim counts when requested;
- pathway frequencies;
- Keyword-category frequencies.

## MMseqs2 table

The `.mmseqs.m8` file contains:

```text
query
target
pident
qcov
alnlen
evalue
tcov
qlen
tlen
```

---

## Reference databases

### Prepared database

The database downloader provides the fastest startup because the FASTA,
annotation file and main SQLite caches are already available.

Use:

```bash
sma3s \
    -i proteins.faa \
    -d db_sma3s/bacteria/uniprot_bacteria.fasta
```

### Building from `.dat.gz`

```bash
sma3s \
    -i proteins.faa \
    -d uniprot_sprot.dat.gz \
    -num_threads 16
```

Sma3s:

1. creates or reuses `uniprot_sprot.dat`;
2. validates that the decompressed file resembles UniProt text format;
3. creates or reuses `uniprot_sprot.fasta`;
4. creates or reuses `uniprot_sprot.annot`;
5. creates the necessary SQLite caches.

Decompression is written to a temporary file and finalized atomically.

### Taxonomically filtered references

Filtered references use deterministic suffixes and are retained for reuse.

Examples:

```text
reference.Salmonella.fasta
reference.without_Vibrio.fasta
reference.order_Enterobacterales.family_Enterobacteriaceae.fasta
```

Inclusion filters are combined with logical AND. Matching any exclusion filter
removes an entry.

---

## Performance and parallelization

### Recommended starting point

```bash
sma3s \
    -i proteins.faa \
    -d reference.fasta \
    -num_threads 16 \
    --annotation-workers 4 \
    --tmpdir /local/scratch/$USER/sma3s
```

### CPU allocation

- `-num_threads` controls MMseqs2 threads.
- `--annotation-workers` controls Python processes.
- `--cluster-threads` controls threads per A3 clustering task.
- `--decompression-threads` controls `rapidgzip`.

Avoid assigning all CPUs independently to every setting. The total concurrent
load can otherwise exceed the allocated resources.

### Memory

Use:

```bash
--split-memory-limit 80G
```

to reduce MMseqs2 out-of-memory risk. Lower limits can increase runtime because
the search is split into more passes.

### Storage

Use fast local scratch for:

```bash
--tmpdir /scratch/$USER/sma3s
```

This is especially important for large MMseqs2 searches and SQLite-intensive
analyses.

### Hit count

Reducing:

```bash
--max-seqs 100
```

reduces disk usage and A2/A3 workload, but also reduces the evidence available
to the multi-hit annotator.

### Reciprocal FASTA mode

- `index`: random access through FASTA offsets;
- `stream`: one sequential scan of the reference FASTA;
- `auto`: select `stream` for large candidate sets.

Sequential mode is often faster on network filesystems.

---

## Caches and restart behavior

Sma3s can create:

```text
reference.fasta.sma3s_fasta_index.sqlite
reference.annot.<flags>.sma3s_annot.sqlite
<prefix>.mmseqs.m8.hits.sqlite
<prefix>.mmseqs.m8.reciprocal.sqlite
```

Reuse rules:

- a non-empty `.mmseqs.m8` is reused unless `--force-search` is supplied;
- FASTA indexes are rebuilt when file size or modification time changes;
- annotation caches are specific to annotation-cleaning options;
- derived FASTA and `.annot` files are reused when both are present;
- decompressed `.dat` files are reused when current relative to `.dat.gz`.

### Clean run outputs

```bash
sma3s \
    -i proteins.faa \
    -d reference.dat.gz \
    --force-clean-outputs \
    --clean-only
```

### Rebuild derived database files

```bash
sma3s \
    -i proteins.faa \
    -d reference.dat.gz \
    --force-clean-db \
    --clean-only
```

For `.dat.gz`, this removes the decompressed `.dat`, derived FASTA, `.annot`
and their caches, but never the original compressed archive.

### Remove both

```bash
sma3s \
    -i proteins.faa \
    -d reference.dat.gz \
    --force-clean \
    --clean-only
```

---

## Troubleshooting

### Check the complete installation

```bash
sma3s --check-install
```

Include this output in bug reports.

### MMseqs2 is not found

Activate the Conda environment:

```bash
conda activate sma3s
```

Or provide the full executable path:

```bash
sma3s --check-install --mmseqs-bin /path/to/mmseqs
```

### Parallel gzip is unavailable

The annotation still works through standard gzip, but decompression is
single-threaded. Ensure `rapidgzip` is included in the Conda environment.

### The prepared annotation cache is ignored

Verify that its filename contains:

```text
q0_go0_goslim0
```

and not:

```text
q0_go0_qoslim0
```

### A FASTA reference reports a missing `.annot`

For:

```text
reference.fasta
```

the following file must exist:

```text
reference.annot
```

### The hit filter reports missing columns

Re-run the annotation with:

```bash
--report-hit-metrics
```

The filter cannot infer identity and coverage from the standard annotation
columns.

### Species are reported as `NA`

Create and provide the UniProt species index:

```bash
sma3s-hit-filter build-index \
    --uniprot-dat reference.dat.gz \
    --species-db reference.species.sqlite
```

Then use:

```bash
sma3s-hit-filter filter \
    -i annotation_hitmetrics.tsv \
    --species-db reference.species.sqlite
```

### GO enrichment reports a missing `GO` column

The annotation was probably generated with `-go`, which creates separate GO
namespace columns.

Either:

1. annotate without `-go`; or
2. select a specific column with `--go-column`.

### Query IDs are missing from the enrichment universe

The IDs in the query list must exactly match the annotation table's ID column.
By default this is `#ID`.

### Treemap output is unavailable

Ensure `squarify` is installed in the environment.

### MMseqs2 runs out of memory

Try:

```bash
--split-memory-limit 80G
--max-seqs 100
--tmpdir /local/scratch/$USER/sma3s
```

### Existing results are unexpectedly reused

Use:

```bash
--force-search
```

or clean the appropriate output/cache category before rerunning.

### Results differ from Sma3s v2

Sma3s v3 preserves the annotation decision logic, but MMseqs2 and BLAST use
different search heuristics. Borderline matches may therefore differ.

---

## Citation

A dedicated Sma3s v3 citation should be added here when the corresponding
publication is available.

For the original Sma3s method, cite:

> Muñoz-Mérida A., Viguera E., Claros M. G., Trelles O., Pérez-Pulido A. J.
> (2014). Sma3s: a three-step modular annotator for large sequence datasets.
> *DNA Research*, 21(4), 341–353.
> https://doi.org/10.1093/dnares/dsu001

For the extended Sma3s implementation and applications, cite:

> Casimiro-Soriguer C. S., Muñoz-Mérida A., Pérez-Pulido A. J.
> (2017). Sma3s: A universal tool for easy functional annotation of proteomes
> and transcriptomes. *Proteomics*, 17, 1700071.
> https://doi.org/10.1002/pmic.201700071

Sma3s v3 uses MMseqs2:

> Steinegger M., Söding J. (2017). MMseqs2 enables sensitive protein sequence
> searching for the analysis of massive data sets.
> *Nature Biotechnology*, 35, 1026–1028.
> https://doi.org/10.1038/nbt.3988

Users should also cite the UniProt or UniRef release used as the reference
database and the Gene Ontology release used for enrichment analyses.

---

## License

Sma3s is distributed under the GNU General Public License v3.0. See
[`LICENSE`](LICENSE).

---

## Issues

Bug reports and feature requests should be submitted through the GitHub issue
tracker.

Please include:

- the complete command;
- the output of `sma3s --check-install`;
- the terminal output or traceback;
- the operating system and Conda environment;
- the reference type (`.dat.gz`, `.dat`, or FASTA/`.annot`);
- whether files are stored locally or on a network filesystem;
- the smallest shareable dataset reproducing the issue.

Do not upload confidential sequence data to a public issue.
