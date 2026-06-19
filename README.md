# MetaSBT

![PyPI](https://img.shields.io/pypi/v/metasbt)
![Python](https://img.shields.io/pypi/pyversions/metasbt)
![Conda](https://img.shields.io/conda/dn/bioconda/metasbt?label=metasbt%20in%20Conda)
[![DOI](https://img.shields.io/badge/bioRxiv-10.1101/2025.08.25.672238-blue.svg)](https://doi.org/10.1101/2025.08.25.672238)

MetaSBT is a scalable framework for indexing microbial genomes and characterizing
metagenome-assembled genomes (MAGs) using **Delta-SBT**, a taxonomy-aware Sequence Bloom Tree
architecture backed by a native Rust engine.

## Key Features

- **Delta-SBT Architecture** — Replaces traditional Bloom filters with
  [FracMinHash](https://doi.org/10.1101/2023.01.12.523840) sketches compressed as
  [Roaring Bitmaps](https://roaringbitmap.org/). Tree nodes store strictly disjoint
  "Delta" k-mers, reducing disk footprint by ~95% compared to classic SBTs.

- **Rust Backend** — Core algorithms (sketching, tree building, accumulator search)
  are implemented in Rust via [PyO3](https://pyo3.rs/) for maximum performance.

- **Dual-Payload Sketches** — Each genome `.bf` file contains both a DNA FracMinHash
  sketch and a 3-frame amino acid translation sketch. The tree is built in two sweeps:
  DNA mode for species/genus/family, AA mode for order/class/phylum/kingdom, enabling
  accurate profiling across all taxonomic ranks.

- **Accumulator Search** — A BFS traversal that exploits the distributive property
  of disjoint Delta filters (`|Q ∩ (Core ∪ Delta)| = |Q ∩ Core| + |Q ∩ Delta|`),
  avoiding costly full-filter decompression.

## Installation

```bash
pip install metasbt
```

Or via Conda:

```bash
conda install -c bioconda metasbt
```

The Rust extension is built automatically during installation via `setuptools-rust`.
No separate Rust toolchain is required for end users.

## Quick Start

```bash
# Index a set of reference genomes
metasbt index --genomes references.tsv

# Update the database (build trees, compute profiles)
metasbt update

# Characterize MAGs against the indexed database
metasbt characterize --genomes mags.txt

# Query a genome against the database
metasbt query --genome genome.fna

# Search for similar genomes
metasbt search --genome genome.fna
```

Input files are TSV/CSV with genome identifiers and paths:

```
genome_id	/path/to/genome.fna
```

## Public Databases

We maintain a collection of pre-built MetaSBT databases for common microbial reference
sets. Browse and download them from [MetaSBT-DBs](https://github.com/cumbof/MetaSBT-DBs).

```bash
metasbt db --list
metasbt db --download <database_name>
```

## Documentation

Full documentation, including configuration options, database creation guides, and
tutorials, is available on the [official wiki](https://github.com/cumbof/MetaSBT/wiki).

## Architecture Overview

```
┌──────────────────────────────────────────────────────────────────┐
│                        Delta-SBT Tree                            │
│                                                                  │
│  Kingdom ─── AA Core (Intersection of phylum AA sketches)       │
│    │                                                             │
│  Phylum ──── AA Core                                             │
│    │                                                             │
│  Class ───── AA Core                                             │
│    │                                                             │
│  Order ───── AA Core                                             │
│    │                                                             │
│  Family ──── DNA Core (Intersection of genus DNA sketches)       │
│    │                                                             │
│  Genus ───── DNA Core                                            │
│    │                                                             │
│  Species ─── DNA Core                                            │
│    │                                                             │
│  Genomes ─── Delta Sketches (disjoint from parent Core)          │
└──────────────────────────────────────────────────────────────────┘
```

Each node contains both DNA and AA bitmaps in a single file. The DNA slot is built
via Core intersection (species → kingdom) in sweep 1; the AA slot is populated via
union propagation and then Core-intersected for order+ levels in sweep 2.

## Citing MetaSBT

```bibtex
@article{Cumbo2025.08.25.672238,
	author    = {Cumbo, Fabio and Blankenberg, Daniel},
	title     = {Characterization of microbial dark matter at scale with MetaSBT and taxonomy-aware Sequence Bloom Trees},
    journal   = {bioRxiv},
    year      = {2025},
    publisher = {Cold Spring Harbor Laboratory},
	doi       = {10.1101/2025.08.25.672238}
}
```

## Contributing

Bug reports and feature requests are managed via GitHub
[Issues](https://github.com/cumbof/MetaSBT/issues). Code contributions are welcome
via [Pull Requests](https://github.com/cumbof/MetaSBT/pulls).

Before contributing:
1. Check for existing issues or PRs related to your topic;
2. Open a discussion first for significant changes;
3. Include tests and documentation with your PR.

## Support

Open an [Issue](https://github.com/cumbof/MetaSBT/issues) or start a
[Discussion](https://github.com/cumbof/MetaSBT/discussions).

Copyright © 2025 [Fabio Cumbo](https://github.com/cumbof),
[Daniel Blankenberg](https://github.com/blankenberg).
See [LICENSE](https://github.com/cumbof/MetaSBT/blob/main/LICENSE).
