# Changelog

## Version 0.1.4

[@cumbof/metasbt@0.1.4](https://github.com/cumbof/MetaSBT/releases/tag/0.1.4)

New object-oriented implementation of MetaSBT.
Clusters are consistent with the definition of Average Nucleotide Identity (ANI).
Clusters' boundaries are defined as the minimum and maximum ANI distance between all the genomes under a specific cluster.

### New features

It provides the following subroutines:
- `db`: List and retrieve public MetaSBT databases;
- `index`: Index a set of reference genomes and build the first baseline of a MetaSBT database;
- `kraken`: Export a MetaSBT database into a custom kraken database;
- `pack`: Build a compressed tarball with a MetaSBT database and report its sha256;
- `profile`: Profile an input genome and report the closest cluster at all the seven taxonomic levels and the closest genome in a MetaSBT database;
- `sketch`: Sketch the input genomes;
- `summarize`: Summarize the content of a MetaSBT database and report some statistics;
- `test`: Check for dependencies and run unit tests. This must be used by code maintainers only;
- `unpack`: Unpack a local MetaSBT tarball database;
- `update`: Update a MetaSBT database with new metagenome-assembled genomes.

The MetaSBT [core](https://github.com/cumbof/MetaSBT/blob/main/metasbt/core.py) provides an interface to the `Database` and `Entry` class abstractions.

### Fixes

None