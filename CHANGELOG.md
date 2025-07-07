# Changelog

## Version 0.1.5

[@cumbof/metasbt@0.1.4.post1](https://github.com/cumbof/MetaSBT/releases/tag/0.1.5)

### New features

- `unpack` can automatically rename an unpacked database with the specified `--database` input argument;
- `update` exposes two new arguments `--uncertainty` and `--pruning-threshold` to tune the profiling performances.

### Fixes

- `db` correctly downloads the selected database version now;
- `unpack` is now trimming the whole database structure out up to the database folder so that `unpack` would eventually work as expected;
- `unpack` automatically fixes the paths to the bloom filter sketches onced a database is unpacked to a new location, usually different from the one where the database was located at the time of packing it;
- `update` correctly generates a new database also in case of no new unknown genomes;

## Version 0.1.4.post1

[@cumbof/metasbt@0.1.4.post1](https://github.com/cumbof/MetaSBT/releases/tag/0.1.4.post1)

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