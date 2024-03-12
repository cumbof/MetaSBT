# Changelog

## Version 0.1.3

[@cumbof/metasbt@0.1.3](https://github.com/cumbof/MetaSBT/releases/tag/0.1.3)

### New features

- New option `--use-representatives` available with the `index` module to use only three representative genomes at the species level;
- New option `--resume` available with the `index` and `update` modules able to resume the index and update processes in case of unexpected errors;
- Both the `index` and `update` modules now display a worning message in case the configuration file under `--resume` has been previously generated with a different version of MetaSBT;
- Both the `index` and `update` modules now integrate `CheckV` and `EukCC` for assessing the quality of viruses and eukaryotes;
- `CheckM` has been updated to `CheckM2`.

### Fixes

- It correctly checks now for new framework versions when starting a new `metasbt` instance;
- Fixed genome quality filtering on completeness and contamination during the `update`;
- Improving docstring adopting the [numpydoc](https://numpydoc.readthedocs.io/en/latest/) documentation format.

## Version 0.1.2

[@cumbof/metasbt@0.1.2](https://github.com/cumbof/MetaSBT/releases/tag/0.1.2)

First public stable release of MetaSBT.

It is composed of the following modules:

- `index`: build a MetaSBT database by building a series of Sequence Bloom Trees at different taxonomic levels;
- `boundaries`: define taxonomy-specific boundaries as the minimum and maximum number of kmers in common between all the genomes under a specific cluster;
- `profile`: taxonomically profile a genome by querying a MetaSBT database at different taxonomic levels;
- `report`: build a report table describing the content of a MetaSBT database;
- `update`: update a MetaSBT database with new genomes;
- `tar`: pack a MetaSBT database into a ready-to-be-distributed tarball;
- `install`: install a MetaSBT database tarball locally under a specific location of the file system.

The framework also comes with a set of utilities:

- `bf_sketch.py`: build minimal bloom filter sketches with cluster-specific marker kmers;
- `esearch_txid.sh`: retrieve GCAs from NCBI GenBank given a specific taxonomic ID;
- `get_ncbi_genomes.py`: retrieve reference genomes and metagenome-assembled genomes under a specific superkingdom and kingdom from NCBI GenBank;
- `howdesbt_index.sh`: index genomes with HowDeSBT;
- `uniform_inputs.sh`: uniform input genome files extension.
