# meta-index
A pipeline for automatically indexing genomes and accurately characterizing metagenome-assembled genomes with sequence bloom trees

## Building a database

The index subroutine allows to automatically retrieve reference genomes and metagenome-assembled genomes from NCBI GenBank and organise them in folders that reflect their taxonomic classification. It finally makes use of kmtricks to rapidly index all the genomes at the species level and create a sequence bloom tree for each of the species. Lower taxonomic levels are indexed by building new sequence bloom trees considering only the root nodes of the upper taxonomic levels.

The following command will trigger the generation of the database with all the available bacterial genomes from isolate sequencing in NCBI GenBank:
```bash
$ bash meta-index index --work-dir=~/myindex \
                        --kmer-len=31 \
                        --filter-size=10000 \
                        --kingdom=Bacteria \
                        --dereplicate \
                        --nproc=4 \
                        --xargs-nproc=2
```

## Updating the database
TBA

## Profiling genomes
TBA

## Improving kraken database quality

TBA

## Credits

Please credit our work in your manuscript by citing:

TBA
