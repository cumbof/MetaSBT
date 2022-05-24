# meta-index
A pipeline for automatically indexing genomes and accurately characterizing metagenome-assembled genomes with sequence bloom trees

## Building a database

The `index` subroutine allows to automatically retrieve genomes from isolate sequencing from the NCBI GenBank and organise them in folders that reflect their taxonomic classification. It finally makes use of `kmtricks` to rapidly index all the genomes at the species level and create a sequence bloom tree for each of the species. Lower taxonomic levels are indexed with `howdesbt` by building new sequence bloom trees considering only the root nodes of the upper taxonomic levels.

The following command will trigger the generation of the database with all the available bacterial genomes from isolate sequencing in NCBI GenBank:
```bash
$ meta-index index --db-dir=~/myindex \
                   --kmer-len=31 \
                   --filter-size=10000 \
                   --kingdom=Bacteria \
                   --dereplicate \
                   --checkm-completeness=50.0 \
                   --checkm-contamination=5.0 \
                   --nproc=4 \
                   --xargs-nproc=2 \
                   --cleanup
```

## Updating the database

This module can be used to add new reference genomes and metagenome-assembled genomes (MAGs) to the database generated with the `index` module. 

In case of new MAGs, it first try to profile them by comparing the input genomes with those present in the database. An input genome is assigned to the closest genome in the database if their set of kmers result similar enough. In case the profiler will not be able to characterize the input MAGs, the `update` subroutine will exploit `kmtricks` to build a kmer matrix for the set of input genomes and automatically group them according to the taxonomic level specific boundaries identified with the `boundaries` utility, and it finally generate new clusters of potentially novel and yet-to-be-named species.

In case of new reference genomes from isolate sequencing, the `update` module simply add the new genomes to the corresponding species by rebuilding the species trees and all the trees at the lower taxonomic levels (please note that the taxonomic labels of the input reference genomes are known). If the taxonomy of an input reference genome is not in the database, the module will compare the input genome with all the cluster of species with no references before generating a new branch in the index database. If the input genome results very close to a cluster of unknown genomes, its taxonomic label will be inherited by all the genomes in the unknown cluster.

The `update` module can be launched with the following command:
```bash
$ meta-index update --input-list=~/mygenomes.txt \
                    --db-dir=~/myindex \
                    --tmp-dir=~/tmp \
                    --kmer-len=31 \
                    --filter-size=10000 \
                    --dereplicate \
                    --checkm-completeness=50.0 \
                    --checkm-contamination=5.0 \
                    --type=references \
                    --extension=fna.gz \
                    --nproc=8 \
                    --cleanup
```

## Profiling genomes

The `profile` module allows to characterize an input genome according to the taxonomic label of the closest genome in the database. It allows to process only one genome in input at a time:
```bash
$ meta-index profile --input=~/mygenome.fna \
                     --tree=~/myindex/tree.txt \
                     --threshold=0.7 \
                     --output-dir=~/myqueries \
                     --output-prefix=query \
                     --expand
```

## Improving kraken database quality

TBA

## Credits

Please credit our work in your manuscript by citing:

TBA
