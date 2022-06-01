# meta-index
A pipeline for automatically indexing genomes and accurately characterizing metagenome-assembled genomes with sequence bloom trees

## Contents

- [Installing `meta-index`](https://github.com/BlankenbergLab/meta-index#installing-meta-index)
- [Building a database](https://github.com/BlankenbergLab/meta-index#building-a-database)
- [Defining boundaries](https://github.com/BlankenbergLab/meta-index#defining-boundaries)
- [Updating the database](https://github.com/BlankenbergLab/meta-index#updating-the-database)
- [Profiling genomes](https://github.com/BlankenbergLab/meta-index#profiling-genomes)
- [Unlocking unknown species profiling with `kraken`](https://github.com/BlankenbergLab/meta-index#unlocking-unknown-species-profiling-with-kraken)
- [Credits](https://github.com/BlankenbergLab/meta-index#credits)
- [Contributing](https://github.com/BlankenbergLab/meta-index#support)

## Installing `meta-index`

The pipeline is available as a `conda` package through the following command:
```bash
conda install -c bioconda meta-index
```

You may need to add the `bioconda` channel first by running:
```bash
conda config --add channels bioconda
```

The `meta-index` pipeline is also available by simply cloning this repository and making all the scripts executables:
```bash
# Clone the meta-index repository
mkdir -p ~/git && cd ~/git
git clone https://github.com/BlankenbergLab/meta-index

# Make the scripts executable
chmod +x meta-index
cd modules && chmod +x *

# Add meta-index to the PATH env variable
PATH=$PATH:~/git/meta-index
```

In this last case, remember to check that the following dependencies are installed and available on your system:
- bc
- [checkm](https://github.com/Ecogenomics/CheckM) (version >=1.1.3)
- gzip
- [howdesbt](https://github.com/medvedevgroup/HowDeSBT) (version >=2.00.02)
- [kmtricks](https://github.com/tlemane/kmtricks) (version >=1.2.1)
- [ncbitax2lin](https://github.com/zyxue/ncbitax2lin)
- [ntcard](https://github.com/bcgsc/ntCard)
- python (version >=3.7)
- pip
- wget

You can check whether all these dependencies are available by running the following command in your terminal:
```bash
meta-index --resolve-dependencies
```

Please note that `meta-index` makes use of some advanced `howdesbt` commands that are not available by default when installing HowDeSBT. They must be enabled by compiling the software with the alternative version of the [Makefile](https://github.com/medvedevgroup/HowDeSBT/blob/master/Makefile_full) available in the root folder of the HowDeSBT repository on GitHub.

For what concerns CheckM, we strongly suggest to install it through `pip` or `conda`, but it will require in any case a couple of extra steps to correctly link the software to its database. This must be necessarily executed manually as reported on the official [CheckM Wiki](https://github.com/Ecogenomics/CheckM/wiki/Installation).

First, you need to download the last available database from the following repository [https://data.ace.uq.edu.au/public/CheckM_databases/](https://data.ace.uq.edu.au/public/CheckM_databases/), decompress it on a dedicated location, and finally inform CheckM about where its database is located by typing:
```bash
checkm data setRoot <checkm_data_dir>
```

Once everything is installed, `meta-index` will be available on your environment. You can check whether it has been correctly installed by typing the following command in your terminal:
```bash
meta-index --version
```

We strongly suggest to permanently add the `meta-index` folder to the PATH environment variable by adding the following line to your `~/.profile` or `~/.bash_profile` (if `bash` is your default shell):
```bash
echo "PATH=$PATH:~/git/meta-index" >> ~/.bash_profile
```

You may finally need to reload your profile to make these changes effective:
```bash
source ~/.bash_profile
```

### Warning (for macOS users)

Unfortunately, the CheckM conda package is still not fully compatible with macOS because of its software dependency `pplacer`. For this reason, macOS users are strongly encouraged to build a Docker container by running the following command from the `meta-index` root directory in which the Dockerfile is located:
```bash
docker build . -t meta-index
```

Once the container is built, you can finally open an interactive shell on the Ubuntu based container already configured to run `meta-index`:
```bash
docker run -it meta-index
```

In case you are able to install all the `meta-index` dependencies without the help of Docker, please remember that the `meta-index` modules makes use of the `date` command in order to take track of the amount of time the whole pipeline requires to run. Unfortunately, on macOS systems, the `date` command uses a different syntax to produce a timestamp accurate to the nanoseconds compared to the GNU equivalent command. For this reason, in order to make the script fully compatible also with macOS systems, users should follow the next additional steps.

Use homebrew to install `coreutils` from your terminal:
```bash
brew install coreutils
```

The GNU equivalent date will be named `gdate`. At this point, we strongly encourage users to not replace the built-in version of `date` with the GNU equivalent in order to avoid issues with other software components that rely on the original version of `date`.

We instead suggest to add an alias with the following command:
```bash
alias date='gdate'
```

In order to make the alias effective every time you open the shell, this last command should be added to your `.bash_profile` (in case `bash` is your default shell). 
```bash
echo "alias date='gdate'" >> ~/.bash_profile
source ~/.bash_profile
```

## Building a database

The `index` subroutine allows to automatically retrieve genomes from isolate sequencing from the NCBI GenBank and organise them in folders that reflect their taxonomic classification. It finally makes use of `kmtricks` to rapidly index all the genomes at the species level and create a sequence bloom tree for each of the species. Lower taxonomic levels are indexed with `howdesbt` by building new sequence bloom trees considering only the root nodes of the upper taxonomic levels.

The following command will trigger the generation of the database with all the available bacterial genomes from isolate sequencing in NCBI GenBank:
```bash
$ meta-index index --db-dir=~/myindex \
                   --kmer-len=31 \
                   --estimate-filter-size \
                   --kingdom=Bacteria \
                   --dereplicate \
                   --checkm-completeness=50.0 \
                   --checkm-contamination=5.0 \
                   --nproc=4 \
                   --xargs-nproc=2 \
                   --tmp-dir=~/tmp \
                   --cleanup
```

## Defining boundaries

The `boundaries` module is crucial for the definition of taxonomy-specific boundaries. It explots `kmtricks` to build a kmer table for each of the taxonomic levels in the database with information about the presence/absence of a kmer in genomes that belong to a particular lineage. It finally build a new table with the minimum and maximum amount of kmers in common between all the genomes in a particular taxonomic level. The average value between all the minimum and maximum common kmers is then used to establish whether a new MAG must be assigned to a new cluster with the `update` module.

The following command will trigger the definition of the kmer boundaries for each taxonomic level in the database:
```bash
$ meta-index boundaries --db-dir=~/myindex \
                        --kingdom=Bacteria \
                        --min-genomes=100 \
                        --output=~/boundaries.txt \
                        --tmp-dir=~/tmp \
                        --nproc=4 \
                        --cleanup
```

## Updating the database

This module can be used to add new reference genomes and metagenome-assembled genomes (MAGs) to the database generated with the `index` module. 

In case of new MAGs, it first try to profile them by comparing the input genomes with those present in the database. An input genome is assigned to the closest genome in the database if their set of kmers result similar enough. In case the profiler will not be able to characterize the input MAGs, the `update` subroutine will exploit `kmtricks` to build a kmer matrix for the set of input genomes and automatically group them according to the taxonomic level specific boundaries identified with the `boundaries` utility, and it finally generate new clusters of potentially novel and yet-to-be-named species.

In case of new reference genomes from isolate sequencing, the `update` module simply add the new genomes to the corresponding species by rebuilding the species trees and all the trees at the lower taxonomic levels (please note that the taxonomic labels of the input reference genomes are known). If the taxonomy of an input reference genome is not in the database, the module will compare the input genome with all the cluster of species with no references before generating a new branch in the index database. If the input genome results very close to a cluster of unknown genomes, its taxonomic label will be inherited by all the genomes in the unknown cluster.

The `update` module can be launched with the following command:
```bash
$ meta-index update --input-list=~/mygenomes.txt \
                    --taxa=~/taxonomies.tsv \
                    --db-dir=~/myindex \
                    --kingdom=Bacteria \
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

Please note that the `--filter-size` value must be the same used while running the `index` module.

Remember also that `meta-index` requires that all your input genomes have the same format and extension before running the pipeline. You can easily uniform your genome files extension by typing the following command in your terminal:
```bash
INPUTS_DIR=~/mygenomes
CURRENT_EXTENSION="fa"
NEW_EXTENSION="fna"
find ${INPUTS_DIR} \
    -type f -iname "*.${CURRENT_EXTENSION}" -follow | xargs -n 1 -i sh -c \
        'INPUT={}; \
         mv "$INPUT" "${INPUT%.'"${CURRENT_EXTENSION}"'}.'"${NEW_EXTENSION}"'";'
```

## Profiling genomes

The `profile` module allows to characterize an input genome according to the taxonomic label of the closest genome in the database. It allows to process only one genome in input at a time:
```bash
$ meta-index profile --input-file=~/mymag.fna \
                     --input-id=mymag \
                     --tree=~/myindex/k__Bacteria/index.detbrief.sbt \
                     --threshold=0.7 \
                     --expand \
                     --output-dir=~/profiles \
                     --output-prefix=mymag
```

## Unlocking unknown species profiling with `kraken`

TBA

## Credits

Please credit our work in your manuscript by citing:

TBA

## Contributing

Long-term discussion and bug reports are maintained via GitHub Issues, while code review is managed via GitHub Pull Requests.

Please, (i) be sure that there are no existing issues/PR concerning the same bug or improvement before opening a new issue/PR; (ii) write a clear and concise description of what the bug/PR is about; (iii) specifying the list of steps to reproduce the behavior in addition to versions and other technical details is highly recommended.

## Support

If you need support, please open an [Issue](https://github.com/BlankenbergLab/meta-index/issues) or a new [Discussion](https://github.com/BlankenbergLab/meta-index/discussions). We will be happy to answer your questions and help you troubleshoot any kind of issue concerning our framework.

Copyright © 2022 Fabio Cumbo, Daniel Blankenberg. See LICENSE for additional details.
