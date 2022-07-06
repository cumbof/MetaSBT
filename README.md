# MetaSBT
A scalable framework for automatically indexing microbial genomes and accurately characterizing metagenome-assembled genomes with Sequence Bloom Trees

> :warning: _This is still under development and some features may produce errors or unexpected results. You may also want to look at the [Contributing](https://github.com/cumbof/MetaSBT#contributing) and [Support](https://github.com/cumbof/MetaSBT#support) sections in case you will encounter a bug_

## Contents

- [Getting started](https://github.com/cumbof/MetaSBT#getting-started)
  - [Installing `MetaSBT`](https://github.com/cumbof/MetaSBT#installing-metasbt)
  - [Dependencies](https://github.com/cumbof/MetaSBT#dependencies)
- [Available features](https://github.com/cumbof/MetaSBT#available-features)
  - [`index`: building a reference database](https://github.com/cumbof/MetaSBT#1-index-building-a-reference-database)
  - [`boundaries`: defining clusters boundaries](https://github.com/cumbof/MetaSBT#2-boundaries-defining-clusters-boundaries)
  - [`profile`: characterizing genomes and metagenome-assembled genomes](https://github.com/cumbof/MetaSBT#3-profile-characterizing-genomes-and-metagenome-assembled-genomes)
  - [`update`: updating the database](https://github.com/cumbof/MetaSBT#4-update-updating-the-database)
  - [`report`: building the database snapshot report](https://github.com/cumbof/MetaSBT#5-report-building-the-database-snapshot-report)
- [Credits](https://github.com/cumbof/MetaSBT#credits)
- [Contributing](https://github.com/cumbof/MetaSBT#contributing)
- [Support](https://github.com/cumbof/MetaSBT#support)

## Getting started

`MetaSBT` is a Python3 framework built for the analysis and characterization of the microbial dark matter. It is available on the [Python Package Index](https://pypi.org/) as well as Conda on the Bioconda[^1] channel and relies on a set of external software dependencies described on the following sections.

### Installing `MetaSBT`

The pipeline is available as a Python3 package that can be install with the following command:
```bash
pip install metasbt
```

It is also available as a `conda` package:
```bash
conda install -c bioconda metasbt
```

You may need to add the `bioconda` channel first by running:
```bash
conda config --add channels bioconda
```

The `MetaSBT` pipeline is also available by simply cloning this repository and making all the scripts executables:
```bash
# Clone the MetaSBT repository
mkdir -p ~/git && cd ~/git
git clone https://github.com/cumbof/MetaSBT.git

# Make the scripts executable
chmod -R +x MetaSBT/*.py

# Add MetaSBT to the PATH env variable
PATH=$PATH:~/git/MetaSBT
```

Once everything is installed, `MetaSBT` will be available on your environment. You can check whether it has been correctly installed by typing the following command in your terminal:
```bash
metasbt --version
```

You can check whether all the dependencies listed above are available on your system by running the following command:
```bash
metasbt --resolve-dependencies
```

You can also access the complete list of available arguments by specifying the `--help` option:
```bash
metasbt --help
```

Please note that the same option is also available for all the `MetaSBT` modules (e.g.: `metasbt profile --help` will print the list of arguments available for the `profile` module). The list of available modules is available by typing:
```bash
metasbt --modules
```

We strongly suggest to permanently add the `MetaSBT` folder to the PATH environment variable by adding the following line to your `~/.profile` or `~/.bash_profile` (if `bash` is your default shell):
```bash
echo "PATH=$PATH:~/git/MetaSBT" >> ~/.bash_profile
```

You may finally need to reload your profile to make these changes effective:
```bash
source ~/.bash_profile
```

### Dependencies

Please note that cloning this repository requires [Git](https://git-scm.com/) to be installed on your system.

In this last case, remember to check that the following dependencies are installed and available on your system:
- [checkm](https://github.com/Ecogenomics/CheckM) (version >=1.2.0)[^2]
- [howdesbt](https://github.com/medvedevgroup/HowDeSBT) (version >=2.00.07)[^3]
- [kmtricks](https://github.com/tlemane/kmtricks) (version >=1.2.1)[^4]
- [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) (version >=0.3.1)
- [ncbitax2lin](https://github.com/zyxue/ncbitax2lin) (version >=2.3.2)
- [ntcard](https://github.com/bcgsc/ntCard) (version >=1.2.2)[^5]
- [pip](https://pip.pypa.io/) (version >=21.2.4)
- [python](http://www.python.org/) (version >=3.7)
- [wget](https://www.gnu.org/software/wget/) (version >=1.21.3)

Please note that `MetaSBT` makes use of some advanced `howdesbt` sub-commands that are not available by default when installing HowDeSBT. They must be enabled by compiling the software with the alternative version of the [Makefile](https://github.com/medvedevgroup/HowDeSBT/blob/master/Makefile_full) available in the root folder of the HowDeSBT repository on GitHub.

For what concerns CheckM, we strongly suggest to install it through `pip` or `conda`, but it will require in any case a couple of extra steps to correctly link the software to its database. This must be necessarily executed manually as reported on the official [CheckM Wiki](https://github.com/Ecogenomics/CheckM/wiki/Installation).

First, you need to download the last available database from the following repository [https://data.ace.uq.edu.au/public/CheckM_databases/](https://data.ace.uq.edu.au/public/CheckM_databases/), decompress it on a dedicated location, and finally inform CheckM about where its database is located by typing:
```bash
checkm data setRoot <checkm_data_dir>
```

Please note that `MetaSBT` is available for Linux and macOS only.

## Available features

The `MetaSBT` framework provides a set of subroutines (here called modules) for (i) **building a taxonomically-organized database** of microbial genomes, (ii) **defining boundaries** for each of the clusters in the database at different taxonomic levels to explore and define the intra-cluster genetic diversity, (iii) **profiling genomes and metagenome-assembled genomes** by querying the database and reporting the closest lineage and reference genome, (iv) **updating a database with new reference genomes and metagenome-assembled genomes** by also defining new and yet-to-be-named clusters, (v) and finally **producing a report** of the current status of a database.

It follows a detailed explanation of all the available modules.

### 1. `index`: building a reference database

The `index` subroutine allows to automatically retrieve genomes from isolate sequencing from the NCBI GenBank and organise them in folders that reflect their taxonomic classification. It finally makes use of `kmtricks` to rapidly index all the genomes at the species level and create a sequence bloom tree for each of the species. Lower taxonomic levels are indexed with `howdesbt` by building new sequence bloom trees considering only the root nodes of the upper taxonomic levels.

The following command will trigger the generation of the database with all the available bacterial genomes from isolate sequencing in NCBI GenBank:
```bash
metasbt index --db-dir ~/myindex \
              --kmer-len 31 \
              --estimate-filter-size \
              --increase-filter-size 5.0 \
              --kingdom Bacteria \
              --dereplicate \
              --similarity 100.0 \
              --completeness 50.0 \
              --contamination 5.0 \
              --parallel 4 \
              --nproc 2 \
              --pplacer-threads 2 \
              --tmp-dir ~/tmp \
              --cleanup \
              --verbose
```

In case you would like to create a database from a specific set of genomes available on your file system, you may use the `--input-list` instead of the `--kingdom` option. It should point to a table file with the list of file paths to the input genomes on the first column and an optional second column with the full taxonomic labels of the input genomes in a tab-separated values format.

You may also want to get rid of the taxonomic organization of genomes by specifying the `--flat-structure` option that will consider all the input genomes together for the generation of a single sequence bloom tree.

#### Available options

| Option                   | Default | Mandatory | Description  |
|:-------------------------|:--------|:---------:|:-------------|
| `--completeness`         | `0.0`   |           | Minimum completeness percentage allowed for input genomes |
| `--contamination`        | `100.0` |           | Maximum contamination percentage allowed for input genomes |
| `--cleanup`              | `False` |           | Remove temporary data at the end of the pipeline |
| `--db-dir`               |         | ⚑         | Database folder path |
| `--dereplicate`          | `False` |           | Dereplicate input genomes |
| `--estimate-filter-size` | `False` |           | Estimate the bloom filter size with ntCard. It automatically override the `--filter-size` option |
| `--filter-size`          |         | ⚑         | Bloom filter size |
| `--flat-structure`       |         |           | Organize genomes without any taxonomic organization. This will lead to the creation of a single sequence bloom tree |
| `--help`                 |         |           | Print the list of arguments and exit |
| `--increase-filter-size` | `0.0`   |           | Increase the estimated filter size by the specified percentage. This is used in conjunction with the `--estimate_filter_size` argument only. It is highly recommended to increase the filter size by a good percentage in case you are planning to update the index with new genomes |
| `--input-list`           |         |           | Path to the input table with a list of genome file paths and an optional column with their taxonomic labels. Please note that the input genome files must be gz compressed with fna extension (i.e.: *.fna.gz) |
| `--kingdom`              |         |           | Consider genomes whose lineage belongs to a specific kingdom |
| `--kmer-len`             |         | ⚑         | This is the length of the kmers used for building bloom filters |
| `--limit-genomes`        | `Inf`   |           | Limit the number of genomes per species. This will remove the exceeding number of genomes randomly to cut the overall number of genomes per species to this number |
| `--log`                  |         |           | Path to the log file |
| `--max-genomes`          | `Inf`   |           | Consider species with this number of genomes at most |
| `--min-genomes`          | `1`     |           | Consider species with a minimum number of genomes |
| `--nproc`                | `1`     |           | This argument refers to the number of processors used for parallelizing the pipeline when possible |
| `--parallel`             | `1`     |           | Maximum number of processors to process each NCBI tax ID in parallel |
| `--pplacer-threads`      | `1`     |           | Maximum number of threads for pplacer. This is required to maximise the CheckM performances |
| `--similarity`           | `100.0` |           | Dereplicate genomes if they have a percentage of common kmers greater than or equals to the specified one. This is used exclusively in conjunction with the `--dereplicate` argument |
| `--tmp-dir`              |         | ⚑         | Path to the temporary folder |
| `--verbose`              | `False` |           | Print results on screen |
| `--version`              |         |           | Print current module version and exit |

#### Warning

Please note that the `--flat-structure` option is not compatible with the `update` module for updating the database with new genomes. Thus, in case you will need to update the sequence bloom tree, there are no other options than building the database from scratch with the new set of genomes.

### 2. `boundaries`: defining clusters boundaries

The `boundaries` module is crucial for the definition of taxonomy-specific boundaries. It explots `kmtricks` to build a kmers table for each of the taxonomic levels in the database with information about the presence/absence of a kmer in genomes that belong to a particular lineage. It finally build a new table with the minimum and maximum amount of kmers in common between all the genomes in a particular taxonomic level. These boundaries are then used by the `update` module in order to establish whether a new genome must be assigned to the closest cluster identified by the `profile` module.

The following command will trigger the definition of the kmer boundaries for each taxonomic level in the database:
```bash
metasbt boundaries --db-dir ~/myindex \
                   --min-genomes 50 \
                   --output ~/boundaries.tsv \
                   --tmp-dir ~/tmp \
                   --nproc 4 \
                   --cleanup \
                   --verbose
```

Please note that the `boundaries` module considers clusters with reference genomes only. These clusters can be considered for establishing boundaries depending on a minimum number of reference genomes that can be set with the `--min-genomes` argument.

#### Available options

| Option                   | Default | Mandatory | Description  |
|:-------------------------|:--------|:---------:|:-------------|
| `--cleanup`              | `False` |           | Remove temporary data at the end of the pipeline |
| `--db-dir`               |         | ⚑         | Database directory with the taxonomically organised sequence bloom trees |
| `--flat-structure`       | `False` |           | Genomes in the database have been organized without a taxonomic structure |
| `--help`                 |         |           | Print the list of arguments and exit |
| `--kingdom`              |         |           | Consider genomes whose lineage belongs to a specific kingdom |
| `--log`                  |         |           | Path to the log file |
| `--min-genomes`          | `3`     |           | Consider clusters with at least this number of genomes |
| `--nproc`                | `1`     |           | This argument refers to the number of processors used for parallelizing the pipeline when possible |
| `--output`               |         | ⚑         | Output file with kmer boundaries for each of the taxonomic labels in the database |
| `--tmp-dir`              |         | ⚑         | Path to the temporary folder |
| `--verbose`              | `False` |           | Print results on screen |
| `--version`              |         |           | Print current module version and exit |

### 3. `profile`: characterizing genomes and metagenome-assembled genomes

The `profile` module allows to characterize an input genome according to the closest lineage in the database. It allows to process only one genome in input at a time:
```bash
metasbt profile --input-file ~/mymag.fna \
                --input-id mymag \
                --input-type genome \
                --tree ~/myindex/k__Bacteria/index.detbrief.sbt \
                --threshold 0.7 \
                --expand \
                --stop-at family \
                --output-dir ~/profiles \
                --output-prefix mymag \
                --verbose
```

Please note that in the example above we explicitly set the `--stop-at` argument to `family`. This argument works in conjunction with the `--expand` option only, and it will prevent epanding the query to all the taxonomic levels lower than the specified one. Also note that the `--expand` argument expands the input query up to the species level by default, by also reporting the closest genome, without the need to use the `--stop-at` argument.

#### Available options

| Option                   | Default | Mandatory | Description  |
|:-------------------------|:--------|:---------:|:-------------|
| `--expand`               | `False` |           | Expand the input query on all the taxonomic levels |
| `--input-file`           |         | ⚑         | Path to the input query |
| `--input-id`             |         |           | Unique identifier of the input query |
| `--input-type`           |         | ⚑         | Accepted input types are genomes or files with one nucleotide sequence per line. In case of genomes, if they contain multiple sequences, results are merged together |
| `--log`                  |         |           | Path to the log file |
| `--output-dir`           |         | ⚑         | Output folders with queries results |
| `--output-prefix`        |         |           | Prefix of the output files with query matches |
| `--stop-at`              |         |           | Stop expanding queries at a specific taxonomic level. Please note that this argument works in conjunction with --expand only |
| `--threshold`            | `0.0`   |           | Fraction of query kmers that must be present in a leaf to be considered a match |
| `--tree`                 |         | ⚑         | This is the tree definition file |
| `--verbose`              | `False` |           | Print results on screen |
| `--version`              |         |           | Print current module version and exit |

### 4. `update`: updating the database

> :warning: _This module is still under development. We strongly suggest to avoid mentioning results produced with this specific module in scientific manuscripts until a stable version will be released_

This module can be used to add new reference genomes and metagenome-assembled genomes (MAGs) to the database generated with the `index` module. 

In case of new MAGs, it first try to profile them by comparing the input genomes with those present in the database. An input genome is assigned to the closest genome in the database if their set of kmers result similar enough. In case the profiler will not be able to characterize the input genome, the `update` subroutine will exploit `kmtricks` to build a kmer matrix for the set of input genomes and automatically group them according to the taxonomy-specific boundaries identified with the `boundaries` utility, and it finally generate new clusters of potentially novel and yet-to-be-named species.

In case of new reference genomes from isolate sequencing, the `update` module simply add the new genomes to the corresponding species by rebuilding the species trees and all the trees at the lower taxonomic levels (please note that the taxonomic labels of the input reference genomes are known). If the taxonomy of an input reference genome is not in the database, the module will compare the input genome with all the cluster of species with no references before generating a new branch in the index database. If the input genome results very close to a cluster of unknown genomes, its taxonomic label will be inherited by all the genomes in the unknown cluster.

The `update` module can be launched with the following command:
```bash
metasbt update --input-list ~/mygenomes.txt \
               --taxa ~/taxonomies.tsv \
               --boundaries ~/boundaries.tsv \
               --db-dir ~/myindex \
               --tmp-dir ~/tmp \
               --dereplicate \
               --similarity 100.0 \
               --completeness 50.0 \
               --contamination 5.0 \
               --type references \
               --extension fna.gz \
               --parallel 4 \
               --nproc 2 \
               --pplacer-threads 2 \
               --cleanup \
               --verbose
```

Please note that `MetaSBT` requires that all your input genomes have the same format and extension before running the pipeline. You can easily uniform your genome files extension by running the `uniform_inputs.sh` script under the `scripts` folder of this repository:
```bash
sh ./scripts/uniform_inputs.sh ~/mygenomes fa fna
```

The first argument is the path to the folder with the set of input genome files, while the second and third arguments are the current and new file extensions respectively. In this particular case, all the `*.fa` files are converted to `*.fna` files.

In case of files with multiple extensions, please run the same script as many times as you need to uniform all the input genome files extensions, until the number of files with the new extension matches the actual number of files in the input folder:
```bash
sh ./scripts/uniform_inputs.sh ~/mygenomes fa fna
sh ./scripts/uniform_inputs.sh ~/mygenomes fasta fna
sh ./scripts/uniform_inputs.sh ~/mygenomes fa.gz fna.gz
```

#### Available options

| Option                   | Default | Mandatory | Description  |
|:-------------------------|:--------|:---------:|:-------------|
| `--boundaries`           |         | ⚑         | Path to the output table produced by the `boundaries` module |
| `--boundary-uncertainty` | `0.0`   |           | Define the percentage of kmers to enlarge and reduce boundaries |
| `--completeness`         | `0.0`   |           | Minimum completeness percentage allowed for input genomes |
| `--contamination`        | `100.0` |           | Maximum contamination percentage allowed for input genomes |
| `--cleanup`              | `False` |           | Remove temporary data at the end of the pipeline |
| `--db-dir`               |         | ⚑         | Database folder path |
| `--dereplicate`          | `False` |           | Dereplicate input genomes |
| `--extension`            |         |           | Specify the input genome files extension. All the input genomes must have the same file extension before running this module |
| `--input-list`           |         | ⚑         | This file contains the list of paths to the new genomes that will be added to the database |
| `--log`                  |         |           | Path to the log file |
| `--nproc`                | `1`     |           | This argument refers to the number of processors used for parallelizing the pipeline when possible |
| `--parallel`             | `1`     |           | Maximum number of processors to process each NCBI tax ID in parallel |
| `--pplacer-threads`      | `1`     |           | Maximum number of threads for pplacer. This is required to maximise the CheckM performances |
| `--similarity`           | `100.0` |           | Dereplicate genomes if they have a percentage of common kmers greater than or equals to the specified one. This is used exclusively in conjunction with the `--dereplicate` argument |
| `--taxa`                 |         | ⚑         | Input file with the mapping between input genome IDs and their taxonomic label. This is used in case of reference genomes only |
| `--tmp-dir`              |         | ⚑         | Path to the temporary folder |
| `--type`                 |         | ⚑         | Define the nature of the input genomes |
| `--verbose`              | `False` |           | Print results on screen |
| `--version`              |         |           | Print current module version and exit |

### 5. `report`: building the database snapshot report

Once the database is built and updated with new MAGs and reference genomes, you can easily extract relevant information about all the species in your database by running the following command:
```bash
metasbt report --db-dir ~/myindex \
               --output-file ~/report.tsv
```

The output file is a table that will contain the number of MAGs and reference genomes, in addition to the mean completeness, contamination, and strain heterogeneity percentages for each lineage in the database. Please note that lineages with no reference genomes correspond to newly defined clusters and potentially new and still-to-be-named species.

#### Available options

| Option                   | Mandatory | Description  |
|:-------------------------|:---------:|:-------------|
| `--db-dir`               | ⚑         | Database folder path |
| `--output-file`          | ⚑         | Path to the output table |
| `--version`              |           | Print current module version and exit |

## Credits

Please credit our work in your manuscript by citing:

> _Manuscript in preparation_

## Contributing

Long-term discussion and bug reports are maintained via GitHub [Discussions](https://github.com/cumbof/MetaSBT/discussions) and [Issues](https://github.com/cumbof/MetaSBT/issues), while code review is managed via GitHub [Pull Requests](https://github.com/cumbof/MetaSBT/pulls).

Please, (i) be sure that there are no existing issues/PR concerning the same bug or improvement before opening a new issue/PR; (ii) write a clear and concise description of what the bug/PR is about; (iii) specifying the list of steps to reproduce the behavior in addition to versions and other technical details is highly recommended.

## Support

If you need support, please open an [Issue](https://github.com/cumbof/MetaSBT/issues) or a new [Discussion](https://github.com/cumbof/MetaSBT/discussions). We will be happy to answer your questions and help you troubleshoot any kind of issue concerning our framework.

Copyright © 2022 [Fabio Cumbo](https://github.com/cumbof), [Daniel Blankenberg](https://github.com/blankenberg). See [LICENSE](https://github.com/cumbof/MetaSBT/blob/main/LICENSE) for additional details.

[^1]: Grüning, Björn, et al. "[Bioconda: sustainable and comprehensive software distribution for the life sciences.](https://doi.org/10.1038/s41592-018-0046-7)" Nature methods 15.7 (2018): 475-476.
[^2]: Parks, Donovan H., et al. "[CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes.](https://doi.org/10.1101/gr.186072.114)" Genome research 25.7 (2015): 1043-1055.
[^3]: Harris, Robert S., and Paul Medvedev. "[Improved representation of sequence bloom trees.](https://doi.org/10.1093/bioinformatics/btz662)" Bioinformatics 36.3 (2020): 721-727.
[^4]: Lemane, Téo, et al. "[kmtricks: Efficient and flexible construction of Bloom filters for large sequencing data collections.](https://doi.org/10.1093/bioadv/vbac029)" Bioinformatics Advances (2022).
[^5]: Mohamadi, Hamid, Hamza Khan, and Inanc Birol. "[ntCard: a streaming algorithm for cardinality estimation in genomics data.](https://doi.org/10.1093/bioinformatics/btw832)" Bioinformatics 33.9 (2017): 1324-1330.