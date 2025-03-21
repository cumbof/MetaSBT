#!/usr/bin/env python3
"""A scalable framework for automatically indexing microbial genomes and accurately
characterizing metagenome-assembled genomes with Sequence Bloom Trees.
"""

import argparse
import gzip
import json
import multiprocessing
import os
import requests
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.request

from datetime import datetime
from packaging import version
from pathlib import Path
from tabulate import tabulate
from typing import Any, List

from metasbt.core import __date__, __version__, DEPENDENCIES, Database, Entry

# Define the URL to the code repository
CODE_REPOSITORY_URL = "https://github.com/cumbof/MetaSBT"

# Define the URL to the public databases repository
DB_REPOSITORY_URL = "https://github.com/cumbof/MetaSBT-DBs"

# Define the URL to the software releases page
CODE_RELEASES_API_URL = "https://api.github.com/repos/cumbof/MetaSBT/releases/latest"

# Define the path to the table with the list of public databases
DB_LIST_URL = "https://raw.githubusercontent.com/cumbof/MetaSBT-DBs/main/databases.tsv"

# Define the list external software dependencies
DEPENDENCIES.extend([
    "kraken2-build",  # Required by the kraken subroutine
    "tar"  # Required by the pack and unpack subroutines
])


class MetaSBT(object):
    """The MetaSBT class object with the definition of its commands and subcommands.
    """

    def __init__(self) -> "MetaSBT":
        """Initialize the MetaSBT object for parsing input commands, subcommands, and arguments.

        Raises
        ------
        Exception
            In case of missing software dependencies.

        Returns
        -------
        MetaSBT
            A MetaSBT object.
        """

        usage = """
        metasbt <command> [<args>]

        The metasbt commands are:
        db          List and retrieve public MetaSBT databases;
        index       Index a set of reference genomes and build the first baseline of a MetaSBT database;
        kraken      Export a MetaSBT database into a custom kraken database;
        pack        Build a compressed tarball with a MetaSBT database and report its sha256;
        profile     Profile an input genome and report the closest cluster at all the seven taxonomic levels
                    and the closest genome in a MetaSBT database;
        sketch      Sketch the input genomes;
        summarize   Summarize the content of a MetaSBT database and report some statistics;
        test        Check for dependencies and run unit tests;
                    This must be used by code maintainers only;
        unpack      Unpack a local MetaSBT tarball database;
        update      Update a MetaSBT database with new metagenome-assembled genomes.
        """

        self.start_at = time.time()

        # Keep track of the Database object
        # This is required because `update` calls `profile` and both initialize a Database object
        self.database = None

        missing_dependencies = list()

        # Check for external software dependencies
        for dependency in DEPENDENCIES:
            if shutil.which(dependency) is None:
                missing_dependencies.append(dependency)

        if missing_dependencies:
            raise Exception(f"Missing software dependencies: {', '.join(missing_dependencies)}")

        try:
            # Check whether a new software version is available
            response = requests.get(CODE_RELEASES_API_URL)

            if response.status_code == 200:
                response = response.json()

                if "tag_name" in response:
                    if version.parse(__version__) < version.parse(response["tag_name"]):
                        print("A new release of MetaSBT is available!")
                        print(f"{CODE_REPOSITORY_URL}\n")

        except requests.ConnectionError:
            # This could happen in case of no internet connection
            pass

        parser = argparse.ArgumentParser(
            prog="MetaSBT",
            usage=usage,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # This is one of the available command
        # Look at the usage
        parser.add_argument("command", help="metasbt command")

        parser.add_argument(
            "-v", "--version",
            action="version",
            version=__version__,
            help="Show version and exit."
        )

        # Read the first argument
        # This is the name of the command
        args = parser.parse_args([sys.argv[1]])

        # Define the command arguments
        subargs = sys.argv[2:]

        # Use dispatch pattern to invoke method with same name of the argument
        getattr(self, args.command)(subargs)

    def __print_credits(self) -> None:
        """Print credits.
        """

        message = f"""
        Thanks for using MetaSBT!

        Please credit our software and databases by citing:

        Fabio Cumbo, Daniel Blankenberg ({__date__.split()[-1]})
        cumbof/MetaSBT: MetaSBT (Version {__version__}) [Computer software]
        {CODE_REPOSITORY_URL}
        {DB_REPOSITORY_URL}
        """

        print(message)

    def db(self, argv: List[Any]) -> None:
        """List and retrieve public MetaSBT databases.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        Exception
            - In case of an unexpected while retrieving the list of public databases from the MetaSBT-DBs repository;
            - If there are no public databases available;
            - If the provided database does not exist;
            - If the provided database version does not exist;
            - In case it is unable to download the tarball.
        """

        parser = argparse.ArgumentParser(
            prog="db",
            description="List and retrieve public MetaSBT databases.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--list",
            action="store_true",
            default=False,
            help="List official public MetaSBT databases.",
        )
        parser.add_argument(
            "--download",
            required="--list" not in argv,
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--version",
            type=str,
            help="The database version. It automatically select the most recent one if a version is not provided."
        )
        parser.add_argument(
            "--folder",
            type=os.path.abspath,
            default=os.getcwd(),
            help="Store the selected database under this folder."
        )

        # Load arguments
        args = parser.parse_args(argv)

        databases = dict()

        # Print the list of official public MetaSBT databases from the MetaSBT-DBs repository
        with tempfile.TemporaryDirectory() as tmpdir:
            db_table_filepath = os.path.join(tmpdir, "databases.tsv")

            try:
                # Download the databases table
                urllib.request.urlretrieve(DB_LIST_URL, db_table_filepath)

            except:
                raise Exception(f"Unable to retrieve the list of databases from {DB_LIST_URL}")

            with open(db_table_filepath) as db_table:
                header = list()

                for line in db_table:
                    line = line.strip()

                    if line:
                        if line.startswith("#"):
                            header = line[1:].strip().split("\t")

                        else:
                            line_split = line.split("\t")

                            # Read databases information
                            if line_split[header.index("id")] not in databases:
                                databases[line_split[header.index("id")]] = list()

                            # id, version, references, mags, size, tarball, sha256, and info
                            db_version = {value: line_split[pos] for pos, value in enumerate(header)}

                            databases[line_split[header.index("id")]].append(db_version)

        if not databases:
            raise Exception("No public databases available!")

        if args.list:
            # Print the list of databases and exit
            table = [["ID", "Version", "References", "MAGs", "Size", "Info"]]

            for db in databases:
                for db_version in sorted(databases[db], key=lambda v: int(v["version"])):
                    table.append(
                        [
                            db,
                            db_version["version"],
                            db_version["references"],
                            db_version["mags"],
                            db_version["size"],
                            db_version["info"]
                        ]
                    )

            print(tabulate(table, headers="firstrow", tablefmt="fancy_grid"))

        elif args.download:
            if args.download not in databases:
                raise Exception("The provided database does not exist!")

            # Retrieve versions of the same database
            versions = [db_version["version"] for db_version in databases[args.database]]

            # Select the latest version by default
            selected_version = str(sorted([int(v) for v in versions])[-1])

            if args.version:
                if args.version not in versions:
                    raise Exception("The provided database version does not exist!")

                selected_version = args.version

            # Retrieve the URL to the tarball
            database_url = [db_version["tarball"] for db_version in databases[args.database]][0]

            try:
                # Download the database tarball
                urllib.request.urlretrieve(database_url, os.path.join(args.folder, os.path.basename(database_url)))

            except:
                raise Exception(f"Unable to retrieve the database tarball from {database_url}")

    def index(self, argv: List[Any]) -> None:
        """Build the first baseline of a MetaSBT database by indexing a set of reference genomes.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        Exception
            If a MetaSBT database with the same name of the provided one already exists.
        """

        parser = argparse.ArgumentParser(
            prog="index",
            description=(
                "Index a set of reference genomes. "
                "This is used to build a first baseline of a MetaSBT database. "
                "Genomes must be known with a fully defined taxonomic label, from the kingdom up to the species level."
            ),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # General arguments
        general_group = parser.add_argument_group("General arguments")

        general_group.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        general_group.add_argument(
            "--database",
            required=True,
            type=str,
            help="The database name."
        )
        general_group.add_argument(
            "--references",
            required=True,
            type=os.path.abspath,
            help=(
                "Path to the tab-separated-values file with the list of reference genomes. "
                "It must contain two columns. The first one with the path to the actual reference genomes. "
                "The second one with their fully defined taxonomic label in the following form: "
                "k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species"
            )
        )
        general_group.add_argument(
            "--dereplicate",
            required=False,
            default=0.0,
            type=float,
            help=(
                "Dereplicate genomes based of their ANI distance according the specified threshold. "
                "The dereplication process is triggered in case of a threshold >0.0."
            )
        )
        general_group.add_argument(
            "--nproc",
            required=False,
            default=os.cpu_count(),
            type=int,
            help="Process the input genomes in parallel."
        )
        general_group.add_argument(
            "--pack",
            action="store_true",
            default=False,
            help="Pack the database into a compressed tarball.",
        )

        # Group of arguments for estimating the bloom filter size
        filter_size_group = parser.add_argument_group("Estimate a proper bloom filter size")

        filter_size_group.add_argument(
            "--filter-size",
            type=int,
            required=False,
            dest="filter_size",
            help=(
                "This is the size of the bloom filters. "
                "It automatically estimates a proper bloom filter size if not provided."
            )
        )
        filter_size_group.add_argument(
            "--increase-filter-size",
            type=float,
            default=50.0,
            dest="increase_filter_size",
            help=(
                "Increase the estimated filter size by the specified percentage. "
                "It is highly recommended to increase the filter size by a good percentage in case you are planning to update the index with new genomes."
            )
        )
        filter_size_group.add_argument(
            "--min-kmer-occurrences",
            type=int,
            default=2,
            dest="min_kmer_occurrences",
            help=(
                "Minimum number of occurrences of kmers to be considered for estimating the bloom filter size "
                "and for building the bloom filter files."
            )
        )

        # Group of arguments for estimating the optimal kmer size
        kitsune_group = parser.add_argument_group("Estimate the optimal kmer size")

        kitsune_group.add_argument(
            "--kmer-size",
            type=int,
            required=False,
            dest="kmer_size",
            help=(
                "The kmer size. "
                "It automatically estimates a proper bloom filter size if not provided."
            )
        )
        kitsune_group.add_argument(
            "--limit-kmer-size",
            type=int,
            default=32,
            dest="limit_kmer_size",
            help="Limit the estimation of the optimal kmer size with kitsune to this size at most.",
        )

        # Group of arguments for assessing the quality of genomes
        quality_group = parser.add_argument_group("Assess the quality of input genomes")

        quality_group.add_argument(
            "--completeness",
            type=float,
            required=False,
            default=0.0,
            help="Percentage threshold on genomes completeness.",
        )
        quality_group.add_argument(
            "--contamination",
            type=float,
            required=False,
            default=100.0,
            help="Percentage threshold on genomes contamination.",
        )

        # Load arguments
        args = parser.parse_args(argv)

        # Create the working dirctory if it does not exist
        os.makedirs(args.workdir, exist_ok=True)

        # Define the path to the database folder
        db_dir = os.path.join(args.workdir, args.database)

        if os.path.isdir(db_dir):
            # The `index` command must be used to initialize a database
            # If a database with the same name of the provided one already exists, consider running the `update` command
            raise Exception(f"A database named '{args.database}' already exist under {args.workdir}")

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        # Initialize the database
        # Use a flat structure by default
        self.database = Database(args.database, db_dir, tmp_dir, flat=True, nproc=args.nproc)

        # Load the list of reference genomes and their taxonomic labels
        references = dict()

        with open(args.references) as input_table:
            for line in input_table:
                line = line.strip()

                if line:
                    # Header and comment lines always start with the pound character
                    if not line.startswith("#"):
                        line_split = line.split("\t")

                        # The path to the genome file is under the first column
                        genome_filepath = line_split[0]

                        # The taxonomic label is under the second column
                        taxonomy = line_split[1]

                        # Index the references dict by the genome filepath
                        references[genome_filepath] = taxonomy

        # Define the set of paths to the reference genomes
        genomes = set(references.keys())

        # Define the database metadata
        # Eventually, estimate the optimal kmer size and a proper bloom filter size
        self.database.set_configs(
            genomes,
            min_kmer_occurrence=args.min_kmer_occurrences,
            kmer_size=args.kmer_size,
            kmer_max=args.limit_kmer_size,
            filter_size=args.filter_size,
            filter_expand_by=args.increase_filter_size
        )

        if args.completeness > 0.0 or args.contamination < 100.0:
            # Retrieve the kingdom
            # Assume the input genomes are all under the same kingdom
            kingdom = references[list(references.keys())[0]].split("|")[0][3:]

            # Asses the quality of the input genomes
            quality = Database.qc(genomes, kingdom, nproc=args.nproc, tmp=tmp_dir)

            # Filter out genomes according to the completeness and contamination thresholds
            genomes = {genome for genome in genomes if quality[genome]["completeness"] >= args.completeness and quality[genome]["contamination"] <= args.contamination}

        if args.dereplicate > 0.0:
            # Dereplicate genomes based on their ANI distance
            genomes = self.database.dereplicate(genomes, threshold=args.dereplicate)

        # Reshape the references dict
        # Consider genomes that passed the dereplication process only
        references = {genome: references[genome] for genome in genomes}

        # Add references to the database
        for genome in references:
            self.database.add(genome, reference=True, taxonomy=references[genome])

        # Finally, index references
        self.database.update()

        if args.pack:
            # Pack the database into a compressed tarball
            self.pack(["--workdir", args.workdir, "--database", args.database])

        # Print credits
        self.__print_credits()

    def kraken(self, argv: List[Any]) -> None:
        """Export a MetaSBT database into a custom kraken database.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        ValueError
            If a node is missing from the NCBI nodes.dmp file.
        Exception
            In case of an error occurred while running `kraken2-build`.
        """

        parser = argparse.ArgumentParser(
            prog="kraken",
            description="Export a MetaSBT database into a custom kraken database.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--genomes",
            type=os.path.abspath,
            required=True,
            help=(
                "Path to the file with the list of paths to the genomes. "
                "Genomes must be in the MetaSBT database in order to be processed."
            ),
        )
        parser.add_argument(
            "--ncbi-names",
            type=os.path.abspath,
            required=True,
            dest="ncbi_names",
            help="Path to the NCBI names.dmp file.",
        )
        parser.add_argument(
            "--ncbi-nodes",
            type=os.path.abspath,
            required=True,
            dest="ncbi_nodes",
            help="Path to the NCBI nodes.dmp file.",
        )

        # Load arguments
        args = parser.parse_args(argv)

        # Load the list of paths to the genome files
        # This is a dict with the paths to the input genomes indexed by their file name
        # Assume the input files are not compressed
        genomes = {os.path.splitext(os.path.basename(line.strip()))[0]: line.strip() for line in open(args.genomes).readlines() if line.strip()}

        ncbi_names = dict()

        # Keep track of the highgest NCBI TaxID
        latest_node_id = 0

        # Load NCBI names
        with open(args.ncbi_names) as ncbi_names_table:
            for line in ncbi_names_table:
                line = line.strip()

                if line:
                    line_split = line.split("|")

                    if line_split[3].strip() == "scientific name":
                        # NCBI TaxID
                        ncbi_taxid = int(line_split[0].strip())

                        # NCBI Tax scientific name
                        # TODO scientific names must be fixed in the same way that MetaSBT does
                        ncbi_tax_scientific_name = line_split[1].strip()

                        ncbi_names[ncbi_tax_scientific_name] = ncbi_taxid

                        if ncbi_taxid > latest_node_id:
                            latest_node_id = ncbi_taxid

        ncbi_nodes = dict()

        # Load NCBI tax relationships
        with open(args.ncbi_nodes) as ncbi_nodes_table:
            for line in ncbi_nodes_table:
                line = line.strip()

                if line:
                    line_split = line.split("|")

                    # NCBI TaxID of the current node
                    node_id = int(line_split[0].strip())

                    # NCBI TaxID of the parent node
                    parent_id = int(line_split[1].strip())

                    # NCBI Tax rank
                    node_rank = line_split[2].strip()

                    ncbi_nodes[node_id] = {
                        "parent": parent_id,
                        "rank": node_rank
                    }

        # Define the path to the clusters.tsv file
        metasbt_report_filepath = os.path.join(args.workdir, args.database, "clusters.tsv")

        metasbt_names = dict()

        metasbt_nodes = dict()

        # Take track of the full taxonomic labels and genomes in MetaSBT
        metasbt_lineages = dict()

        # Load MetaSBT clusters
        with open(metasbt_report_filepath) as metasbt_report:
            header = list()

            for line in metasbt_report:
                line = line.strip()

                if line:
                    if line.startswith("#"):
                        header = line[1:].strip().split("\t")

                    else:
                        line_split = line.split("\t")

                        lineage = line_split[header.index("taxonomy")].split("|")

                        if lineage[-1].startswith("s__"):
                            lineage_genomes = line_split[header.index("references_list")].split(",") + line_split[header.index("mags_list")].split(",")

                            metasbt_lineages["|".join(lineage)] = lineage_genomes

                        for rank_position, rank in enumerate(lineage):
                            cluster_name = rank[3:]

                            # The parent node name of a superkingdom in NCBI is root
                            parent_name = "root"

                            if rank_position > 0:
                                parent_name = lineage[rank_position-1][3:]

                            if cluster_name.startswith("MSBT"):
                                latest_node_id += 1

                                metasbt_names[cluster_name] = latest_node_id

                            if metasbt_names[cluster_name] not in metasbt_nodes:
                                parent_id = None

                                if parent_name in ncbi_names:
                                    parent_id = ncbi_names[parent_name]

                                elif parent_name in metasbt_names:
                                    parent_id = metasbt_names[parent_name]

                                else:
                                    raise ValueError('Node "{}" is missing!'.format(parent_name))

                                metasbt_nodes[metasbt_names[cluster_name]] = {
                                    "parent": parent_id,
                                    "rank": TAXONOMIC_RANKING[rank_position]
                                }

        # Finally, write the MetaSBT names and nodes additions
        with open(args.ncbi_names, "a+") as ncbi_names_table:
            for node_name, node_id in sorted(metasbt_names.items(), key=lambda node_name: metasbt_names[node_name]):
                # Add the tax id
                metasbt_node = NAME_TEMPLATE.replace("taxid", node_id)

                # Add the node name
                metasbt_node = metasbt_node.replace("var", node_name)

                ncbi_names_table.write("{}\n".format(metasbt_node))

        with open(args.ncbi_nodes, "a+") as ncbi_nodes_table:
            for node_id, node_data in sorted(metasbt_nodes.items()):
                # Add the tax id
                metasbt_node = NODE_TEMPLATE.replace("taxid", node_id)

                # Add the parent node name
                metasbt_node = metasbt_node.replace("parent", node_data["parent"])

                # Add the rank
                metasbt_node = metasbt_node.replace("rank", node_data["rank"])

                ncbi_nodes_table.write("{}\n".format(metasbt_node))

        for lineage in metasbt_lineages:
            # Intersect the list of genomes under this lineage with the input genomes
            selected_genomes = set(genomes.keys()).intersection(set(metasbt_lineages[lineage]))

            if selected_genomes:
                # Lineages are always defined up to the species level here
                lineage = lineage.split("|")

                # Retrieve the list of genome files for the selection
                genome_files = [genomes[genome] for genome in selected_genomes]

                metasbt_cluster_name = lineage[-1][3:]

                metasbt_cluster_id = metasbt_names[metasbt_cluster_name] if metasbt_cluster_name in metasbt_names else ncbi_names[metasbt_cluster_name]

                for genome_file in genome_files:
                    with tempfile.NamedTemporaryFile() as tmp_genome, open(genome_file) as genome:
                        with open(tmp_genome.name, "wt"), open(genome_file) as genome:
                            for record in SeqIO.parse(genome, "fasta"):
                                # Rebuild the record with the new ID
                                record = SeqRecord(
                                    record.seq, 
                                    id="{}|kraken:taxid|{} {}".format(record.id, metasbt_cluster_id, record.id), 
                                    name=record.name, 
                                    description=record.description
                                )

                                SeqIO.write(record, tmp_genome_file, "fasta")

                        # Run kraken2-build over the temporary genome file
                        # kraken2-build --add-to-library $file --db $DBNAME
                        command_line = [
                            "kraken2-build",
                            "--add-to-library",
                            tmp_genome.name,
                            "--db",
                            args.database
                        ]

                        try:
                            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                        except subprocess.CalledProcessError as e:
                            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                            raise Exception(error_message).with_traceback(e.__traceback__)

    def pack(self, argv: List[Any]) -> None:
        """Build a compressed tarball with a MetaSBT database and report its sha256.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        FileNotFoundError
            If the specified database does not exist.
        Exception
            - In case of an error occurred while producing the compressed tarball;
            - In case of an error occurred while computing the sha256 hash of the compressed tarball.
        """

        parser = argparse.ArgumentParser(
            prog="pack",
            description="Pack a MetaSBT database into a compressed tarball.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )

        # Load arguments
        args = parser.parse_args(argv)

        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Get the current time
        now = datetime.now()

        # Define a timestamp
        timestamp = f"{now.year}{now.month}{now.day}"

        # Search for other tarballs under the same working directory
        tarballs = Path(args.workdir).glob(f"MetaSBT-{args.database}-*.tar.gz")

        # Define the path to the output tarball
        output_filepath = os.path.join(args.workdir, f"MetaSBT-{args.database}-{timestamp}-{len(list(tarballs))+1}.tar.gz")

        # Build the compressed tarball
        command_line = [
            "tar",
            "-czvf",
            output_filepath,
            "-C",
            args.workdir,
            db_dir
        ]

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message).with_traceback(e.__traceback__)

        # Compute the sha256 hash
        sha256 = subprocess.run(["sha256sum", output_filepath], capture_output=True, text=True)

        if result.returncode == 0:
            print(sha256.stdout)

        raise Exception(f"An error has occurred while computing the sha256 hash of {output_filepath}")

    def profile(self, argv: List[Any]) -> None:
        """Profile a set of genomes against a specific MetaSBT database.
        Profile tables are stored under a dedicated folder in the workdir temporary directory.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        FileNotFoundError
            If the specified database does not exist.
        """

        parser = argparse.ArgumentParser(
            prog="profile",
            description=(
                "Profile a set of genomes. "
                "This is used to report the closest kingdom, phylum, class, order, family, genus, species, "
                "and the closest genome in a specific database. "
            ),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--genome",
            required="--genomes" not in argv,
            type=os.path.abspath,
            help="Path to the input genome."
        )
        parser.add_argument(
            "--genomes",
            required="--genome" not in argv,
            type=os.path.abspath,
            help="Path to the file with a list of paths to the input genomes."
        )
        parser.add_argument(
            "--nproc",
            required=False,
            type=int,
            default=os.cpu_count(),
            help="Process the input genomes in parallel."
        )

        # Load arguments
        args = parser.parse_args(argv)
        
        # Define the path to the database folder
        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        if self.database is None:
            # Initialize the database
            # Use a flat structure by default
            self.database = Database(args.database, db_dir, tmp_dir, flat=True, nproc=args.nproc)

        # Load the list of paths to the input genomes
        genomes = [args.genome] if args.genome else [line.strip() for line in open(args.genomes).readlines() if line.strip()]

        # Genomes must be sketched first
        # Note that this function and `sketch()` have the same set of arguments
        # Sketches are stored under the dedicated folder in the database directory
        sketches = self.sketch(argv)

        # Profile genomes in parallel
        with multiprocessing.Pool(processes=args.nproc) as pool, tqdm.tqdm(total=len(genomes)) as progress_bar:
            def progress(*args):
                progress_bar.update()

            # TODO Uncertainty and pruning threshold should not be hardcoded
            jobs = [
                pool.apply_async(
                    Database._profile, 
                    args=(
                        self.database,
                        genome_filepath,
                        sketch_filepath, 
                        20.0,  # uncertainty
                        0.0,  # pruning threshold
                    ),
                    callback=progress
                ) for genome_filepath, sketch_filepath in zip(genomes, sketches)
            ]

            for job in jobs:
                _, _ = job.get()

    def sketch(self, argv: List[Any]) -> List[os.path.abspath]:
        """Sketch the input genomes.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        FileNotFoundError
            If the specified database does not exist.

        Returns
        -------
        list
            A list with paths to the sketch files.
        """

        parser = argparse.ArgumentParser(
            prog="sketch",
            description="Sketch the input genomes.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--genome",
            required="--genomes" not in argv,
            type=os.path.abspath,
            help="Path to the input genome."
        )
        parser.add_argument(
            "--genomes",
            required="--genome" not in argv,
            type=os.path.abspath,
            help="Path to the file with a list of paths to the input genomes."
        )
        parser.add_argument(
            "--nproc",
            required=False,
            type=int,
            default=os.cpu_count(),
            help="Process the input genomes in parallel."
        )

        # Load arguments
        args = parser.parse_args(argv)
        
        # Define the path to the database folder
        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        if self.database is None:
            # Initialize the database
            # Use a flat structure by default
            self.database = Database(args.database, db_dir, tmp_dir, flat=True, nproc=args.nproc)

        # Load the list of paths to the input genomes
        genomes = [args.genome] if args.genome else [line.strip() for line in open(args.genomes).readlines() if line.strip()]

        sketches = list()

        # Genomes must be sketched first
        for genome_filepath in genomes:
            # Retrieve the genome name
            genome_name = os.path.splitext(os.path.basename(genome_filepath))[0]

            # Define the genome object as MetaSBT Entry
            genome_object = Entry(self.database, genome_name, genome_name, "genome")

            # Build the bloom filter representation of the genome
            sketch_filepath = genome_object.sketch(genome_filepath)

            # Keep track of the path to the sketch file
            sketches.append(sketches)

        return sketches

    def summarize(self, argv: List[Any]) -> None:
        """Summarize the content of a MetaSBT database and report some statistics.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        FileNotFoundError
            In case one of the following files does not exist in the database: `clusters.tsv`, `genomes.tsv`, and `metadata.json`.
        """

        parser = argparse.ArgumentParser(
            prog="summarize",
            description="Summarize the content of a MetaSBT database.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )

        # Load arguments
        args = parser.parse_args(argv)

        db_dir = os.path.join(args.workdir, args.database)

        # Search for the following files
        search_for = [os.path.join(db_dir, "metadata.json"), os.path.join(db_dir, "genomes.tsv"), os.path.join(db_dir, "clusters.tsv")]

        for filepath in search_for:
            if not os.path.isfile(filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

        table = list()

        # Load metadata.json
        with open(search_for[0]) as metadata_json_file:
            metadata = json.loads("".join(metadata_json_file.readlines()))

            # Add metadata to the summary table
            table.extend([[key, metadata[key]] for key in metadata])

        # Count the number of clusters at all the seven taxonomic levels
        clusters = {"kingdom": 0, "phylum": 0, "class": 0, "order": 0, "family": 0, "genus": 0, "species": 0}

        # Also count the number of known clusters at all the seven taxonomic levels
        knowns = {level: 0 for level in clusters}

        # Count the number of reference genomes and MAGs
        genomes = {"references": 0, "mags": 0}

        # Keep track of the bloom filter density at the root
        density = 0.0

        with open(search_for[2]) as clusters_table:
            header = list()

            for line in clusters_table:
                if line.strip():
                    if line.startswith("#"):
                        header = line[1:].strip().split("\t")

                    else:
                        line_split = line.split("\t")

                        # Strip the new line character out of the last element 
                        line_split[-1] = line_split[-1].strip()

                        # Increment the clusters counter based on the taxonomic level
                        clusters[line_split[header.index("level")]] += 1

                        if eval(line_split[header.index("known")]):
                            # Increment the known clusters counter based on the taxonomic level
                            knowns[line_split[header.index("level")]] += 1 

                        if line_split[header.index("level")] == "kingdom":
                            # Update the bloom filter density
                            density = float(line_split[header.index("density")])

                        if line_split[header.index("level")] == "species":
                            # Increment the number of reference genomes and MAGs
                            genomes["references"] += int(line_split[header.index("references_count")])

                            genomes["mags"] += int(line_split[header.index("mags_count")])

        # Add the number of genomes to the summary table
        table.append(["references", genomes["references"]])

        table.append(["mags", genomes["mags"]])

        # Finally, add stats to the summary table
        # Add the number of known clusters over the total number of clusters per taxonomic level
        for level in clusters:
            table.append([level, f"{knowns[level]}/{clusters[level]}"])

        # Add the root density
        table.append(["density", density])

        print(tabulate(table, tablefmt="fancy_grid"))

    def test(self, argv: List[Any]) -> None:
        """Check for software dependencies and run unit tests.

        Parameters
        ----------
        argv : list
            The list of arguments.
        """

        parser = argparse.ArgumentParser(
            prog="test",
            description="Check for software dependencies and run unit tests.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "test",
            default="all",
            type=str,
            choices=["all", "db", "index", "kraken", "pack", "profile", "sketch", "summarize", "unpack", "update"],
            help="The database name."
        )

        # Load arguments
        args = parser.parse_args(argv)

        class Test(unittest.TestCase):
            """Unit tests.
            """

            common_root = "k__MSBT1|p__MSBT2|c__MSBT3|o__MSBT4|f__MSBT5"

            # Reference genomes and their taxonomic labels
            references = {
                "GCA_023531975.1_ASM2353197v1": f"{common_root}|g__MSBT6|s__Orf_virus",
                "GCA_023536435.1_ASM2353643v1": f"{common_root}|g__MSBT6|s__Orf_virus",
                "GCA_024425995.1_ASM2442599v1": f"{common_root}|g__MSBT6|s__Orf_virus",
                "GCA_025133395.1_ASM2513339v1": f"{common_root}|g__MSBT7|s__Monkeypox_virus",
                "GCA_025133395.1_ASM2513339v1": f"{common_root}|g__MSBT7|s__Monkeypox_virus",
                "GCA_025627565.1_ASM2562756v1": f"{common_root}|g__MSBT7|s__Monkeypox_virus",
                "GCA_001745695.1_ViralProj344115": f"{common_root}|g__MSBT7|s__Skunkpox_virus",
                "GCA_001744115.1_ViralProj344208": f"{common_root}|g__MSBT7|s__Skunkpox_virus",
            }

            # Metagenome-assembled genomes
            mags = [
                "GCA_023701625.1_ASM2370162v1",  # s__Monkeypox_virus
                "GCA_007575705.1_ASM757570v1",  # s__Monkeypox_virus
                "GCA_006451315.1_ASM645131v1",  # s__Lumpy_skin_disease_virus
                "GCA_009651155.1_ASM965115v1",  # s__Orf_virus
            ]

            @classmethod
            def setUpClass(cls):
                """ Set up the unit class by retrieving a bunch of genomes to use within all the unit tests.
                """

                # Create the temporary working directory
                cls.working_dir = tempfile.TemporaryDirectory()

                # Define the database folder
                cls.db_dir = os.path.join(cls.working_dir.name, "MetaSBT")

                # Define the temporary folder
                cls.tmp_dir = os.path.join(cls.working_dir.name, "tmp")

                # Retrieve a bunch of genomes from NCBI GenBank to use during the tests
                genomes = set(cls.mags)

                genomes.update(set([genome for taxonomy in cls.references for genome in cls.references[taxonomy]]))

                for genome in genomes:
                    # Define the URL to the fna.gz file
                    fna_gz_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/{4}/{4}_genomic.fna.gz".format(
                        genome[0:3],
                        genome[4:7],
                        genome[7:10],
                        genome[10:13],
                        genome
                    )

                    # Define the path to the local fna.gz file
                    fna_gz_path = os.path.join(cls.working_dir.name, "{}_genomic.fna.gz".format(genome))

                    # Download the genome into the temporary working directory
                    urllib.request.urlretrieve(fna_gz_url, fna_gz_path)

                    # Unzip the genome
                    with gzip.open(fna_gz_path, "rb") as fna_gz_in:
                        with open(os.path.splitext(fna_gz_path)[0], "wb") as fna_out:
                            shutil.copyfileobj(fna_gz_in, fna_out)

            @classmethod
            def tearDownClass(cls):
                """ Get rid of the temporary working directory with the test genomes.
                """

                # Delete the temporary working directory
                shutil.rmtree(cls.__working_dir.name, ignore_errors=False)

            def db(self):
                """Test the `db()` function.
                """

                # TODO

            def index(self):
                """Test the `index()` function.
                """

                # TODO

            def kraken(self):
                """Test the `kraken()` function.
                """

                # TODO

            def pack(self):
                """Test the `pack()` function.
                """

                # TODO

            def profile(self):
                """Test the `profile()` function.
                """

                # TODO

            def sketch(self):
                """Test the `sketch()` function.
                """

                # TODO

            def summarize(self):
                """Test the `summarize()` function.
                """

                # TODO

            def unpack(self):
                """Test the `unpack()` function.
                """

                # TODO

            def update(self):
                """Test the `update()` function.
                """

                # TODO

        def run_test(test_id: str="all") -> None:
            """Run unit tests.

            Parameters
            ----------
            test_id : str, default "all"
                The name of the test or the wildcard "all" to run all the unit tests.
            """

            # We should manually call the setUpClass
            Test.setUpClass()

            if test_id == "all":
                # Initialize the test loader in case of "all"
                loader = unittest.TestLoader()

                # Load all the defined tests in Test
                suite = loader.loadTestsFromTestCase(Test)

            else:
                # Use a suite for all the other cases
                suite = unittest.TestSuite()

                # Load the specific test
                suite.addTest(Test(test_id))

            runner = unittest.TextTestRunner()

            # Finally, run the unit tests
            runner.run(suite)

            # Remove the test data
            Test.tearDownClass()

        run_test(args.test)

    def unpack(self, argv: List[Any]) -> None:
        """Unpack a local MetaSBT tarball database.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        Exception
            - If a database with the same of the tarball one already exists;
            - In case of an error occurred while unpacking the tarball.
        """

        parser = argparse.ArgumentParser(
            prog="unpack",
            description="Unack a local MetaSBT tarball database.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--tarball",
            type=os.path.abspath,
            help="Path to the MetaSBT tarball database."
        )

        # Load arguments
        args = parser.parse_args(argv)

        # Retrieve the database name from the tarball file name
        database = os.path.basename(args.tarball).split("-")[1]

        db_dir = os.path.join(args.workdir, database)

        if os.path.isdir(db_dir):
            raise Exception(f"A database with the same name already exists: {database}")

        # Build the compressed tarball
        command_line = [
            "tar",
            "-xzvf",
            args.tarball,
            "-C",
            args.workdir
        ]

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message).with_traceback(e.__traceback__)

    def update(self, argv: List[Any]) -> None:
        """Update a specific MetaSBT database with new metagenome-assembled genomes.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        FileNotFoundError
            If the specified database does not exist.
        """

        parser = argparse.ArgumentParser(
            prog="update",
            description="Update a MetaSBT database with new genomes.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
            "--workdir",
            required=True,
            type=os.path.abspath,
            help="Path to the working directory."
        )
        parser.add_argument(
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--genome",
            required="--genomes" not in argv,
            type=os.path.abspath,
            help="Path to the input genome."
        )
        parser.add_argument(
            "--genomes",
            required="--genome" not in argv,
            type=os.path.abspath,
            help="Path to the file with a list of paths to the input genomes."
        )
        parser.add_argument(
            "--dereplicate",
            required=False,
            default=0.0,
            type=float,
            dest="dereplicate",
            help=(
                "Dereplicate genomes based of their ANI distance according the specified threshold. "
                "The dereplication process is triggered in case of a threshold >0.0."
            )
        )
        parser.add_argument(
            "--completeness",
            type=float,
            required=False,
            default=0.0,
            help="Percentage threshold on genomes completeness.",
        )
        parser.add_argument(
            "--contamination",
            type=float,
            required=False,
            default=100.0,
            help="Percentage threshold on genomes contamination.",
        )
        parser.add_argument(
            "--nproc",
            required=False,
            type=int,
            default=os.cpu_count(),
            help="Process the input genomes in parallel."
        )
        parser.add_argument(
            "--pack",
            action="store_true",
            default=False,
            help="Pack the database into a compressed tarball.",
        )

        # Load arguments
        args = parser.parse_args(argv)
        
        # Define the path to the database folder
        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        # Initialize the database
        # Use a flat structure by default
        self.database = Database(args.database, db_dir, tmp_dir, flat=True, nproc=args.nproc)

        # Load the set of paths to the input genomes
        # TODO There could be reference genomes as well
        genomes = {args.genome} if args.genome else {line.strip() for line in open(args.genomes).readlines() if line.strip()}

        if args.completeness > 0.0 or args.contamination < 100.0:
            # Retrieve the kingdom from the root node of the target database
            # Assume the input genomes are all under the same kingdom
            kingdom = list(self.database.clusters["kingdom"].keys())[0].split("|")[0][3:]

            # Asses the quality of the input genomes
            quality = Database.qc(genomes, kingdom, nproc=args.nproc, tmp=tmp_dir)

            # Filter out genomes according to the completeness and contamination thresholds
            genomes = {genome for genome in genomes if quality[genome]["completeness"] >= args.completeness and quality[genome]["contamination"] <= args.contamination}

        if args.dereplicate > 0.0:
            # Dereplicate genomes based on their ANI distance
            # Reshape the set of genomes
            # Consider genomes that passed the dereplication process only (input-vs-input)
            genomes = self.database.dereplicate(genomes, threshold=args.dereplicate)

        # Profile the input genomes first
        # This also produce the bloom filter representation of the input genomes
        # Note that this function and `profile()` have the same set of arguments
        # Genomes profiles are stored under the dedicated folder in the workdir temporary directory
        self.profile(argv)

        if args.dereplicate > 0.0:
            # Dereplicate the input genomes again versus the genomes in the database
            genomes = self.database.dereplicate(genomes, threshold=args.dereplicate, compare_with="database")

        species_assignments = dict()

        with multiprocessing.Pool(processes=args.nproc) as pool, tqdm.tqdm(total=len(genomes)) as progress_bar:
            def progress(*args):
                progress_bar.update()

            # Parse the genome profiles and check whether they are close enough to already defined clusters at the species level
            jobs = [
                pool.apply_async(
                    Database._is_known,
                    args=(
                        self.database,
                        genome_filepath,
                    ),
                    callback=progress
                ) for genome_filepath in genomes
            ]

            for job in jobs:
                # Retrieve the assigned species
                # It could be None if the input genome is not close enough to any species clusters in the database
                genome_filepath, taxonomy = job.get()

                # Keep track of the species assignments
                species_assignments[genome_filepath] = taxonomy

        for genome_filepath in species_assignments:
            # Immediately add genomes to their species assignment
            self.database.add(genome_filepath, reference=False, taxonomy=species_assignments[genome_filepath])

        # Cluster all the unassigned genomes together and define new clusters at different taxonomic levels
        characterized, unassigned = self.database.characterize()

        # Finally, index the new genomes
        self.database.update()

        if args.pack:
            # Pack the database into a compressed tarball
            self.pack(["--workdir", args.workdir, "--database", args.database])

        # Print credits
        self.__print_credits()


if __name__ == "__main__":
    # Initialize the MetaSBT object and run it
    framework = MetaSBT()

    # Print the running time and exit
    print(f"Total elapsed time: {time.time-framework.start_at}s")