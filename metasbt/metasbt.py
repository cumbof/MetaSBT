#!/usr/bin/env python3
"""A scalable framework for automatically indexing microbial genomes and accurately
characterizing metagenome-assembled genomes with Sequence Bloom Trees.
"""

import argparse
import errno
import gzip
import json
import os
import re
import requests
import shutil
import subprocess
import sys
import tempfile
import time
import tqdm
import unittest
import urllib.request

from datetime import datetime
from packaging import version
from pathlib import Path
from tabulate import tabulate
from typing import Any, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from metasbt import __date__, __version__
from metasbt.core import DEPENDENCIES, Database, Entry

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
    "sha256sum",  # Required by the pack subroutine
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
        export      Export a MetaSBT database taxonomy as a Newick phylogenetic tree;
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
            if dependency == "deltatree":
                # deltatree is a Python module bound to Rust, not a shell executable.
                try:
                    import deltatree
                except ImportError:
                    missing_dependencies.append(dependency)
            elif shutil.which(dependency) is None:
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

        if len(sys.argv) == 1:
            # Print usage in case of no arguments
            sys.argv.append("--help")

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
            versions = [db_version["version"] for db_version in databases[args.download]]

            # Select the latest version by default
            selected_version = str(sorted([int(v) for v in versions])[-1])

            if args.version:
                if args.version not in versions:
                    raise Exception("The provided database version does not exist!")

                selected_version = args.version

            # Retrieve the URL to the tarball
            database_url = [db_version["tarball"] for db_version in databases[args.download] if db_version["version"] == selected_version][0]

            database_filepath = os.path.join(args.folder, f"MetaSBT-{args.download}-{selected_version}.tar.gz")

            try:
                # Download the database tarball
                urllib.request.urlretrieve(database_url, database_filepath)

            except:
                raise Exception(f"Unable to retrieve the database tarball from {database_url}")

            # Compute the sha256 hash
            sha256 = subprocess.run(["sha256sum", database_filepath], capture_output=True, text=True)

            if sha256.returncode == 0:
                # Perform a sanity check on the tarball
                # Compare the sha256 hashes of the downloaded tarball with the hash in the databases table
                if sha256.stdout.strip().split()[0] != [db_version["sha256"] for db_version in databases[args.download] if db_version["version"] == selected_version][0]:
                    raise Exception("sha256 mismatch! Consider downloading the database again")

            else:
                raise Exception(f"An error has occurred while computing the sha256 hash of {database_filepath}")

    def export(self, argv: List[Any]) -> None:
        """Export a MetaSBT database taxonomy as a Newick phylogenetic tree.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        FileNotFoundError
            If the database directory does not exist.
        ValueError
            If the specified taxonomic level is not valid.
        """

        parser = argparse.ArgumentParser(
            prog="export",
            description="Export a MetaSBT database taxonomy as a Newick phylogenetic tree.",
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
            required=True,
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--level",
            default="species",
            choices=["phylum", "class", "order", "family", "genus", "species"],
            help="Taxonomic resolution of the tree leaves."
        )
        parser.add_argument(
            "--output",
            type=os.path.abspath,
            default=None,
            help=(
                "Path to the output Newick file.  "
                "Defaults to {workdir}/{database}-{YYYYMMDD}.nwk."
            )
        )

        args = parser.parse_args(argv)

        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        tmp_dir = os.path.join(args.workdir, "tmp")

        if self.database is None:
            self.database = Database(args.database, db_dir, tmp_dir)

        newick_str = self.database.to_newick(level=args.level)

        if args.output:
            output_path = args.output
        else:
            timestamp = datetime.now().strftime("%Y%m%d")
            output_path = os.path.join(args.workdir, f"{args.database}-{timestamp}.nwk")

        with open(output_path, "w") as nwk_file:
            nwk_file.write(newick_str + "\n")

        print(f"Newick tree written to: {output_path}")

    @staticmethod
    def _reset_for_resume(db_dir: str) -> None:
        """Strip a database folder down to its genome sketches, ready for a `--resume` rebuild.

        Everything under `db_dir` is removed except the `sketches` folder and the `dereplicated.tsv`
        cache, and the metadata is rewritten keeping only the two entries the sketches depend on (the
        kmer size and the scaled factor) plus a zeroed cluster counter. Anything else the interrupted
        run may have written (cluster folders, the report, the learned species radii) is derived
        state that the rebuild recomputes, and reusing a partially written copy of it would corrupt
        the new build. The dereplication cache is kept because it is an expensive, input-determined
        result: `dereplicate` reuses it only when the input set and threshold still match, and
        recomputes (overwriting it) otherwise, so keeping a stale copy is always safe.

        Parameters
        ----------
        db_dir : str
            Path to the database folder.

        Raises
        ------
        Exception
            If the folder has no sketches to resume from, or no metadata recording how they were
            built (without the kmer size and scaled factor the sketches cannot be safely reused).
        """

        sketches_dir = os.path.join(db_dir, "sketches")

        if not os.path.isdir(sketches_dir) or not os.listdir(sketches_dir):
            raise Exception(
                f"There are no genome sketches to resume from under {db_dir}. "
                "Remove the folder and run `index` without --resume."
            )

        metadata_filepath = os.path.join(db_dir, "metadata.json")

        if not os.path.isfile(metadata_filepath):
            raise Exception(
                f"No metadata.json under {db_dir}: there is no record of the kmer size and scaled "
                "factor the existing sketches were built with, so they cannot be safely reused."
            )

        with open(metadata_filepath) as metadata_file:
            metadata = json.load(metadata_file)

        if "kmer_size" not in metadata or "scaled_factor" not in metadata:
            raise Exception(
                f"The metadata under {db_dir} does not record both the kmer size and the scaled "
                "factor, so the existing sketches cannot be safely reused."
            )

        sketches_count = len(os.listdir(sketches_dir))

        for entry in os.listdir(db_dir):
            if entry in ("sketches", "metadata.json", "dereplicated.tsv"):
                continue

            entry_path = os.path.join(db_dir, entry)

            if os.path.isdir(entry_path) and not os.path.islink(entry_path):
                shutil.rmtree(entry_path)

            else:
                os.unlink(entry_path)

        with open(metadata_filepath, "w+") as metadata_file:
            json.dump(
                {
                    "kmer_size": metadata["kmer_size"],
                    "scaled_factor": metadata["scaled_factor"],
                    "clusters_count": 0,
                },
                metadata_file,
            )

        print(
            f"Resuming from {sketches_count} genome sketches under {db_dir} "
            f"(kmer size {metadata['kmer_size']}, scaled factor {metadata['scaled_factor']})"
        )

    def index(self, argv: List[Any]) -> None:
        """Build the first baseline of a MetaSBT database by indexing a set of reference genomes.

        Run against an already complete database (without --resume), `index` instead adds the input
        genomes as a brand-new kingdom, keeping the existing clusters and genomes untouched: the input
        is sketched, dereplicated and clustered on its own (kingdoms never cross-compare) and spliced
        in beside the existing kingdoms. It refuses input whose kingdom is already in the database
        (that is `update`'s job) and refuses a half-built database (which must be resumed).

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        Exception
            - If the database exists but is incomplete and --resume was not given;
            - If the input genomes belong to a kingdom already present in the database.
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
        general_group.add_argument(
            "--resume",
            action="store_true",
            default=False,
            help=(
                "Rebuild an existing database, reusing the genome sketches already under it. "
                "Everything else is rebuilt from scratch, and the kmer size and scaled factor are "
                "taken from the existing metadata so that the sketches stay valid."
            ),
        )

        general_group.add_argument(
            "--scaled-factor",
            type=int,
            default=1000,
            dest="scaled_factor",
            help="FracMinHash scaled factor (1 in N k-mers are sampled). Lower values yield larger sketches."
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

        db_exists = os.path.isdir(db_dir)

        # A database is "complete" only once `update()` has dumped clusters.tsv (see Database.__init__).
        # A complete database can take a brand-new kingdom incrementally (new-kingdom mode, below); an
        # incomplete one (a crashed `index`) can only be resumed.
        db_complete = db_exists and os.path.isfile(os.path.join(db_dir, "clusters.tsv"))

        if db_exists and not args.resume and not db_complete:
            # A half-built database: the only safe move is to resume it. Running afresh would either
            # collide with the partial state or, with new-kingdom mode, splice onto clusters that were
            # never finished.
            raise Exception(
                f"A database named '{args.database}' already exist under {args.workdir} but it is not "
                "complete. Use --resume to rebuild it reusing the genome sketches already under it."
            )

        if args.resume and db_exists:
            # Drop everything but the sketches so the rebuild starts from a clean slate. A crashed
            # `index` can leave half-populated clusters and a metadata file carrying counters and
            # learned radii from the interrupted run; keeping any of that would silently mix two
            # builds. The sketches are the one artifact that is safe to carry over: they depend on
            # nothing but the genome, the kmer size, and the scaled factor.
            self.__class__._reset_for_resume(db_dir)

        # New-kingdom mode: add genomes from a kingdom NOT yet in a complete database, keeping its
        # existing clusters and genomes untouched. The input genomes are sketched, dereplicated and
        # clustered on their own (kingdoms never cross-compare), then `add()`/`update()` splice the new
        # kingdom subtree in beside the existing ones. `--resume` (a full rebuild) takes precedence.
        # The kingdom-overlap guard is enforced below, once the input taxonomy has been read.
        incremental = db_complete and not args.resume

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        # Initialize the database
        self.database = Database(args.database, db_dir, tmp_dir, nproc=args.nproc)

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

        # The kingdoms the input genomes belong to (the first taxonomic level of every label)
        input_kingdoms = {Database._format_taxonomy(taxonomy).split("|")[0] for taxonomy in references.values()}

        if incremental:
            # `index` only ever adds a kingdom that is not yet in the database. Adding genomes to a
            # kingdom already present is `update`'s job (it would profile them against the existing
            # tree); a full rebuild is `--resume`. Refuse anything that overlaps an existing kingdom.
            existing_kingdoms = set(self.database.clusters["kingdom"].keys())
            overlap = existing_kingdoms.intersection(input_kingdoms)

            if overlap:
                raise Exception(
                    f"The database '{args.database}' already contains genomes from "
                    f"{', '.join(sorted(overlap))}. The `index` command only adds genomes from a "
                    "kingdom not yet in the database; use `update` to add genomes to an existing "
                    "kingdom, or --resume to rebuild the whole database from scratch."
                )

        if (args.resume or incremental) and Database._validate_metadata(self.database.metadata):
            # The existing sketches were built with the kmer size and scaled factor already recorded
            # in the metadata: keep them. Re-running set_configs could estimate a different kmer
            # size, and a sketch compared under the wrong kmer size yields a meaningless distance
            # without raising anything (the .bf file stays a structurally valid pair of bitmaps).
            if args.kmer_size and args.kmer_size != self.database.metadata["kmer_size"]:
                raise Exception(
                    f"--kmer-size {args.kmer_size} conflicts with the kmer size the existing sketches "
                    f"were built with ({self.database.metadata['kmer_size']}). Drop --resume to rebuild "
                    "the sketches from scratch."
                )

            # `--scaled-factor` always carries a value, so only an explicitly provided one (i.e. one
            # differing from the parser default) can contradict the stored configuration
            scaled_factor_given = args.scaled_factor != parser.get_default("scaled_factor")

            if scaled_factor_given and args.scaled_factor != self.database.metadata["scaled_factor"]:
                raise Exception(
                    f"--scaled-factor {args.scaled_factor} conflicts with the scaled factor the existing "
                    f"sketches were built with ({self.database.metadata['scaled_factor']}). Drop --resume "
                    "to rebuild the sketches from scratch."
                )

        else:
            # Define the database metadata
            # Eventually, estimate the optimal kmer size
            self.database.set_configs(
                genomes,
                kmer_size=args.kmer_size,
                kmer_max=args.limit_kmer_size,
                scaled_factor=args.scaled_factor
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
            # Dereplicate genomes based on their ANI distance. The references carry taxonomic labels,
            # so the comparison is partitioned by genus when the threshold allows it (see dereplicate).
            # The result is cached under the database folder (and preserved across --resume) so an
            # interrupted build that already dereplicated does not repeat it on the next run. A
            # new-kingdom run caches under a per-kingdom name so it never clobbers the baseline
            # kingdom's dereplication record.
            cache_name = "dereplicated.tsv"

            if incremental:
                kingdom_tag = "_".join(sorted(kingdom[3:] for kingdom in input_kingdoms))
                cache_name = f"dereplicated_{kingdom_tag}.tsv"

            genomes = self.database.dereplicate(
                genomes,
                threshold=args.dereplicate,
                taxonomy=references,
                cache_filepath=os.path.join(db_dir, cache_name),
            )

        # Reshape the references dict
        # Consider genomes that passed the dereplication process only
        references = {genome: references[genome] for genome in genomes}

        # Cluster the references into ANI-coherent species clusters without a fixed threshold
        # and relabel each cluster by majority vote of its members' taxonomy
        references = self.database.cluster_references(references)

        # Add references to the database under their refined lineage
        for genome in references:
            self.database.add(genome, reference=True, taxonomy=references[genome])

        # Finally, index references
        self.database.update()

        if args.pack:
            # Pack the database into a compressed tarball
            self.pack(argv, parse_known_args=True)

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
            "--kmer-size",
            type=int,
            default=27,
            dest="kmer_size",
            help="The kmer size in bp."
        )
        parser.add_argument(
            "--minimizer-length",
            type=int,
            default=21,
            dest="minimizer_length",
            help="The minimizer length in bp."
        )
        parser.add_argument(
            "--minimizer-spaces",
            type=int,
            default=5,
            dest="minimizer_spaces",
            help="Number of characters in minimizer that are ignored in comparisons."
        )
        parser.add_argument(
            "--threads",
            type=int,
            default=1,
            help="Number of threads for kraken2-build.",
        )

        # Load arguments
        args = parser.parse_args(argv)

        # Define the ordered list of taxonomic levels
        TAXONOMIC_RANKING = [
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        ]

        # Template for node entries
        NODE_TEMPLATE = "taxid\t|\tparent\t|\trank\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|"

        # Template for name entries
        NAME_TEMPLATE = "taxid\t|\tvar\t|\t\t|\tscientific name\t|"

        # Load the list of paths to the genome files
        # This is a dict with the paths to the input genomes indexed by their file name
        # Gzip-compressed fasta files (e.g. .fna.gz) are supported
        with open(args.genomes) as fh:
            genomes = {Database._basename(line.strip()): line.strip() for line in fh if line.strip()}

        # Keep track of the current working directory
        curr_workdir = os.getcwd()

        # Set --workdir as the current working directory
        os.chdir(args.workdir)

        # Initialize the kraken2 database
        # This is going to download the NCBI taxonomy into a `taxonomy` folder under the `args.database` directory
        command_line = [
            "kraken2-build",
            "--download-taxonomy",
            "--skip-maps",
            "--db",
            args.database
        ]

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            # Set the current working directory back to the original one
            os.chdir(curr_workdir)

            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message) from e

        # Define the path to the names.dmp and nodes.dmp files
        ncbi_names_filepath = os.path.join(args.workdir, args.database, "taxonomy", "names.dmp")

        ncbi_nodes_filepath = os.path.join(args.workdir, args.database, "taxonomy", "nodes.dmp")

        # Keep track of NCBI names
        ncbi_names = dict()

        # Keep track of the highgest NCBI TaxID
        latest_node_id = 0

        # Load NCBI names
        with open(ncbi_names_filepath) as ncbi_names_table:
            for line in ncbi_names_table:
                line = line.strip()

                if line:
                    line_split = line.split("|")

                    if line_split[3].strip().lower() in ["scientific name", "equivalent name", "synonym"]:
                        # NCBI TaxID
                        ncbi_taxid = int(line_split[0].strip())

                        # NCBI Tax scientific names, equivalent names, and synonyms
                        # Names must be fixed in the same way that MetaSBT does
                        ncbi_tax_name = re.sub(r"_+", "_", re.sub(r"\W+", "_", line_split[1].strip())).strip("_")

                        ncbi_names[ncbi_tax_name] = ncbi_taxid

                        if ncbi_taxid > latest_node_id:
                            latest_node_id = ncbi_taxid

        # Keep track of NCBI node IDs
        ncbi_nodes = dict()

        # Load NCBI tax relationships
        with open(ncbi_nodes_filepath) as ncbi_nodes_table:
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
        metasbt_clusters_filepath = os.path.join(args.workdir, args.database, "clusters.tsv")

        metasbt_names = dict()

        metasbt_nodes = dict()

        # Take track of the full taxonomic labels and genomes in MetaSBT
        metasbt_lineages = dict()

        # Load MetaSBT clusters
        with open(metasbt_clusters_filepath) as metasbt_report:
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

                                if cluster_name not in ncbi_names and cluster_name not in metasbt_names:
                                    # This could be one of the following clusters:
                                    # (1) MSBT cluster;
                                    # (2) known species clustered into multiple clades;
                                    # (3) known species whose name has been deprecated in the current version of names.dmp;
                                    # We consider these clusters as new species
                                    latest_node_id += 1

                                    metasbt_names[cluster_name] = latest_node_id

                                if cluster_name in metasbt_names and metasbt_names[cluster_name] not in metasbt_nodes:
                                    parent_id = None

                                    if parent_name in ncbi_names:
                                        parent_id = ncbi_names[parent_name]

                                    elif parent_name in metasbt_names:
                                        parent_id = metasbt_names[parent_name]

                                    else:
                                        # Set the current working directory back to the original one
                                        os.chdir(curr_workdir)

                                        raise ValueError('Node "{}" is missing!'.format(parent_name))

                                    metasbt_nodes[metasbt_names[cluster_name]] = {
                                        "parent": parent_id,
                                        "rank": TAXONOMIC_RANKING[rank_position]
                                    }

        # Finally, write the MetaSBT names and nodes additions
        with open(ncbi_names_filepath, "a+") as ncbi_names_table:
            for node_name in sorted(metasbt_names.keys(), key=lambda node_name: metasbt_names[node_name]):
                # Retrieve the node id
                node_id = metasbt_names[node_name]

                # Add the tax id
                metasbt_node = NAME_TEMPLATE.replace("taxid", str(node_id))

                # Add the node name
                metasbt_node = metasbt_node.replace("var", node_name)

                ncbi_names_table.write("{}\n".format(metasbt_node))

        with open(ncbi_nodes_filepath, "a+") as ncbi_nodes_table:
            for node_id, node_data in sorted(metasbt_nodes.items()):
                # Add the tax id
                metasbt_node = NODE_TEMPLATE.replace("taxid", str(node_id))

                # Add the parent node name
                metasbt_node = metasbt_node.replace("parent", str(node_data["parent"]))

                # Add the rank
                metasbt_node = metasbt_node.replace("rank", node_data["rank"])

                ncbi_nodes_table.write("{}\n".format(metasbt_node))

        for lineage in metasbt_lineages:
            print(f"Processing {lineage}")

            # Intersect the list of genomes under this lineage with the input genomes
            selected_genomes = set(genomes.keys()).intersection(set(metasbt_lineages[lineage]))

            if selected_genomes:
                # Lineages are always defined up to the species level here
                lineage = lineage.split("|")

                # Retrieve the list of genome files for the selection
                genome_files = [genomes[genome] for genome in selected_genomes]

                metasbt_cluster_name = lineage[-1][3:]

                if metasbt_cluster_name in metasbt_names:
                    metasbt_cluster_id = metasbt_names[metasbt_cluster_name]

                elif metasbt_cluster_name in ncbi_names:
                    metasbt_cluster_id = ncbi_names[metasbt_cluster_name]

                else:
                    # This should never happen
                    continue

                for genome_file in genome_files:
                    with tempfile.NamedTemporaryFile() as tmp_genome:
                        with open(tmp_genome.name, "wt"), open(genome_file) as genome:
                            for record in SeqIO.parse(genome, "fasta"):
                                # Rebuild the record with the new ID
                                record = SeqRecord(
                                    record.seq, 
                                    id="{}|kraken:taxid|{} {}".format(record.id, metasbt_cluster_id, record.id), 
                                    name=record.name, 
                                    description=record.description
                                )

                                SeqIO.write(record, tmp_genome.name, "fasta")

                        # Add the temporary genome file to the library
                        command_line = [
                            "kraken2-build",
                            "--add-to-library",
                            tmp_genome.name,
                            "--db",
                            args.database,
                            "--no-masking",
                            "--threads",
                            str(args.threads)
                        ]

                        try:
                            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                        except subprocess.CalledProcessError as e:
                            # Set the current working directory back to the original one
                            os.chdir(curr_workdir)

                            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                            raise Exception(error_message) from e

        # Finally, build the database
        command_line = [
            "kraken2-build",
            "--build",
            "--db",
            args.database,
            "--kmer-len",
            str(args.kmer_size),
            "--minimizer-len",
            str(args.minimizer_length),
            "--minimizer-spaces",
            str(args.minimizer_spaces),
            "--threads",
            str(args.threads)
        ]

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            # Set the current working directory back to the original one
            os.chdir(curr_workdir)

            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message) from e

        # Set the current working directory back to the original one
        os.chdir(curr_workdir)

    def pack(self, argv: List[Any], parse_known_args=False) -> None:
        """Build a compressed tarball with a MetaSBT database and report its sha256.

        Parameters
        ----------
        argv : list
            The list of arguments.
        parse_known_args : bool, default False
            Parse known arguments only without raising any exceptions for unrecognized arguments.

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
        if parse_known_args:
            args, _ = parser.parse_known_args(argv)

        else:
            args = parser.parse_args(argv)

        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Get the current time
        now = datetime.now()

        # Define a timestamp
        timestamp = f"{now.year}{str(now.month).zfill(2)}{str(now.day).zfill(2)}"

        # Search for other tarballs under the same working directory
        tarballs = list(Path(args.workdir).glob(f"MetaSBT-{args.database}-{timestamp}*.tar.gz"))

        # Define the output file name
        output_filename = f"MetaSBT-{args.database}-{timestamp}"

        if tarballs:
            # Add an incremental number in case of multiple versions of the same database
            output_filename += f"-{len(tarballs)+1}"

        # Define the path to the output tarball
        output_filepath = os.path.join(args.workdir, f"{output_filename}.tar.gz")

        # Build the compressed tarball
        # Use args.database (relative) as the source so that -C args.workdir produces
        # a tarball whose root entry is just the database folder name, not an absolute path.
        command_line = [
            "tar",
            "-czvf",
            output_filepath,
            "-C",
            args.workdir,
            args.database
        ]

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message) from e

        # Compute the sha256 hash
        sha256 = subprocess.run(["sha256sum", output_filepath], capture_output=True, text=True)

        if sha256.returncode == 0:
            print(sha256.stdout)

        else:
            raise Exception(f"An error has occurred while computing the sha256 hash of {output_filepath}")

    def profile(self, argv: List[Any], parse_known_args=False, print_summary=True) -> None:
        """Profile a set of genomes against a specific MetaSBT database.
        Profile tables are stored under a dedicated folder in the workdir temporary directory.

        Parameters
        ----------
        argv : list
            The list of arguments.
        parse_known_args : bool, default False
            Parse known arguments only without raising any exceptions for unrecognized arguments.
        print_summary : bool, default True
            Print the per-genome human-readable profile table. Disabled when `update` calls
            this for its side effects (sketching + profiling) so it does not flood the log.

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
            "--uncertainty",
            required=False,
            type=float,
            default=20.0,
            help=(
                "Percentage by which the beam is expanded around the closest hit when keeping "
                "multiple best hits: a sibling is kept within best * (1 + uncertainty/100). On an "
                "exact match (best distance 0) the expansion is anchored to the nearest non-zero "
                "competitor instead, so close alternatives are still admitted."
            )
        )
        parser.add_argument(
            "--pruning-threshold",
            dest="pruning_threshold",
            required=False,
            type=float,
            default=0.0,
            help="Threshold for pruning the Sequence Bloom Tree."
        )
        parser.add_argument(
            "--nproc",
            required=False,
            type=int,
            default=os.cpu_count(),
            help="Process the input genomes in parallel."
        )

        # Load arguments
        if parse_known_args:
            args, _ = parser.parse_known_args(argv)

        else:
            args = parser.parse_args(argv)
        
        # Define the path to the database folder
        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        if self.database is None:
            # Initialize the database
                self.database = Database(args.database, db_dir, tmp_dir, nproc=args.nproc)

        # Load the list of paths to the input genomes
        if args.genome:
            genomes = [args.genome]
        else:
            with open(args.genomes) as fh:
                genomes = [line.strip() for line in fh if line.strip()]

        # Genomes must be sketched first
        # Note that this function and `sketch()` have the same set of arguments
        # Sketches are stored under the dedicated folder in the database directory
        sketches = self.sketch(argv, parse_known_args=True)

        # Profile genomes in parallel
        results = self.database.profile_genomes(genomes, sketches, uncertainty=args.uncertainty, pruning_threshold=args.pruning_threshold, mode="split")

        if not print_summary:
            # Called by `update` for its side effects only; skip the per-genome dump.
            return

        # Print human-readable summary with confidence scores
        col_w = 70
        header = f"{'Level':<12}  {'Closest':<{col_w}}  {'ANI':>8}  {'Confidence':>10}"
        sep = "-" * len(header)

        for genome_filepath in genomes:
            profile = results.get(genome_filepath, {})
            genome_name = Database._basename(genome_filepath)
            print(f"\n## {genome_name}")
            print(header)
            print(sep)

            confidences = profile.get("confidence", {})

            for level in Database.LEVELS + ["genome"]:
                if level not in profile or not profile[level]:
                    continue
                label, ani = min(profile[level].items(), key=lambda x: x[1])
                confidence = confidences.get(level)
                conf_str = f"{confidence:.3f}" if confidence is not None else "N/A"
                warning = "  [LOW CONFIDENCE]" if confidence is not None and confidence < 0.5 else ""
                print(f"{level:<12}  {label:<{col_w}}  {ani:>8.4f}  {conf_str:>10}{warning}")

    def sketch(self, argv: List[Any], parse_known_args=False) -> List[str]:
        """Sketch the input genomes.

        Parameters
        ----------
        argv : list
            The list of arguments.
        parse_known_args : bool, default False
            Parse known arguments only without raising any exceptions for unrecognized arguments.

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
        if parse_known_args:
            args, _ = parser.parse_known_args(argv)

        else:
            args = parser.parse_args(argv)

        # Define the path to the database folder
        db_dir = os.path.join(args.workdir, args.database)

        if not os.path.isdir(db_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), db_dir)

        # Define the path to the temporary folder
        tmp_dir = os.path.join(args.workdir, "tmp")

        if self.database is None:
            # Initialize the database
                self.database = Database(args.database, db_dir, tmp_dir, nproc=args.nproc)

        # Load the list of paths to the input genomes
        if args.genome:
            genomes = [args.genome]
        else:
            with open(args.genomes) as fh:
                genomes = [line.strip() for line in fh if line.strip()]

        sketch_map = self.database.sketch_genomes(genomes)

        # Preserve original order
        sketches = [sketch_map[g] for g in genomes]

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
            metadata = json.load(metadata_json_file)

            # Add metadata to the summary table
            table.extend([[key, metadata[key]] for key in metadata])

        # Count the number of clusters at all the seven taxonomic levels
        clusters = {"kingdom": 0, "phylum": 0, "class": 0, "order": 0, "family": 0, "genus": 0, "species": 0}

        # Also count the number of known clusters at all the seven taxonomic levels
        knowns = {level: 0 for level in clusters}

        # Count the number of reference genomes and MAGs
        genomes = {"references": 0, "mags": 0}

        # Keep track of the root cardinality (from the kingdom cluster)
        root_cardinality = 0

        with open(search_for[2]) as clusters_table:
            header = list()

            for line in clusters_table:
                if line.strip():
                    if line.startswith("#"):
                        # The last commented line is the header
                        header = line[1:].strip().split("\t")

                    else:
                        line_split = line.split("\t")

                        # Strip the new line character out of the last element 
                        line_split[-1] = line_split[-1].strip()

                        # Increment the clusters counter based on the taxonomic level
                        clusters[line_split[header.index("level")]] += 1

                        if line_split[header.index("known")].strip() == "True":
                            # Increment the known clusters counter based on the taxonomic level
                            knowns[line_split[header.index("level")]] += 1 

                        if line_split[header.index("level")] == "kingdom":
                            # Update the root cardinality
                            root_cardinality = int(line_split[header.index("cardinality")])

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

        # Add the root cardinality
        table.append(["cardinality", root_cardinality])

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
            "--references",
            required=True,
            type=os.path.abspath,
            help="Path to the file with the list of paths to the reference genomes and their taxonomies."
        )
        parser.add_argument(
            "--mags",
            required=True,
            type=os.path.abspath,
            help="Path to the file with the list of paths to the metagenome-assembled genomes."
        )

        # Load arguments
        args = parser.parse_args(argv)

        class Test(unittest.TestCase):
            """Unit tests.
            """

            # These are used in the setUpClass function in order to load the set of reference genomes and MAGs
            references_filepath = None

            mags_filepath = None

            @classmethod
            def setUpClass(cls):
                """Set up the unit class by retrieving a bunch of genomes to use within all the unit tests.
                """

                # Create the temporary working directory
                cls.working_dir = tempfile.TemporaryDirectory()

                # Define the database name
                cls.db_name = "Test"

                # Define the database folder
                cls.db_dir = os.path.join(cls.working_dir.name, cls.db_name)

                # Define the temporary folder
                cls.tmp_dir = os.path.join(cls.working_dir.name, "tmp")

                os.makedirs(cls.tmp_dir)

                # Define the path to the folder with genomes
                cls.genomes_dir = os.path.join(cls.working_dir.name, "genomes")

                os.makedirs(cls.genomes_dir)

                genomes = set()

                if cls.references_filepath:
                    genomes.update({line.strip().split("\t")[0] for line in open(cls.references_filepath).readlines() if line.strip() and not line.startswith("#")})

                if cls.mags_filepath:
                    genomes.update({line.strip() for line in open(cls.mags_filepath).readlines() if line.strip() and not line.startswith("#")})

                for genome_url in genomes:
                    try:
                        # Define the path to the local fna.gz file
                        fna_gz_path = os.path.join(cls.genomes_dir, os.path.basename(genome_url))

                        # Download the genome into the temporary working directory
                        urllib.request.urlretrieve(genome_url, fna_gz_path)

                        # Unzip the genome
                        with gzip.open(fna_gz_path, "rb") as fna_gz_in:
                            with open(os.path.splitext(fna_gz_path)[0], "wb") as fna_out:
                                shutil.copyfileobj(fna_gz_in, fna_out)

                    except urllib.error.HTTPError as e:
                        error_message = f"Unable to retrieve {genome_url}\n\n"

                        raise Exception(error_message) from e

            @classmethod
            def tearDownClass(cls):
                """ Get rid of the temporary working directory with the test genomes.
                """

                # Delete the temporary working directory
                shutil.rmtree(cls.working_dir.name, ignore_errors=False)

            def test_pipeline(self):
                """Test the MetaSBT pipeline.
                """

                # Define the path to the file with the list of reference genomes and their taxonomic labels
                references_filepath = os.path.join(self.__class__.working_dir.name, "references.tsv")

                with open(self.__class__.references_filepath) as references_file1, open(references_filepath, "w+") as references_file2:
                    for line in references_file1:
                        line = line.strip()

                        if line:
                            if not line.startswith("#"):
                                line_split = line.split("\t")

                                # Genomes are .fna files
                                reference_genome_filepath = os.path.join(self.__class__.genomes_dir, os.path.splitext(os.path.basename(line_split[0]))[0])

                                # Define the new references file
                                references_file2.write(f"{reference_genome_filepath}\t{line_split[1]}\n")

                # Build the baseline with the `index` command
                # Define the index command line
                # Test genomes are viruses, so we can set --kmer-size a priori here
                sys.argv = [
                    "metasbt",
                    "index",
                    "--workdir",
                    self.__class__.working_dir.name,
                    "--database",
                    self.__class__.db_name,
                    "--references",
                    references_filepath,
                    "--kmer-size",
                    "9",
                    "--nproc",
                    "1"
                ]

                # Run the command line and index the reference genomes
                MetaSBT()

                with self.subTest():
                    # Define the path to clusters.tsv and genomes.tsv
                    clusters_filepath = os.path.join(self.__class__.db_dir, "clusters.tsv")

                    genomes_filepath = os.path.join(self.__class__.db_dir, "genomes.tsv")

                    # clusters.tsv and genomes.tsv should exist in the database folder
                    self.assertTrue(os.path.isfile(clusters_filepath) and os.path.isfile(genomes_filepath))

                # Define the path to the file with the list of metagenome-assembled genomes
                mags_filepath = os.path.join(self.__class__.working_dir.name, "mags.txt")

                with open(self.__class__.mags_filepath) as mags_file1, open(mags_filepath, "w+") as mags_file2:
                    for line in mags_file1:
                        line = line.strip()

                        if line:
                            if not line.startswith("#"):
                                # Genomes are .fna files
                                mags_genome_filepath = os.path.join(self.__class__.genomes_dir, os.path.splitext(os.path.basename(line))[0])

                                # Define the new mags file
                                mags_file2.write(f"{mags_genome_filepath}\n")

                # Update the database with the `update` command
                # This is going to pack the database into a compressed tarball
                # Define the update command line
                sys.argv = [
                    "metasbt", 
                    "update", 
                    "--workdir", 
                    self.__class__.working_dir.name, 
                    "--database", 
                    self.__class__.db_name, 
                    "--genomes", 
                    mags_filepath, 
                    "--pack",
                    "--nproc",
                    "1"
                ]

                # Run the command line and update the database
                MetaSBT()

                with self.subTest():
                    # Get the current time
                    now = datetime.now()

                    # Define a timestamp
                    timestamp = f"{now.year}{str(now.month).zfill(2)}{str(now.day).zfill(2)}"

                    # Define the path to the compressed tarball
                    tarball_filepath = os.path.join(self.__class__.working_dir.name, f"MetaSBT-{self.__class__.db_name}-{timestamp}.tar.gz")

                    # The tarball should exist here
                    self.assertTrue(os.path.isfile(tarball_filepath))

        def run_test(references: str, mags: str) -> None:
            """Run unit tests.

            Parameters
            ----------
            references : str
                Path to the file with the list of reference genomes and their taxonomic labels.
            mags : str
                Path to the file with the list of metagenome-assembled genomes.
            """

            # Manually set the path to the files with the list of reference genomes and MAGs 
            # Warning: an active internet connection is required to retrieve genomes
            Test.references_filepath = references

            Test.mags_filepath = mags

            # Initialize the test loader
            loader = unittest.TestLoader()

            # Load all the defined tests in Test
            suite = loader.loadTestsFromTestCase(Test)

            # Increment the verbosity level and stop running tests if a test fails
            runner = unittest.TextTestRunner(verbosity=2, failfast=True)

            # Finally, run the unit tests
            runner.run(suite)

        run_test(args.references, args.mags)

    def unpack(self, argv: List[Any]) -> None:
        """Unpack a local MetaSBT tarball database.

        Parameters
        ----------
        argv : list
            The list of arguments.

        Raises
        ------
        Exception
            - If it is unable to retrieve the database name from the tarball;
            - If a database with the same of the tarball one already exists;
            - In case of an error occurred while unpacking the tarball.
        FileNotFoundError
            If the input tarball does not exist.
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
            "--database",
            default="MetaSBT",
            type=str,
            help="The database name."
        )
        parser.add_argument(
            "--tarball",
            type=os.path.abspath,
            help="Path to the MetaSBT tarball database."
        )

        # Load arguments
        args = parser.parse_args(argv)

        # Create the working directory if it does not exist
        os.makedirs(args.workdir, exist_ok=True)

        if args.database:
            # --database is an optional argument
            db_dir = os.path.join(args.workdir, args.database)

            if os.path.isdir(db_dir):
                raise Exception(f"A database with the same name already exists: {args.database}")

        if not os.path.isfile(args.tarball):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.tarball)

        # This reads the first entry
        # The first entry in a MetaSBT database tarball is always the database root folder
        proc = subprocess.run(["tar", "-tzf", args.tarball], capture_output=True, text=True, check=True)
        db_dir_path_in_tarball = proc.stdout.strip().split("\n")[0]

        if db_dir_path_in_tarball.endswith(os.sep):
            # Trim the last char out
            db_dir_path_in_tarball = db_dir_path_in_tarball[:-1]

        # Keeps track of the root folder (the original database name)
        tarball_database = db_dir_path_in_tarball.split(os.sep)[-1]

        if not tarball_database.strip():
            raise Exception(f"Unable to retrieve the database name from the tarball")

        tarball_db_dir = os.path.join(args.workdir, tarball_database)

        if os.path.isdir(tarball_db_dir):
            raise Exception(f"Cannot extract the database in {args.workdir}. A database with the same name already exists: {tarball_database}")

        # Build the compressed tarball
        command_line = [
            "tar",
            "-xzvf",
            args.tarball,
            "-C",
            args.workdir
        ]

        # Count how many folder levels are before the database folder in the tarball
        dir_levels = len(db_dir_path_in_tarball.split(os.sep))-1

        if dir_levels > 0:
            # Add --strip-components to the command line
            # This must be specified before -C
            command_line.insert(3, f"--strip-components={dir_levels}")

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message) from e

        if args.database:
            # Rename the extracted database folder to args.database
            os.rename(tarball_db_dir, db_dir)

            tarball_db_dir = db_dir

        # Fix absolute paths to sketches
        # Search for all the txt files under the clusters folder
        # txt files contain the list of paths to the genome sketches that are used to build the SBTs
        sketch_lists = Path(os.path.join(db_dir, "clusters")).glob("**/*.txt")

        for sketch_list_filepath in sketch_lists:
            with open(sketch_list_filepath) as fh:
                sketch_filepaths = [line.strip() for line in fh if line.strip()]

            # Retrieve the partial path to the sketches starting from the clusters or sketches folder
            # Search for the last occurrence of "clusters" or "sketches" if it appears in multiple positions
            sketches_or_clusters_dir_index = max([pos for pos, v in enumerate(sketch_filepaths[0].split(os.sep)) if v == "clusters" or v == "sketches"])

            # Rebase the sketch filepaths to the new path on the local filesystem
            rebased_sketch_filepaths = [os.path.join(db_dir, os.sep.join(line.strip().split(os.sep)[sketches_or_clusters_dir_index:])) for line in sketch_filepaths]

            with open(sketch_list_filepath, "w+") as sketch_list_file:
                for sketch_filepath in rebased_sketch_filepaths:
                    sketch_list_file.write(f"{sketch_filepath}\n")

            # Note: index.delta files (Delta-SBT) store only binary Roaring Bitmap data
            # and contain no embedded absolute paths, so no rebasing is needed for them.

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
            required=True,
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
            help=(
                "Path to the file with the list of input genomes. "
                "With '--type mags' it contains one column with the paths to the genome files. "
                "With '--type references' it is a two-columns tab-separated-values file, the first "
                "column with the paths to the genome files and the second with their fully defined "
                "taxonomic label (k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species)."
            )
        )
        parser.add_argument(
            "--type",
            required=False,
            type=str.lower,
            default="mags",
            choices=["mags", "references"],
            dest="genomes_type",
            help=(
                "The type of the input genomes. "
                "'mags' are metagenome-assembled genomes with no known taxonomy: they are profiled "
                "against the database and either assigned to the closest species or clustered into new "
                "MSBT clusters. 'references' come with a fully defined taxonomic label and are added "
                "under their own lineage exactly like with the 'index' command."
            )
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
        parser.add_argument(
            "--uncertainty",
            required=False,
            type=float,
            default=20.0,
            help=(
                "Percentage by which the beam is expanded around the closest hit when profiling "
                "input genomes: a sibling is kept within best * (1 + uncertainty/100). On an exact "
                "match (best distance 0) the expansion is anchored to the nearest non-zero competitor "
                "instead, so close alternatives are still admitted."
            )
        )
        parser.add_argument(
            "--pruning-threshold",
            dest="pruning_threshold",
            required=False,
            type=float,
            default=0.0,
            help="Threshold for pruning the Sequence Bloom Tree while profiling input genomes."
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
        self.database = Database(args.database, db_dir, tmp_dir, nproc=args.nproc)

        if args.genomes_type == "references":
            # Reference genomes come with a fully defined taxonomic label, so they must be provided
            # through the two-columns `--genomes` file (path + taxonomy); a single `--genome` would
            # carry no taxonomy
            if not args.genomes:
                raise ValueError("Reference genomes must be provided through the two-columns --genomes file!")

            # Load the list of reference genomes and their taxonomic labels
            references = dict()

            with open(args.genomes) as input_table:
                for line in input_table:
                    line = line.strip()

                    if line and not line.startswith("#"):
                        line_split = line.split("\t")

                        if len(line_split) < 2:
                            raise ValueError(
                                "The --genomes file must contain two columns (path and taxonomy) with --type references!"
                            )

                        # The path to the genome file is under the first column, the taxonomy under the second
                        references[os.path.abspath(line_split[0])] = line_split[1]

            genomes = set(references.keys())

            if args.completeness > 0.0 or args.contamination < 100.0:
                # Retrieve the kingdom from the taxonomic label of the input references
                kingdom = references[list(references.keys())[0]].split("|")[0][3:]

                quality = Database.qc(genomes, kingdom, nproc=args.nproc, tmp=tmp_dir)

                genomes = {genome for genome in genomes if quality[genome]["completeness"] >= args.completeness and quality[genome]["contamination"] <= args.contamination}

            if args.dereplicate > 0.0:
                # Dereplicate references based on their ANI distance (input-vs-input), partitioned by
                # genus when the threshold allows it since the references carry taxonomic labels
                genomes = self.database.dereplicate(genomes, threshold=args.dereplicate, taxonomy=references)

            # Cluster the surviving references into ANI-coherent species clusters (reusing the
            # species radius learned at index time) and relabel each cluster by majority vote;
            # a cluster whose label matches an existing species merges into it through add()
            references = self.database.cluster_references({genome: references[genome] for genome in genomes})

            # Add the references under their refined lineage
            for genome in references:
                self.database.add(genome, reference=True, taxonomy=references[genome])

        else:
            # Load the set of paths to the input MAGs
            if args.genome:
                genomes = {args.genome}
            else:
                with open(args.genomes) as fh:
                    genomes = {line.strip() for line in fh if line.strip() and not line.startswith("#")}

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
            self.profile(argv, parse_known_args=True, print_summary=False)

            if args.dereplicate > 0.0:
                # Dereplicate the input genomes again versus the genomes in the database
                genomes = self.database.dereplicate(genomes, threshold=args.dereplicate, compare_with="database")

            species_assignments = self.database.is_known_batch(genomes)

            for genome_filepath in species_assignments:
                # Immediately add genomes to their species assignment
                self.database.add(genome_filepath, reference=False, taxonomy=species_assignments[genome_filepath])

            try:
                # Cluster all the unassigned genomes together and define new clusters at different taxonomic levels
                characterized, unassigned = self.database.characterize()

            except Exception as e:
                # Only ignore the exception if it is the expected "no unknown genomes" state
                if "There are no unknown genomes" not in str(e):
                    raise e

        # Finally, index the new genomes
        self.database.update()

        if args.pack:
            # Pack the database into a compressed tarball
            self.pack(argv, parse_known_args=True)

        # Print credits
        self.__print_credits()


def run() -> None:
    """Run MetaSBT.
    """

    # Initialize the MetaSBT object and run it
    framework = MetaSBT()

    # Print the running time and exit
    print(f"Total elapsed time: {time.time()-framework.start_at}s")


if __name__ == "__main__":
    run()
