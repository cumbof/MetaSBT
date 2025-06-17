"""Object-Oriented implementation of MetaSBT with Database and Entry abstractions.
"""

import copy
import datetime
import errno
import json
import math
import multiprocessing as mp
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

from collections import Counter, OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import fastcluster
import scipy.cluster.hierarchy as hier
import tqdm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metasbt import __date__, __version__

# Define the list of external software dependencies
DEPENDENCIES = [
    "busco",
    "checkm",
    "checkv",
    "howdesbt",  # This must be installed with the `Makefile_full` configuration
    "kitsune",
    "ntcard"
]


class Database(object):
    """Database object."""

    # Define the list of taxonomic levels
    LEVELS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    def __init__(
        self,
        name: str,
        folder: os.path.abspath,
        tmp: os.path.abspath,
        flat: bool=True,
        nproc: int=os.cpu_count(),
    ) -> "Database":
        """Initialize a Database object.

        Parameters
        ----------
        name : str
            The database id.
        folder : os.path.abspath
            Path to the database folder.
        tmp : os.path.abspath
            Path to the temporary folder.
        flat : bool, default True
            Skip the clustering of bloom filters and the definition of the Sequence Bloom Trees.
            A tree definition file is manually defined with all the leaves as children to the same root node.
            `flat=False` is advantageous if using a specific threshold for pruning the trees while profiling genomes.
        nproc : int, default os.cpu_count()
            Used to compute in parallel when possible.
            It automatically uses all the available CPU cores if its provided value is <1 or >os.cpu_count().

        Raises
        ------
        Exception
            - In case it is unable to retrieve the database report;
            - If a sketch does not exist under the sketches folder.
        FileNotFoundError
            - If a database already exists under `folder` but it does not contain a "clusters" folder;
            - If a database already exists under `folder` but it does not contain a "sketches" folder;
            - If a database already exists under `folder` but it does not have any metadata;
            - If a cluster folder does not exist in the database.
        ValueError
            - In case `name` is None;
            - In case `folder` is None.

        Returns
        -------
        Database
            A Database object.
        """

        # Register version and date
        self.version = f"{__date__} ({__version__})"

        if not name:
            raise ValueError("Invalid Database name!")

        if not folder:
            raise ValueError("Invalid Database path!")

        # Initialize a new database if `self.root` does not exist
        self.name = name

        # Register folder
        self.root = folder

        # Create the temporary folder
        os.makedirs(tmp, exist_ok=True)

        # Register tmp
        self.tmp = tmp

        # Register nproc
        self.nproc = nproc

        if self.nproc < 1 or self.nproc > os.cpu_count():
            self.nproc = os.cpu_count()

        # Define a list to keep track of clusters that have been created or modified during `self.add()` and `self.characterize()`
        self.__clusters: List[str] = list()

        # Define a list to keep track of the genomes that do not match with any species clusters in the database
        self.__unknowns: List["Entry"] = list()

        # Define a dictionary to keep track of the clusters in the database indexed by their id.
        # The hierarchical organization of clusters is maintained through clusters attributes `parent` and `children`
        # e.g.: {"kingdom": {"k__Viruses": <Entry>}, "phylum": {...}, ..., "species": {...}}
        self.clusters: Dict[str, Dict[str, "Entry"]] = {level: dict() for level in self.__class__.LEVELS}

        # Also define a dictionary to keep track of the genomes in the database indexed by their id.
        self.genomes: Dict[str, "Entry"] = dict()

        if not os.path.isdir(self.root):
            self.flat = flat

            os.makedirs(self.root)

            # Create a folder for keeping track of clusters
            os.makedirs(os.path.join(self.root, "clusters"))

            # Also create a folder for keeping track of genome sketches
            os.makedirs(os.path.join(self.root, "sketches"))

            # Init database metadata
            # This is supposed to keep track of k-mer size, bloom filter size, minimum occurrence of k-mers, number of clusters
            self.metadata: Dict[str, Any] = dict()

            # There are no clusters yet
            self.metadata["clusters_count"] = 0

            self.metadata["flat"] = self.flat

            # Init the database report with the list of clusters, theirs stats, and boundaries
            self.report: Dict[str, Any] = dict()

        else:
            # A database already exists under `folder`
            # There should be a "clusters" folder under the database
            clusters_folder = os.path.join(self.root, "clusters")

            if not os.path.isdir(clusters_folder):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), clusters_folder)

            # There should also be a "sketches" folder under the database
            sketches_folder = os.path.join(self.root, "sketches")

            if not os.path.isdir(sketches_folder):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), sketches_folder)

            # There should be a `metadata.json` file with info about the database
            metadata_json_filepath = os.path.join(self.root, "metadata.json")

            if not os.path.isfile(metadata_json_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), metadata_json_filepath)

            # Load metadata
            with open(metadata_json_filepath) as metadata_json_file:
                self.metadata = json.loads("".join(metadata_json_file.readlines()))

            if not self.__class__._validate_metadata(self.metadata):
                raise Exception("Database metadata did not pass the validation!")

            self.flat = self.metadata["flat"]

            # There should also be a report file under the database folder
            report_filepath = os.path.join(self.root, "clusters.tsv")

            if not os.path.isfile(report_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), report_filepath)

            # Load the report table with the list of clusters in the database
            self.report = self.__class__._load_report(report_filepath)

            if not self.report:
                raise Exception("Unable to retrieve the database report!")

            # Iterate over the species clusters only
            for cluster_id in self.report:
                # Retrieve the assigned taxonomic label
                taxonomy = self.report[cluster_id]["taxonomy"]

                # Get the taxonomic level
                cluster_level = self.report[cluster_id]["level"]

                # Get the last taxonomic level as cluster name
                cluster_name = taxonomy.split("|")[-1]

                # Retrieve the batch number from the cluster id
                batch_number = self.__class__._get_cluster_batch(cluster_id)

                # Define the cluster folder in the database
                cluster_folder = os.path.join(clusters_folder, batch_number, cluster_id)

                # Define the name of the parent cluster from the taxonomic label
                parent = None if cluster_level == "kingdom" else taxonomy.split("|")[-2]

                # Retrieve the set of reference genomes and mags
                cluster_references = self.report[cluster_id]["references"]

                cluster_mags = self.report[cluster_id]["mags"]

                # Define the set of entry children names
                children = cluster_references.union(cluster_mags)

                # Define the Entry object representing the cluster
                cluster_obj = Entry(self, cluster_id, cluster_name, cluster_level, folder=cluster_folder, parent=parent, children=children)

                if not cluster_obj.sketch_filepath:
                    # These clusters should have been sketched already
                    raise Exception(f"No sketch found for cluster {cluster_id}")

                # Register the cluster object
                self.clusters[cluster_level][cluster_name] = cluster_obj

                # Register references and mags in case of species clusters only
                if cluster_level == "species":
                    for genome_id in children:
                        # Define the entry representation of genomes
                        genome_obj = Entry(self, genome_id, genome_id, "genome", parent=cluster_name, taxonomy=taxonomy, known=genome_id in cluster_references)

                        if not genome_obj.sketch_filepath:
                            # These genomes should have been sketched already
                            raise Exception(f"No sketch found for genome {genome_id}")

                        # Register the genome object
                        self.genomes[genome_id] = genome_obj

    @staticmethod
    def _get_cluster_batch(cluster_id: str) -> str:
        """Given a MSBT cluster id, return the cluster batch number.

        Parameters
        ----------
        cluster_id : str
            The MSBT cluster id.

        Raises
        ------
        ValueError
            If `cluster_id` does not match with the regular expression ^MSBT[0-9]*$.

        Returns
        -------
        str
            The cluster batch number.
        """

        if not re.search("^MSBT[0-9]*$", cluster_id):
            raise ValueError("The provided argument is not a valid MSBT cluster id!")

        # Retrieve the batch number from the cluster id
        # Batches contain 1,000 clusters at most
        numerical_id = int(cluster_id[4:])

        _, batch_id = math.modf(numerical_id/1000)

        # The batch id is mapped over a mask of six 0s
        batch_mask = str(int(batch_id)).zfill(6)

        return batch_mask

    def __str__(self) -> None:
        """Print the Database object properties.

        Returns
        -------
        str
            A description of the Database object.
        """

        # Retrieve the number of clusters per taxonomic level
        kingdoms, phyla, classes, orders, families, genera, species = self.size()

        # Retrieve the total number of reference genomes
        references = sum([1 for genome in self.genomes if self.genomes[genome].is_known()])

        return f"""
            Class:                        metasbt.objects.Database
            Version:                      {self.version}
            Name:                         {self.name}
            Kingdoms:                     {kingdoms}
            Phyla:                        {phyla}
            Classes:                      {classes}
            Orders:                       {orders}
            Families:                     {families}
            Genera:                       {genera}
            Species:                      {species}
            Reference Genomes:            {references}
            Metagenome-Assembled Genomes: {len(self.genomes)-references}
        """

    def __contains__(self, genome: str) -> bool:
        """Check whether a genome is in the database by its name.

        Parameters
        ----------
        genome : str
            A genome name.

        Returns
        -------
        bool
            True if `genome` is in the database, False otherwise.
        """

        return genome in self.genomes

    def taxonomy_exists(self, taxonomy: str) -> bool:
        """Check whether a taxonomic label exists in the database.

        Parameters
        ----------
        taxonomy : str
            The taxonomic label.

        Returns
        -------
        bool
            True if the full `taxonomy` is defined in the database, False otherwise.
        """

        if not self.__class__._is_defined(taxonomy):
            return False

        for level, level_id in zip(self.__class__.LEVELS, taxonomy.split("|")):
            if level not in self.clusters:
                return False

            if level_id not in self.clusters[level]:
                return False

        return True

    def size(self) -> Tuple[int, int, int, int, int, int, int]:
        """Return the size of the database in terms of number of kingdoms,
        phyla, classes, orders, families, genera, and species clusters.

        Returns
        -------
        tuple
            A tuple with the number of clusters for each of the seven taxonomic levels.
        """

        return tuple([len(self.clusters[level]) for level in self.__class__.LEVELS])

    def set_configs(
        self,
        filepaths: List[os.path.abspath],
        min_kmer_occurrence: int=None,
        kmer_size: int=None,
        kmer_max: int=None,
        filter_size: int=None,
        filter_expand_by: float=None,
    ) -> None:
        """Set global configuration parameters.
        Search for the best configuration if no input is provided.

        Parameters
        ----------
        filepaths : list
            List of paths to the input uncompressed genome files.
        min_kmer_occurrence : int, default None
            Minimum number of kmer occurrences for establishing a proper bloom filter size.
        kmer_size : int, default None
            Kmers size.
        kmer_max : int, default None
            Max kmer size for kitsune, to be used if `kmer_size=None`.
        filter_size : int, default None
            Bloom filter size.
        filter_expand_by : float, default None
            Expand the bloom filter size by a predefined percentage.
            This is used in case of `filter_size=None` only.
            If `filter_size` is not None, this is ignored.

        Raises
        ------
        Exception
            If a metadata.json file already exists for the current database.
        ValueError:
            - if `filepaths` is empty;
            - if `min_kmer_occurrence` is None;
            - if both `kmer_size` and `kmer_max` are None.
        """

        metadata_json_filepath = os.path.join(self.root, "metadata.json")

        if os.path.isfile(metadata_json_filepath):
            raise Exception("A metadata.json file is already defined for this database!")

        if not filepaths:
            raise ValueError("One or more genome files must be provided in input!")

        if not min_kmer_occurrence:
            raise ValueError("The minimum number of kmer occurrences must be provided in input!")

        self.metadata["min_kmer_occurrence"] = min_kmer_occurrence

        input_list_filepath = os.path.join(self.tmp, "genomes.txt")

        if (not kmer_size and kmer_max) or not filter_size:
            # Dump the list of paths to the input genomes
            with open(input_list_filepath, "w+") as input_list:
                for filepath in filepaths:
                    input_list.write(f"{filepath}\n")

        if not kmer_size:
            if not kmer_max:
                raise ValueError("Max kmer size must be provided!")

            # Kitsune output file path
            kitsune_out_filepath = os.path.join(self.tmp, "kitsune.txt")

            try:
                # Search for the best kmer size based on the input set of genomes with kitsune
                command_line = [
                    "kitsune",
                    "kopt",
                    "--filenames",
                    input_list_filepath,
                    "--k-max",
                    str(kmer_max),
                    "--canonical",
                    "--fast",
                    "--nproc",
                    str(self.nproc),
                    "--threads",
                    "1",
                    "--in-memory",
                    "--output",
                    kitsune_out_filepath
                ]

                subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                raise Exception(error_message).with_traceback(e.__traceback__)

            # Get kitsune output message
            kitsune_out_content = open(kitsune_out_filepath).read().strip()

            # Retrieve the optimal kmer size
            kmer_size = int(kitsune_out_content.split(" ")[-1])

        self.metadata["kmer_size"] = kmer_size

        if not filter_size:
            # Search for the best bloom filter size based on the input set of genomes with ntcard
            try:
                # Estimate the bloom filter size with ntcard
                command_line = [
                    "ntcard",
                    f"--kmer={self.metadata['kmer_size']}",
                    f"--threads={self.nproc}",
                    f"--pref={os.path.join(self.tmp, 'genomes')}",
                    f"@{input_list_filepath}"
                ]

                subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                raise Exception(error_message).with_traceback(e.__traceback__)

            # Total number of kmers in reads
            F1 = 0

            # Number of distinct kmers
            F0 = 0

            # List with number of kmers occurring less than `min_kmer_occurrence`
            fs = list()

            with open(os.path.join(self.tmp, f"genomes_k{self.metadata['kmer_size']}.hist")) as hist_file:
                for line in hist_file:
                    line = line.strip()

                    if line:
                        line_split = line.split()

                        if line_split[0] == "F1":
                            F1 = int(line_split[-1])

                        elif line_split[0] == "F0":
                            F0 = int(line_split[-1])

                        elif isinstance(line_split[0], int) and int(line_split[0]) < min_kmer_occurrence:
                            fs.append(int(line_split[-1]))

                        else:
                            break

            # Use F1 as the estimated bloom filter size if F0 == 0
            # It could happen for very small genomes
            filter_size = F1 if F0 == 0 else F0-sum(fs)

            if filter_expand_by:
                # Expand the filter size
                # Transform `filter_expand_by` from percentage to absolute value based on `filter_size`
                filter_expand_by = int((filter_size*filter_expand_by)/100.0)

                filter_size += filter_expand_by

        self.metadata["filter_size"] = filter_size

        # Store the global parameters into the metadata.json file under the database root folder
        self._dump_metadata()

    def cluster(self, genomes: Dict[str, str], threshold: float=0.05) -> Dict[str, str]:
        """Cluster group of genomes based on a specific threshold on the ANI distance.

        Parameters
        ----------
        genomes : dict
            A dictionary with genomes and their taxonomic label.
            Genomes could be defined with paths to their fasta file or paths to their bloom filter sketch representation.
        threshold : float, default 0.05
            Threshold on the ANI distance used to cut the dendrogram and reshape clusters.

        Raises
        ------
        FileNotFoundError
            In case an input file does not exist.
        ValueError
            In case the extension of an input file is not supported.
            See `_is_supported` for additional information.

        Returns
        -------
        dict
            A new dictionary with genomes and their new taxonomic label.
        """

        # Check whether the genomes exist
        for genome_filepath in genomes:
            if not os.path.isfile(genome_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), genome_filepath)

            elif not self.__class__._is_supported(genome_filepath):
                raise ValueError(f"This file type is not supported: {genome_filepath}")

        assignments = dict()

        sketches = list()

        for genome_filepath in genomes:
            if os.path.splitext(genome_filepath)[1] == ".bf":
                sketches.append(genome_filepath)

            else:
                # In case of fasta files
                filename = os.path.splitext(os.path.basename(genome_filepath))[0]

                genome_obj = Entry(self, filename, filename, "genome")

                # Build their bloom filter sketch representation
                genome_sketch_filepath = genome_obj.sketch(genome_filepath)

                sketches.append(genome_sketch_filepath)

        # Build a condensed ANI distance matrix
        condensed_distance_matrix = list()

        # Rescale nproc
        nproc = self.nproc if len(sketches) > self.nproc else len(sketches)

        if nproc > 1:
            # Keep track of the ANI distances
            dists = dict()

            with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(sketches)) as pbar:
                # Wrapper around the update function of tqdm
                def progress(*args):
                    pbar.update()

                jobs = [
                    pool.apply_async(
                        self.__class__.dist, 
                        args=(
                            sketch_filepath, 
                            sketches[pos+1:],
                            self.metadata["kmer_size"],
                            self.tmp,
                            False,
                        ),
                        callback=progress
                    ) for pos, sketch_filepath in enumerate(sketches)
                ]

                for job in jobs:
                    sketch_filepath, sketch_dists = job.get()

                    # `sketch_dists` is an OrderedDict so it maintains the same order of elements in `sketches`
                    dists[sketch_filepath] = list(sketch_dists.values())

            for sketch_filepath in sketches:
                condensed_distance_matrix.extend(dists[sketch_filepath])

        else:
            # Avoid using multiprocessing if `nproc` is 1
            for pos, sketch_filepath in enumerate(sketches):
                _, sketch_dists = self.__class__.dist(sketch_filepath, sketches[pos+1:], self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                condensed_distance_matrix.extend(list(sketch_dists.values()))

        # Build a dendrogram based on the ANI distances between unknown genomes
        # Method: average-linkage
        dendro = fastcluster.linkage(condensed_distance_matrix, method="average")

        # Finally, cut the dendrogram on the input threshold
        if len(sketches) > 1:
            # The dendrogram exists in case of >1 sketches
            clusters = hier.fcluster(dendro, threshold, criterion="distance")

        else:
            # There is only one sketch here
            clusters = [1]

        if len(set(clusters)) == 1:
            # Define the cluster label by majority voting
            # Count the number of occurrences for each of the taxonomic labels first
            counts = Counter(genomes.values())

            # Get the most occurring one
            max_count = max(counts.values())

            most_frequent = [label for label, count in counts.items() if count == max_count]

            # Get the first one in alphabetical order
            # This is required in case of equally occurring taxonomic labels
            label = sorted(most_frequent)[0]

            # Apply the assignment
            assignments[label] = set(genomes.keys())

        else:
            # Group genomes according to the new clustering
            clusters_map = {cluster: set() for cluster in set(clusters)}

            for genome_filepath, cluster in zip(genomes.keys(), clusters):
                clusters_map[cluster].add(genome_filepath)

            # Keep track of the label occurrences in different clusters
            labels_count = dict()

            for cluster in clusters_map:
                # Assign a taxonomic label to the new clusters based on the majority voting
                # Count the number of occurrences for each of the taxonomic labels first
                counts = Counter([genomes[genome_filepath]] for genome_filepath in clusters_map[cluster])

                # Get the most occurring one
                max_count = max(counts.values())

                most_frequent = [label for label, count in counts.items() if count == max_count]

                # Get the first one in alphabetical order
                # This is required in case of equally occurring taxonomic labels
                label = sorted(most_frequent)[0]

                # Check whether a cluster with this taxonomic label already exists
                if label not in labels_count:
                    labels_count[label] = 1

                else:
                    if label in assignment:
                        # There are more than 1 clusters with the same label
                        # Rename the first one as clade 1
                        clade_1 = f"{label}__clade_1"

                        assignments[clade_1] = assignments.pop(label)

                    labels_count[label] += 1

                    # Add an incremental number
                    label += f"__clade_{labels_count[label]}"

                # Apply the assignment
                assignments[label] = clusters_map[cluster]

        return {genome: label for label in assignments for genome in assignments[label]}

    @staticmethod
    def _is_known(instance: "Database", filepath: os.path.abspath) -> Tuple[os.path.abspath, Optional[str]]:
        """Just a wrapper aroung the `is_known` function to make it callable in multiprocessing.
        It is safe to run in multiprocessing because it does not modify any instance attributes.

        Parameters
        ----------
        instance : Database
            A database instance.
        filepath : os.path.abspath
            Path to the uncompressed genome file.

        Returns
        -------
        str
            The assigned taxonomic label.
            See the `is_known` function description for additional information.
        """

        return (filepath, instance.is_known(filepath))

    def is_known(self, filepath: os.path.abspath) -> Optional[str]:
        """Check whether an input MAG could be characterized to any species in the database.
        This must be run for MAGs only and always before `add()`.

        Parameters
        ----------
        filepath : os.path.abspath
            Path to the uncompressed genome file.

        Returns
        -------
        str
            The assigned taxonomy defined up to the species level or None is case there are
            not species closed enough to the input genome.
        """
    
        if not self.__class__._validate_metadata(self.metadata):
            raise Exception("No database metadata found!")

        if not self.__class__._is_supported(filepath):
            raise ValueError("Input file format is not supported!")

        filename = os.path.splitext(os.path.basename(filepath))[0]

        if filename in self.genomes:
            # In case the input genome already exists in the database return its taxonomic label
            genome_obj = self.genomes[filename]

            return genome_obj.get_full_taxonomy()

        # Check whether the database report with info about clusters exists
        # This is probably not needed because we already check whether clusters are defined in the database
        if not self.report:
            raise Exception("No database report defined!")

        # Define an Entry object for the current genome
        genome_obj = Entry(self, filename, filename, "genome", taxonomy=None, known=False)

        # Build the genome sketch
        # Use k-mer size, minimum k-mer occurrence, and filter size from database metadata
        genome_sketch_filepath = genome_obj.sketch(filepath)

        # Search for the closest species and genome in the database
        genome_obj.profile = self.profile(filepath, genome_sketch_filepath)

        characterized = False

        taxonomy = None

        # The taxonomy of the closest species cluster could differ from the taxonomy of the closest genome
        # We need to test it against the species assigned to the closest genome first
        # The, if it is not assigned, we need to test it against the closest species cluster
        # Otherwise, it is not assigned
        for level in ["genome", "species"]:
            matches = genome_obj.profile[level]

            # Closest clusters and genomes come with their ANI distances
            # The smaller the better
            taxonomies = sorted(matches.keys(), key=lambda match: matches[match])

            # Get the best match only
            # We could comment this line to test it against all the other closest cluster under this level
            taxonomies = [taxonomies[0]]

            # Sort matches based on the best score
            if level == "genome":
                # Trim the last level t out of the taxonomic label
                taxonomies = ["|".join(match.split("|")[:-1]) for match in taxonomies]

            species_ids = list()

            centroid_sketches = list()

            # Check whether the current genome is close enough to its closest cluster
            for closest_species_taxonomy in taxonomies:
                closest_species_name = closest_species_taxonomy.split("|")[-1]

                # Retrieve the closest species identifier
                closest_species_id = self.clusters["species"][closest_species_name].identifier

                # Get the name of the closest species centroid
                closest_species_centroid = self.report[closest_species_id]["centroid"]

                # Retrieve the centroid of the closest species from the report and locate its sketch
                closest_species_centroid_sketch = os.path.join(self.root, "sketches", f"{closest_species_centroid}.bf")

                species_ids.append(closest_species_id)

                centroid_sketches.append(closest_species_centroid_sketch)

            # Compute the distance between the input genome and the closest cluster centroids
            _, dists = self.__class__.dist(genome_sketch_filepath, centroid_sketches, self.metadata["kmer_size"], tmp=self.tmp, resume=False)

            for pos, closest_species_taxonomy in enumerate(taxonomies):
                closest_species_id = species_ids[pos]

                distance_from_centroid = dists[centroid_sketches[pos]]

                # Retrieve the boundaries of the closest species cluster
                _, max_boundary = self._estimate_boundaries(self.report[closest_species_id]["taxonomy"])

                # Compare the distance from the centroid with the closest cluster boundaries for assignment
                if distance_from_centroid <= max_boundary:
                    taxonomy = closest_species_taxonomy

                    # Enable `characterized` to assign the current genome to the fully defined taxonomic label
                    characterized = True

                    break

            if characterized:
                # We have found an assignment
                break

        return taxonomy

    def add(
        self,
        filepath: os.path.abspath,
        reference: bool=False,
        taxonomy: Optional[str]=None,
    ) -> None:
        """Assign a genome to a specific taxonomy.
        This must be run also for genomes for which their taxonomy is not known (i.e., `taxonomy=None`).

        Parameters
        ----------
        filepath : os.path.abspath
            Path to the uncompressed genome file.
        reference : bool, default False
            The type of the input genome. True if `filepath` is a reference genome, False otherwise.
        taxonomy : str, optional
            The full taxonomic label in case `reference` is True.
            it is supposed to be defined from the kingdom up to the species level with a total of seven taxonomic levels.

        Raises
        ------
        Exception
            - If no database metadata is found;
            - If `reference` is False and there are no clusters in the database;
            - If a genome with the same name of `filepath` already exists in the database;
        ValueError
            - If the input file format is not supported;
            - If `reference` is True but `taxonomy` is not provided;
            - If `reference` is True but `taxonomy` does not pass the validation.
        """

        if not self.__class__._validate_metadata(self.metadata):
            raise Exception("No database metadata found!")

        if not self.__class__._is_supported(filepath):
            raise ValueError("Input file format is not supported!")

        if reference and not taxonomy:
            raise ValueError("A full taxonomic label must be defined for reference genomes!")

        if not reference and not self.clusters:
            raise Exception("The database does not contain any clusters!")

        filename = os.path.splitext(os.path.basename(filepath))[0]

        if filename in self.genomes:
            raise Exception("A genome with the same name already exists in the database!")

        # Define an Entry object for the current genome
        genome_obj = Entry(self, filename, filename, "genome", taxonomy=taxonomy, known=reference)

        # Retrieve the genome sketch
        # Use the already defined sketch if it already exists
        genome_sketch_filepath = genome_obj.sketch(filepath)

        if taxonomy:
            # Check whether the assigned taxonomy is fully defined
            if not self.__class__._is_defined(taxonomy):
                raise ValueError("The assigned taxonomy is malformed or not defined!")

            # Format the taxonomic label
            taxonomy = self.__class__._format_taxonomy(taxonomy)

            taxonomy_split = taxonomy.split("|")

            new_clusters = False

            # Start from the species all the way up to the kingdom
            for taxonomic_position, taxonomic_level in reversed(list(enumerate(taxonomy_split))):
                cluster_obj = None

                if taxonomic_level in self.clusters[self.__class__.LEVELS[taxonomic_position]]:
                    # Retrieve the cluster from the database
                    cluster_obj = self.clusters[self.__class__.LEVELS[taxonomic_position]][taxonomic_level]

                else:
                    if re.search("^MSBT[0-9]*$", taxonomic_level[3:]):
                        # The name of the taxonomic level could match with the internal MSBT cluster naming convention
                        # In this case, do not bump the clusters counter
                        cluster_id = taxonomic_level[3:]

                    else:
                        # Increment the clusters counter
                        # Always start from 1 (MSBT0 is reserved)
                        self.metadata["clusters_count"] += 1

                        # Define the internal cluster identifier
                        cluster_id = f"MSBT{self.metadata['clusters_count']}"

                    # Retrieve the batch number from the cluster id
                    batch_number = self.__class__._get_cluster_batch(cluster_id)

                    # Define the cluster folder
                    cluster_folder = os.path.join(self.root, "clusters", batch_number, cluster_id)

                    os.makedirs(cluster_folder)

                    # Define the parent level
                    parent = None if taxonomic_position == 0 else taxonomy_split[taxonomic_position-1]

                    # Initialize the Entry object representing the cluster
                    cluster_obj = Entry(self, cluster_id, taxonomic_level, self.__class__.LEVELS[taxonomic_position], folder=cluster_folder, parent=parent)

                    new_clusters = True

                if self.__class__.LEVELS[taxonomic_position] == "species":
                    # Set the parent species cluster
                    genome_obj.parent = taxonomic_level

                    # Register the input genome
                    self.genomes[filename] = genome_obj

                    # This cluster is a species and its children are genomes
                    cluster_obj.add_children(filename)

                else:
                    cluster_obj.add_children(taxonomy_split[taxonomic_position+1])

                # Register the cluster
                # Note: genomes cannot be indexed here!
                #       once all the input genomes have been processed, the `self.update()` function should be run as the last step
                #       in order to go through all the new and modified custers and build the Sequence Bloom Trees.
                self.clusters[self.__class__.LEVELS[taxonomic_position]][taxonomic_level] = cluster_obj

            if new_clusters:
                self._dump_metadata()

            self.__clusters.append(taxonomy)

        else:
            # Load the genome profiles
            genome_obj.profile = self.profile(filepath, genome_sketch_filepath)

            # This is a MAG and there are no species clusters closed enough in the database to characterize it
            # Keep track of the input genome object for further processing
            self.__unknowns.append(genome_obj)

    def characterize(self) -> Tuple[Dict[str, Set[str]], Set[str]]:
        """Process all the unassigned genomes together to define new clusters at different taxonomic levels (step 2).

        Raises
        ------
        Exception
            - If no database metadata is found;
            - If there are no unknown genomes to be processed.

        Returns
        -------
        tuple
            A tuple containing a dictionary with the genomes organized by taxonomic assignment, and
            the set of genomes that too far from everything in the database.
        """

        if not self.__class__._validate_metadata(self.metadata):
            raise Exception("No database metadata found!")

        if not self.__unknowns:
            raise Exception("There are no unknown genomes to be characterized!")

        # Retrieve the paths to the genome sketches
        sketches = [genome.sketch_filepath for genome in self.__unknowns]

        if len(sketches) > 1:
            # Build a condensed ANI distance matrix
            condensed_distance_matrix = list()

            print(f"Computing pair-wise distances between {len(sketches)} uncharacterized genomes")

            # Rescale nproc
            nproc = self.nproc if len(sketches) > self.nproc else len(sketches)

            if nproc > 1:
                # Keep track of the ANI distances
                dists = dict()

                with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(sketches)) as pbar:
                    # Wrapper around the update function of tqdm
                    def progress(*args):
                        pbar.update()

                    jobs = [
                        pool.apply_async(
                            self.__class__.dist, 
                            args=(
                                sketch_filepath, 
                                sketches[pos+1:],
                                self.metadata["kmer_size"],
                                self.tmp,
                                False,
                            ),
                            callback=progress
                        ) for pos, sketch_filepath in enumerate(sketches)
                    ]

                    for job in jobs:
                        sketch_filepath, sketch_dists = job.get()

                        # `sketch_dists` is an OrderedDict so it maintains the same order of elements in `sketches`
                        dists[sketch_filepath] = list(sketch_dists.values())

                for sketch_filepath in sketches:
                    condensed_distance_matrix.extend(dists[sketch_filepath])

            else:
                # Avoid using multiprocessing if `nproc` is 1
                for pos, sketch_filepath in enumerate(sketches):
                    print(f"Computing distance: {os.path.splitext(os.path.basename(sketch_filepath))[0]} [{pos+1}/{len(sketches)}]")

                    _, sketch_dists = self.__class__.dist(sketch_filepath, sketches[pos+1:], self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                    condensed_distance_matrix.extend(list(sketch_dists.values()))

            # Build a dendrogram based on the ANI distances between unknown genomes
            # Method: average-linkage
            dendro = fastcluster.linkage(condensed_distance_matrix, method="average")

        # Assignments map
        assignments = dict()

        # This contains genomes that have not been characterized not even at the kingdom level
        unassigned = set()

        # These genomes have been already tested against their closest species
        processed = set()

        print(f"Clustering {len(sketches)} genomes. This may take a while..")

        # Start processing from the genus level up to the kingdom
        # We can skip the species level here since we already went through it during `add()`
        for level in reversed(self.__class__.LEVELS[:-1]):
            if len(processed) == len(self.__unknowns):
                break

            print(f"Clustering at the {level} level")

            for pos, genome_obj in enumerate(self.__unknowns):
                if genome_obj.sketch_filepath not in processed:
                    matches = genome_obj.profile[level]

                    # Closest clusters come with their ANI distances
                    # The smaller the better
                    taxonomies = sorted(matches.keys(), key=lambda match: matches[match])

                    # Get the best match only
                    # We could comment this line to test it against all the other closest cluster under this level
                    taxonomies = [taxonomies[0]]

                    centroid_sketches = list()

                    # Check whether the current genome is close enough to its closest cluster
                    for closest_cluster_taxonomy in taxonomies:
                        current_level = level

                        closest_cluster_name = closest_cluster_taxonomy.split("|")[-1]

                        # Retrieve the closest cluster identifier
                        closest_cluster_id = self.clusters[current_level][closest_cluster_name].identifier

                        closest_cluster_centroid = self.report[closest_cluster_id]["centroid"]

                        # Centroids are always genomes
                        closest_cluster_centroid_sketch = self.genomes[closest_cluster_centroid].sketch_filepath

                        centroid_sketches.append(closest_cluster_centroid_sketch)

                    # Compute the distance between the unknown genome and the closest cluster centroids
                    _, dists = self.__class__.dist(genome_obj.sketch_filepath, centroid_sketches, self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                    for pos, closest_cluster_taxonomy in enumerate(taxonomies):
                        distance_from_centroid = dists[centroid_sketches[pos]]

                        try:
                            # Retrieve the closest cluster boundaries
                            _, max_boundary = self._estimate_boundaries(self.report[closest_cluster_id]["taxonomy"])

                        except Exception:
                            # This happens in case the distance is outside the kingdom boundaries
                            # The input genome will be marked as unassigned
                            continue

                        # Compare the distance from the centroid with the closest cluster boundaries for assignment
                        if distance_from_centroid <= max_boundary:
                            if len(sketches) > 1:
                                # The dendrogram exists in case of >1 sketches
                                # The current genome must be assigned to the closest cluster
                                # We should cut the dendrogram using the closest cluster boundaries to check for other assignments
                                clusters = hier.fcluster(dendro, max_boundary, criterion="distance")

                            else:
                                # There is only one sketch here
                                clusters = [1]

                            # Search for the cluster id assigned to the current genome
                            cluster_id = clusters[pos]

                            # Collect all the genomes that have been assigned to the same cluster id
                            for sketch_filepath, assigned_cluster in zip(sketches, clusters):
                                if assigned_cluster == cluster_id and sketch_filepath not in processed:
                                    # Get the genome name from the sketch filepath
                                    sketch_name = os.path.splitext(os.path.basename(sketch_filepath))[0]

                                    # All genomes that fall under this cluster are assigned to an already defined cluster in the database
                                    assignments[sketch_name] = closest_cluster_taxonomy

                                    # Mark the assigned genomes as processed
                                    processed.add(sketch_filepath)

                                    print(f"\t[{len(processed)}/{len(self.__unknowns)}] {level}={closest_cluster_taxonomy}\t{sketch_name}")

                            # We have found a match under the current level
                            break

                if level == "kingdom" and genome_obj.sketch_filepath not in processed:
                    # This genome has not been characterized, not even at the kingdom level
                    # We should simply report these genomes, we are not going to define new kingdoms here
                    unassigned.add(genome_obj.name)

        print(f"{len(unassigned)} genomes have not been characterized")

        # Get full and partial assignments first
        characterized = dict()

        partially_characterized = dict()

        partially_characterized_genomes = set()

        for sketch_name in assignments:
            taxonomy = assignments[sketch_name]

            if (taxonomy.count("|") + 1) == len(self.__class__.LEVELS):
                # The assigned taxonomic label is defined at the species level
                if taxonomy not in characterized:
                    characterized[taxonomy] = set()

                characterized[taxonomy].add(sketch_name)

            else:
                # The assigned taxonomic label is partially defined
                if taxonomy not in partially_characterized:
                    partially_characterized[taxonomy] = set()

                partially_characterized[taxonomy].add(sketch_name)

                partially_characterized_genomes.add(sketch_name)

        print(f"{len(partially_characterized_genomes)} genomes have been partially characterized")

        if partially_characterized:
            print("Processing partially characterized genomes. This will create new clusters..")

        processed = 0

        new_clusters = False

        # We can now finish characterizing the partially characterized genomes
        # Genomes in `partially_characterized` are all characterized at the kingdom level at least
        while len(partially_characterized) > 0:
            for taxonomy in list(partially_characterized.keys()):
                print(f"Expanding {taxonomy}")

                # We previously defined the assignments starting from the species level up to the kingdom
                # We should now start from the kingdom level down to the species level to finish characterizing genomes
                # If a genome has not been characterized, it means that it is not close enough to its closest cluster according to its boundaries
                # We can now create new clusters and estimate their boundaries, then cut the dendrogram to see how many genomes fall in them
                tmp_taxonomy = f"{taxonomy}|{self.__class__.LEVELS[taxonomy.count('|')+1][0]}__tmp"

                # Estimate the boundaries for the temporary cluster
                # This is a totally new cluster. There is no need to search for a centroid here
                _, max_boundary = self._estimate_boundaries(tmp_taxonomy)

                if len(sketches) > 1:
                    # The dendrogram exists in case of >1 sketches
                    # Cluster genomes according to the temporary cluster's boundaries
                    clusters = hier.fcluster(dendro, max_boundary, criterion="distance")

                else:
                    # There is only one sketch here
                    clusters = [1]

                # Keep track of the mapping between fcluster ids and database cluster ids
                clusters_map = dict()

                for sketch_filepath, assigned_cluster in zip(sketches, clusters):
                    # Get the sketch file name
                    sketch_name = os.path.splitext(os.path.basename(sketch_filepath))[0]

                    # Consider genomes under the current partially characterized taxonomy only
                    if sketch_name in partially_characterized[taxonomy]:
                        if assigned_cluster not in clusters_map:
                            # Create a new MSBT cluster
                            # Increment the cluster counter
                            self.metadata["clusters_count"] += 1

                            clusters_map[assigned_cluster] = self.metadata["clusters_count"]

                            new_clusters = True

                        # Define the assigned taxonomic label
                        assigned_taxonomy = f"{taxonomy}|{self.__class__.LEVELS[taxonomy.count('|')+1][0]}__MSBT{clusters_map[assigned_cluster]}"

                        if self.__class__._is_defined(assigned_taxonomy):
                            # This is fully characterized now
                            if assigned_taxonomy not in characterized:
                                characterized[assigned_taxonomy] = set()

                            characterized[assigned_taxonomy].add(sketch_name)

                            # This is a mag and we already established its taxonomy
                            # We should force assign it to avoid going through the characterization again
                            # The first argument should be the genome file path, but we can safely use its sketch file path here
                            self.add(sketch_filepath, reference=False, taxonomy=assigned_taxonomy)

                            # Just a counter of fully characterized genomes
                            processed += 1

                            print(f"\t[{processed}/{len(partially_characterized_genomes)}] species={assigned_taxonomy}\t{sketch_name}")

                        else:
                            # This is not fully characterized yet
                            if assigned_taxonomy not in partially_characterized:
                                partially_characterized[assigned_taxonomy] = set()

                            partially_characterized[assigned_taxonomy].add(sketch_name)

                # The partial taxonomy must be removed from the bucket of partially characterized labels
                del partially_characterized[taxonomy]

        # In case of an update with only reference genomes with no force assignment `force=False`
        # these new clusters should actually be transformed to known clusters based on a majority voting
        # mechanism on their genomes' taxonomic labels
        # TODO

        # Reset the list of unknowns
        self.__unknowns = list()

        if new_clusters:
            self._dump_metadata()

        return characterized, unassigned

    def _dump_metadata(self) -> None:
        """Dump database metadata to a JSON file under the database root folder.
        """

        metadata_json_filepath = os.path.join(self.root, "metadata.json")

        # Update metadata.json with the new clusters counter
        with open(metadata_json_filepath, "w+") as metadata_json_file:
            json.dump(self.metadata, metadata_json_file)

    def _estimate_boundaries(self, taxonomy: str) -> Tuple[float, float]:
        """Estimate the boundaries of a given taxonomic entry.

        Parameters
        ----------
        taxonomy : str
            A taxonomic label.

        Raises
        ------
        Exception
            If it is unable to estimate the cluster boundaries for the input `taxonomy`.
        ValueError
            If a the input `taxonomy` is not provided.

        Returns
        -------
        tuple
            A tuple with cluster boundaries.
        """

        if not taxonomy:
            raise ValueError("Taxonomy is not provided!")

        cluster_level = self.__class__.LEVELS[taxonomy.count("|")]

        if cluster_level == "species":
            # Force the radius of species clusters to 5% of genetic distance
            return (0.0, 0.05)

        cluster_name = taxonomy.split("|")[-1]

        # Check whether the input taxonomy is already defined in the database
        if self.taxonomy_exists(taxonomy):
            # Retrieve the cluster identifier
            cluster_id = self.clusters[cluster_level][cluster_name].identifier

            min_ani, max_ani = self.report[cluster_id]["boundaries"]

            if min_ani is not None and max_ani is not None:
                return (min_ani, max_ani)

        if cluster_level == "kingdom":
            raise Exception(f"Unable to retrieve boundaries for {taxonomy}")

        # Keep track of the minimum and maximum ANIs to define boundaries
        min_bounds = list()
        max_bounds = list()

        # Keep track of the input cluster level
        current_cluster_level = cluster_level

        # Retrieve the parent taxonomy and its cluster name
        current_cluster_parent_taxonomy = "|".join(taxonomy.split("|")[:-1])
        
        current_cluster_parent = current_cluster_parent_taxonomy.split("|")[-1]

        while not min_bounds and not max_bounds:
            parent_level = self.__class__.LEVELS[self.__class__.LEVELS.index(current_cluster_level)-1]

            if current_cluster_parent in self.clusters[parent_level]:
                parent_obj = self.clusters[parent_level][current_cluster_parent]

                # Retrieve all the children up to the original input cluster level starting from the current cluster
                parent_children = parent_obj.get_children(up_to=cluster_level)

                if parent_children:
                    for child in parent_children:
                        parent_children_obj = self.clusters[cluster_level][child]

                        if parent_children_obj.identifier in self.report:
                            min_ani, max_ani = self.report[parent_children_obj.identifier]["boundaries"]

                            if min_ani is not None and max_ani is not None:
                                min_bounds.append(min_ani)
                                max_bounds.append(max_ani)

            if parent_level == "kingdom":
                break

            # Update the current cluster level, parent taxonomy, and parent cluster name
            current_cluster_level = parent_level

            current_cluster_parent_taxonomy = "|".join(current_cluster_parent_taxonomy.split("|")[:-1])

            current_cluster_parent = current_cluster_parent_taxonomy.split("|")[-1]

        if not min_bounds and not max_bounds:
            #return 0.0, 1.0
            raise Exception(f"Unable to retrieve boundaries for {taxonomy}")

        return (round(statistics.mean(min_bounds), 5), round(statistics.mean(max_bounds), 5))

    def update(self) -> None:
        """Process all the clusters created or modified during the `self.add()` run, and build Sequence Bloom Trees (step 3).

        Raises
        ------
        Exception
            If there are no clusters that have been created or modified after `self.add()`.
        """

        if not self.__clusters:
            raise Exception("There is nothing to update in the database!")

        if self.__unknowns:
            raise Exception("There are still unknown genomes in cache! Please run `self.characterize()` first")

        processed = set()

        # Process clusters based on their taxonomic level
        # All the species first, then genera, and so on up to the kingdom
        # Note that clusters are fuly defined here, but could contain some level of unknownness
        for level in reversed(self.__class__.LEVELS):
            pos = self.__class__.LEVELS.index(level)

            for taxonomy in self.__clusters:
                # Keep levels up to `level` and trim out lower taxonomic levels
                partial_taxonomy = "|".join(taxonomy.split("|")[:pos+1])

                if partial_taxonomy not in processed:
                    print(f"Update: {partial_taxonomy}")

                    # Get the last cluster id
                    cluster = partial_taxonomy.split("|")[-1]

                    # Retrieve the cluster object
                    cluster_obj = self.clusters[level][cluster]

                    # Finally, build the Sequence Bloom Tree
                    cluster_obj.index()

                    processed.add(partial_taxonomy)

        # Retrieve the set of kingdoms in the database
        kingdoms = set(self.clusters["kingdom"].keys())

        # MSBT0 is reserved
        db_id = "MSBT0"

        # Retrieve the batch number from MSBT0
        batch_number = self.__class__._get_cluster_batch(db_id)

        # Define the db folder
        db_folder = os.path.join(self.root, "clusters", batch_number, db_id)

        # Define a new cluster to index the kingdom entries
        # The id MSBT0 is reserved for this specific entry
        # It must stay out of the clusters folder
        db_obj = Entry(self, "MSBT0", "db", None, folder=db_folder, parent=None, children=kingdoms)

        # Index the kingdom entries
        db_obj.index()

        # Dump the report
        self._dump_report()

        # Dump the list of genomes and their assignments
        self._dump_genomes()

        # Load the report
        self.report = self.__class__._load_report(os.path.join(self.root, "clusters.tsv"))

        # Reset the list of clusters
        self.__clusters = list()

    def _dump_genomes(self) -> None:
        """Dump the list of reference genomes and mags with their assignments.
        """

        # Search for already existing genomes.txt files
        genomes_files = len(list(Path(self.root).glob("genomes*.tsv")))

        # Define the path to the genomes.tsv file
        genomes_filepath = os.path.join(self.root, "genomes.tsv")

        if os.path.isfile(genomes_filepath):
            # Archive the latest genomes table if it exists
            shutil.move(genomes_filepath, os.path.join(self.root, f"genomes_{genomes_files+1}.tsv"))

        # Dump the genomes table
        with open(genomes_filepath, "w+") as genomes_file:
            # The first line should contain the timestamp
            today = datetime.date.today()

            genomes_file.write(f"# MetaSBT genomes {today.year}-{today.month}-{today.day}\n")

            # Define the header
            header = [
                "genome_id",
                "type",
                "assignment",
                "internal",
            ]

            genomes_file.write("# {}\n".format("\t".join(header)))

            # Dump information about all the genomes in the database
            for genome in self.genomes:
                genomes_info = [
                    genome,
                    "reference" if self.genomes[genome].is_known() else "mag",
                    self.genomes[genome].get_full_taxonomy(),
                    self.genomes[genome].get_full_taxonomy(internal=True),
                ]

                genomes_file.write("{}\n".format("\t".join(genomes_info)))

    def _dump_report(self) -> None:
        """Dump the list of clusters at all the seven taxonomic levels.
        """

        # Search for already existing report clusters.tsv files
        report_files = len(list(Path(self.root).glob("clusters*.tsv")))

        # Define the path to the clusters.tsv file
        report_filepath = os.path.join(self.root, "clusters.tsv")

        if os.path.isfile(report_filepath):
            # Archive the latest report file if it exists
            shutil.move(report_filepath, os.path.join(self.root, f"clusters_{report_files+1}.tsv"))

        # Dump the report table
        with open(report_filepath, "w+") as report_file:
            # The first line should contain a timestamp
            today = datetime.date.today()

            report_file.write(f"# MetaSBT report {today.year}-{today.month}-{today.day}\n")

            # Define the header
            header = [
                "cluster_id",
                "level",
                "density",
                "references_count",
                "mags_count",
                "references_list",
                "mags_list",
                "centroid",
                "known",
                "taxonomy",
                "internal",
                "min_ani",
                "max_ani",
            ]

            report_file.write("# {}\n".format("\t".join(header)))

            # Dump information about clusters at all the seven taxonomic levels
            # starting from the species all the way up to the kingdom
            for level in reversed(self.__class__.LEVELS):
                for cluster_name in self.clusters[level]:
                    # Retrieve the cluster object
                    cluster_obj = self.clusters[level][cluster_name]

                    # Retrieve the list of reference genomes and mags
                    if level == "species":
                        cluster_references = sorted([genome for genome in cluster_obj.children if self.genomes[genome].is_known()])

                    else:
                        cluster_references = sorted([cluster for cluster in cluster_obj.children if self.clusters[self.__class__.LEVELS[self.__class__.LEVELS.index(level)+1]][cluster].is_known()])

                    cluster_mags = list(cluster_obj.children.difference(cluster_references))

                    # The bloom filter density must be compute every time
                    cluster_density = cluster_obj.get_density()

                    processed = False

                    if cluster_obj.identifier in self.report:
                        # Check if there is any difference between the current and previous set of references and mags
                        same_references = len(self.report[cluster_obj.identifier]["references"].difference(cluster_references)) == 0

                        same_mags = len(self.report[cluster_obj.identifier]["mags"].difference(cluster_mags)) == 0

                        if same_references and same_mags:
                            cluster_min_ani, cluster_max_ani = self.report[cluster_obj.identifier]["boundaries"]

                            if cluster_min_ani is None:
                                cluster_min_ani = ""

                            if cluster_max_ani is None:
                                cluster_max_ani = ""

                            # Retrieve information from the previous report if there are no changes
                            cluster_info = [
                                cluster_obj.identifier,
                                cluster_obj.level,
                                str(cluster_density),
                                str(len(cluster_references)),
                                str(len(cluster_mags)),
                                ",".join(cluster_references),
                                ",".join(cluster_mags),
                                self.report[cluster_obj.identifier]["centroid"],
                                str(self.report[cluster_obj.identifier]["known"]),
                                self.report[cluster_obj.identifier]["taxonomy"],
                                self.report[cluster_obj.identifier]["internal"],
                                str(cluster_min_ani),
                                str(cluster_max_ani),
                            ]

                            processed = True

                    if not processed:
                        cluster_is_known = cluster_obj.is_known()

                        cluster_taxonomy = cluster_obj.get_full_taxonomy()

                        cluster_internal = cluster_obj.get_full_taxonomy(internal=True)

                        # Retrieve the cluster boundaries
                        cluster_min_ani, cluster_max_ani, cluster_centroid = cluster_obj.get_boundaries()

                        if cluster_min_ani is None:
                            cluster_min_ani = ""

                        if cluster_max_ani is None:
                            cluster_max_ani = ""

                        cluster_info = [
                            cluster_obj.identifier,
                            cluster_obj.level,
                            str(cluster_density),
                            str(len(cluster_references)),
                            str(len(cluster_mags)),
                            ",".join(cluster_references),
                            ",".join(cluster_mags),
                            cluster_centroid,
                            str(cluster_is_known),
                            cluster_taxonomy,
                            cluster_internal,
                            str(cluster_min_ani),
                            str(cluster_max_ani),
                        ]

                        # Update the report immediately to avoid recomputing the species clusters boundaries at the higher levels
                        self.report[cluster_obj.identifier] = dict()

                        # Retrieve the cluster level
                        self.report[cluster_obj.identifier]["level"] = cluster_obj.level

                        # Retrieve the cluster density
                        self.report[cluster_obj.identifier]["density"] = cluster_density
    
                        # Retrieve the set of reference genomes
                        self.report[cluster_obj.identifier]["references"] = set(cluster_references)

                        # Retrieve the set of metagenome-assembled genomes
                        self.report[cluster_obj.identifier]["mags"] = set(cluster_mags)

                        # Retrieve the cluster centroid
                        self.report[cluster_obj.identifier]["centroid"] = cluster_centroid

                        self.report[cluster_obj.identifier]["known"] = cluster_is_known

                        # Retrieve the taxonomic assignment
                        self.report[cluster_obj.identifier]["taxonomy"] = cluster_taxonomy

                        self.report[cluster_obj.identifier]["internal"] = cluster_internal

                        # Retrieve the cluster boundaries
                        cluster_min_ani = float(cluster_min_ani) if str(cluster_min_ani).strip() else None

                        cluster_max_ani = float(cluster_max_ani) if str(cluster_max_ani).strip() else None

                        self.report[cluster_obj.identifier]["boundaries"] = (cluster_min_ani, cluster_max_ani)

                    report_file.write("{}\n".format("\t".join(cluster_info)))

    @staticmethod
    def dist(
        sketch_filepath: os.path.abspath, 
        sketches: List[os.path.abspath], 
        kmer_size: int,
        tmp: os.path.abspath=None,
        resume: bool=False,
    ) -> Tuple[os.path.abspath, Dict[str, float]]:
        """Compute the Average Nucleotide Identity (ANI) between genome sketches.

        TODO: If the input genomes have been sketched starting from their translated aminoacid sequences,
              we could use this function to compute the Average Aminoacid Identity (AAI) that we could use
              as a metric to better cluster genomes at higher taxonomic levels.

        Parameters
        ----------
        sketch_filepath : os.path.abspath
            The path to the first genome sketch.
        sketches : list
            List with paths to the sketch files.
        kmer_size : int
            The kmer size.
        tmp : os.path.abspath, default None
            Path to the temporary folder.
            Use the current working directory if None.
        resume : bool, default False
            If True, load a distance table if it already exists.
            Otherwise, overwrite the results.

        Raises
        ------
        FileNotFoundError
            If `sketch_filepath` does not exist.

        Returns
        -------
        tuple
            A tuple with the path to the input `sketch_filepath` and a dictionary with the ANI distances between 
            the input genome and the sketches in `sketches` indexed by sketches names in the same order of `sketches`.
        """

        if not sketches:
            return (sketch_filepath, dict())

        if not tmp:
            # Use the current working directory as temporary folder if not specified
            tmp = os.getcwd()

        # Store the distance tables into a dedicated subfolder
        tmp = os.path.join(tmp, "distances")

        if not os.path.isdir(tmp):
            os.makedirs(tmp, exist_ok=True)

        # Testing whether all the sketches in `sketches` exist could be expensive
        # Check whether `sketch_filepath` exists only, and assume all the sketches in `sketches` exist
        if not os.path.isfile(sketch_filepath):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), sketch1_filepath)

        # Keep track of the ANI distances between sketches
        # The order of keys must be the same of the elements in `scketches`
        distances = OrderedDict()

        # Retrieve the sketch file name
        sketch_filename = os.path.splitext(os.path.basename(sketch_filepath))[0]

        # Define the path to the file with the list of sketches (filepaths)
        sketches_list_filepath = os.path.join(tmp, f"{sketch_filename}__list.txt")

        # Define the path to the distance table (intersect)
        intersect_filepath = os.path.join(tmp, f"{sketch_filename}__intersect.txt")

        # Define the path to the distance table (union)
        union_filepath = os.path.join(tmp, f"{sketch_filename}__union.txt")

        if not resume or (resume and not os.path.isfile(intersect_filepath) and not os.path.isfile(union_filepath)):
            # Always overwrite already existing results
            with open(sketches_list_filepath, "w+") as sketches_list:
                for sketch in sketches:
                    sketches_list.write(f"{sketch}\n")

            try:
                # Compute the intersection of kmers between `sketch_filepath` and all the other sketches
                command_line = [
                    "howdesbt",
                    "bfdistance",
                    f"--list={sketches_list_filepath}",
                    f"--focus={sketch_filepath}",
                    "--show:intersect",
                ]

                with open(intersect_filepath, "w+") as intersect_file:
                    subprocess.check_call(command_line, stdout=intersect_file, stderr=intersect_file)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                raise Exception(error_message).with_traceback(e.__traceback__)

            if not os.path.isfile(intersect_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), intersect_filepath)

            try:
                # Compute the union of kmers between `sketch_filepath` and all the other sketches
                command_line = [
                    "howdesbt",
                    "bfdistance",
                    f"--list={sketches_list_filepath}",
                    f"--focus={sketch_filepath}",
                    "--show:union",
                ]

                with open(union_filepath, "w+") as union_file:
                    subprocess.check_call(command_line, stdout=union_file, stderr=union_file)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                raise Exception(error_message).with_traceback(e.__traceback__)

            if not os.path.isfile(union_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), union_filepath)

        with open(intersect_filepath) as intersect_file, open(union_filepath) as union_file:
            # Skip the first line of both the intersect and union files
            next(intersect_file)
            next(union_file)

            # Iterate through the content of bfdistance results and compute the ANI distance
            # `intersect_file` and `union_file` have the same number of rows sorted in the same order of the elements in `sketches`
            for intersect_line, union_line in zip(intersect_file, union_file):
                # Split the line and remove extra spaces first
                intersect_line = " ".join(intersect_line.split()).split(" ")

                union_line = " ".join(union_line.split()).split(" ")

                # Get rid of non informative fields (intersection) and (union)
                if intersect_line[-1] == "(intersection)":
                    intersect_line = intersect_line[:-1]

                if union_line[-1] == "(union)":
                    union_line = union_line[:-1]

                # Retrieve the sketch filepath and the bfdistance results
                # The sketch filepath is the same in the intersect and union file tables
                target_filepath = intersect_line[0].split(":")[0]

                intersect = int(intersect_line[-1])

                union = int(union_line[-1])

                # Compute the Jaccard index
                jaccard_index = round(intersect/union, 5)

                if jaccard_index == 0.0:
                    # There is nothing in common here
                    # Return the max ANI distance
                    ani_distance = 1.0

                else:
                    # Compute the ANI as a distance measure
                    ani_distance = 1 - (1 + (1/kmer_size) * math.log((2*jaccard_index) / (1+jaccard_index)))

                    if ani_distance > 1.0:
                        # A very small jaccard index could lead to a very high result in the ANI distance computation
                        ani_distance = 1.0

                distances[target_filepath] = ani_distance

        return sketch_filepath, distances

    @staticmethod
    def _profile(
        instance: "Database",
        genome_filepath: os.path.abspath,
        sketch_filepath: os.path.abspath,
        uncertainty: float=50.0,
        pruning_threshold: float=0.0,
    ) -> Tuple[os.path.abspath, Dict[str, Dict[str, float]]]:
        """Just a wrapper around the `profile` function to make it callable in multiprocessing.
        It is safe to run in multiprocessing because it does not modify any instance attributes.

        Parameters
        ----------
        instance : Database
            A database instance.
        genome_filepath : os.path.abspath
            Path to the input genome file in fasta format.
        sketch_filepath : os.path.abspath
            Path to the sketch representation of the input genome.
        uncertainty : float, default 50.0
            Percentage of uncertainty used to expand the selection of best matches.
        pruning_threshold : float, deafult 0.0
            Percentage of number of kmer hits under which HowDeSBT prunes the Sequence Bloom Trees.
            This is applied at the kingdom level only in order to avoid selecting all the species clusters in the database.
            It must be between 0.0 and 1.0.

        Returns
        -------
        tuple
            A tuple with `genome_filepath` and a dictionary with a genome's profiles.
            See the `profile` function description for additional information.
        """

        return (genome_filepath, instance.profile(genome_filepath, sketch_filepath, uncertainty=uncertainty, pruning_threshold=pruning_threshold))

    def profile(
        self,
        genome_filepath: os.path.abspath,
        sketch_filepath: os.path.abspath,
        uncertainty: float=50.0,
        pruning_threshold: float=0.0,
    ) -> Dict[str, Dict[str, float]]:
        """Profile the input genome by querying the root node to establish the closest kingdom, and expanding 
        the subsequent queries up to the species level in order to establish the closest clusters at all the 
        seven taxonomic levels and the closest genome in the database.

        Parameters
        ----------
        genome_filepath : os.path.abspath
            Path to the input genome file in fasta format.
        sketch_filepath : os.path.abspath
            Path to the sketch representation of the input genome.
        uncertainty : float, default 50.0
            Percentage of uncertainty used to expand the selection of best matches.
        pruning_threshold : float, default 0.0
            Percentage of number of kmer hits under which HowDeSBT prunes the Sequence Bloom Trees.
            This is applied at the kingdom level only in order to avoid selecting all the species clusters in the database.
            It must be between 0.0 and 1.0.

        Raises
        ------
        Exception
            If an error occurs while running HowDeSBT.
        FileNotFoundError
            - If a tree definition file does not exist.
            - If the output of HowDeSBT does not exist.

        Returns
        -------
        dict
            A dictionary with the closest cluster and its ANI distance, indexed by the name of the taxonomic level.
        """

        # Define the list of taxonomic levels
        # Prepend the db level with the index of all the kingdoms in the database
        # Also append the genome level as a break point
        levels = ["db"] + self.__class__.LEVELS + ["genome"]

        # Keep track of the profiles
        # Mapping between the closest cluster and its distance indexexd by the taxonomic level
        # It may map multiple closest clusters
        profiles = {level: dict() for level in levels[1:]}

        clusters = list()

        # Get the input file name
        genome_filename = os.path.splitext(os.path.basename(genome_filepath))[0]

        profiles_dir = os.path.join(self.tmp, "profiles")

        if not os.path.isdir(profiles_dir):
            os.makedirs(profiles_dir, exist_ok=True)

        # Define the path to the output of `howdesbt query`
        # It will then contain the final profiles
        query_result_filepath = os.path.join(profiles_dir, f"{genome_filename}.txt")

        if os.path.isfile(query_result_filepath):
            try:
                # Read the output file if it exists and return the genomes' profiles
                with open(query_result_filepath) as profile_file:
                    for line in profile_file:
                        line = line.strip()

                        if line:
                            if not line.startswith("#"):
                                line_split = line.split("\t")

                                # Retrieve the taxonomic level
                                level = line_split[0]

                                # Retrieve the taxonomic label
                                label = line_split[1]

                                # Retrieve the ANI distance between the input genome and the closest cluster
                                ani = float(line_split[2])

                                profiles[level][label] = ani

                return profiles

            except Exception:
                # The query result file could be corrupted or could contain an intermediate result
                # because of previous unexpectedly interrupted jobs;
                # Remove `query_result_filepath` and query the database from scratch
                os.unlink(query_result_filepath)

        # 1
        # In case the input genome contains multiple contigs, we should collapse all of them into a single one
        # This is because of how HowDeSBT treats input queries: different reads are different queries
        # We should concatenate everything into a single sequence with the N character as a workaround
        # 
        # 2
        # We could also use the bloom filter representation of the input genome
        # However, we will not be able to apply the query command and we should then query every single node in a tree individually
        # This is a more expensive solution
        # 
        # We are going to proceed with the solution 1
        # We should first check how many contigs are in the input genome
        records = list(SeqIO.parse(genome_filepath, format="fasta"))

        collapsed = False

        if len(records) > 1:
            # This is a bottleneck!
            # Collapse records with the N character
            collapsed_sequence = Seq("N".join([str(record.seq) for record in records]))

            # Create a new SeqRecord entry with the collapsed sequences
            # Use the id, name, and description of the first record
            collapsed_record = SeqRecord(collapsed_sequence, id=records[0].id, name=records[0].name, description=records[0].description)

            # Finally dump the merged sequences into a temporary fasta file
            # We need this file to be persistent, and we will eventually delete it as the final step in this function
            with tempfile.NamedTemporaryFile(mode="w+t", delete=False) as temp_fasta_file:
                SeqIO.write(collapsed_record, temp_fasta_file, format="fasta")

            # We need this flag to get rid of the temporary fasta file
            collapsed = True

            # Replace the input genome file path with the path to the temporary fasta file
            genome_filepath = temp_fasta_file.name

        # This is to keep track of the ANI distances between the input genome and the centroids on the closest clusters
        # We need to store these information after the first query to avoid recomputing the distances again for all the other levels
        dists = dict()

        # This triggers the computation of the ANI distances versus all the species centroids under a particular level
        # It is set to False after the first iteration to avoid computing the ANI distances again
        first_iter = True

        # Iterate over the taxonomic levels
        for pos, level in enumerate(levels):
            if level == "genome":
                # We cannot query genomes
                break

            # Define the next taxonomic level
            next_level = levels[pos+1]

            if level == "db":
                # Start querying the database super Entry
                # identifier=MSBT0
                # name=db
                clusters = ["db"]

            best_level_matches = dict()

            for cluster in clusters:
                if level == "db":
                    cluster_id = "MSBT0"

                else:
                    # Retrieve the cluster id
                    cluster_id = self.clusters[level][cluster].identifier

                tree_filepath = os.path.join(self.root, "clusters", Database._get_cluster_batch(cluster_id), cluster_id, "tree", "index.detbrief.sbt")

                if not os.path.isfile(tree_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), tree_filepath)

                # We could use a minimal threshold at the kingdom level only to avoid selecting all the species clusters
                threshold = pruning_threshold if level == "kingdom" else 0.0

                command_line = [
                    "howdesbt",
                    "query",
                    "--sort",
                    "--distinctkmers",
                    f"--tree={tree_filepath}",
                    f"--threshold={threshold}",
                    genome_filepath,
                ]

                try:
                    with open(query_result_filepath, "w+") as query_result_file:
                        subprocess.check_call(command_line, stdout=query_result_file, stderr=query_result_file)

                except subprocess.CalledProcessError as e:
                    error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                    raise Exception(error_message).with_traceback(e.__traceback__)

                if not os.path.isfile(query_result_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), query_result_filepath)

                matches = dict()

                with open(query_result_filepath) as query_result_file:
                    for line in query_result_file:
                        line = line.strip()

                        if line:
                            if not line.startswith("#") and not line.startswith("*"):
                                line_split = line.split(" ")

                                node = line_split[0]

                                if node not in matches:
                                    matches[node] = {"common": 0, "total": 0}

                                hits = line_split[1].split("/")

                                # Keep track of the number of common kmers between the input query and the target bloom filter
                                matches[node]["common"] = int(hits[0])

                                # Keep track of the total number of kmers in the query
                                matches[node]["total"] = int(hits[1])

                if not matches:
                    # In case of `threshold` >0.0, `matches` could be empty
                    # Process the next cluster
                    continue

                matches = {node: matches[node]["common"]/matches[node]["total"] for node in matches}

                # Search for the best match
                best_match = sorted(matches.keys(), key=lambda match: matches[match])[-1]

                # Get the best score
                best_score = matches[best_match]

                # Define a threshold on the best score
                score_threshold = float((best_score*uncertainty)/100.0)

                best_matches = {match: matches[match] for match in matches if matches[match] >= best_score-score_threshold}

                # Keep track of the best matches under the same taxonomic level
                best_level_matches.update(best_matches)

            if level == "db":
                # We want to treat the db level differently
                # We compute the ANI distance between the input query and the kingdom bloom filters, without retrieving the species clusters
                # This allows to use the `--threshold` on the next iteration, when `level` is "kingdom", so that we can select the best subset of species clusters
                best_level_matches_arr = list(best_level_matches.keys())

                # Retrieve the sketch representation of the best matches
                # Best matches are kingdoms here
                best_level_match_sketches = [self.clusters[next_level][best_level_match].sketch_filepath for best_level_match in best_level_matches_arr]

                # Get the distances to the kingdom bloom filters
                _, kingdom_dists = self.__class__.dist(sketch_filepath, best_level_match_sketches, self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                for pos, best_level_match in enumerate(best_level_matches_arr):
                    # Retrieve the best match full taxonomic label
                    best_level_match_taxonomy = self.clusters[next_level][best_level_match].get_full_taxonomy()

                    best_level_match_sketch = best_level_match_sketches[pos]

                    # Keep track of the profile
                    profiles[next_level][best_level_match_taxonomy] = round(kingdom_dists[best_level_match_sketch], 5)

            elif level == "species":
                # At the species level, we have to search for the closest genomes
                # We don't want to report the ANI distances versus all the genomes in a particular species, so we select the very best one only
                # Best level matches are sorted based on the number of common kmers over the total number of kmers in the query
                # The higher the better
                top_level_match = sorted(best_level_matches.keys(), key=lambda match: best_level_matches[match])[-1]

                # Retrieve the Entry object of the best match and its sketch representation
                best_level_match_sketch = self.genomes[top_level_match].sketch_filepath

                # Compute the ANI distance between the input genome and all the best matches under the current level
                _, dists = self.__class__.dist(sketch_filepath, [best_level_match_sketch], self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                # Keys in the dists dictionary are paths to bloom filter sketches, so we should retrieve the file name from them
                best_level_match = os.path.splitext(os.path.basename(best_level_match_sketch))[0]

                # Retrieve the Entry object of the best match in terms of ANI distance
                best_level_match_obj = self.genomes[best_level_match]

                # Retrieve the best match full taxonomic label
                # Also add the closest genome to the last taxonomic level
                best_level_match_taxonomy = f"{best_level_match_obj.get_full_taxonomy()}|t__{best_level_match_obj.name}"

                # Keep track of the profile
                profiles[next_level][best_level_match_taxonomy] = round(dists[best_level_match_obj.sketch_filepath], 5)

            else:
                for best_level_match in best_level_matches:
                    # Retrieve the Entry object of the best match in terms of number of hits
                    best_level_match_obj = self.clusters[next_level][best_level_match]

                    # Retrieve the best match full taxonomic label
                    best_level_match_taxonomy = best_level_match_obj.get_full_taxonomy()

                    if next_level == "species":
                        # There are no species levels under the species level
                        # The species cluster is the best level match
                        species_entries = {best_level_match_obj.name}

                    else:
                        # Retrieve all the species under this specific taxonomic level
                        species_entries = best_level_match_obj.get_children(up_to="species")

                    # Retrieve the species centroids
                    # These are genomes
                    species_centroid = [self.report[self.clusters["species"][species].identifier]["centroid"] for species in species_entries]

                    species_centroid_sketches = [self.genomes[centroid].sketch_filepath for centroid in species_centroid]

                    if first_iter:
                        # Compute the ANI distance between the input genoms and all the species clusters under this specific taxonomic level
                        # At the kingdom level, this is going to select all the species clusters in the database
                        # These distances are shared with the lower levels to avoid computing them again up to seven times
                        _, best_level_match_dists = self.__class__.dist(sketch_filepath, species_centroid_sketches, self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                        dists.update(best_level_match_dists)

                    # Get the distances to the species centroids only
                    species_centroid_dists = [dists[species_centroid_sketch] for species_centroid_sketch in species_centroid_sketches]

                    # Compute the average distance from the centroids
                    best_level_match_distance = statistics.mean(species_centroid_dists)

                    # Keep track of the profile
                    profiles[next_level][best_level_match_taxonomy] = round(best_level_match_distance, 5)

                if first_iter:
                    # Avoid computing ANI distances again
                    first_iter = False

            # Keep querying the closest clusters at the immediate lower taxonomic level
            clusters = list(best_level_matches.keys())

        with open(query_result_filepath, "w+") as profiles_table:
            profiles_table.write("# level\tclosest\tani\n")

            for level in self.__class__.LEVELS + ["genome"]:
                for match in profiles[level]:
                    profiles_table.write(f"{level}\t{match}\t{profiles[level][match]}\n")

        if collapsed:
            # Get rid of the temporary fasta file with the collapsed records
            os.unlink(genome_filepath)

        return profiles

    @classmethod
    def qc(
        cls,
        genomes: Set[os.path.abspath], 
        kingdom: str, 
        nproc: int=os.cpu_count(), 
        tmp: os.path.abspath=os.getcwd()
    ) -> Dict[os.path.abspath, Dict[str, float]]:
        """Perform a quality control based on the genomes' kingdom.
        Genomes from different kingdoms are not allowed to be in the same set.

        Warning: here we suppose that the quality control pipelines and their databases are
                 all available and configured.

        Parameters
        ----------
        genomes : set
            Set with paths to the input genomes.
            Compressed files are not allowed here.
        kingdom : str
            The genomes' kingdom.
            All genomes in the input set of genomes must belong to the same kingdom.
            Possible values: "Viruses", "Bacteria", "Archaea", "Fungi".
        nproc : int, default os.cpu_count()
            Maximum number of CPUs for multiprocessing.
        tmp : os.path.abspath, default os.getcwd()
            Path to the temporary folder.

        Raises
        ------
        Exception
            - If there are no genomes in `genomes`;
            - If the input genomes have an unsupported format;
            - In case of an unexpected error while running a quality control pipeline.
        ValueError
            If the provided `kingdom` is not among the accepted values.
        FileNotFoundException
            - If the path to the tmp folder does not exist;
            - If the output quality summary table does not exist.

        Returns
        -------
        dict
            A series of dictionaries as values of another dictionary indexed by the paths
            to the input genomes, with the quality control statistics of completeness and
            contamination.
        """

        if not genomes:
            raise Exception("There are no genomes to quality control!")

        for genome in genomes:
            # Iterate over the paths to the genome files and check whether their format is supported
            if not cls._is_supported(genome):
                raise Exception(f"Genome file format not supported {genome}")

        if not os.path.isdir(tmp):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), tmp)

        if nproc <= 0 or nproc > os.cpu_count():
            # Use all the available CPUs in case of negative values
            # Otherwise, downscale nproc to the maximum number of CPUs available
            nproc = os.cpu_count()

        # Keep track of the quality control statistics of completeness and contamination
        # for each of the input genomes
        quality = {genome: {"completeness": 0.0, "contamination": 0.0} for genome in genomes}

        if kingdom == "Viruses":
            tmp = os.path.join(tmp, "checkv")

            os.makedirs(tmp, exist_ok=True)

            # CheckV requires a single input file with all the viral sequences in it
            # We should merge all the genomes into a single fna file first
            merged_filepath = os.path.join(tmp, "merged.fna")

            if os.path.isfile(merged_filepath):
                # Get rid of the merged.fna file if it already exists
                os.unlink(merged_filepath)

            # We need to map the contid IDs to the genome files
            contig_to_genomes = dict()

            with open(merged_filepath, "w+") as merged_file:
                for genome in genomes:
                    with open(genome) as genome_file:
                        lines = genome_file.readlines()

                        # Assuming the first line is not empty
                        # Extract the contig ID from the first line
                        contig_id = lines[0].strip()[1:]

                        contig_to_genomes[contig] = genome

                        # Write the genome content into the merged file
                        merged_file.write("".join(lines))

            try:
                # Define the command line to run CheckV
                command_line = [
                    "checkv",
                    "end_to_end",
                    merged_filepath,
                    tmp,
                    "-t",
                    str(nproc)
                ]

                # Execute the command lines
                subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            # The output quality summary table is in the tmp folder
            output_filepath = os.path.join(tmp, "quality_summary.tsv")

            if not os.path.isfile(output_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_filepath)

            # There could be multiple contigs for each of the genomes
            # The final completeness and contamination is the average on all contigs
            contigs = {genome: {"completeness": list(), "contamination": list()} for genome in genomes}

            with open(output_filepath) as output:
                # Retrieve the header line
                header = output_table.readline().strip().split("\t")

                for line in output:
                    line = line.strip()

                    if line:
                        line_split = line.split("\t")

                        # Retrieve the contg ID
                        contig = line_split[header.index("contig")]

                        # Retrieve the genomes from the `contig_to_genome` dictionary
                        genome = contig_to_genome[contig]

                        # Retrieve the completeness
                        completeness = line_split[header.index("completeness")]

                        if completeness != "NA":
                            # Completeness could be NA if CheckV is unable to determine it
                            completeness = 0.0

                        completeness = float(completeness)

                        # Retrieve the contamination
                        contamination = float(line_split[header.index("contamination")])

                        # Keep track of the completeness and contamination stats for this contig
                        contigs[genome]["completeness"].append(completeness)

                        contigs[genome]["contamination"].append(contamination)

            for genome in contigs:
                # Compute the average completeness and contamination in case of multiple contigs
                completeness = statistics.average(contigs[genome]["completeness"])

                contamination = statistics.average(contigs[genome]["contamination"])

                # These are the final completeness and contamination stats
                quality[genome]["completeness"] = completeness

                quality[genome]["contamination"] = contamination

        elif kingdom == "Bacteria" or kingdom == "Archaea":
            tmp = os.path.join(tmp, "checkm")

            # CheckM requires all input genomes to be in a single directory
            genomes_dir = os.path.join(tmp, "genomes")

            os.makedirs(genomes_dir, exist_ok=True)

            for genome in genomes:
                # CheckM requires that the input genomes must all have the same file extension
                # We can force it to be fna
                symlink_filepath = os.path.join(genomes_dir, f"{os.path.splitext(os.path.basename(genome))[0]}.fna")

                if not os.path.isfile(symlink_filepath):
                    os.symlink(genome, symlink_filepath)

            # Also define the path to the output folder
            output_dir = os.path.join(tmp, "output")

            os.makedirs(output_dir, exist_ok=True)

            try:
                # Run CheckM lineage workflow
                command_line = [
                    "checkm",
                    "lineage_wf",
                    "-x",
                    "fna",
                    "-t",
                    str(nproc),
                    genomes_dir,
                    output_dir
                ]

                # Execute the command lines
                subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            # The quality summary file is generated in output_dir
            output_filepath = os.path.join(output_dir, "qa_summary.tsv")

            if not os.path.isfile(output_filepath):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_filepath)

            with open(output_filepath) as output:
                # Retrieve the header line
                header = output.readline().strip().split("\t")

                # Retrieve the column index for bin id, completeness, and contamination
                genome_col = header.index("Bin Id")

                completeness_col = header.index("Completeness")

                contamination_col = header.index("Contamination")

                for line in output:
                    line = line.strip()

                    if line:
                        parts = line.split("\t")

                        # CheckM reports the genome name, so we should recover the whole path
                        bin_id = parts[genome_col]

                        completeness = float(parts[completeness_col])

                        contamination = float(parts[contamination_col])

                        # Reconstruct full path from symlink name
                        original_genome = None

                        for genome in genomes:
                            if os.path.splitext(os.path.basename(genome))[0] == bin_id:
                                original_genome = genome

                                break

                        if original_genome:
                            # Finally, keep track of the genome quality stats
                            quality[original_genome]["completeness"] = completeness

                            quality[original_genome]["contamination"] = contamination

        elif kingdom == "Fungi":
            tmp = os.path.join(tmp, "busco")

            # Define the latest BUSCO database version for Fungi
            # TODO This is hardcoded!
            #      We should automatically point to the latest BUSCO database for Fungi
            busco_db = "fungi_odb10"

            # BUSCO runs over one genome at a time
            for genome in genomes:
                # Define the genome name
                # This is used as the folder name with the busco results in the tmp directory
                genome_name = os.path.splitext(os.path.basename(genome))[0]

                try:
                    # Define the command line to run BUSCO
                    command_line = [
                        "busco",
                        "--in",
                        genome,
                        "--out_path",
                        tmp,
                        "--out",
                        genome_name,
                        "--mode",
                        "genome",
                        "--lineage_dataset",
                        busco_db,
                        "--cpu",
                        str(nproc),
                        "--force",
                        "--metaeuk",
                        "--skip_bbtools",
                        "--opt-out-run-stats",
                        "--quiet"
                    ]

                    # Execute the command lines
                    subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                except subprocess.CalledProcessError as e:
                    error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                # Search for the BUSCO result as JSON file
                busco_json_filepath = os.path.join(tmp, genome_name, f"short_summary.specific.{busco_db}.{genome_name}.json")

                if not os.path.isfile(busco_json_filepath):
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), busco_json_filepath)

                with open(busco_json_filepath) as busco_json:
                    # Load the JSON file
                    busco_results = json.load(busco_json)

                # Compute the completeness and contamination
                completeness = 100.0*(1.0-busco_results["results"]["Missing BUSCOs"]/busco_results["results"]["n_markers"])

                contamination = 100.0*busco_results["results"]["Multi copy BUSCOs"]/busco_results["results"]["Complete BUSCOs"]

                # Finally, keep track of the genome quality stats
                quality[genome]["completeness"] = completeness

                quality[genome]["contamination"] = contamination

        else:
            raise ValueError(f"Invalid kingdom {kingdom}!")

        return quality

    def dereplicate(
        self, 
        genomes: Set[os.path.abspath], 
        threshold: float=0.01, 
        compare_with: str="self"
    ) -> List[os.path.abspath]:
        """Dereplicate a set of genomes versus themselves or versus the genomes in the database.
        The dereplication process is based on their ANI distance according to a specific threshold.

        Parameters
        ----------
        genomes : set
            Set with paths to the genome fasta files.
        threshold : float, default 0.01
            Dereplicate genomes according to this threshold.
            The dereplication process is triggered if `threshold` is >0.0.
        compare_with : str, default "self"
            Dereplicate the input genomes versus themselves or versus the database, using "self" or "database" respectively.
            Possible values: "self", "database".

        Raises
        ------
        Exception
            If there are no genomes in `genomes`.
        ValueError
            - If `threshold` is <0.0 or >=1.0;
            - In case `compare_with` is difference than "self" and "database".

        Returns
        -------
        set
            A set with paths to the genomes that passed the dereplication processed.
            All the other genomes are discarded.
        """

        if not genomes:
            raise Exception("There are no genomes to dereplicate!")

        if threshold < 0.0 or threshold >= 1.0:
            # ANI distances are always between 0.0 and 1.0
            # An ANI distance equals to 1.0 would discard all the input genomes
            raise ValueError("The ANI threshold must be >0.0 and <1.0!")

        if compare_with not in ["self", "database"]:
            # This argument is used to select the genomes against with the dereplication process is performed
            # It can be run against the input genomes themselves using "self", or against the genomes in the database with "database"
            raise ValueError("The dereplication can be performed against the input itself or the genomes in the database!")

        sketches = list()

        # Map genome names to the input file paths
        names = dict()

        for genome_filepath in genomes:
            # Define the input file name
            filename = os.path.splitext(os.path.basename(genome_filepath))[0]

            names[filename] = genome_filepath

            # Assume the input genomes are all in fasta format
            genome_obj = Entry(self, filename, filename, "genome")

            # Build their bloom filter sketch representation
            genome_sketch_filepath = genome_obj.sketch(genome_filepath)

            sketches.append(genome_sketch_filepath)
    
        # Rescale nproc
        nproc = self.nproc if len(sketches) > self.nproc else len(sketches)

        # Keep track of replicas
        # A genome is considered a replica if there exists another genome at a genetic distance < threshold
        # key=replica, value={dereplicated: distance}
        replicas = dict()

        # Dereplicate the input genomes versus themselves
        if compare_with == "self":
            if nproc > 1:
                with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(sketches)) as pbar:
                    # Wrapper around the update function of tqdm
                    def progress(*args):
                        pbar.update()

                    jobs = [
                        pool.apply_async(
                            self.__class__.dist, 
                            args=(
                                sketch_filepath, 
                                sketches[pos+1:],
                                self.metadata["kmer_size"],
                                self.tmp,
                                False,
                            ),
                            callback=progress
                        ) for pos, sketch_filepath in enumerate(sketches)
                    ]

                    for job in jobs:
                        sketch_filepath, sketch_dists = job.get()

                        # Select sketch replicas
                        sketch_clones = {os.path.splitext(os.path.basename(sketch_target))[0]: sketch_dist for sketch_target in sketch_dists if sketch_dist <= threshold}

                        if sketch_clones:
                            # Define the input file name
                            sketch_filename = os.path.splitext(os.path.basename(sketch_filepath))[0]

                            # Keep track of replicas
                            replicas[sketch_filename] = sketch_clones

            else:
                # Avoid using multiprocessing if `nproc` is 1
                for pos, sketch_filepath in enumerate(sketches):
                    _, sketch_dists = self.__class__.dist(sketch_filepath, sketches[pos+1:], self.metadata["kmer_size"], tmp=self.tmp, resume=False)

                    # Select sketch replicas
                    sketch_clones = {os.path.splitext(os.path.basename(sketch_target))[0]: sketch_dist for sketch_target in sketch_dists if sketch_dist <= threshold}

                    if sketch_clones:
                        # Define the input file name
                        sketch_filename = os.path.splitext(os.path.basename(sketch_filepath))[0]

                        # Keep track of replicas
                        replicas[sketch_filename] = sketch_clones

        elif compare_with == "database":
            if nproc > 1:
                # Profile genomes in parallel
                with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(sketches)) as pbar:
                    def progress(*args):
                        pbar.update()

                    # TODO Uncertainty and pruning threshold shouldn't be hardcoded
                    # The uncertainty can be very low here since we are searching for replicas
                    jobs = [
                        pool.apply_async(
                            self.__class__._profile, 
                            args=(
                                self,
                                genome_filepath,
                                sketch_filepath, 
                                1.0,  # uncertainty
                                0.0,  # pruning threshold
                            ),
                            callback=progress
                        ) for genome_filepath, sketch_filepath in zip(genomes, sketches)
                    ]

                    for job in jobs:
                        genome_filepath, genome_profile = job.get()

                        # The profile function report the first closest genome only
                        # The name of the closest genome is reported under the t__ level
                        closest_genome = list(genome_profile["genome"].keys())[0].split("|")[-1][3:]

                        # Retrieve the ANI distance with the closest genome
                        closest_genome_distance = genome_profile["genome"][closest_genome]

                        if closest_genome_distance <= threshold:
                            # Define the input file name
                            genome_filename = os.path.splitext(os.path.basename(genome_filepath))[0]

                            if closest_genome not in replicas:
                                replicas[closest_genome] = dict()

                            # Keep track of the replica in the database
                            # It doesn't matter if `closest_genome` is a full path here
                            replicas[closest_genome][genome_filename] = closest_genome_distance

            else:
                # Avoid using multiprocessing if `nproc` is 1
                for genome_filepath, sketch_filepath in zip(genomes, sketches):
                    # Again, the uncertainty can be very low here since we are searching for replicas
                    _, genome_profile = self.__class__._profile(self, genome_filepath, sketch_filepath, 1.0, 0.0)

                    # The profile function report the first closest genome only
                    # The name of the closest genome is reported under the t__ level
                    closest_genome = list(genome_profile["genome"].keys())[0].split("|")[-1][3:]

                    # Retrieve the ANI distance with the closest genome
                    closest_genome_distance = genome_profile["genome"][closest_genome]

                    if closest_genome_distance <= threshold:
                        # Define the input file name
                        genome_filename = os.path.splitext(os.path.basename(genome_filepath))[0]

                        if closest_genome not in replicas:
                            replicas[closest_genome] = dict()

                        # Keep track of the replica in the database
                        # It doesn't matter if `closest_genome` is a full path here
                        replicas[closest_genome][genome_filename] = closest_genome_distance

        # Define the set of excluded genomes based on their ANI distance
        excluded = set()

        for replica in replicas:
            # Retrieve the input path to the fasta file or the bloom filter file
            # In case of input versus database, `sketch_filename` does not exist in the input set of genomes
            input_filepath = names.get(replica, None)

            if input_filepath not in excluded:
                # Exclude all the replicas to the current genome
                excluded = excluded.union({names[sketch_filename] for sketch_filename in replicas[input_filepath].keys()})

        # Compute the difference between the input set of genomes and the excluded ones
        return genomes.difference(excluded)

    def merge(self, clusters: Set[str], into: str=None, update: bool=True) -> None:
        """Merge two or more species clusters.

        Parameters
        ----------
        clusters : Set
            Set of MetaSBT cluster IDs.
        into : str, default None
            Merge clusters in `clusters` into the cluster specified under `into`.
            Note that this MetaSBT cluster ID must be also present in `clusters`.
            If None, clusters are going to be merged into the oldest cluster in `clusters`.
        update : bool, default True
            The merging of clusters leads to the removal of the merged ones.
            By default, each call to this function will trigger the database update.
            If False, remember to call the `update()` function manually once you are done merging clusters.

        Raises
        ------
        Exception
            - If one of the provided clusters is not in the database at the species level;
            - If `into` is not in `clusters`.
        """

        for cluster in clusters:
            if not (cluster in self.report and self.report[cluster]["level"] == "species"):
                # Check whether the provided clusters are in the database at the species level
                raise Exception(f"Cluster {cluster} is not a species!")

        if not into:
            # Select the oldest cluster in `clusters` if `into` is not specified
            into = f"MSBT{min([int(cluster[4:]) for cluster in clusters])}"

        if into not in clusters:
            raise Exception(f"The {into} cluster is not in the list of clusters!")

        # Retrieve the `into` cluster name
        into_name = self.report[into]["taxonomy"].split("|")[-1]

        # Retrieve the `into` species cluster object
        into_obj = self.clusters["species"][into_name]

        for cluster in clusters:
            if cluster != into:
                # Retrieve the cluster name
                cluster_name = self.report[cluster]["taxonomy"].split("|")[-1]

                # Define the cluster Entry object
                cluster_obj = self.clusters["species"][cluster_name]

                # Get cluster genomes
                cluster_genomes = cluster_obj.get_children()

                for genome in cluster_genomes:
                    # Redirect the genome parent to the `into` cluster
                    self.genomes[genome].parent = into_name

                    # Add genomes to the `into` cluster
                    into_obj.add_children(genome)

                # Empty the set of children for the current cluster
                cluster_obj.children = set()

        # We are going to remove `into` from `cluster`
        # Make a copy of `clusters` to avoid affecting the input set
        clusters = copy.deepcopy(clusters)

        # Exclude the `into` cluster from the list of clusters that must be removed
        clusters.remove(into)

        # Add the `into` full taxonomic label to the list of clusters that must be indexed
        self.__clusters.append(into_obj.get_full_taxonomy())

        # We could use the `remove()` function to get rid of the merged cluster excluding the `into` one
        # This takes care of removing branches eventually, and indexing clusters
        # Do not remove the sketch representation of genomes `keep_sketches=True`
        # Genomes objects are not going to be removed because we redirected their parents to `into` earlier
        self.remove(species=clusters, keep_sketches=True, update=update)

    def remove(
        self, 
        genomes: Set[str]=None, 
        species: Set[str]=None, 
        keep_sketches: bool=False, 
        update: bool=True
    ) -> None:
        """Remove genomes and clusters.

        Parameters
        ----------
        genomes : set, default None
            Set of genome names to be removed from the database.
            Note that this could potentially lead to the removal of clusters if they remain empty.
        species : Set, default None
            Set of MetaSBT species cluster IDs to be removed from the database.
            The genomes assigned to these clusters are also removed from the database.
        keep_sketches : bool, default False
            Keep the sketch representation of genomes.
        update : bool, default True
            By default, each call to this function will trigger the database update.
            If False, remember to call the `update()` function manually once you are done merging clusters.

        Raises
        ------
        Exception
            If no genomes and no clusters are provided.
        """

        def remove_genomes(genomes: Set[str], keep_sketches: bool=False) -> Set[str]:
            """Remove a set of genomes from the database and report the empty clusters if any.

            Parameters
            ----------
            genomes : set
                Set of genome names to be removed from the database.
            keep_sketches : bool, default False
                Keep the sketch representation of genomes.

            Returns
            -------
            set
                A set of MetaSBT species cluster IDs that become empty due to the removal of genomes.
            """

            # Keep track of the clusters that become empy due to the removal of genomes
            empty_clusters = set()

            if genomes:
                # Keep track of clusters that must be indexed
                index_clusters = set()

                for genome in genomes:
                    # Check whether these genomes are in the database
                    if genome not in self.genomes:
                        # Just skip a genome if it's not in the database
                        # Move to the next genome
                        continue

                    # Retrieve the genome Entry object
                    genome_obj = self.genomes[genome]

                    if not keep_sketches:
                        # Remove the sketch
                        os.unlink(genome_obj.sketch_filepath)

                    # Retrieve the species where the genome resides in the database
                    species_obj = self.clusters["species"][genome_obj.parent]

                    # Remove the genome from the set of children of its species
                    species_obj.children.remove(genome)

                    # We want to keep the genome object if we keep the sketch representation of the genome
                    if not keep_sketches:
                        # Remove the genome from the instance dictionary of genomes
                        del self.genomes[genome]

                    # Check whether the species has no more children
                    if not species_obj.children:
                        # This cluster should be removed
                        empty_clusters.add(species_obj.identifier)

                        if genome_obj.parent in index_clusters:
                            index_clusters.remove(genome_obj.parent)

                    else:
                        # We should fix the SBT definition file of the species cluster
                        # Index the species cluster
                        index_clusters.add(genome_obj.parent)

                if index_clusters:
                    # Index clusters to fix the SBT definition file
                    for sp in index_clusters:
                        # Define the species cluster Entry object
                        sp_obj = self.clusters["species"][sp]

                        # Index the whole branch and the upper taxonomic levels
                        self.__clusters.append(sp_obj.get_full_taxonomy())

            return empty_clusters

        def remove_clusters(clusters: Set[str], level: str="species", keep_sketches: bool=False) -> None:
            """Remove a set of clusters under a specific taxonomic level and report
            their parents if they become empty.

            Parameters
            ----------
            clusters : set
                Set of MetaSBT cluster IDs to be removed from the database.
            level : str, default "species"
                The taxonomic level.
            keep_sketches : bool, default False
                Keep the sketch representation of genomes.
            """

            if clusters:
                # Keep track of clusters that must be indexed
                # key=cluster, value=level
                index_clusters = dict()

                stack = [(cluster, level) for cluster in clusters]

                while stack:
                    cluster, level = stack.pop()

                    # Check whether this cluster exists in the database
                    if not (cluster in self.report and self.report[cluster]["level"] == level):
                        # Just skip a cluster if it's not in the database
                        # This cluster could actually be in the database but it's not under this level
                        # Move to the next cluster
                        continue

                    # Retrieve the cluster name
                    cluster_name = self.report[cluster]["taxonomy"].split("|")[-1]

                    # Retrieve the cluster Entry object
                    cluster_obj = self.clusters[level][cluster_name]

                    # Retrieve the set of genomes assigned under this cluster
                    genomes = cluster_obj.get_children()

                    if genomes:
                        # Remove the genomes under this cluster first
                        # This is not necessarily a species cluster
                        _ = remove_genomes(genomes, keep_sketches=keep_sketches)
                    
                    # What if this cluster was the only one assigned to its higher level cluster?
                    # We should remove clusters all the way up to the kingdom level
                    # Retrieve the parent taxonomic level first
                    parent_level = self.LEVELS[self.LEVELS.index(level)-1] if level != "kingdom" else None

                    if parent_level:
                        # Also retrieve the parent object
                        parent_obj = self.clusters[parent_level][cluster_obj.parent]

                        # Remove the cluster from the set of its parent children
                        parent_obj.children.remove(cluster_name)

                        if not parent_obj.children:
                            # If there are no other children
                            # Remove the parent cluster
                            stack.append((parent_obj.identifier, parent_level))

                            if parent_obj.name in index_clusters:
                                del index_clusters[parent_obj.name]

                        else:
                            # We should fix the parent SBT definition file
                            # Index the parent cluster
                            index_clusters[parent_obj.name] = parent_level

                    # We should now remove the cluster folder
                    # Retrieve the cluster batch id
                    cluster_batch = self.__class__._get_cluster_batch(cluster_obj.identifier)

                    # Define the path to the cluster folder
                    cluster_dir = os.path.join(self.root, "clusters", cluster_batch, cluster_obj.identifier)

                    # Remove the cluster folder
                    shutil.rmtree(cluster_dir, ignore_errors=True)

                    # It is safe to remove the cluster now
                    del self.clusters[level][cluster_name]

                    # Also remove the cluster from the report
                    del self.report[cluster]

                if index_clusters:
                    # Index clusters to fix the SBT definition file
                    for cluster in index_clusters:
                        # Define the cluster Entry object
                        cluster_obj = self.clusters[index_clusters[cluster]][cluster]

                        # Retrieve the species cluster with the smaller number of genomes under `cluster` and use it to index the whole branch
                        # Use the cluster with the smaller number of genomes to speed up the indexing
                        # Retrieve the set of species under `cluster` first
                        species = list(cluster_obj.get_children(up_to="species"))

                        # Sort species according to their number of genomes
                        species = sorted(species, key=lambda sp: len(self.clusters["species"][sp].children))

                        # Remove clusters that appear in the input set of clusters
                        species = [sp for sp in species if self.clusters["species"][sp].identifier not in clusters]

                        # There must be at least one cluster in `species`
                        # Retrieve the species cluster Entry object
                        species_obj = self.clusters["species"][species[0]]

                        # Retrieve the species full taxonomic label
                        species_taxonomy = species_obj.get_full_taxonomy()

                        # Index the species cluster and its whole branch
                        self.__clusters.append(species_taxonomy)

        if not genomes and not species:
            raise Exception("Please provide at least one genome or one cluster!")

        # The input parameter `species` is called that way to instruct users to using species clusters only
        # However, the nested function `remove_clusters` is able to remove clusters at all the taxonomic levels
        # Just rename `species` to `clusters`
        clusters = species

        # Remove genomes first
        # This could return a set with MetaSBT species cluster IDs that become empty due to the removal of genomes
        empty_clusters = remove_genomes(genomes, keep_sketches=keep_sketches)

        if empty_clusters:
            if not clusters:
                # Initialize the set of species clusters if None
                clusters = set()

            # Add the empty species clusters to the set of clusters that must be removed
            clusters = clusters.union(empty_clusters)

        remove_clusters(clusters, level="species", keep_sketches=keep_sketches)

        # Fix `self.__clusters`
        # Taxonomies added in self._clusters during a specific iteration could be removed during the next iteration
        # We should remove redundant partial branches
        branches = set()

        for taxonomy in self.__clusters:
            # Split the taxonomy
            taxonomy_split = taxonomy.split("|")

            # Iterate over the taxonomic levels
            for taxonomy_level, level_name in zip(taxonomy_split, self.LEVELS):
                # Check whether this level still exists in the database
                # If not, we should cut it and build a new taxonomy up to an existing species (the smallest one under the partial taxonomy)
                if taxonomy_level not in self.clusters[level_name]:
                    # Retrieve the name of the previous cluster
                    # This cluster exists in the database
                    cluster_name = taxonomy_split[taxonomy_split.index(taxonomy_level)-1]

                    # Consider the previous level and retrieve the cluster Entry object
                    cluster_obj = self.clusters[self.LEVELS[self.LEVELS.index(level_name)-1]][cluster_name]

                    # Retrieve all the species under `cluster_obj`
                    species = list(cluster_obj.get_children(up_to="species"))

                    # Sort species according to their number of genomes
                    species = sorted(species, key=lambda sp: len(self.clusters["species"][sp].children))

                    # There must be at least one cluster in `species`
                    # Retrieve the species cluster Entry object
                    species_obj = self.clusters["species"][species[0]]

                    # Retrieve the species full taxonomic label
                    taxonomy = species_obj.get_full_taxonomy()

                    break

            # Register the taxonomy
            branches.add(taxonomy)

        self.__clusters = list(branches)

        if self.__clusters and update:
            # The removal of genomes and/or clusters would trigger the update
            # This will index the affected clusters again and rebuild the MetaSBT report files
            self.update()

    @classmethod
    def _is_defined(cls, taxonomy: str) -> bool:
        """Check whether the input taxonomic label is fully defined from the kingdom up to the species level.

        Parameters
        ----------
        taxonomy : str
            The taxonomic label.

        Raises
        ------
        ValueError
            If the provided taxonomy is None;

        Returns
        -------
        bool
            True if `taxonomy` is defined from the kingdom up to the species level for a total of seven taxonomic levels.
            False otherwise.
        """

        if not taxonomy:
            raise ValueError("No taxonomy provided!")

        if taxonomy.count("|") != 6:
            # Valid taxonomic labels must be defined at all the seven taxonomic levels
            # From the kingdom down to the species level
            return False

        else:
            for level, level_id in zip(cls.LEVELS, taxonomy.split("|")):
                if len(level_id) <= 3 or level_id[:3] != f"{level[0]}__":
                    return False

        return True

    @classmethod
    def _format_taxonomy(cls, taxonomy: str) -> str:
        """Format a taxonomic label.

        Parameters
        ----------
        taxonomy : str
            The input taxonomic label.

        Raises
        ------
        ValueError
            If the provided taxonomy is malformed.
            See `_is_defined` function for additional details.

        Returns
        -------
        str
            The new formatted taxonomic label.
        """

        if not cls._is_defined(taxonomy):
            raise ValueError(f"Malformed taxonomic label: {taxonomy}")

        def format_level(current_level: str, prev_level: str):
            """Define a taxonomic level name.

            Parameters
            ----------
            current_level : str
                Current level name.
            prev_level : str
                Previous level name in case of unclassified.

            Returns
            -------
            str
                The new level name.
            """

            # Remove special characters from current and previous level names
            current_level = re.sub(r"_+", "_", re.sub(r"\W+", "_", current_level)).strip("_")

            prev_level = re.sub(r"_+", "_", re.sub(r"\W+", "_", prev_level)).strip("_")

            # Build the new level name
            level_prefix = current_level.strip()

            level_suffix = ""

            if not level_prefix:
                level_prefix = prev_level

                # Fill empty taxa levels with unclassified
                level_suffix = "_unclassified"

            return f"{level_prefix}{level_suffix}"

        taxonomic_levels = {cls.LEVELS[pos]: taxonomic_level[3:] for pos, taxonomic_level in enumerate(taxonomy.split("|"))}

        # Build the new taxonomic label
        taxonomy = "|".join([f"{level[0]}__{taxonomic_levels[level]}" for level in cls.LEVELS])

        return taxonomy

    @staticmethod
    def _is_supported(filepath: os.path.abspath) -> bool:
        """Check whether an input file is supported base on its extension.

        Parameters
        ----------
        filepath : os.path.abspath
            The input file path.

        Raises
        ------
        FileNotFoundError
            If the input file `filepath` does not exist.

        Returns
        -------
        bool
            True if the input file is supported, False otherwise.
        """

        if not os.path.isfile(filepath):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

        # Define the set of supported extensions
        supported = {".fa", ".fasta", ".fna", ".bf"}

        return os.path.splitext(filepath)[1] in supported

    @staticmethod
    def _load_report(filepath: os.path.abspath) -> Dict[str, Any]:
        """Load the report table with the list of clusters in the database.

        Parameters
        ----------
        filepath : os.path.abspath
            Path to the report table.

        Raises
        ------
        FileNotFoundError
            If the input file `filepath` does not exist.

        Returns
        -------
        dict
            A dictionary with the list of clusters and their information.
        """

        if not os.path.isfile(filepath):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)

        report_table = dict()

        with open(filepath) as table:
            for line in table:
                if line.strip():
                    if line.startswith("#"):
                        # The last commented line is the header
                        header = line.strip().split("\t")

                        # Get rid of the pound char
                        header[0] = header[0][1:].strip()

                    else:
                        line_split = line.split("\t")

                        # Retrieve the cluster id
                        cluster_id = line_split[header.index("cluster_id")]

                        # The report table contains information about clusters at all the seven taxonomic levels,
                        # from the kingdom (at least one cluster) up to the species level
                        report_table[cluster_id] = dict()

                        # Retrieve the cluster level
                        report_table[cluster_id]["level"] = line_split[header.index("level")]

                        # Retrieve the cluster density
                        report_table[cluster_id]["density"] = float(line_split[header.index("density")])
    
                        # Retrieve the set of reference genomes
                        report_table[cluster_id]["references"] = set()

                        # There could be no references in a cluster
                        if line_split[header.index("references_list")].strip():
                            report_table[cluster_id]["references"] = set(line_split[header.index("references_list")].split(","))

                        # Retrieve the set of metagenome-assembled genomes
                        report_table[cluster_id]["mags"] = set()

                        if line_split[header.index("mags_list")].strip():
                            report_table[cluster_id]["mags"] = set(line_split[header.index("mags_list")].split(","))

                        # Retrieve the cluster centroid
                        report_table[cluster_id]["centroid"] = line_split[header.index("centroid")]

                        report_table[cluster_id]["known"] = line_split[header.index("known")]

                        # Retrieve the taxonomic assignment
                        report_table[cluster_id]["taxonomy"] = line_split[header.index("taxonomy")]

                        report_table[cluster_id]["internal"] = line_split[header.index("internal")]

                        # Retrieve the cluster boundaries
                        min_ani = None

                        if line_split[header.index("min_ani")].strip():
                            min_ani = float(line_split[header.index("min_ani")])

                        max_ani = None

                        if line_split[header.index("max_ani")].strip():
                            # `max_ani` is the last field in the clusters table
                            # strip it to get rid of the new line character
                            max_ani = float(line_split[header.index("max_ani")].strip())

                        report_table[cluster_id]["boundaries"] = (min_ani, max_ani)

        return report_table

    @staticmethod
    def _validate_metadata(metadata: Dict[str, Any]) -> bool:
        """Validate database metadata.

        Parameters
        ----------
        metadata : dict
            The metadata dictionary that must be validated.

        Returns
        -------
        bool
            True if metadata pass the validation, False otherwise.
        """

        # Define the set of attributes that a database is required to have
        required_attributes = {"clusters_count", "filter_size", "kmer_size", "min_kmer_occurrence", "flat"}

        return len(set(metadata.keys()).intersection(required_attributes)) == len(required_attributes)


class Entry(object):
    """Entry object representing an a specific taxonomic level and genomes."""

    def __init__(
        self,
        database: "Database",
        identifier: str,
        name: str,
        level: str,
        folder: Optional[os.path.abspath]=None,
        parent: Optional[str]=None,
        children: Optional[Set[str]]=None,
        taxonomy: Optional[str]=None,
        known: bool=False,
    ) -> "Entry":
        """Initialize an Entry object.

        Parameters
        ----------
        identifier : str
            The entry ID.
            This can be the same as `name` in case of `level="genome"`.
            Otherwise, it is usually an incremental numerical ID with the "MSBT" prefix.
        name : str
            The entry name.
        level : str
            Taxonomic level.
            Possible values: {'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'genome'}
        folder : os.path.abspath, optional, default None
            Path to the entry root folder.
        parent : str, optional, default None
            The name of the parent entry.
        children : set, optional, default None
            The set of children names.
        taxonomy : str, optional, default None
            The taxonomic label, used in case `level="genome"`
        known : bool, default False
            Clusters should always be initialized as `known=False`.
            This is used particularly in case of genomes. Thus, if it is known, the genome is a reference and must have a specific taxonomic label.

        Raises
        ------
        ValueError
            - If `identifier` is None;
            - If `name` is None;
            - If `level` is None or it is not one of the following values: kingdom, phylum, class, order, family, genus, species, and genome;
            - If `level` is `genome` and `children` is not empty;
            - If `level` is `genome` and `identifier` does not match with `name`;
            - If `level` is not `genome` and `identifier` does not match with the regular expression `^MSBT[0-9]*$`;
            - If the current entry is a genome, it is a reference, but `taxonomy` is not provided;
            - If the current entry is a genome, it is a reference, and the provided `taxonomy` is malformed.

        Returns
        -------
        Entry
            A new Entry object.
        """

        self.database = database

        if not identifier:
            raise ValueError("Invalid Entry ID!")

        self.identifier = identifier

        if not name:
            raise ValueError("Invalid Entry name!")

        self.name = name

        if not level:
            # It could be MSBT0
            pass

        elif level not in self.database.__class__.LEVELS + ["genome"]:
            raise ValueError("Invalid taxonomic level!")

        if level == "genome" and children:
            raise ValueError("A genome cannot contain children!")

        self.level = level

        if self.level == "genome" and not self.identifier == self.name:
            raise ValueError("A genome entry ID should match with the entry name!")

        elif self.level in self.database.__class__.LEVELS and not re.search("^MSBT[0-9]*$", self.identifier):
            raise ValueError("A cluster entry ID should match with the following regular expression: ^MSBT[0-9]*$")

        if not folder:
            # Use the current working directory if folder is not provided
            folder = os.path.join(os.getcwd(), name)

        self.folder = folder

        if not self.level or self.level in self.database.__class__.LEVELS:
            # Do not create a folder in case of genomes.
            os.makedirs(os.path.join(self.folder, "tree"), exist_ok=True)

        # The immediate higher taxonomic level
        self.parent = parent

        # The immediate lower taxonomic levels or genomes
        # This contains the genome names in case the entry represents a species
        self.children = children if children else set()

        # An Entry node is known if it contains at least one reference genome by definition
        self.__known = known

        # Set the taxonomic label
        # This is always None except in case of genomes, especially references
        self.taxonomy = taxonomy

        if self.level == "genome" and self.__known and not self.taxonomy:
            raise ValueError("The entry is a reference but no taxonomy is provided!")

        if self.level == "genome" and self.__known and (self.taxonomy and not self.database.__class__._is_defined(self.taxonomy)):
            raise ValueError("The provided taxonomy is malformed!")

        # Used to temporary keep track of the genome profile indexed by the taxonomic level
        # Only in case of `level="genome"`
        self.profile = dict()

        self.sketch_filepath = None

        if not self.level or self.level in self.database.__class__.LEVELS:
            sketch_filepath = os.path.join(self.database.root, "clusters", self.database.__class__._get_cluster_batch(self.identifier), self.identifier, f"{self.name}.bf")

        else:
            # This is a genome
            # Search whether a sketch representation of this entry already exists in the database
            sketch_filepath = os.path.join(self.database.root, "sketches", f"{self.name}.bf")

        if os.path.isfile(sketch_filepath):
            self.sketch_filepath = sketch_filepath

    def __len__(self) -> int:
        """Get the Entry size.

        Returns
        -------
        int
            The Entry size in terms of number of children.
        """

        return len(self.children)

    def __str__(self) -> None:
        """Print the Entry object properties.

        Returns
        -------
        str
            A description of the Entry object.
        """

        return f"""
            Class:    metasbt.objects.Entry
            Version:  {database.version}
            Name:     {self.name}
            Level:    {self.level}
            Known:    {self.is_known()}
            Parent:   {self.parent}
            Children: {len(self.children)}
        """

    def get_children(self, up_to: str=None) -> Set[str]:
        """Get the set of genomes under a specific cluster or the set of clusters under a specific taxonomic level.

        Parameters
        ----------
        up_to : str, default None
            Up to a specific taxonomic level.

        Raises
        ------
        Exception
            If the current entry represents a genome.

        Returns
        -------
        set
            The set of genomes under a specific cluster if `up_to` is None,
            otherwis the set of clusters under a specific taxonomic level.
        """

        levels = self.database.__class__.LEVELS + ["genome"]

        if not up_to:
            # Return genomes by default
            up_to = "genome"

        else:
            if up_to not in levels:
                raise Exception(f"{up_to} is not a valid taxonomic level!")

            elif levels.index(up_to) <= levels.index(self.level):
                raise Exception(f"There are no {up_to} levels under the {self.level} level!")

        # Keep track of children
        children = set()

        stack = [self]

        while stack:
            current = stack.pop()

            # Genomes do not have children
            if current.level == "genome":
                raise Exception("The current entry represents a genome!")

            # If we reached the species level, return children
            if current.level == "species":
                children.update(current.children)

                # We reached the end of the branch
                continue

            next_level = levels[levels.index(current.level)+1]

            # If the next level is the desired level, return children directly
            if up_to and next_level == up_to:
                children.update(current.children)

                # We reached the `up_to` stop condition
                continue

            for child in current.children:
                # Keep processing new entries
                stack.append(self.database.clusters[next_level][child])

        return children

    def get_density(self) -> float:
        """Retrieve the density of the bloom filter representation of a cluster.

        Returns
        -------
        float
            The cluster density.
        """

        density = 0.0

        with tempfile.NamedTemporaryFile() as dump_bloom_filter:
            command_line = [
                "howdesbt",
                "dumpbf",
                os.path.join(self.folder, f"{self.name}.bf"),
                "--show:density"
            ]

            try:
                with open(dump_bloom_filter.name, "w+") as dump_bloom_filter_file:
                    subprocess.check_call(command_line, stdout=dump_bloom_filter_file, stderr=dump_bloom_filter_file)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                raise Exception(error_message).with_traceback(e.__traceback__)

            with open(dump_bloom_filter.name) as dump_bloom_filter_file:
                density = float(dump_bloom_filter_file.readline().strip().split(" ")[-1])

        return density

    def get_full_taxonomy(self, current_entry: "Entry"=None, taxonomy: str=None, internal: bool=False) -> str:
        """Recursively define the full taxonomic label based on parents.

        Parameters
        ----------
        current_entry : Entry, default None
            Current entry object.
        taxonomy : str, default None
            Taxonomic label.
        internal : bool, default False
            Use the internal MSBT IDs if True.

        Returns
        -------
        str
            The full taxonomic label.
        """

        entry = current_entry if current_entry else self

        if not entry.level:
            return ""

        if entry.level == "genome":
            # Move to the species level
            species_obj = self.database.clusters["species"][entry.parent]

            return entry.get_full_taxonomy(current_entry=species_obj, internal=internal)

        entry_name = f"{entry.level[0]}__{entry.identifier}" if internal else entry.name

        taxonomy = f"{entry_name}|{taxonomy}" if taxonomy else entry_name

        if entry.level == "kingdom":
            # There is nothing after this level
            # Stop the recursion
            return taxonomy

        levels = self.database.__class__.LEVELS + ["genome"]

        parent_level = levels[levels.index(entry.level)-1]

        parent_obj = self.database.clusters[parent_level][entry.parent]

        return entry.get_full_taxonomy(current_entry=parent_obj, taxonomy=taxonomy, internal=internal)

    def add_children(self, child: str) -> None:
        """Add a children.

        Parameters
        ----------
        child : str
            The child name.

        Raises
        ------
        Exception
            In case the current entry is a genome.
        """

        if self.level == "genome":
            raise Exception("Genomes cannot have children!")

        self.children.add(child)

    def is_known(self) -> bool:
        """Check whether this Entry is known or unknown.
        By definition, a entry is known if it contains at least one reference genome.

        Returns
        -------
        bool
            True if it is known, False otherwise.
        """

        if not self.level:
            return False

        if self.level == "genome":
            return self.__known

        if self.__known:
            # If it has already been marked as known, just report it
            # An unknown entry could become known, but a known entry could not become unknown
            return True

        entries = self.children

        if self.level != "species":
            for level in self.database.__class__.LEVELS[self.database.__class__.LEVELS.index(self.level)+1:]:
                next_entries = set()

                for entry in entries:
                    next_entries = next_entries.union(self.database.clusters[level][entry].children)

                entries = next_entries

        for genome in entries:
            self.__known = genome in self.database.genomes and self.database.genomes[genome].is_known()

            if self.__known:
                return True

        return False

    def index(self) -> str:
        """Build the Sequence Bloom Tree.
        It overwrites the index if it already exists.

        Raises
        ------
        Exception
            - In case no database metadata is found;
            - In case the current entry is a genome;
            - In case the entry does not contain any children (empty);
            - In case of an unexpected error while performing clustering with HowDeSBT;
            - In case of an unexpected error while building the Sequence Bloom Tree with HowDeSBT;
            - In case of an unexpected error while building the entry representative with HowDeSBT.
        """

        if not self.database.__class__._validate_metadata(self.database.metadata):
            raise Exception("No database metadata found!")

        if self.level == "genome":
            raise Exception("Genome entries cannot be indexed!")

        if not self.children:
            raise Exception("This entry is empty!")

        # Delete the index folder if it already exists
        if os.path.isdir(self.folder):
            shutil.rmtree(os.path.join(self.folder, "tree"), ignore_errors=True)

            # Recreate the cluster folder
            os.makedirs(os.path.join(self.folder, "tree"), exist_ok=True)

        # Colect paths to the children sketches
        sketches = set()

        if self.level == "species":
            # This entry is a species
            # We should search for Sequence objects here
            for child in self.children:
                # Assume the sketch file exists
                sketches.add(self.database.genomes[child].sketch_filepath)

        else:
            # Children are located at the next higher taxonomic level
            next_level = "kingdom" if not self.level else self.database.__class__.LEVELS[self.database.__class__.LEVELS.index(self.level)+1]

            # This entry represent a taxonomic level, different than a species
            # Thus, children are Entry objects
            for child in self.children:
                # Retrieve the Entry object from the database
                entry_obj = self.database.clusters[next_level][child]

                # The bloom filter representation of the entry is located in its folder
                entry_sketch = os.path.join(entry_obj.folder, f"{entry_obj.name}.bf")

                # Assume the sketch file exists
                sketches.add(entry_sketch)

        # Dump the list of sketches to file
        sketches_list_filepath = os.path.join(self.folder, f"{self.name}.txt")

        with open(sketches_list_filepath, "w+") as file:
            for sketch_filepath in sketches:
                file.write(f"{sketch_filepath}\n")

        # Perform the clustering first
        union_tree_filepath = os.path.join(self.folder, "tree", "union.sbt")

        if len(sketches) == 1:
            # It does not make sense to perform a clustering with one sketch only
            single_sketch_filepath = list(sketches)[0]

            shutil.copy(single_sketch_filepath, os.path.join(self.folder, f"{self.name}.bf"))

            # Manually define the union.sbt file with a single node only
            with open(union_tree_filepath, "w+") as file:
                file.write(f"{single_sketch_filepath}\n")

        if not self.database.flat:
            # Keep track of the original working directory
            current_working_directory = os.getcwd()

            # Move to the index folder
            # This will force howdesbt to build the compressed nodes into the index folder
            os.chdir(os.path.join(self.folder, "tree"))

            if len(sketches) > 1:
                command_line = [
                    "howdesbt",
                    "cluster",
                    f"--list={sketches_list_filepath}",
                    f"--bits={self.database.metadata['filter_size']}",
                    f"--tree={union_tree_filepath}",
                    f"--nodename={os.path.join(self.folder, 'tree', 'node{number}')}",
                    "--keepallnodes"
                ]

                try:
                    subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                except subprocess.CalledProcessError as e:
                    error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                    # Move back to the original working directory
                    os.chdir(current_working_directory)

                    raise Exception(error_message).with_traceback(e.__traceback__)

            # Build the Sequence Bloom Tree
            command_line = [
                "howdesbt",
                "build",
                "--howde",
                f"--tree={union_tree_filepath}".format(union_tree_filepath),
                f"--outtree={os.path.join(self.folder, 'tree', 'index.detbrief.sbt')}"
            ]

            try:
                subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                # Move back to the original working directory
                os.chdir(current_working_directory)

                raise Exception(error_message).with_traceback(e.__traceback__)

            # Move back to the original working directory
            os.chdir(current_working_directory)

        if os.path.isfile(union_tree_filepath):
            # Get rid of the union.sbt file
            os.unlink(union_tree_filepath)

        if len(sketches) > 1:
            # Build the bloom filter representation of the entry
            # This is not required for entries with 1 genome only, since that genome is the representative
            command_line = [
                "howdesbt",
                "bfoperate",
                f"--list={sketches_list_filepath}",
                "--or",
                f"--out={os.path.join(self.folder, '{}.bf'.format(self.name))}"
            ]

            try:
                subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            except subprocess.CalledProcessError as e:
                error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

                raise Exception(error_message).with_traceback(e.__traceback__)

        if self.database.flat:
            # Use a flat structure by manually defining the index sbt file
            # All the sketches are at the same level and all under a single root node
            with open(os.path.join(self.folder, "tree", "index.detbrief.sbt"), "w+") as flat_tree:
                # The root node is the bloom filter representation of the cluster
                flat_tree.write(f"{os.path.join(self.folder, '{}.bf'.format(self.name))}\n")

                for sketch_filepath in sketches:
                    flat_tree.write(f"*{sketch_filepath}\n")

        # Set the file path to the sketch representation of the cluster
        self.sketch_filepath = os.path.join(self.folder, f"{self.name}.bf")

    def sketch(self, filepath: os.path.abspath) -> os.path.abspath:
        """Build a sketch representation of the input genome.

        Parameters
        ----------
        filepath : os.path.abspath
            Path to the input genome.

        Raises
        ------
        Exception
            - In case no database metadata is found;
            - In case the current entry is not a genome.
        ValueError
            If the input file format is not supported.

        Returns
        -------
        os.path.abspath
            The path to the output sketch file.
        """

        if not self.database.__class__._validate_metadata(self.database.metadata):
            raise Exception("No database metadata found!")

        if self.level in self.database.__class__.LEVELS:
            # The list of taxonomic levels in `database` do not contain "genome"
            raise Exception("This function can be used in case of genome entries only!")

        if not self.database.__class__._is_supported(filepath):
            raise ValueError("Input file format is not supported!")

        # Use the entry name, not the file name
        sketch_filepath = os.path.join(self.database.root, "sketches", f"{self.name}.bf")

        if os.path.isfile(sketch_filepath):
            self.sketch_filepath = sketch_filepath

            return self.sketch_filepath

        command_line = [
            "howdesbt",
            "makebf",
            f"--k={self.database.metadata['kmer_size']}",
            f"--min={self.database.metadata['min_kmer_occurrence']}",
            f"--bits={self.database.metadata['filter_size']}",
            "--hashes=1",
            "--seed=0,0",
            filepath,
            f"--out={sketch_filepath}",
            f"--threads={self.database.nproc}"
        ]

        try:
            subprocess.check_call(command_line, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            error_message = f"An error has occurred while running\n{' '.join(command_line)}\n\n"

            raise Exception(error_message).with_traceback(e.__traceback__)

        self.sketch_filepath = sketch_filepath

        return self.sketch_filepath

    def get_boundaries(self, limit_number: int=0, limit_percentage: float=100.0) -> str:
        """Search for the cluster boundaries as the minimum and maximum ANI distance
        from the centroid versus all the other genomes in the same cluster.

        WARNING: boundaries can be computed with a minimum of 3 entries.

        Parameters
        ----------
        limit_number : int, default 0
            Limit the pair-wise computation of distances to `limit_number` children at most.
            If `limit_number <= 0`, consider all the whole set of children.
        limit_percentage : float, default 100.0
            Consider a limited percentage on the number of children to limit the pair-wise computation of distances.
            If `limit_percentage >= 100.0`, consider the whole set of childre.
            Note that `limit_percentage` has precedence over `limit_number`.

        Raises
        ------
        Exception
            If the current entry represents a genome.

        Returns
        -------
        tuple
            A tuple with he minimum and maximum ANI distance, and the cluster centroid.
        """

        if self.level == "genome":
            raise Exception("Cannot compute boundaries over a genome entry!")

        if self.level == "species":
            # Use `get_children()` if the current cluster is a species
            # This returns all the genomes under the species cluster
            # Convert the set of children to list to preserve their order
            children = list(self.get_children())

        else:
            # Otherwise, retrieve the species clusters under this specific level 
            species_clusters = list(self.get_children(up_to="species"))

            children = list()

            for species_cluster in species_clusters:
                species_cluster_identifier = self.database.clusters["species"][species_cluster].identifier

                if species_cluster_identifier in self.database.report:
                    children.append(self.database.report[species_cluster_identifier]["centroid"])

                else:
                    # Compute the species cluster centroid
                    # This function is usually called in `_dump_report`, but the report is updated as soon as the boundaries of a species cluster are computed
                    # Thus, this should never happen
                    _, _, species_cluster_centroid = self.database.clusters["species"][species_cluster].get_boundaries(limit_number=limit_number, limit_percentage=limit_percentage)

                    children.append(species_cluster_centroid)

        # Prioritize `limit_percentage` over `limit_number`
        if limit_percentage < 100.0:
            limit = math.ceil(len(children)*limit_percentage/100.0)

        elif limit_number > 0:
            limit = limit_number if limit_number <= len(children) else len(children)

        else:
            limit = len(children)

        if 0 < limit < len(children):
            # Always set the random seed to preserve reproducibility
            random.seed(0)

            # Random sample children up to `limit`
            children = random.sample(children, limit)

        # Boundaries can be computed with a minimum of 3 children
        if len(children) < 3:
            # Pick the first child in lexicographical order as the centroid
            centroid = sorted(children)[0]

            # These boundaries must be estimated
            return (None, None, centroid)

        # Children are genomes
        search_in = self.database.genomes

        # Keep track of the pair-wise distances between children
        pairwise_dists = {child: list() for child in children}

        # Rescale nproc
        nproc = self.database.nproc if len(children) > self.database.nproc else len(children)

        if nproc > 1:
            with mp.Pool(processes=nproc) as pool:
                jobs = [
                    pool.apply_async(
                        self.database.__class__.dist, 
                        args=(
                            search_in[source].sketch_filepath, 
                            [search_in[target].sketch_filepath for target in children[pos+1:]],
                            self.database.metadata["kmer_size"],
                            self.database.tmp,
                            False,
                        )
                    ) for pos, source in enumerate(children)
                ]

                for job in jobs:
                    source_sketch, sketch_dists = job.get()

                    source = os.path.splitext(os.path.basename(source_sketch))[0]

                    for target_sketch in sketch_dists:
                        target = os.path.splitext(os.path.basename(target_sketch))[0]

                        pairwise_dists[source].append(sketch_dists[target_sketch])

                        pairwise_dists[target].append(sketch_dists[target_sketch])

        else:
            for pos, source in enumerate(children):
                if pos < len(children)-1:
                    # Retrieve the source sketch filepath
                    source_sketch = search_in[source].sketch_filepath

                    # Retrieve the target sketch filepaths
                    target_sketches = [search_in[target].sketch_filepath for target in children[pos+1:]]

                    # Compute the ANI distance between source and targets
                    _, dists = self.database.__class__.dist(source_sketch, target_sketches, self.database.metadata["kmer_size"], tmp=self.tmp, resume=False)

                    # Keep track of the ANI distances
                    for target_sketch in dists:
                        # Retrieve the target name
                        target = os.path.splitext(os.path.basename(target_sketch))[0]

                        pairwise_dists[source].append(dists[target_sketch])

                        pairwise_dists[target].append(dists[target_sketch])

        # Search for the centroid
        avg_ani = sys.float_info.max

        centroid = None

        for child in pairwise_dists:
            # Minimum average ANI
            # Centroid: the child that minimizes the distances versus all the other children
            child_avg_ani = statistics.mean(pairwise_dists[child])

            if child_avg_ani < avg_ani:
                avg_ani = child_avg_ani

                centroid = child

        # The centroid is always defined here
        min_ani = round(min(pairwise_dists[centroid]), 5)

        max_ani = round(max(pairwise_dists[centroid]), 5)

        if self.level != "kingdom" and max_ani == 1.0:
            # In case of max ANI distance = 1.0
            # This is too large and it is going to mess up with the assignment
            # Every genome that matches with this node in terms of number of kmers, will eventually be assigned to this node
            # Do not report the actual boudaries so that they can eventually be estimated differently
            # WARNING: We cannot estimate boundaries in case of a kingdom node
            return (None, None, centroid)

        return (round(min_ani, 5), round(max_ani, 5), centroid)
