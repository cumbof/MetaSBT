#!/usr/bin/env python

__author__ = ("Fabio Cumbo (fabio.cumbo@gmail.com)")
__version__ = "0.1.0"
__date__ = "Jun 13, 2022"

import sys, os, time, errno, re, shutil
import argparse as ap
from pathlib import Path
from itertools import partial
from collections import Counter
from utils import checkm, cluster, filter_genomes, howdesbt, init_logger, it_exists, kmtricks_matrix, number, println, run

# Define the module name
TOOL_ID = "update"

# Define the software root directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

def read_params():
    p = ap.ArgumentParser(description="Update a database with new genomes",
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument( "--boundaries",
                    type = str,
                    required = True,
                    help = ("Path to the output table produced by the boundaries module. "
                            "It is required in case of MAGs as input genomes only") )
    p.add_argument( "--boundary-uncertainty",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 0.0,
                    dest = "boundary_uncertainty",
                    help = "Define the percentage of kmers to enlarge and reduce boundaries" )
    p.add_argument( "--completeness",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 0.0,
                    help = "Input genomes must have a minimum completeness percentage before being processed and added to the database" )
    p.add_argument( "--contamination",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 100.0,
                    help = "Input genomes must have a maximum contamination percentage before being processed and added to the database" )
    p.add_argument( "--cleanup",
                    action = "store_true",
                    default = False,
                    help = "Remove temporary data at the end of the pipeline" )
    p.add_argument( "--db-dir",
                    type = str,
                    required = True,
                    dest = "db_dir",
                    help = "This is the database directory with the taxonomically organised sequence bloom trees" )
    p.add_argument( "--dereplicate",
                    action = "store_true",
                    default = False,
                    help = "Enable the dereplication of genomes" )
    p.add_argument( "--extension",
                    type = str,
                    required = True,
                    choices=["fa", "fa.gz", "fasta", "fasta.gz", "fna", "fna.gz"],
                    help = ("Specify the input genome files extension. "
                            "All the input genomes must have the same file extension before running this module") )
    p.add_argument( "--filter-size",
                    type = number(int, minv=10000),
                    dest = "filter_size",
                    help = "This is the size of the bloom filters" )
    p.add_argument( "--input-list",
                    type = str,
                    required = True,
                    dest = "input_list",
                    help = "This file contains the list of paths to the new genomes that will be added to the database" )
    p.add_argument( "--kmer-len",
                    type = number(int, minv=6),
                    dest = "kmer_len",
                    help = "This is the length of the kmers used for building bloom filters" )
    p.add_argument( "--log",
                    type = str,
                    help = "Path to the log file" )
    p.add_argument( "--nproc",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "This argument refers to the number of processors used for parallelizing the pipeline when possible" )
    p.add_argument( "--pplacer_threads",
                    type = number(int, minv=1, maxv=os.cpu_count()),
                    default = 1,
                    help = "Maximum number of threads for pplacer. This is required to maximise the CheckM performances" )
    p.add_argument( "--similarity",
                    type = number(float, minv=0.0, maxv=100.0),
                    default = 100.0,
                    help = ("Dereplicate genomes if they have a percentage of common kmers greater than or equals to the specified one. "
                            "This is used exclusively in conjunction with the --dereplicate argument") )
    p.add_argument( "--taxa",
                    type = str,
                    help = ("Input file with the mapping between input genome IDs and their taxonomic label. "
                            "This is used in case of reference genomes only \"--type=references\"") )
    p.add_argument( "--tmp-dir",
                    type = str,
                    required = True,
                    dest = "tmp_dir",
                    help = "Path to the folder for storing temporary data" )
    p.add_argument( "--type",
                    type = str,
                    required = True,
                    choices=["MAGs", "references"],
                    help = "Define the nature of the input genomes" )
    p.add_argument( "--verbose",
                    action = "store_true",
                    default = False,
                    help = "Print results on screen" )
    p.add_argument( "-v",
                    "--version",
                    action = "version",
                    version = "{} version {} ({})".format(TOOL_ID, __version__, __date__),
                    help = "Print the current {} version and exit".format(TOOL_ID) )
    return p.parse_args()

def update(input_list, input_type, extension, db_dir, kingdom, tmp_dir, boundaries,
           boundary_uncertainty=None, taxa_map=None, completeness=0.0, contamination=100.0, 
           dereplicate=False, similarity=100.0, logger=None, verbose=False, nproc=1, pplacer_threads=1):
    """
    Update a database with a new set of metagenome-assembled genomes and reference genomes.
    Also create new clusters in case the input genomes result too far from everything in the database

    :param input_list:              Path to the input file with the list of input genome paths
    :param input_type:              Nature of the input genomes (MAGs or references)
    :param extension:               File extension of the input genome files
    :param db_dir:                  Path to the database root folder
    :param kingdom:                 Retrieve genomes that belong to a specific kingdom
    :param tmp_dir:                 Path to the temporary folder
    :param boundaries:              Path to the boundaries table file
    :param boundary_uncertainty:    Percentage of kmers to enlarge and reduce boundaries
    :param taxa_map:                Path to the file with the mapping between the input genome names and their taxonomic labels
                                    Used only in case of input reference genomes
    :param completeness:            Threshold on the CheckM completeness
    :param contamination:           Threshold on the CheckM contamination
    :param dereplicate:             Enable the dereplication step to get rid of replicated genomes
    :param similarity:              Get rid of genomes according to this threshold in case the dereplication step is enabled
    :param logger:                  Logger object
    :param verbose:                 Print messages on screen
    :param nproc:                   Make the process parallel when possible
    :param pplacer_threads:         Maximum number of threads to make pplacer parallel with CheckM
    """

    # Define a partial println function to avoid specifying logger and verbose
    # every time the println function is invoked
    printline = partial(println, logger=logger, verbose=verbose)

    # Check whether the selected kingdom exists in the database
    if not os.path.exists(os.path.join(db_dir, "k__{}".format(kingdom))):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), os.path.join(db_dir, "k__{}".format(kingdom)))

    # Check whether the manifest file exists in the database
    manifest_filepath = os.path.join(db_dir, "manifest.txt")
    if not os.path.exists(manifest_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), manifest_filepath)

    # Retrieve both the kmer length and the filter size from the manifest file
    kmer_len = 0
    filter_size = 0
    with open(manifest_filepath) as manifest:
        for line in manifest:
            line = line.strip()
            if line:
                line_split = line.split(" ")
                if line_split[0] == "--kmer-len":
                    kmer_len = int(line_split[1])
                elif line_split[0] == "--filter-size":
                    filter_size = int(line_split[1])

    # Check whether the kmer length and the filter size have been successfully retrieved
    if kmer_len == 0 or filter_size == 0:
        raise Exception("Unable to retrieve data from the manifest file:\n{}".format(manifest_filepath))
    
    # Load the list of genome paths
    genomes_paths = [path.strip() for path in open(input_list).readlines() if path.strip()]

    # Load the file with the mapping between genome names and taxonomic labels
    # Only in case of input reference genomes
    taxonomies = dict()
    if input_type == "references":
        taxonomies = {line.strip().split("\t")[0]: line.strip().split("\t")[1] for line in open(taxa_map).readlines() \
                        if line.strip() and not line.strip().startswith("#")}
        
        # Check whether a taxonomic label is defined for each of the input genome
        genome_names = list()
        for genome_path in genomes_path:
            genome_name = os.path.splitext(os.path.basename(genome_path))[0]
            if genome_path.endswith(".gz"):
                genome_name = os.path.splitext(genome_name)[0]
            gneome_names.append(genome_name)
        
        if len(set(genome_names).intersection(set(taxonomies.keys()))) < len(set(genome_names)):
            raise Exception("Unable to find all the taxonomic labels for the input reference genomes")

    printline("Processing {} input genomes ({})".format(len(genomes_paths), input_type))

    # Take track of the CheckM data
    # Mapping between genome names and lines in the CheckM output tables
    checkm_header = ""
    checkm_data = dict()

    # Check completeness and contamination requirements
    if genomes_paths and (completeness > 0.0 or contamination < 100.0):
        printline("Quality controlling {} genomes".format(len(genomes_pahts)))
        printline("Minimum completeness threshold: {}".format(completeness))
        printline("Maximum contamination threshold: {}".format(contamination))

        checkm_tmp_dir = os.path.join(tmp_dir, "checkm")
        os.makedirs(checkm_tmp_dir, exist_ok=True)
        
        # Run CheckM to quality control genomes
        checkm_tables = checkm(genomes_paths, checkm_tmp_dir, 
                               file_extension=extension, nproc=nproc, pplacer_threads=pplacer_threads)
        
        # Load CheckM tables
        for checkm_table in checkm_tables:
            with open(checkm_table) as table:
                line_count = 0
                for line in table:
                    line = line.strip()
                    if line:
                        if line_count == 0:
                            checkm_header = line
                        else:
                            line_split = line.split("\t")
                            checkm_data[line_split[0]] = line
                        line_count += 1

        # Filter genomes according to the input --completeness and --contamination thresholds
        genome_ids = filter_checkm_tables(checkm_tables, completeness=0.0, contamination=100.0)
        # Build the new list of genome paths that passed the filter
        filtered_genomes_paths = list()
        for genome_path in genomes_paths:
            genome_id = os.path.splitext(os.path.basename(genome_path))[0]
            if genome_path.endswith(".gz"):
                genome_id = os.path.splitext(genome_id)[0]
            if genome_id in genome_ids:
                filtered_genomes_paths.append(genome_path)

        printline("{} genomes have been excluded".format(len(genomes_paths)-len(filtered_genomes_paths)))

        # Define the new list of genomes
        genomes_paths = filtered_genomes_paths
    
    # Check whether the input genomes must be dereplicated
    if len(genomes_paths) > 1 and dereplicate:
        printline("Dereplicating {} genomes".format(len(genomes_pahts)))

        # Define the kmtricks temporary folder
        kmtricks_tmp_dir = os.path.join(tmp_dir, "kmtricks")
        os.makedirs(kmtricks_tmp_dir, exist_ok=True)

        # Create a fof file with the list of genomes
        genomes_fof_filepath = os.path.join(kmtricks_tmp_dir, "genomes.fof")
        ordered_genome_names = list()
        with open(genomes_fof_filepath, "w+") as genomes_fof:
            for genome_path in genomes_paths:
                genome_name = os.path.splitext(os.path.basename(genome_path))[0]
                if genome_path.endswith(".gz"):
                    genome_name = os.path.splitext(genome_name)[0]
                genomes_fof.write("{} : {}\n".format(genome_name, genome_path))
                ordered_genome_names.append(genome_name)
        
        # Run kmtricks to produce the kmer matrix
        output_table = os.path.join(kmtricks_tmp_dir, "matrix.txt")
        kmtricks_matrix(genomes_fof_filepath, kmtricks_tmp_dir, nproc, output_table)

        # Add the header line to the kmer matrix
        with open(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), "w+") as file1:
            file1.write("# kmer\t{}\n".format("\t".join(ordered_genome_names)))
            with open(output_table) as file2:
                for line in file2:
                    line = line.strip()
                    if line:
                        file1.write("{}\n".format(line))
        
        # Get rid of the matrix withour header line
        os.unlink(output_table)
        shutil.move(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), output_table)

        # Filter genomes according to their percentage of common kmers defined with --similarity
        dereplicated_genomes_filepath = os.path.join(tmp_dir, "kmtricks", "dereplicated.txt")
        filter_genomes(output_table, dereplicated_genomes_filepath, similarity=similarity)

        # Load the list of dereplicated genomes
        dereplicated_genomes = [genome_id.strip() for genome_id in open(dereplicated_genomes_filepath).readlines() if genome_id.strip()]

        # Build the new list of genome paths that passed the filter
        filtered_genomes_paths = list()
        for genome_path in genomes_paths:
            genome_id = os.path.splitext(os.path.basename(genome_path))[0]
            if genome_path.endswith(".gz"):
                genome_id = os.path.splitext(genome_id)[0]
            if genome_id in dereplicated_genomes:
                filtered_genomes_paths.append(genome_path)
        
        printline("{} genomes have been excluded".format(len(genomes_paths)-len(filtered_genomes_paths)))

        # Define the new list of genomes
        genomes_paths = filtered_genomes_paths

    # Check whether at least one genome survived both the quality control and dereplication steps
    if not genomes_paths:
        printline("No input genomes available")

    # Create a temporary folder in case input genomes are gzip compressed
    tmp_genomes_dir = os.path.join(tmp_dir, "genomes")
    os.makedirs(tmp_genomes_dir, exist_ok=True)

    # In case one or more reference genomes are assigned to an unknown cluster
    # Keep track of the unknown taxonomy and change the lineage with
    # the most occurring taxa among the assigned reference genomes
    to_known_taxa = dict()

    # Keep track of the lineages that must be rebuilt
    rebuild = list()

    # Keep track of the unassigned genomes
    # They will be assigned to new groups
    unassigned = list()

    # Start processing input genomes
    for genome_path in genomes_paths:
        filepath = genome_path
        genome_name = os.path.splitext(os.path.basename(filepath))[0]
        genome_ext = os.path.splitext(os.path.basename(filepath))[1][1:]
        # in case the current genome is gzip compressed
        if genome_path.endswith(".gz"):
            filepath = os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0])
            genome_name = os.path.splitext(os.path.basename(filepath))[0]
            genome_ext = os.path.splitext(os.path.basename(filepath))[1][1:]
            # Unzip the current genome to the tmp folder
            run(["gunzip", "-c", genome_path], stdout=filepath)
        
        printline("Profiling {}".format(genome_name))

        # Run the profiler to establish the closest genome and the closest group 
        # for each taxonomic level in the tree
        # The profile module is in the same folder of the update module
        run([sys.executable, os.path.join(SCRIPT_DIR, "profile.py"), "--input-file", filepath,
                                                                     "--input-id", filepath,
                                                                     "--tree", os.path.join(db_dir, "index", "index.detbrief.sbt"),
                                                                     "--expand",
                                                                     "--output-dir", os.path.join(tmp_dir, "profiling"),
                                                                     "--output-prefix", genome_name],
            silence=True)
        
        # Define the path to the profiler output file
        profile_path = os.path.join(tmp_dir, "profiling", "{}__profiles.tsv".join(genome_name))

        # Check whether the profile exists
        if os.path.exists(profile_path):
            # Load the profile
            profile_data = dict()
            with open(profile_path) as profile:
                for line in profile:
                    line = line.strip()
                    if line:
                        line_split = line.split("\t")
                        profile_data[line_split[1]] = {
                            "taxonomy": line_split[2],
                            "common_kmers": int(line_split[3]),
                            "score": float(line_split[4])
                        }
            
            # Discard the genome if it is too similar with the closest genome according to the similarity score
            # Also, do not discard the input genome if it is a reference genome and the closest genome in the database is
            # a MAG with a high number of overlapped kmers according to the similarity score
            closest_genome = profile_data["genome"]

            # Reconstruct the full lineage of the closest taxa
            closest_taxa = list()
            closest_common_kmers = 0
            closest_score = 0.0
            for level in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
                closest_taxa.append(profile_data[level]["taxonomy"])
                if level == "species":
                    closest_common_kmers = profile_data[level]["common_kmers"]
                    closest_score = profile_data[level]["score"]
            closest_taxa = "|".join(closest_taxa)

            printline("Closest lineage: {} (score {})".format(closest_taxa, closest_score))
            printline("Closest genome: {} (score {})".format(closest_genome["taxonomy"], closest_genome["score"]))

            # Define the folder path of the closest taxonomy under the database
            closest_taxadir = os.path.join(db_dir, closest_taxa.replace("|", os.sep))

            # Load the set of reference genomes that belongs to the closest species
            references_filepath = os.path.join(closest_taxadir, "references.txt")
            references = list()
            if os.path.exists(references_filepath):
                references = [ref.strip() for ref in open(references_filepath).readlines() if ref.strip()]

            # Also load the set of MAGs that belongs to the closest species
            mags_filepath = os.path.join(closest_taxadir, "mags.txt")
            mags = list()
            if os.path.exists(mags_filepath):
                mags = [mag.strip() for mag in open(mags_filepath).readlines() if mag.strip()]

            # Check whether the input genome must be discarded
            skip_genome = False
            if dereplicate and closest_genome["score"]*100.0 >= similarity:
                printline("Dereplicating genome")
                if input_type == "MAGs":
                    # In case the input genome is a MAG
                    if closest_genome["taxonomy"] in references:
                        # In case the closest genome is a reference genome
                        skip_genome = True
                        # Print the reason why the current genome is excluded
                        printline("Discarding genome\nInput genome is a MAG and the closest genome is a reference")
                    elif closest_genome["taxonomy"] in mags:
                        # In case the closest genome is a MAG
                        skip_genome = True
                        # Print the reason why the current genome is excluded
                        printline("Discarding genome\nInput genome and the closest genome are both MAGs")
                
                elif input_type == "references":
                    # In case the input genome is a reference genome
                    if closest_genome["taxonomy"] in references:
                        # In case the closest genome is a reference genome
                        skip_genome = True
                        # Print the reason why the current genome is excluded
                        printline("Discarding genome\nInput genome and the closest genome are both reference genomes")
            
            # In case the input genome survived the dereplication with the closest genome in the database
            if not skip_genome:
                # Retrieve the minimum and maximum common kmers for the closest taxa
                min_bound, max_bound = get_level_boundaries(boundaries, closest_taxa)

                if input_type == "MAGs":
                    # In case the input genome is a MAG
                    if closest_common_kmers <= max_bound and closest_common_kmers >= min_bound:
                        # Assign the current genome to the closest lineage
                        target_genome = oa.path.join(closest_taxadir, "genomes", "{}.{}".format(genome_name, genome_ext))
                        
                        # The input genome is always uncompressed
                        # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                        run(["gzip", target_genome], silence=True)
                        target_genome = "{}.gz".format(target_genome)

                        # Add the current genome to the list of MAGs
                        with open(mags_filepath, "a+") as file:
                            file.write("{}\n".format(genome_name))
                        
                        # Also add its CheckM statistics if available
                        if genome_name in checkm_data:
                            checkm_filepath = os.path.join(closest_taxadir, "checkm.tsv")
                            checkm_exists = os.path.exists(checkm_filepath)
                            with open(checkm_filepath, "a+") as checkm_file:
                                if not checkm_exists:
                                    checkm_file.write("{}\n".format(checkm_header))
                                checkm_file.write("{}\n".format(checkm_data[genome_name]))
                        
                        # Do not remove the index here because there could be other input genome
                        # with the same closest taxonomic label
                        rebuild.append(closest_taxa)

                        printline("{} has been characterised as {}".format(genome_name, closest_taxa))
                
                    else:
                        # Mark the input genome as unassigned
                        # The CheckM statistics for this genome will be reported after the assignment
                        unassigned.append(genome_path)

                        printline("{} is still unassigned".format(genome_name))
                
                elif input_type == "references":
                    # In case the input genome is a reference genome
                    # Retrieve the taxonomic label from the input mapping file
                    taxalabel = taxonomies[genome_name]

                    if closest_common_kmers <= max_bound and closest_common_kmers >= min_bound:
                        # Check whether the closest taxonomy contains any reference genome
                        how_many_references = 0
                        if os.path.exists(os.path.join(closest_taxadir, "references.txt")):
                            with open(os.path.join(closest_taxadir, "references.txt")) as file:
                                for line in file:
                                    line = line.strip()
                                    if line:
                                        how_many_references += 1
                        
                        if how_many_references == 0:
                            # If the closest genome belongs to a new cluster with no reference genomes
                            # Assign the current reference genome to the new cluster and rename its lineage with the taxonomic label of the reference genome
                            # Do not change lineage here because there could be more than one reference genome assigned to the same unknown cluster
                            # Assign a new taxonomy according to a majority rule applied on the reference genomes taxa
                            if closest_taxadir not in to_known_taxa:
                                to_known_taxa[closest_taxadir] = list()
                            to_known_taxa[closest_taxadir].append(genome_name)

                            printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                        elif how_many_references >= 1:
                            # If the closest genome belongs to a cluster with at least one reference genome
                            if taxalabel != closest_taxa:
                                # If the taxonomic labels of the current reference genome and that one of the closest genomedo not match
                                # Report the inconsistency
                                printline("Inconsistency found:\nInput genome: {}\nClosest lineage: {}".format(taxalabel, closest_taxa))
                            
                            # Assign the current genome to the closest lineage
                            target_genome = oa.path.join(closest_taxadir, "genomes", "{}.{}".format(genome_name, genome_ext))
                            
                            # The input genome is always uncompressed
                            # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                            run(["gzip", target_genome], silence=True)
                            target_genome = "{}.gz".format(target_genome)

                            # Add the current genome to the list of reference genomes
                            with open(references_filepath, "a+") as file:
                                file.write("{}\n".format(genome_name))

                            # Do not remove the index here because there could be other input genome
                            # with the same closest taxonomic label
                            rebuild.append(closest_taxa)

                            printline("{} has been characterised as {}".format(genome_name, closest_taxa))
                        
                        # Also add its CheckM statistics if available
                        if genome_name in checkm_data:
                            checkm_filepath = os.path.join(closest_taxadir, "checkm.tsv")
                            checkm_exists = os.path.exists(checkm_filepath)
                            with open(checkm_filepath, "a+") as checkm_file:
                                if not checkm_exists:
                                    checkm_file.write("{}\n".format(checkm_header))
                                checkm_file.write("{}\n".format(checkm_data[genome_name]))
                    
                    else:
                        # If nothing is closed enough to the current genome and its taxonomic label does not exist in the database
                        # Create a new branch with the new taxonomy for the current genome
                        taxdir = os.path.join(db_dir, taxalabel.replace("|", os.sep))
                        os.makedirs(os.path.join(taxdir, "genomes"), exist_ok=True)

                        # Assign the current genome to the closest lineage
                        target_genome = oa.path.join(taxdir, "genomes", "{}.{}".format(genome_name, genome_ext))
                        
                        # The input genome is always uncompressed
                        # It must be gzip compressed before moving it to the genomes folder of the closest taxonomy
                        run(["gzip", target_genome], silence=True)
                        target_genome = "{}.gz".format(target_genome)

                        # Add the current genome to the list of reference genomes
                        with open(os.path.join(taxdir, "references.txt"), "a+") as file:
                            file.write("{}\n".format(genome_name))
                        
                        # Do not remove the index here because there could be other input genome
                        # with the same closest taxonomic label
                        rebuild.append(closest_taxa)

                        printline("{} has been characterised as {}".format(genome_name, closest_taxa))

                        # Also add its CheckM statistics if available
                        if genome_name in checkm_data:
                            checkm_filepath = os.path.join(closest_taxadir, "checkm.tsv")
                            checkm_exists = os.path.exists(checkm_filepath)
                            with open(checkm_filepath, "a+") as checkm_file:
                                if not checkm_exists:
                                    checkm_file.write("{}\n".format(checkm_header))
                                checkm_file.write("{}\n".format(checkm_data[genome_name]))

        # Remove the uncompressed version of the input genome in the temporary folder
        if os.path.exists(os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0])):
            os.unlink(os.path.join(tmp_genomes_dir, os.path.splitext(os.path.basename(genome_path))[0]))

    # In case a reference has been assigned to an unknown cluster
    # Update the cluster taxonomy by applying a majority voting on the reference genomes taxa
    if to_knwon_taxa:
        printline("Characterising {} unknown taxa".format(len(to_known_taxa)))

        # Iterate over the unknown clusters in dictionary
        for unknown_taxonomy in to_known_taxa:
            # Rebuild the folder path to the unknown taxonomy in the database
            unknown_taxonomy_path = os.path.join(db_dir, unknown_taxonomy.replace("|", os.sep))

            # Take track of the reference genomes taxonomic labels
            known_taxa = list()

            # Retrieve the reference genomes taxa from the input mapping file
            for genome in to_unknown_taxa[unknown_taxonomy]:
                known_taxa.append(taxonomies[genome])
            
            # Get the most occurring taxonomic label
            known_taxa_counter = Counter(known_taxa)
            assigned_taxonomy = known_taxa_counter.most_common(1)[0][0]

            # Finally move the genome files
            for genome_path in genomes_paths:
                genome_name = os.path.splitext(os.path.basename(genome_path))[0]
                if genome_path.endswith(".gz"):
                    genome_name = os.path.splitext(genome_name)[0]
                
                # Check whether the genome name is in the current list of genomes for the unknown taxonomy
                if genome_name in to_known_taxa[unknown_taxonomy]:
                    # Copy the input genome file to the genomes folder
                    shutil.copy(genome_path, os.path.join(unknown_taxonomy_path, "genomes"))
                    
                    # In case the input genome is not gzip compressed
                    if not genome_path.endswith(".gz"):
                        run(["gzip", os.path.join(unknown_taxonomy_path, "genomes", os.path.basename(genome_path))], silence=True)
                    
                    # Also update the list of reference genomes
                    with open(os.path.join(unknown_taxonomy_path, "references.txt"), "a+") as file:
                        file.write("{}\n".format(genome_name))
            
            # Split taxonomies into levels
            unknown_levels = unknown_taxonomy.split("|")
            assigned_levels = assigned_taxonomy.aplit("|")

            # Characterise unknown clusters
            # Iterate backwards over level
            for i in range(6, -1, -1):
                if unknown_levels[i] == assigned_levels[i]:
                    break
                else:
                    # Get the partial levels for the new taxonomy
                    new_levels = assigned_levels[:i+1]
                    new_levels_subpath = os.path.join(db_dir, os.sep.join(new_levels))
                    makedirs(new_levels_subpath, exist_ok=True)

                    # Get the partial levels for the old taxonomy
                    old_levels = unknown_levels[:i+1]
                    old_levels_subpath = os.path.join(db_dir, os.sep.join(old_levels))

                    # Start moving folders
                    for folder_path in Path(old_levels_subpath).glob("**/*"):
                        if os.path.isdir(str(folder_path)):
                            shutil.move(str(folder_path), new_levels_subpath)
                    
                    # Also move remaining files
                    for file_path in Path(old_levels_subpath).glob("**/*"):
                        if os.path.isfile(str(file_path)):
                            shutil.move(str(file_path), new_levels_subpath)

                    # Remove the old bloom filter root node
                    bloom_filter_node = os.path.join(new_levels_subpath, "{}.bf".format(unknown_levels[i]))
                    if os.path.exists(bloom_filter_node):
                        os.unlink(bloom_filter_node)
            
            # Add the new renamed cluster to the rebuild list of taxa
            rebuild.append(assigned_taxonomy)

    # Cluster unassigned genomes before rebuilding the updated lineages
    if unassigned:
        printline("Defining new clusters for {} unassigned genomes".format(len(unassigned)))

        # The kmers matrix already exists in case the dereplication step has been performed
        # Otherwise, run kmtricks
        kmtricks_tmp_dir = os.path.join(tmp_dir, "kmtricks")
        output_table = os.path.join(kmtricks_tmp_dir, "matrix.txt")
        if not os.path.exists(output_table):
            # Create a fof file with the list of unassigned genomes
            genomes_fof_filepath = os.path.join(kmtricks_tmp_dir, "genomes.fof")
            ordered_genome_names = list()
            with open(genomes_fof_filepath, "w+") as genomes_fof:
                for genome_path in unassigned:
                    genome_name = os.path.splitext(os.path.basename(genome_path))[0]
                    if genome_path.endswith(".gz"):
                        genome_name = os.path.splitext(genome_name)[0]
                    genomes_fof.write("{} : {}\n".format(genome_name, genome_path))
                    ordered_genome_names.append(genome_name)

            # Build the kmers matrix
            kmtricks_matrix(genomes_fof_filepath, kmtricks_tmp_dir, nproc, output_table)

            # Add the header line to the kmer matrix
            with open(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), "w+") as file1:
                file1.write("# kmer\t{}\n".format("\t".join(ordered_genome_names)))
                with open(output_table) as file2:
                    for line in file2:
                        line = line.strip()
                        if line:
                            file1.write("{}\n".format(line))
            
            # Get rid of the matrix withour header line
            os.unlink(output_table)
            shutil.move(os.path.join(kmtricks_tmp_dir, "matrix_with_header.txt"), output_table)

        # Cluster genomes according to the boundaries defined by the boundaries module
        # Define a cluster for each taxomomic level
        # Look at the genomes profiles and update or build new clusters 
        assignments_filepath = os.path.join(db_dir, "assignments.txt")
        cluster(output_table, boundaries, manifest_filepath, os.path.join(tmp_dir, "profiling"), assignments_filepath)

        # Load the new assignments
        assignments = dict()
        with open(assignments_filepath) as file:
            for line in file:
                line = line.strip()
                if line:
                    if not line.startswith("#"):
                        line_split = line.split("\t")
                        assignments[line_split[0]] = line_split[1]

        # Iterate over the input genomes
        for genome_path in genomes_paths:
            genome_name = os.path.splitext(os.path.basename(genome_path))[0]
            if genome_path.endswith(".gz"):
                genome_name = os.path.splitext(genome_name)[0]
            
            if genome_name in assignments:
                # Create the new cluster folder in the database
                tax_dir = os.path.join(db_dir, line_split[1].replace("|", os.sep))
                tax_genomes_dir = os.path.join(tax_dir, "genomes")
                makedirs(tax_genomes_dir, exist_ok=True)

                # Copy the input genome into the genomes folder of the new cluster
                shutil.copy(genome_path, tax_genomes_dir)

                # In case the input genome is not gzip compressed
                if not genome_path.endswith(".gz"):
                    run(["gzip", os.path.join(tax_genomes_dir, os.path.basename(genome_path))], silence=True)

                # Also update the mags.txt file
                with open(os.path.join(tax_dir, "mags.txt"), "a+") as file:
                    file.write("{}\n".format(genome_name))

                # Also report the CheckM statistics of the genomes in the new clusters
                if genome_name in checkm_data:
                    checkm_filepath = os.path.join(tax_dir, "checkm.tsv")
                    checkm_exists = os.path.exists(checkm_filepath)
                    with open(checkm_filepath, "a+") as checkm_file:
                        if not checkm_exists:
                            checkm_file.write("{}\n".format(checkm_header))
                        checkm_file.write("{}\n".format(checkm_data[genome_name]))

                # Add the full taxonomy to the list of taxonomic labels that must be rebuilt
                rebuild.append(assignments[genome_name])

    # Check whether there is at least one lineage that must be rebuilt
    rebuild = list(set(rebuild))
    if rebuild:
        printline("Updating the database")
        # Extract the taxonomic levels from the list of taxa that must be rebuilt
        # Process all the species first, then all the genera, and so on up to the kingdom level
        for i in range(6, -1, -1):
            taxalist = set()
            for label in rebuild:
                levels = label.split("|")
                # Temporarily skip partial taxonomic labels if the number of levels is lower than i
                # Waiting for the right i
                if len(levels) < i+1:
                    continue

                # Build the partial taxonomic label
                taxalist.add("|".join(levels[:i+1]))
            
            # Rebuild the sequence bloom trees
            for taxonomy in taxalist:
                printline("\t{}".format(taxonomy))
                tax_dir = os.path.join(db_dir, taxonomy.replace("|", os.sep))
                
                # Remove the old index if it exists
                if os.path.exists(os.path.join(tax_dir, "index")):
                    shutil.rmtree(os.path.join(tax_dir, "index"), ignore_errors=True)
                
                # Also remove the copy of the root node if it exists
                if os.path.exists(os.path.join(tax_dir, "{}.bf".format(os.path.basename(tax_dir)))):
                    os.unlink(os.path.join(tax_dir, "{}.bf".format(os.path.basename(tax_dir))))
                
                # Rebuild the index with HowDeSBT
                howdesbt(tax_dir, kmer_len, filter_size, nproc)
    
    else:
        printline("No lineages have been updated")

def main():
    # Load command line parameters
    args = read_params()

    # Initialise the logger
    logger = init_logger(filepath=args.log, verbose=args.verbose)

    # Check whether the database folder exists
    if not it_exists(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)
    
    # Check whether the file with the list of input genome paths exists
    if not it_exists(args.input_list):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.input_list)

    # Check whether the input file with the mapping between genome name and taxonomic labels exists
    # Only in case of input reference genomes
    if args.type == "references" and not it_exists(args.taxa):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.taxa)

    # Also create the temporary folder
    # Do not raise an exception in case it already exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    t0 = time.time()

    update(args.input_list, args.type, args.extension, args.db_dir, args.kingdom, args.tmp_dir, args.boundaries,
           boundary_uncertainty=args.boundary_uncertainty, taxa_map=taxa, completeness=completeness, 
           contamination=contamination, dereplicate=dereplicate, similarity=similarity, logger=logger, 
           verbose=args.verbose, nproc=args.nproc, pplacer_threads=args.pplacer_threads)

    if args.cleanup:
        # Remove the temporary folder
        println("Cleaning up temporary space".format(level), logger=logger, verbose=verbose)
        shutil.rmtree(args.tmp_dir, ignore_errors=True)

    t1 = time.time()
    println("Total elapsed time {}s".format(int(t1 - t0)), 
            logger=logger, verbose=args.verbose)

if __name__ == "__main__":
    main()