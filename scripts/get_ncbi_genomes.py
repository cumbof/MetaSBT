#!/usr/bin/env python3
"""
Retrieve reference genomes and metagenome-assembled genomes for a specific 
superkingdom and kingdom from NCBI GenBank
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 27, 2023"

import argparse as ap
import datetime
import gzip
import multiprocessing as mp
import os
import re
import subprocess
import tarfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import tqdm
from urllib.request import urlretrieve

TOOL_ID = "get_ncbi_genomes"

# Define the url to the NCBI taxdump
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Define the url to the NCBI GenBank Assembly Summary
# https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
# https://www.ncbi.nlm.nih.gov/assembly/help/
ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

# Consider the excluded_from_refseq tags for discriminating reference genomes and MAGs
# https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
REFERENCE_TAGS = [
    "derived from single cell",
    "derived from surveillance project",
    "assembly from type material",
    "assembly from synonym type material",
    "assembly designated as neotype",
    "assembly designated as reftype",
    "assembly from pathotype material",
    "assembly from proxytype material",
    "missing strain identifier",
    "genus undefined",
    "from large multi-isolate project"
]

# In case of a MAG, exclude the genome if at least one of the following tags
# are reported under the excluded_from_refseq column
EXCLUDE_TAGS = [
    "abnormal gene to sequence ratio",
    "chimeric",
    "contaminated",
    "genome length too large",
    "genome length too small",
    "hybrid",
    "low gene count",
    "low quality sequence",
    "many frameshifted proteins",
    "metagenome",
    "misassembled",
    "mixed culture",
    "untrustworthy as type"
]


def read_params():
    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description=(
            "Retrieve reference genomes and metagenome-assembled genomes for a specific "
            "superkingdom and kingdom from NCBI GenBank"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        dest="db_dir",
        help="Path to the root folder of a MetaSBT database",
    )
    p.add_argument(
        "--download",
        action="store_true",
        default=False,
        help="Download genome files"
    )
    p.add_argument(
        "--kingdom",
        type=str,
        help=(
            "Consider genomes whose lineage belongs to a specific kingdom. "
            "It is optional and must be provided in conjunction with --superkingdom"
        ),
    )
    p.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Retrieve genomes in parallel",
    )
    p.add_argument(
        "--out-dir",
        type=os.path.abspath,
        required=True,
        dest="out_dir",
        help="Path to the output folder",
    )
    p.add_argument(
        "--superkingdom",
        type=str,
        required=True,
        choices=["Archaea", "Bacteria", "Eukaryota", "Viruses"],
        help="Consider genomes whose lineage belongs to a specific superkingdom",
    )
    p.add_argument(
        "--taxa-level-id",
        type=str,
        choices=["phylum", "class", "order", "family", "genus", "species"],
        dest="taxa_level_id",
        help="Taxonomic level identifier"
    )
    p.add_argument(
        "--taxa-level-name",
        type=str,
        dest="taxa_level_name",
        help=(
            "Name of the taxonomic level. "
            "Must be used in conjunction with \"--taxa-level-id\""
        )
    )
    p.add_argument(
        "--type",
        type=str,
        choices=["reference", "mag"],
        help="Retrieve reference genomes or metagenome-assembled genomes only"
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def get_genomes_in_db(db_dir: str) -> List[str]:
    """
    Retrieve the list of reference genomes and MAGs in the MetaSBT database

    :param db_dir:  Path to the root folder of the MetaSBT database
    :return:        List with genome IDs in the database
    """

    genomes = list()

    for genome_type in ["references", "mags"]:
        genome_files = Path(db_dir).glob("**/{}.txt".format(genome_type))

        for filepath in genome_files:
            genomes.extend(
                [genome.strip() for genome in open(str(filepath)).readlines() \
                    if genome.strip() and not genome.strip().startswith("#")]
            )

    return genomes


def level_name(current_level: str, prev_level: str) -> str:
    """
    Define a taxonomic level name

    :param current_level:   Current level name
    :param prev_level:      Previous level name in case of unclassified
    :return:                The new level name
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

    return "{}{}".format(level_prefix, level_suffix)


def download_taxdump(taxdump_url: str, folder_path: str) -> Tuple[str, str]:
    """
    Download and extract the NCBI taxdump tarball

    :param taxdump_url:     URL to the NCBI taxdump
    :param folder_path:     Path to the folder in which the taxdump tarball will be unpacked
    :return:                The nodes.dmp and names.dmp file paths
    """

    # Create the taxdump folder in the temporary directory
    taxdump_dir = os.path.join(folder_path, "taxdump")
    os.makedirs(taxdump_dir, exist_ok=True)

    nodes_dmp = os.path.join(taxdump_dir, "nodes.dmp")
    names_dmp = os.path.join(taxdump_dir, "names.dmp")

    if os.path.isfile(nodes_dmp) and os.path.isfile(names_dmp):
        return nodes_dmp, names_dmp

    taxdump = os.path.join(folder_path, os.path.basename(taxdump_url))

    if not os.path.isfile(taxdump):
        try:
            urlretrieve(taxdump_url, taxdump)

        except Exception:
            raise Exception("Unable to retrieve data from remote location\n{}".format(taxdump_url))

    # Decompress the archive
    with tarfile.open(taxdump, "r:gz") as tar:
        tar.extractall(taxdump_dir)
    
    return nodes_dmp, names_dmp


def ncbitax2lin(
    tmpdir: str,
    nodes_dmp: str,
    names_dmp: str,
    superkingdom: Optional[str] = None,
    kingdom: Optional[str] = None
) -> Dict[str, str]:
    """
    Run ncbitax2lin over nodes and names dumps and produce the
    mapping between NCBI tax IDs and ful taxonomic labels

    :param tmpdir:          Path to the tmp directory
    :param nodes_dmp:       Path to the NCBI nodes dump
    :param names_dmp:       Path to the NCBI names dump
    :param superkingdom:    Filter results on this superkingdom only
    :param kingdom:         Filter results on this kingdom only
    :return:                Dictionary with the mapping between NCBI tax IDs and full taxonomic labels
    """
    
    ncbitax2lin_table = os.path.join(tmpdir, "ncbi_lineages.csv.gz")

    if not os.path.isfile(ncbitax2lin_table):
        # Run ncbitax2lin
        subprocess.check_call(
            [
                "ncbitax2lin",
                "--nodes-file",
                nodes_dmp,
                "--names-file",
                names_dmp,
                "--output",
                ncbitax2lin_table,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

    taxa_map = dict()

    with gzip.open(ncbitax2lin_table, "rt") as ncbi_table:
        # Load the first line as header and search for "superkingdom" and "kingdom" columns
        header = ncbi_table.readline().split(",")
        superkingdom_pos = header.index("superkingdom")  # Archaea, Bacteria, Eukaryota, Viruses
        kingdom_pos = header.index("kingdom")

        for line in ncbi_table:
            line = line.strip()
            if line:
                line_split = line.split(",")

                # Check whether the current taxonomy must be processed
                skip = True

                if not superkingdom and line_split[superkingdom_pos].strip():
                    skip = False

                elif line_split[superkingdom_pos] == superkingdom:
                    if not kingdom:
                        skip = False

                    elif line_split[kingdom_pos] == kingdom:
                        skip = False

                if not skip:
                    # Build the current full taxonomic label
                    label = "k__{}|p__{}|c__{}|o__{}|f__{}|g__{}|s__{}".format(
                        superkingdom,  # Superkingdom
                        level_name(line_split[2], superkingdom if superkingdom else line_split[1]),  # Phylum
                        level_name(line_split[3], line_split[2]),  # Class
                        level_name(line_split[4], line_split[3]),  # Order
                        level_name(line_split[5], line_split[4]),  # Family
                        level_name(line_split[6], line_split[5]),  # Genus
                        level_name(line_split[7], line_split[6]),  # Species
                    )

    return taxa_map


def get_assembly_summary(assembly_summary_url: str, tmpdir: str) -> Dict[str, List[Dict[str, str]]]:
    """
    Download and load the last available NCBI GenBank Assembly Report table

    :param assembly_summary_url:    URL to the NCBI GenBank Assembly Report table
    :param tmpdir:                  Path to the tmp folder
    :return:                        Dictionary with genomes info
    """

    assembly_summary = dict()

    assembly_summary_filepath = os.path.join(tmpdir, os.path.basename(assembly_summary_url))

    if not os.path.isfile(assembly_summary_filepath):
        # Download the NCBI GenBank Assembly Summary table
        # and load the list of genomes grouped by species taxid
        urlretrieve(ASSEMBLY_SUMMARY_URL, assembly_summary_filepath)

    with open(assembly_summary_filepath) as asf:
        # Skip the first line, it is just a comment
        next(asf)

        # Load the header line
        header = next(asf)[1:].strip().split("\t")

        for line in asf:
            if line.strip():
                line_split = line.split("\t")

                species_taxid = line_split[header.index("species_taxid")]
                ftp_path = line_split[header.index("ftp_path")]
                genome_url = os.path.join(ftp_path, "{}_genomic.fna.gz".format(os.path.basename(ftp_path)))

                if species_taxid:
                    if species_taxid not in assembly_summary:
                        assembly_summary[species_taxid] = list()

                    species_info = dict()

                    for h in header:
                        species_info[h] = line_split[header.index(h)].strip()

                    genome_type = "na"

                    if not species_info["excluded_from_refseq"].strip() or \
                        all([ex.strip() in REFERENCE_TAGS for ex in species_info["excluded_from_refseq"].split(";")]):
                        genome_type = "reference"
                    
                    else:
                        excluded = False
                        for ex in species_info["excluded_from_refseq"].split(";"):
                            if ex in EXCLUDE_TAGS:
                                excluded = True
                                break
                        
                        if not excluded:
                            genome_type = "mag"
                    
                    species_info["ftp_filepath"] = genome_url

                    local_filename = os.path.splitext(os.path.splitext(os.path.basename(genome_url))[0])[0]
                    species_info["local_filename"] = os.path.basename(local_filename)

                    species_info["genome_type"] = genome_type

                    assembly_summary[species_taxid].append(species_info)
    
    return assembly_summary


def get_genomes_in_ncbi(
    superkingdom: str,
    tmpdir: str,
    kingdom: Optional[str] = None,
    taxa_level_id: Optional[str] = None,
    taxa_level_name: Optional[str] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Retrieve links and taxonomic information about reference genomes and MAGs in NCBI GenBank

    :param superkingdom:    Archaea, Bacteria, Eukaryota, or Viruses
    :param tmpdir:          Path to the temporary folder
    :param kingdom:         A specific kingdom related to the superkingdom. Optional
    :param taxa_level_id:   Taxonomic level identifier (phylum, class, order, family, genus, or species)
    :param taxa_level_name: Name of the taxonomic level as appear in NCBI
    :return:                A dictionary with the genome IDs as keys and URL and taxonomic info as values
                            and optionally the name of the specified cluster
    """

    target_cluster = None

    if taxa_level_id and taxa_level_name:
        # Search for genomes belonging to a specific cluster
        target_cluster = "{}__{}".format(
            taxa_level_id.lower()[0],
            re.sub(r"_+", "_", re.sub(r"\W+", "_", taxa_level_name)).strip("_")
        )

    # Download the NCBI nodes and names dumps
    nodes_dmp, names_dmp = download_taxdump(TAXDUMP_URL, tmpdir)

    # Produce a mapping between NCBI tax IDs and full taxonomic labels
    taxa_map = ncbitax2lin(tmpdir, nodes_dmp, names_dmp, superkingdom=superkingdom, kingdom=kingdom)

    # Download and load the most recent NCBI GenBank Assembly Report table
    assembly_summary = get_assembly_summary(ASSEMBLY_SUMMARY_URL, tmpdir)

    ncbi_genomes = dict()

    # Get genome info from the assembly reporta table
    for species_taxid in assembly_summary:
        if species_taxid in taxa_map:
            taxonomy = taxa_map[species_taxid]

            if not target_cluster or "|{}|".format(target_cluster) in "{}|".format(taxonomy):
                for species_info in assembly_summary[species_taxid]:
                    ncbi_genomes[species_info["local_filename"]] = {
                        "type": species_info["genome_type"],
                        "taxonomy": taxonomy,
                        "excluded_from_refseq": species_info["excluded_from_refseq"] if species_info["excluded_from_refseq"].strip() else "na",
                        "url": species_info["ftp_filepath"]
                    }

    return ncbi_genomes, target_cluster


def urlretrieve_wrapper(url: str, filepath: str, retry: int = 5) -> Tuple[str, bool]:
    """
    Wrapper around urlretrieve

    :param url:         Input URL
    :param filepath:    Output file path
    :return:            The output file path
    """

    exists_and_passed_integrity = False

    while retry > 0 and not exists_and_passed_integrity:
        try:
            urlretrieve(url, filepath)
            
            # Check file integrity
            subprocess.check_call(
                [
                    "gzip",
                    "-t",
                    filepath,
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            exists_and_passed_integrity = True

        except Exception:
            if os.path.isfile(filepath):
                os.unlink(filepath)

            retry -= 1

    return filepath, exists_and_passed_integrity


def main() -> None:
    args = read_params()

    if (args.taxa_level_id and not args.taxa_level_name) or (args.taxa_level_name and not args.taxa_level_id):
        raise ValueError("\"--taxa-level-id\" must always be used in conjunction with \"--taxa-level-name\" and the other way around")

    os.makedirs(args.out_dir, exist_ok=True)

    if args.download:
        out_folder = "genomes" if not args.type else "{}s".format(args.type)

        genomes_dir = os.path.join(args.out_dir, out_folder)

        os.makedirs(genomes_dir, exist_ok=True)

    tmp_dir = os.path.join(args.out_dir, "tmp")

    os.makedirs(tmp_dir, exist_ok=True)

    # Retrieve genomes from NCBI
    ncbi_genomes, target_cluster = get_genomes_in_ncbi(
        args.superkingdom,
        tmp_dir,
        kingdom=args.kingdom,
        taxa_level_id=args.taxa_level_id,
        taxa_level_name=args.taxa_level_name,
    )

    if args.db_dir and os.path.isdir(args.db_dir):
        # Get the list of genome IDs in the database
        db_genomes = get_genomes_in_db(args.db_dir)

        # Remove genomes already in the database
        in_db = set(ncbi_genomes.keys()).intersection(db_genomes)

        for genome in in_db:
            del ncbi_genomes[genome]

    if ncbi_genomes:
        out_file_name = "genomes" if not args.type else "{}s".format(args.type)
        out_file_path = os.path.join(args.out_dir, "{}.tsv".format(out_file_name))

        exclude_genomes = list()

        if not os.path.isfile(out_file_path):
            with open(out_file_path, "w+") as genomes_table:
                genomes_table.write("# {} v{} ({})\n".format(TOOL_ID, __version__, __date__))
                genomes_table.write("# timestamp {}\n".format(datetime.datetime.utcnow()))
                genomes_table.write("# id\ttype\ttaxonomy\texcluded_from_refseq\turl\n")
        
        else:
            exclude_genomes = [
                line.strip().split("\t")[0] for line in open(out_file_path).readlines() if line.strip() and not line.strip().startswith("#")
            ]

        # Get the genomes of the same type as the input --type
        # Exclude genomes if they already appear in an existing output table
        # Exclude the unclassified genomes or consider them in case the input --type is "mag"
        genomes = list(
            filter(
                lambda genome: (ncbi_genomes[genome]["type"] == args.type or not args.type) and genome not in exclude_genomes and \
                               ("unclassified" not in ncbi_genomes[genome]["taxonomy"] or ("unclassified" in ncbi_genomes[genome]["taxonomy"] and args.type == "mag")),
                ncbi_genomes.keys()
            )
        )

        downloaded = list()

        print(
            "Retrieving {} genomes (Superkingdom \"{}\"; Kingdom \"{}\"; Cluster \"\"; Type \"{}\")".format(
                len(genomes),
                args.superkingdom,
                args.kingdom,
                target_cluster,
                args.type
            )
        )

        with mp.Pool(processes=args.nproc) as pool, tqdm.tqdm(total=len(genomes)) as pbar:
            # Wrapper around the update function of tqdm
            def progress(*args):
                pbar.update()

            # Process input genomes
            jobs = [
                pool.apply_async(
                    urlretrieve_wrapper,
                    args=(
                        ncbi_genomes[genome]["url"],
                        os.path.join(genomes_dir, os.path.basename(ncbi_genomes[genome]["url"]))
                    ),
                    callback=progress,
                )
                for genome in genomes
            ]

            # Get results from jobs
            for job in jobs:
                filepath, exists = job.get()

                if exists:
                    downloaded.append(filepath)

        if downloaded:
            with open(out_file_path, "a+") as genomes_table:
                for filepath in downloaded:
                    genome = os.path.splitext(os.path.splitext(os.path.basename(filepath))[0])[0]

                    genomes_table.write(
                        "{}\t{}\t{}\t{}\t{}\n".format(
                            genome,
                            ncbi_genomes[genome]["type"],
                            ncbi_genomes[genome]["taxonomy"],
                            ncbi_genomes[genome]["excluded_from_refseq"],
                            ncbi_genomes[genome]["url"]
                        )
                    )


if __name__ == "__main__":
    main()
