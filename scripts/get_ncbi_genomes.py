#!/usr/bin/env python3
"""
Retrieve reference genomes and metagenome-assembled genomes for a specific 
superkingdom and kingdom from NCBI GenBank
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Apr 5, 2023"

import argparse as ap
import datetime
import gzip
import os
import re
import subprocess
import tarfile
from pathlib import Path
from typing import Dict, List, Optional

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

    # Search for references
    references = Path(db_dir).glob("**/references.txt")

    for filepath in references:
        genomes.extend(
            [ref.strip() for ref in open(str(filepath)).readlines() \
                if ref.strip() and not ref.strip().startswith("#")]
        )

    # Search for MAGs
    mags = Path(db_dir).glob("**/mags.txt")

    for filepath in mags:
        genomes.extend(
            [mag.strip() for mag in open(str(filepath)).readlines() \
                if mag.strip() and not mag.strip().startswith("#")]
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


def get_genomes_in_ncbi(
    superkingdom: str,
    tmpdir: str,
    kingdom: Optional[str] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Retrieve links and taxonomic information about reference genomes and MAGs in NCBI GenBank

    :param superkingdom:    Archaea, Bacteria, Eukaryota, or Viruses
    :param tmpdir:          Path to the temporary folder
    :param kingdom:         A specific kingdom related to the superkingdom. Optional
    :return:                A dictionary with the genome IDs as keys and URL and taxonomic info as values
    """

    taxdump_dir = os.path.join(tmpdir, "taxdump")

    # Retrieve nodes and names dumps
    nodes_dmp = os.path.join(taxdump_dir, "nodes.dmp")
    names_dmp = os.path.join(taxdump_dir, "names.dmp")

    if not os.path.isfile(nodes_dmp) or not os.path.isfile(names_dmp):
        os.makedirs(taxdump_dir, exist_ok=True)

        taxdump_filepath = os.path.join(taxdump_dir, os.path.basename(TAXDUMP_URL))
        
        if not os.path.isfile(taxdump_filepath):
            urlretrieve(TAXDUMP_URL, taxdump_filepath)

        # Decompress the archive
        with tarfile.open(taxdump_filepath, "r:gz") as tar:
            tar.extractall(taxdump_dir)

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

                    # Exclude unclassified taxonomic labels
                    if "unclassified" not in label:
                        tax_id = int(line_split[0])

                        taxa_map[tax_id] = label

    assembly_summary = dict()

    assembly_summary_filepath = os.path.join(tmpdir, os.path.basename(ASSEMBLY_SUMMARY_URL))

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

    ncbi_genomes = dict()

    for species_taxid in assembly_summary:
        if species_taxid in taxa_map:
            taxonomy = taxa_map[species_taxid]

            for species_info in assembly_summary[species_taxid]:
                ncbi_genomes[species_info["local_filename"]] = {
                    "type": species_info["genome_type"],
                    "taxonomy": taxonomy,
                    "excluded_from_refseq": species_info["excluded_from_refseq"] if species_info["excluded_from_refseq"].strip() else "na",
                    "url": species_info["ftp_filepath"]
                }

    return ncbi_genomes


def main() -> None:
    args = read_params()
    
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)
    
    if args.download:
        genomes_dir = os.path.join(args.out_dir, "genomes")

        if not os.path.isdir(genomes_dir):
            os.makedirs(genomes_dir, exist_ok=True)
    
    tmp_dir = os.path.join(args.out_dir, "tmp")

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    # Retrieve genomes from NCBI
    ncbi_genomes = get_genomes_in_ncbi(
        args.superkingdom,
        tmp_dir,
        kingdom=args.kingdom
    )

    if args.db_dir:
        if os.path.isdir(args.db_dir):
            # Get the list of genome IDs in the database
            db_genomes = get_genomes_in_db(args.db_dir)

            # Remove genomes already in the database
            in_db = set(ncbi_genomes.keys()).intersection(db_genomes)

            for genome in in_db:
                del ncbi_genomes[genome]

    if ncbi_genomes:
        with open(os.path.join(args.out_dir, "genomes.tsv"), "w+") as genomes_table:
            genomes_table.write("# {} v{} ({})\n".format(TOOL_ID, __version__, __date__))
            genomes_table.write("# timestamp {}\n".format(datetime.datetime.utcnow()))
            genomes_table.write("# id\ttype\ttaxonomy\texcluded_from_refseq\turl\n")

            # Download genomes
            for genome in ncbi_genomes:
                process = True

                if args.type:
                    if args.type != ncbi_genomes[genome]["type"]:
                        process = False

                if args.download and process:
                    try:
                        urlretrieve(ncbi_genomes[genome]["url"], os.path.join(genomes_dir, os.path.basename(ncbi_genomes[genome]["url"])))
                    
                    except Exception:
                        process = False

                if process:
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
