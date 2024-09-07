#!/usr/bin/env python3
"""A Python script to build a custom kraken database with MetaSBT clustered sequences.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Sep 7, 2024"

import gzip
import multiprocessing as mp
import os
import tempfile
from pathlib import Path

from Bio import SeqIO

TOOL_ID = "kraken"

# Define the list of dependencies
DEPENDENCIES = [
    "kraken2-build",
]

# Assume MetaSBT is installed
# Define the modules root directory
MODULES_DIR = os.path.join(os.path.dirname(metasbt.__file__), "modules")

# Define the path to utils.py
UTILS_FILEPATH = os.path.join(MODULES_DIR, "utils.py")

# Check whether utils.py exists
if not os.path.isfile(UTILS_FILEPATH):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), UTILS_FILEPATH)

# Add the modules dir to the system path
sys.path.append(MODULES_DIR)

# Finally import utils
import utils

NODE_TEMPLATE = "taxid\t|\tparent\t|\trank\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|"
NAME_TEMPLATE = "taxid\t|\tvar\t|\t\t|\tscientific name\t|"

TAXONOMIC_RANKING = [
    "superkingdom", 
    "phylum", 
    "class", 
    "order", 
    "family", 
    "genus", 
    "species"
]


def read_params():
    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="A Python script to build a custom kraken database with MetaSBT clustered sequences.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--db-dir",
        type=os.path.abspath,
        required=True,
        dest="db_dir",
        help="Path to the MetaSBT database root folder",
    )
    p.add_argument(
        "--metasbt-report",
        type=os.path.abspath,
        required=True,
        dest="metasbt_report",
        help="Path to MetaSBT report file with the list of clusters and their taxonomic labels",
    )
    p.add_argument(
        "--ncbi-names",
        type=os.path.abspath,
        required=True,
        dest="ncbi_names",
        help="Path to the NCBI names.dmp file",
    )
    p.add_argument(
        "--ncbi-nodes",
        type=os.path.abspath,
        required=True,
        dest="ncbi_nodes",
        help="Path to the NCBI nodes.dmp file",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()                


def process_cluster(lineage, ncbi_names, metasbt_names, db_name):
    # Lineages are always defined up to the species level
    lineage = lineage.split("|")

    # Genomes are always Gzip compressed
    genomes_files = Path(os.path.join(args.db_dir, os.sep.join(lineage), "genomes")).glob("*.gz")

    metasbt_cluster_name = lineage[-1][3:]

    metasbt_cluster_id = metasbt_names[metasbt_cluster_name] if metasbt_cluster_name in metasbt_names else ncbi_names[metasbt_cluster_name]

    for genome_file in genomes_files:
        with tempfile.NamedTemporaryFile() as tmp_genome:
            with open(tmp_genome.name, "wt") as tmp_genome_file:
                with gzip.open(genome_file, "rt") as genome:
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
            kraken2build = [
                "kraken2-build",
                "--add-to-library",
                tmp_genome.name,
                "--db",
                db_name
            ]

            run(kraken2build, silence=True)


def main() -> None:
    # Load command line parameters
    args = read_params()

    if not os.path.isdir(args.db_dir):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.db_dir)

    if not os.path.isfile(args.metasbt_report):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.metasbt_report)

    if not os.path.isfile(args.ncbi_names):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.ncbi_names)

    if not os.path.isfile(args.ncbi_nodes):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.ncbi_nodes)

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
                    ncbi_tax_scientific_name = line_split[1].strip()

                    # TODO scientific names must be fixed in the same way that MetaSBT does

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

    metasbt_names = dict()

    metasbt_nodes = dict()

    # Take track of the full taxonomic labels in MetaSBT
    metasbt_lineages = list()

    # Load MetaSBT clusters
    with open(args.metasbt_report) as metasbt_report:
        for line in metasbt_report:
            line = line.strip()

            if line:
                if not line.startswith("#"):
                    line_split = line.split("\t")

                    lineage = line_split[1].split("|")

                    if lineage[-1].startswith("s__"):
                        metasbt_lineages.append("|".join(lineage))

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
            metasbt_node = NAME_TEMPLATE.replace("taxid", node_id)
            metasbt_node = metasbt_node.replace("var", node_name)

            ncbi_names_table.write("{}\n".format(metasbt_node))

    with open(args.ncbi_nodes, "a+") as ncbi_nodes_table:
        for node_id, node_data in sorted(metasbt_nodes.items()):
            metasbt_node = NODE_TEMPLATE.replace("taxid", node_id)
            metasbt_node = metasbt_node.replace("parent", node_data["parent"])
            metasbt_node = metasbt_node.replace("rank", node_data["rank"])

            ncbi_nodes_table.write("{}\n".format(metasbt_node))

    with mp.Pool(processes=args.nproc) as pool:
        jobs = [pool.apply_async(process_cluster, args=(lineage, ncbi_names, metasbt_names, os.path.basename(args.db_dir),)) for lineage in metasbt_lineages]

        for job in jobs:
            job.get()
