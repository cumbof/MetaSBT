#!/usr/bin/env python3
"""Retrieve reference genomes and metagenome-assembled genomes for a specific superkingdom and kingdom from NCBI GenBank.
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.2.0"
__date__ = "Jul 13, 2026"

import argparse as ap
import datetime
import gzip
import json
import multiprocessing as mp
import os
import re
import subprocess
import tarfile
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple
from urllib.request import Request, urlopen, urlretrieve

import tqdm
import numpy as np

TOOL_ID = "get_ncbi_genomes"

# Define the list of dependencies
DEPENDENCIES = [
    "gzip",
    "ncbitax2lin",
]

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

# NCBI rate-limits concurrent connections and starts returning HTTP 503 errors when too many
# downloads run in parallel. Cap the number of parallel download workers to a server-friendly value
# regardless of the requested --nproc, otherwise the vast majority of the downloads silently fail.
MAX_DOWNLOAD_PROCESSES = 16

# Ordered taxonomic assembly levels, from the most to the least complete
ASSEMBLY_LEVELS = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]

# Define the url to the NCBI Datasets API endpoint reporting the assembly metadata.
# NCBI runs CheckM on the prokaryotic assemblies and publishes the resulting completeness and
# contamination estimates under "checkm_info". They are not part of the Assembly Summary table,
# so they must be retrieved here. Nothing is computed locally: CheckM is never invoked
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/
DATASETS_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report"

# Number of assembly accessions per Datasets API request
CHECKM_BATCH_SIZE = 1000

# The Datasets API allows 3 requests per second without an API key and 10 with a key.
# Wait at least this long between two consecutive requests to stay under the limit
DATASETS_API_DELAY = 0.35
DATASETS_API_DELAY_WITH_KEY = 0.11


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
        "--max-genomes-per-species",
        type=int,
        default=0,
        dest="max_genomes_per_species",
        help=(
            "Limit the number of downloaded genomes per species. "
            "Not limited by default (--max-genomes-per-species 0)"
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
        action="append",
        required=True,
        choices=["Archaea", "Bacteria", "Eukaryota", "Viruses"],
        help="Specify one or more superkingdoms (repeatable).",
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
        help=(
            "Retrieve reference genomes or metagenome-assembled genomes (MAGs) only. "
            "Genomes are categorized as references or MAGs according to their tags under the \"excluded_from_refseq\" "
            "column in the NCBI GenBank Assembly Summary Report table"
        )
    )
    p.add_argument(
        "--reference-genome",
        action="store_true",
        default=False,
        dest="reference_genome",
        help=(
            "Retrieve genomes marked as \"reference genome\" under the \"refseq_category\" column in the "
            "NCBI GenBank Assembly Suppary Report table. Can be used if --type is not provided"
        )
    )
    p.add_argument(
        "--representative-genome",
        action="store_true",
        default=False,
        dest="representative_genome",
        help=(
            "Retrieve genomes marked as \"representative genome\" under the \"refseq_category\" column in the "
            "NCBI GenBank Assembly Suppary Report table. Can be used if --type is not provided"
        )
    )
    p.add_argument(
        "--full-only",
        action="store_true",
        default=False,
        dest="full_only",
        help=(
            "Retrieve fully-represented genomes only (i.e. \"genome_rep\" is \"Full\" in the NCBI GenBank "
            "Assembly Summary table). Note that this does not filter on the assembly level: draft genomes "
            "(Contig, Scaffold) can still be \"Full\". Use \"--assembly-level\" to filter on the assembly level"
        )
    )
    p.add_argument(
        "--assembly-level",
        action="append",
        dest="assembly_level",
        choices=ASSEMBLY_LEVELS,
        help=(
            "Retrieve genomes with a specific assembly level only (repeatable). "
            "E.g. \"--assembly-level 'Complete Genome' --assembly-level Chromosome\" to exclude draft genomes. "
            "All assembly levels are retrieved by default"
        )
    )
    p.add_argument(
        "--min-completeness",
        type=float,
        default=0.0,
        dest="min_completeness",
        help=(
            "Retrieve genomes whose CheckM completeness is greater than or equal to this percentage. "
            "The completeness estimates are the ones precomputed by NCBI and retrieved through the "
            "Datasets API (CheckM is never run locally). Note that the assembly level says nothing about "
            "the completeness of a genome: a fragmented assembly can still be complete, while a "
            "\"Complete Genome\" is not necessarily so. Disabled by default (--min-completeness 0.0)"
        )
    )
    p.add_argument(
        "--max-contamination",
        type=float,
        default=100.0,
        dest="max_contamination",
        help=(
            "Retrieve genomes whose CheckM contamination is lower than or equal to this percentage. "
            "See \"--min-completeness\". Disabled by default (--max-contamination 100.0)"
        )
    )
    p.add_argument(
        "--require-checkm",
        action="store_true",
        default=False,
        dest="require_checkm",
        help=(
            "Discard genomes for which NCBI does not report any CheckM estimate. NCBI computes CheckM "
            "on the prokaryotic assemblies only, so this always discards every eukaryotic and viral genome. "
            "Genomes with no CheckM estimate are retained by default"
        )
    )
    p.add_argument(
        "--include-superseded",
        action="store_true",
        default=False,
        dest="include_superseded",
        help=(
            "Also retrieve genomes whose \"version_status\" is not \"latest\" in the NCBI GenBank Assembly "
            "Summary table (i.e. assemblies that have been replaced by a newer version or suppressed). "
            "Superseded assemblies are discarded by default"
        )
    )
    p.add_argument(
        "--api-key",
        type=str,
        default=os.environ.get("NCBI_API_KEY"),
        dest="api_key",
        help=(
            "NCBI API key, used to raise the Datasets API rate limit from 3 to 10 requests per second "
            "while retrieving the CheckM estimates. Defaults to the NCBI_API_KEY environment variable"
        )
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version='"{}" version {} ({})'.format(TOOL_ID, __version__, __date__),
        help='Print the "{}" version and exit'.format(TOOL_ID),
    )
    return p.parse_args()


def level_name(current_level: str, prev_level: str) -> str:
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

    return "{}{}".format(level_prefix, level_suffix)


def download_taxdump(taxdump_url: str, folder_path: os.path.abspath) -> Tuple[os.path.abspath, os.path.abspath]:
    """Download and extract the NCBI taxdump tarball.

    Parameters
    ----------
    taxdump_url : str
        URL to the NCBI taxdump.
    folder_path : os.path.abspath
        Path to the folder in which the taxdump tarball will be unpacked.

    Raises
    ------
    Exception
        If it is unable to retrieve data from the remote location.

    Returns
    -------
    tuple
        The nodes.dmp and names.dmp file paths.
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
    tmpdir: os.path.abspath,
    nodes_dmp: os.path.abspath,
    names_dmp: os.path.abspath,
    superkingdom: Optional[str]=None,
    kingdom: Optional[str]=None
) -> Dict[str, str]:
    """Run ncbitax2lin over nodes and names dumps and produce the mapping between NCBI tax IDs and ful taxonomic labels.

    Parameters
    ----------
    tmpdir : os.path.abspath
        Path to the tmp directory.
    nodes_dmp : os.path.abspath
        Path to the NCBI nodes dump.
    names_dmp : os.path.abspath
        Path to the NCBI names dump.
    superkingdom : str
        Filter results on this superkingdom only.
    kingdom : str
        Filter results on this kingdom only.

    Returns
    -------
    dict
        Dictionary with the mapping between NCBI tax IDs and full taxonomic labels.
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
        # Load the first line as header and search for "domain" and "kingdom" columns
        header = ncbi_table.readline().split(",")
        kingdom_pos = header.index("kingdom")

        # domain: Archaea, Bacteria, Eukaryota
        # acellular root: Viruses
        superkingdom_pos = header.index("acellular root") if superkingdom == "Viruses" else header.index("domain")

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

                    taxa_map[line_split[0]] = label

    return taxa_map


def get_assembly_summary(
    assembly_summary_url: str,
    tmpdir: os.path.abspath,
    full_only: bool=False,
    include_superseded: bool=False
) -> Dict[str, List[Dict[str, str]]]:
    """Download and load the last available NCBI GenBank Assembly Report table.

    Parameters
    ----------
    assembly_summary_url : str
        URL to the NCBI GenBank Assembly Report table.
    tmpdir : os.path.abspath
        Path to the tmp folder.
    full_only : bool, default False
        Retrieve full genomes only.
    include_superseded : bool, default False
        Also retrieve the assemblies whose "version_status" is not "latest".

    Returns
    -------
    dict
        Dictionary with genomes info.
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

                # The "ftp_path" column in the NCBI Assembly Summary table now ends with a trailing slash.
                # Strip it, otherwise os.path.basename() returns an empty string and both the download URL
                # and the local file name end up malformed (every genome would collapse to the same name).
                ftp_path = line_split[header.index("ftp_path")].strip().rstrip("/")

                if not ftp_path or ftp_path == "na":
                    # No downloadable genome is available for this assembly
                    continue

                genome_url = os.path.join(ftp_path, "{}_genomic.fna.gz".format(os.path.basename(ftp_path)))

                if species_taxid:
                    species_info = dict()

                    for h in header:
                        species_info[h] = line_split[header.index(h)].strip()

                    if full_only and species_info["genome_rep"] == "Partial":
                        # genome_rep could be Full or Partial only
                        # Skip the current iteration in case of non-Full genomes (if full_only=True)
                        continue

                    if not include_superseded and species_info.get("version_status", "latest") != "latest":
                        # The Assembly Summary table also lists the assemblies that have been replaced by a
                        # newer version or suppressed. Their "version_status" is anything but "latest" and
                        # they must not end up in the database next to the version that superseded them
                        continue

                    genome_type = "na"

                    if not species_info["excluded_from_refseq"].strip() or species_info["excluded_from_refseq"].strip() == "na" or \
                        all([ex.strip().lower() in REFERENCE_TAGS for ex in species_info["excluded_from_refseq"].split(";") if ex.strip()]):
                        genome_type = "reference"

                    else:
                        excluded = False
                        for ex in species_info["excluded_from_refseq"].split(";"):
                            if ex.strip().lower() in EXCLUDE_TAGS:
                                excluded = True
                                break

                        if not excluded:
                            genome_type = "mag"

                    species_info["ftp_filepath"] = genome_url

                    local_filename = os.path.splitext(os.path.splitext(os.path.basename(genome_url))[0])[0]
                    species_info["local_filename"] = os.path.basename(local_filename)

                    species_info["genome_type"] = genome_type

                    if species_taxid not in assembly_summary:
                        assembly_summary[species_taxid] = list()

                    assembly_summary[species_taxid].append(species_info)

    return assembly_summary


def query_datasets_api(
    accessions: List[str],
    api_key: Optional[str]=None,
    retry: int=5
) -> Dict[str, Tuple[float, float]]:
    """Query the NCBI Datasets API for the CheckM estimates of a batch of assembly accessions.

    Parameters
    ----------
    accessions : list
        List of assembly accessions (e.g. GCA_000008005.1).
    api_key : str, optional
        NCBI API key.
    retry : int, default 5
        Number of attempts before giving up on a request.

    Raises
    ------
    Exception
        If the Datasets API keeps failing after `retry` attempts.

    Returns
    -------
    dict
        Dictionary with the assembly accessions as keys and the (completeness, contamination)
        tuples as values. Accessions for which NCBI reports no CheckM estimate are not in the result.
    """

    checkm_info = dict()

    page_token = None

    while True:
        payload = {
            "accessions": accessions,
            "returned_content": "COMPLETE",
            "page_size": CHECKM_BATCH_SIZE,
        }

        if page_token:
            payload["page_token"] = page_token

        headers = {"Content-Type": "application/json", "Accept": "application/json"}

        if api_key:
            headers["api-key"] = api_key

        response_data = None

        attempt = 0

        while attempt < retry:
            try:
                request = Request(
                    DATASETS_API_URL,
                    data=json.dumps(payload).encode("utf-8"),
                    headers=headers,
                    method="POST",
                )

                with urlopen(request, timeout=300) as response:
                    response_data = json.loads(response.read().decode("utf-8"))

                break

            except Exception:
                attempt += 1

                if attempt >= retry:
                    raise Exception(
                        "Unable to retrieve the CheckM estimates from the NCBI Datasets API\n{}".format(DATASETS_API_URL)
                    )

                # Same exponential backoff as the genome downloads: NCBI throttles with HTTP 503
                time.sleep(min(2 ** attempt, 30))

        for report in response_data.get("reports", list()) or list():
            accession = report.get("accession")

            # NCBI only runs CheckM on the prokaryotic assemblies, so "checkm_info" is
            # missing for the eukaryotic and viral ones
            report_checkm = report.get("checkm_info") or dict()

            completeness = report_checkm.get("completeness")
            contamination = report_checkm.get("contamination")

            if accession and completeness is not None and contamination is not None:
                checkm_info[accession] = (float(completeness), float(contamination))

        page_token = response_data.get("next_page_token")

        if not page_token:
            break

    return checkm_info


def get_checkm_info(
    accessions: Iterable[str],
    tmpdir: os.path.abspath,
    api_key: Optional[str]=None
) -> Dict[str, Optional[Tuple[float, float]]]:
    """Retrieve the CheckM completeness and contamination estimates precomputed by NCBI.

    The estimates are cached on disk so that a subsequent run does not query the Datasets API again
    for the same assemblies. Accessions with no CheckM estimate are cached as well, otherwise every
    run would keep asking for them.

    Parameters
    ----------
    accessions : iterable
        Assembly accessions (e.g. GCA_000008005.1).
    tmpdir : os.path.abspath
        Path to the tmp folder.
    api_key : str, optional
        NCBI API key.

    Returns
    -------
    dict
        Dictionary with the assembly accessions as keys and the (completeness, contamination) tuples
        as values, or None for the accessions with no CheckM estimate.
    """

    checkm_info: Dict[str, Optional[Tuple[float, float]]] = dict()

    checkm_filepath = os.path.join(tmpdir, "checkm_info.tsv")

    if os.path.isfile(checkm_filepath):
        with open(checkm_filepath) as checkm_table:
            for line in checkm_table:
                line = line.strip()

                if line and not line.startswith("#"):
                    line_split = line.split("\t")

                    if len(line_split) >= 3:
                        checkm_info[line_split[0]] = (
                            None if line_split[1] == "na" else (float(line_split[1]), float(line_split[2]))
                        )

    missing = sorted({accession for accession in accessions if accession and accession not in checkm_info})

    if not missing:
        return checkm_info

    print("Retrieving the CheckM estimates of {} assemblies from the NCBI Datasets API".format(len(missing)))

    delay = DATASETS_API_DELAY_WITH_KEY if api_key else DATASETS_API_DELAY

    with open(checkm_filepath, "a+") as checkm_table:
        if os.path.getsize(checkm_filepath) == 0:
            checkm_table.write("# accession\tcompleteness\tcontamination\n")

        for position in tqdm.tqdm(range(0, len(missing), CHECKM_BATCH_SIZE)):
            batch = missing[position: position + CHECKM_BATCH_SIZE]

            batch_checkm_info = query_datasets_api(batch, api_key=api_key)

            for accession in batch:
                # Cache the accessions with no CheckM estimate as None so that they are not queried again
                info = batch_checkm_info.get(accession)

                checkm_info[accession] = info

                checkm_table.write(
                    "{}\t{}\t{}\n".format(
                        accession,
                        "na" if info is None else info[0],
                        "na" if info is None else info[1],
                    )
                )

            checkm_table.flush()

            time.sleep(delay)

    return checkm_info


def get_genomes_in_ncbi(
    superkingdom: str,
    tmpdir: os.path.abspath,
    full_only: bool=False,
    include_superseded: bool=False,
    kingdom: Optional[str]=None,
    taxa_level_id: Optional[str]=None,
    taxa_level_name: Optional[str]=None,
) -> Dict[str, Dict[str, str]]:
    """Retrieve links and taxonomic information about reference genomes and MAGs in NCBI GenBank.

    Parameters
    ----------
    superkingdom : str
        Archaea, Bacteria, Eukaryota, or Viruses.
    tmpdir : os.path.abspath
        Path to the temporary folder.
    full_only : bool, default False
        Retrieve full genomes only.
    include_superseded : bool, default False
        Also retrieve the assemblies whose "version_status" is not "latest".
    kingdom : str, optional
        A specific kingdom related to the superkingdom. Optional.
    taxa_level_id : str, optional
        Taxonomic level identifier (phylum, class, order, family, genus, or species).
    taxa_level_name : str, optional
        Name of the taxonomic level as appear in NCBI.

    Returns
    -------
    dict
        A dictionary with the genome IDs as keys and URL and taxonomic info as values
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
    assembly_summary = get_assembly_summary(
        ASSEMBLY_SUMMARY_URL, tmpdir, full_only=full_only, include_superseded=include_superseded
    )

    ncbi_genomes = dict()

    # Get genome info from the assembly reporta table
    for species_taxid in assembly_summary:
        if species_taxid in taxa_map:
            taxonomy = taxa_map[species_taxid]

            if not target_cluster or "|{}|".format(target_cluster) in "{}|".format(taxonomy):
                for species_info in assembly_summary[species_taxid]:
                    ncbi_genomes[species_info["local_filename"]] = {
                        "type": species_info["genome_type"],
                        "refseq_category": species_info["refseq_category"],
                        "taxonomy": taxonomy,
                        "excluded_from_refseq": species_info["excluded_from_refseq"] if species_info["excluded_from_refseq"].strip() else "na",
                        "url": species_info["ftp_filepath"],
                        "assembly_level": species_info["assembly_level"],
                        # The CheckM estimates are not part of the Assembly Summary table. They are
                        # retrieved from the Datasets API, which is keyed on the assembly accession
                        "assembly_accession": species_info["assembly_accession"],
                        "completeness": "na",
                        "contamination": "na",
                    }

    return ncbi_genomes, target_cluster


def urlretrieve_wrapper(url: str, filepath: os.path.abspath, retry: int=5) -> Tuple[str, bool]:
    """Just a wrapper around urlretrieve.

    Parameters
    ----------
    url : str
        Input URL.
    filepath : os.path.abspath
        Output file path.

    Returns
    -------
    tuple
        A tuple with the output file path and a boolean (True if it passes the integrity check).
    """

    exists_and_passed_integrity = False

    attempt = 0

    while attempt < retry and not exists_and_passed_integrity:
        try:
            if not os.path.isfile(filepath):
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

            attempt += 1

            if attempt < retry:
                # NCBI returns HTTP 503 errors under heavy concurrency.
                # Back off with an exponentially increasing delay (capped at 30 seconds)
                # so that the retries do not all hit the same throttling window.
                time.sleep(min(2 ** attempt, 30))

    return filepath, exists_and_passed_integrity


def _urlretrieve_wrapper_star(args: Tuple[str, os.path.abspath]) -> Tuple[str, bool]:
    url, filepath = args
    return urlretrieve_wrapper(url, filepath)


def main() -> None:
    args = read_params()

    if (args.taxa_level_id and not args.taxa_level_name) or (args.taxa_level_name and not args.taxa_level_id):
        raise ValueError("\"--taxa-level-id\" must always be used in conjunction with \"--taxa-level-name\" and the other way around")

    if (args.type and args.reference_genome) or (args.type and args.representative_genome):
        raise ValueError("\"--reference-genome\" and \"--representative-genome\" cannot be used in conjunction with \"--type\"")

    if not 0.0 <= args.min_completeness <= 100.0:
        raise ValueError("\"--min-completeness\" must be a percentage in [0.0, 100.0]")

    if not 0.0 <= args.max_contamination <= 100.0:
        raise ValueError("\"--max-contamination\" must be a percentage in [0.0, 100.0]")

    os.makedirs(args.out_dir, exist_ok=True)

    if args.download:
        out_folder = "genomes" if not args.type else "{}s".format(args.type)

        genomes_dir = os.path.join(args.out_dir, out_folder)

        os.makedirs(genomes_dir, exist_ok=True)

    tmp_dir = os.path.join(args.out_dir, "tmp")

    os.makedirs(tmp_dir, exist_ok=True)

    # Retrieve genomes from NCBI
    ncbi_genomes = dict()

    for superkingdom in args.superkingdom:
        superkingdom_genomes, target_cluster = get_genomes_in_ncbi(
            superkingdom,
            tmp_dir,
            full_only=args.full_only,
            include_superseded=args.include_superseded,
            kingdom=args.kingdom,
            taxa_level_id=args.taxa_level_id,
            taxa_level_name=args.taxa_level_name,
        )

        ncbi_genomes.update(superkingdom_genomes)

    if ncbi_genomes:
        out_file_name = "genomes" if not args.type else "{}s".format(args.type)
        out_file_path = os.path.join(args.out_dir, "{}.tsv".format(out_file_name))

        exclude_genomes = list()

        if not os.path.isfile(out_file_path):
            with open(out_file_path, "w+") as genomes_table:
                genomes_table.write("# {} v{} ({})\n".format(TOOL_ID, __version__, __date__))
                genomes_table.write("# timestamp {}\n".format(datetime.datetime.utcnow()))
                genomes_table.write("# id\ttype\ttaxonomy\texcluded_from_refseq\tassembly_level\turl\tcompleteness\tcontamination\n")

        else:
            exclude_genomes = [
                line.strip().split("\t")[0] for line in open(out_file_path).readlines() if line.strip() and not line.strip().startswith("#")
            ]

        species = dict()

        for genome in ncbi_genomes.keys():
            selected = False

            if not args.type and (args.reference_genome or args.representative_genome):
                # Get genomes marked as "reference genome" or "representative genome" in the Assembly Summary table
                if args.reference_genome and ncbi_genomes[genome]["refseq_category"] == "reference genome":
                    selected = True

                elif args.representative_genome and ncbi_genomes[genome]["refseq_category"] == "representative genome":
                    selected = True

                if selected:
                    # Override the genome type
                    ncbi_genomes[genome]["type"] = ncbi_genomes[genome]["refseq_category"]

            elif (ncbi_genomes[genome]["type"] == args.type or not args.type) and genome not in exclude_genomes and \
                ("unclassified" not in ncbi_genomes[genome]["taxonomy"] or ("unclassified" in ncbi_genomes[genome]["taxonomy"] and args.type == "mag")):
                # Get genomes of the same type as the input --type
                # Exclude genomes if they already appear in an existing output table
                # Exclude unclassified genomes or consider them in case the input --type is "mag"
                selected = True

            if selected and args.assembly_level and ncbi_genomes[genome]["assembly_level"] not in args.assembly_level:
                # Discard genomes whose assembly level is not in the requested set
                selected = False

            if selected:
                taxonomy = ncbi_genomes[genome]["taxonomy"]

                if taxonomy not in species:
                    species[taxonomy] = list()

                species[taxonomy].append(genome)

        if args.min_completeness > 0.0 or args.max_contamination < 100.0 or args.require_checkm:
            # The assembly level measures how contiguous an assembly is, not how complete it is.
            # A fragmented assembly loses just the k-mers spanning the contig breaks, which is
            # negligible, while a genome missing a fraction of its k-mers inflates every distance
            # measured against it and widens the boundaries of the cluster it lands in.
            # Filter on the CheckM estimates instead, before capping the number of genomes per
            # species, so that a species does not spend its quota on low-quality assemblies
            checkm_info = get_checkm_info(
                {ncbi_genomes[genome]["assembly_accession"] for sp in species for genome in species[sp]},
                tmp_dir,
                api_key=args.api_key,
            )

            low_quality = 0
            no_checkm = 0

            for sp in list(species.keys()):
                retained = list()

                for genome in species[sp]:
                    info = checkm_info.get(ncbi_genomes[genome]["assembly_accession"])

                    if info is None:
                        # NCBI does not report any CheckM estimate for this assembly
                        # (e.g. every eukaryotic and viral genome)
                        if args.require_checkm:
                            no_checkm += 1

                        else:
                            retained.append(genome)

                        continue

                    completeness, contamination = info

                    ncbi_genomes[genome]["completeness"] = completeness
                    ncbi_genomes[genome]["contamination"] = contamination

                    if completeness < args.min_completeness or contamination > args.max_contamination:
                        low_quality += 1

                    else:
                        retained.append(genome)

                if retained:
                    species[sp] = retained

                else:
                    del species[sp]

            print(
                "{} genomes discarded (CheckM completeness < {} or contamination > {})".format(
                    low_quality, args.min_completeness, args.max_contamination
                )
            )

            if args.require_checkm:
                print("{} genomes discarded (no CheckM estimate available)".format(no_checkm))

        if args.max_genomes_per_species > 0:
            # Limit the number of genomes per species
            for sp in species:
                # Define the priority order for assembly levels
                priority_order = ASSEMBLY_LEVELS

                # Group the initial genome IDs by their assembly level
                grouped_by_level = {level: list() for level in priority_order}

                for gid in species[sp]:
                    level = ncbi_genomes[gid].get("assembly_level")

                    if level in grouped_by_level:
                        grouped_by_level[level].append(gid)

                # Subsampling genomes (with priorities)
                selected_genomes = list()

                # Iterate through the levels in order of priority
                for level in priority_order:
                    # Determine how many more genomes we still need
                    needed = args.max_genomes_per_species - len(selected_genomes)

                    if needed <= 0:
                        break  # Stop if we have already selected enough genomes

                    # Get the available genomes at the current priority level
                    available_at_level = grouped_by_level.get(level, list())

                    # Always use the same seed for reproducibility
                    rng = np.random.default_rng(0)

                    # Shuffle for random selection if we don't take all of them
                    rng.shuffle(available_at_level)

                    # Decide how many genomes to take from this group
                    num_to_take = min(needed, len(available_at_level))

                    # Add the selected genomes to our final list
                    selected_genomes.extend(available_at_level[:num_to_take])

                species[sp] = selected_genomes

        genomes = list()

        for sp in species:
            genomes += species[sp]

        print(
            "{} genomes (Superkingdom \"{}\"; Kingdom \"{}\"; Cluster \"{}\"; Type \"{}\")".format(
                len(genomes),
                ",".join(args.superkingdom),
                args.kingdom,
                target_cluster,
                args.type
            )
        )

        if args.download:
            # Keep track of the genomes that could not be downloaded so that they are
            # reported instead of being silently dropped from the output table
            failed = list()

            def record_genome(genome: str) -> None:
                with open(out_file_path, "a+") as genomes_table:
                    genomes_table.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            genome,
                            ncbi_genomes[genome]["type"],
                            ncbi_genomes[genome]["taxonomy"],
                            ncbi_genomes[genome]["excluded_from_refseq"],
                            ncbi_genomes[genome]["assembly_level"],
                            ncbi_genomes[genome]["url"],
                            ncbi_genomes[genome]["completeness"],
                            ncbi_genomes[genome]["contamination"]
                        )
                    )

            # Cap the number of parallel downloads, otherwise NCBI throttles the connections
            # with HTTP 503 errors and most of the downloads fail
            download_processes = min(args.nproc, MAX_DOWNLOAD_PROCESSES)

            if args.nproc > download_processes:
                print(
                    "Warning: capping the number of parallel downloads to {} (out of the requested {}) "
                    "to avoid NCBI rate-limiting".format(download_processes, args.nproc)
                )

            if download_processes > 1:
                with mp.Pool(processes=download_processes) as pool:
                    args_list = [(ncbi_genomes[genome]["url"], os.path.join(genomes_dir, os.path.basename(ncbi_genomes[genome]["url"]))) for genome in genomes]
                    for filepath, exists in tqdm.tqdm(pool.imap_unordered(_urlretrieve_wrapper_star, args_list), total=len(args_list)):
                        genome = os.path.splitext(os.path.splitext(os.path.basename(filepath))[0])[0]
                        if exists:
                            record_genome(genome)
                        else:
                            failed.append(genome)
            else:
                for genome in tqdm.tqdm(genomes):
                    filepath, exists = urlretrieve_wrapper(ncbi_genomes[genome]["url"], os.path.join(genomes_dir, os.path.basename(ncbi_genomes[genome]["url"])))
                    if exists:
                        record_genome(genome)
                    else:
                        failed.append(genome)

            if failed:
                # Dump the list of genomes that could not be retrieved. The same command can be
                # re-run to retry the missing genomes only (the ones already in the output table are skipped)
                failed_file_path = os.path.join(args.out_dir, "{}_failed.txt".format(out_file_name))

                with open(failed_file_path, "w+") as failed_file:
                    for genome in failed:
                        failed_file.write("{}\t{}\n".format(genome, ncbi_genomes[genome]["url"]))

                print(
                    "{} out of {} genomes could not be downloaded (NCBI throttling or missing files).\n"
                    "The list of failed genomes has been written to {}\n"
                    "Re-run the same command to retry the missing genomes only".format(
                        len(failed), len(genomes), failed_file_path
                    )
                )

            else:
                print("All {} genomes have been successfully downloaded".format(len(genomes)))

        elif genomes:
            with open(out_file_path, "a+") as genomes_table:
                for genome in genomes:
                    genomes_table.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            genome,
                            ncbi_genomes[genome]["type"],
                            ncbi_genomes[genome]["taxonomy"],
                            ncbi_genomes[genome]["excluded_from_refseq"],
                            ncbi_genomes[genome]["assembly_level"],
                            ncbi_genomes[genome]["url"],
                            ncbi_genomes[genome]["completeness"],
                            ncbi_genomes[genome]["contamination"]
                        )
                    )


if __name__ == "__main__":
    main()
