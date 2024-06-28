import errno
import os
import platform
import sys

import setuptools

from metasbt.metadata import __version__

if platform.system() not in ["Darwin", "Linux"]:
    sys.exit("MetaSBT does not work on {} platforms".format(platform.system()))

if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and sys.version_info[1] < 8):
    sys.exit(
        "MetaSBT requires Python 3.8 or higher. Your current Python version is {}.{}.{}\n".format(
            sys.version_info[0], sys.version_info[1], sys.version_info[2]
        )
    )

REQUIREMENTS = os.path.join(os.path.dirname(os.path.realpath(__file__)), "metasbt", "requirements.txt")

if not os.path.isfile(REQUIREMENTS):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), REQUIREMENTS)

setuptools.setup(
    author="Fabio Cumbo",
    author_email="fabio.cumbo@gmail.com",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    description=(
        "A scalable framework for automatically indexing microbial genomes and accurately characterizing "
        "metagenome-assembled genomes with sequence bloom trees"
    ),
    download_url="https://pypi.org/project/MetaSBT/",
    entry_points={"console_scripts": ["metasbt=metasbt.metasbt:main"]},
    install_requires=[
        requirement.strip() for requirement in open(REQUIREMENTS).readlines() if requirement.strip()
    ],
    keywords=[
        "bioinformatics",
        "characterization",
        "classification",
        "genomes",
        "metagenome-assembled genomes",
        "metagenomics",
        "microbiome",
        "profiling",
        "sequence bloom trees",
    ],
    license="MIT",
    license_files=["LICENSE"],
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    name="MetaSBT",
    packages=setuptools.find_packages(),
    package_data={
        "metasbt": ["requirements.txt"],
    },
    platforms=["Linux", "Mac OSX"],
    project_urls={
        "Databases": "https://github.com/cumbof/MetaSBT-DBs",
        "Discussions": "https://github.com/cumbof/MetaSBT/discussions",
        "Issues": "https://github.com/cumbof/MetaSBT/issues",
        "Source": "https://github.com/cumbof/MetaSBT",
        "Wiki": "https://github.com/cumbof/MetaSBT/wiki",
    },
    python_requires=">=3.8",
    scripts=[
        "scripts/bf_sketch.py",
        "scripts/cluster.py",
        "scripts/esearch_txid.sh",
        "scripts/expand_fasta.py",
        "scripts/get_ncbi_genomes.py",
        "scripts/howdesbt_index.sh",
        "scripts/minimizers.sh",
        "scripts/uniform_inputs.sh",
        "scripts/uniform_strands.py",
    ],
    url="http://github.com/cumbof/MetaSBT",
    version=__version__,
    zip_safe=False,
)
