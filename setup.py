import platform
import sys

import setuptools

from metasbt import __version__

if platform.system() not in ["Darwin", "Linux"]:
    sys.exit("MetaSBT does not work on {} platforms".format(platform.system()))

if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and sys.version_info[1] < 9):
    sys.exit(
        "MetaSBT requires Python 3.9 or higher. Your current Python version is {}.{}.{}\n".format(
            sys.version_info[0], sys.version_info[1], sys.version_info[2]
        )
    )

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
        "metagenome-assembled genomes with Sequence Bloom Trees"
    ),
    download_url="https://pypi.org/project/MetaSBT/",
    entry_points={"console_scripts": ["metasbt=metasbt.metasbt:run"]},
    install_requires=[
        "biopython>=1.85",
        "fastcluster<1.3.0",
        "numpy==1.26.4",
        "packaging>=25.0",
        "requests>=2.32.4",
        "scipy>=1.13.1",
        "tabulate>=0.9.0",
        "tqdm>=4.67.1",
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
    platforms=["Linux", "Mac OSX"],
    project_urls={
        "Databases": "https://github.com/cumbof/MetaSBT-DBs",
        "Discussions": "https://github.com/cumbof/MetaSBT/discussions",
        "Issues": "https://github.com/cumbof/MetaSBT/issues",
        "Source": "https://github.com/cumbof/MetaSBT",
        "Wiki": "https://github.com/cumbof/MetaSBT/wiki",
    },
    python_requires=">=3.9",
    scripts=[
        "scripts/get_ncbi_genomes.py",
    ],
    include_package_data=True,
    url="http://github.com/cumbof/MetaSBT",
    version=__version__,
    zip_safe=False,
)
