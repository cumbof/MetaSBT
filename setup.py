import sys, setuptools

if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and  sys.version_info[1] < 8):
    sys.stdout.write("MetaSBT requires Python 3.8 or higher. Your Python your current Python version is {}.{}.{}"
                     .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

setuptools.setup(name="MetaSBT",
                 version="0.1.0-alpha",
                 author="Fabio Cumbo",
                 author_email="fabio.cumbo@gmail.com",
                 url="http://github.com/cumbof/MetaSBT",
                 license="MIT",
                 license_files=[
                    "LICENSE"
                 ],
                 packages=setuptools.find_packages(),
                 package_data = {
                     "metasbt": [
                        "requirements.txt"
                     ],
                 },
                 scripts=[
                     "scripts/esearch_txid.sh",
                     "scripts/howdesbt_index.sh",
                     "scripts/uniform_inputs.sh"
                 ],
                 entry_points={
                     "console_scripts": [
                        "metasbt = metasbt.metasbt:main"
                     ]
                 },
                 description=("A scalable framework for automatically indexing microbial genomes and accurately characterizing "
                              "metagenome-assembled genomes with sequence bloom trees"),
                 long_description=open("README.md").read(),
                 long_description_content_type="text/markdown",
                 install_requires=[
                     "ncbi-genome-download",
                     "ncbitax2lin",
                     "numpy",
                     "requests",
                     "tqdm"
                 ],
                 python_requires='>=3.8',
                 zip_safe=False,
                 classifiers=[
                     "Environment :: Console",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: MIT License",
                     "Programming Language :: Python :: 3",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"
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
                     "sequence bloom trees"
                 ])
