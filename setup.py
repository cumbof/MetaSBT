import sys, setuptools

if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and  sys.version_info[1] < 8):
    sys.stdout.write("meta-index requires Python 3.8 or higher. Your Python your current Python version is {}.{}.{}"
                     .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

setuptools.setup(name="meta-index",
                 version="0.1.0",
                 author="Fabio Cumbo",
                 author_email="fabio.cumbo@gmail.com",
                 url="http://github.com/BlankenbergLab/meta-index",
                 license="LICENSE",
                 packages=setuptools.find_packages(),
                 package_data={
                    "meta-index": [
                        "requirements.txt",
                        "LICENSE",
                        "modules/*",
                 ]},
                 entry_points={
                    "console_scripts": ["meta-index = meta-index:main"]
                 },
                 description=("A pipeline for automatically indexing microbial genomes and accurately characterizing "
                              "metagenome-assembled genomes with sequence bloom trees"),
                 long_description=open("README.md").read(),
                 long_description_content_type="text/markdown",
                 install_requires=[
                    "numpy==1.16.3"
                 ],
                 python_requires='>=3.8',
                 zip_safe=False,
                 classifiers=[
                    "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: MIT License"
                 ])