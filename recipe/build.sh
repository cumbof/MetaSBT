#!/bin/bash

mkdir -p $PREFIX/bin

# Conda package requires the modules folder in addition to the main meta-index controller, HELP, and LICENSE files
# Everything else will be removed
find . -mindepth 1 -maxdepth 1 -type d -not -name "modules" -not -name "requirements" -exec rm -rf {} \;
find . -mindepth 1 -maxdepth 1 -type f -not -name "meta-index" -not -name "LICENSE" -not -name "HELP" -exec rm -f {} \;

# Copy the main meta-index script and make it executable
cp meta-index $PREFIX/bin/meta-index
chmod +x $PREFIX/bin/meta-index

# Copy LICENSE and HELP assets
cp LICENSE $PREFIX/bin/LICENSE
cp HELP $PREFIX/bin/HELP

# Copy modules and requirements folders
cp -r modules $PREFIX/bin/modules
cp -r requirements $PREFIX/bin/requirements

# Make all modules executables
chmod -R +x $PREFIX/bin/modules
