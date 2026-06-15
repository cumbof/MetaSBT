.PHONY: clean deploy install mount sdist test uninstall upload

# Remove unnecessary data
clean:
	rm -rf build dist MetaSBT.egg-info deltatree/target
	find metasbt -type f -iname "*.pyc" -delete
	find metasbt -type d -iname "__pycache__" -delete

# Shortcut for building and uploading package
deploy: sdist upload clean

# Install package from the source distribution
install: dist
	find dist -type f -iname "MetaSBT-*.tar.gz" -exec pip install {} \;

# Shortcut for building and installing package
mount: sdist install clean

# Create the new distribution
sdist: pyproject.toml
	python -m build --sdist

# Run all unit tests
test:
	python metasbt/metasbt.py test all

# Uninstall package
uninstall:
	pip uninstall -y metasbt

# Upload the new distribution to the Python Package Index
upload: sdist
	twine upload dist/*
