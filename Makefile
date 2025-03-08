.PHONY: clean deploy format install mount sdist test uninstall upload

# Remove unnecessary data
clean:
	rm -rf build dist MetaSBT.egg-info
	find metasbt -type f -iname "*.pyc" -delete
	find metasbt -type d -iname "__pycache__" -delete

# Shortcut for building and uploading package
deploy: sdist upload clean

# Format code with black
format:
	black --line-length 120 .

# Install requirements and software
install: requirements.txt dist
	pip install -r requirements.txt
	find dist -type f -iname "MetaSBT-*.tar.gz" -exec pip install {} \;

# Run linting with tox
lint: tox.ini
	tox -e lint

# Shortcut for building and installing package
mount: sdist install clean

# Create the new distribution
sdist: setup.py
	python setup.py sdist

# Run all unit tests
test:
	python metasbt/metasbt.py test all

# Uninstall package
uninstall:
	pip uninstall metasbt

# Upload the new distribution to the Python Package Index
upload: sdist
	twine upload dist/*