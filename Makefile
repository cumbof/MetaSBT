.PHONY: clean deploy install test uninstall

# Remove unnecessary data
clean:
	rm -rf build dist MetaSBT.egg-info deltatree/target
	find metasbt -type f -iname "*.pyc" -delete
	find metasbt -type d -iname "__pycache__" -delete

# Setup environment, build source distribution, install it, and clean up
install:
	pip install ".[dev]"
	python -m build --sdist
	find dist -type f -iname "MetaSBT-*.tar.gz" -exec pip install {} \;
	$(MAKE) clean

# Setup environment and run all unit tests using the test data
test:
	pip install ".[dev]"
	python metasbt/metasbt.py test all --references test/references.tsv --mags test/mags.txt

# Setup environment, build source distribution, upload to PyPI, and clean up
deploy:
	pip install ".[dev]"
	python -m build --sdist
	twine upload dist/*
	$(MAKE) clean

# Uninstall package
uninstall:
	pip uninstall -y metasbt
