import os
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name = "phylo-treetime",
        version = "0.4.1",
        author = "Pavel Sagulenko, Emma Hodcroft, and Richard Neher",
        author_email = "richard.neher@unibas.ch",
        description = ("Maximum-likelihood phylodynamic inference"),
        long_description = long_description,
        long_description_content_type="text/markdown",
        license = "MIT",
        keywords = "Time-stamped phylogenies, phylogeography, virus evolution",
        url = "https://github.com/neherlab/treetime",
        packages=['treetime'],
        install_requires = [
            'biopython >=1.66',
            'matplotlib >=2.0',
            'numpy >=1.10.4',
            'pandas >=0.17.1',
            'scipy >=0.16.1'
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6"
            ],
        scripts=['bin/mugration.py',
        		'bin/homoplasy_scanner.py',
        		'bin/timetree_inference.py',
        		'bin/temporal_signal.py',
        		'bin/ancestral_reconstruction.py'
        ]
    )
