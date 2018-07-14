"""
Construct time-stamped phylogenies from a precomputed trees
Author: Pavel Sagulenko and Richard Neher
"""
import os
from setuptools import setup

setup(
        name = "treetime",
        version = "0.5.0",
        author = "Pavel Sagulenko and Richard Neher",
        author_email = "richard.neher@unibas.ch",
        description = ("Maximum-likelihood phylodynamic inference"),
        license = "MIT",
        keywords = "Time-stamped phylogenies, phylogeography, virus evolution",
        url = "https://github.com/neherlab/treetime",
        packages=['treetime'],
        install_requires = [
            'biopython>=1.66',
            'matplotlib>=2.0',
            'numpy>=1.10.4',
            'pandas>=0.17.1',
            'scipy>=0.16.1'
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6"
            ],
        scripts=['bin/treetime']
    )
