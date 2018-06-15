"""
Construct time-stamped phylogenies from a precomputed trees
Author: Pavel Sagulenko and Richard Neher
"""
import os
from setuptools import setup

setup(
        name = "treetime",
        version = "0.4.0",
        author = "Pavel Sagulenko and Richard Neher",
        author_email = "richard.neher@unibas.ch",
        description = ("Maximum-likelihood phylodynamic inference"),
        license = "MIT",
        keywords = "Time-stamped phylogenies, evolution",
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
            "Topic :: Science",
            "License :: OSI Approved :: MIT License",
            ],
        scripts=['bin/mugration.py',
        		'bin/homoplasy_scanner.py',
        		'bin/timetree_inference.py',
        		'bin/temporal_signal.py',
        		'bin/ancestral_reconstruction.py'
        ]
        )
