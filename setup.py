"""
Construct time-stamped phylogenies from a precomputed trees
Author: Pavel Sagulenko and Richard Neher
"""
import os
from setuptools import setup

setup(
        name = "treetime",
        version = "0.2.4",
        author = "Pavel Sagulenko and Richard Neher",
        author_email = "richard.neher@unibas.ch",
        description = ("Maximum-likelihood phylodynamic inference"),
        license = "MIT",
        keywords = "Time-stamped phylogenies, evolution",
        url = "https://github.com/neherlab/treetime",
        packages=['treetime'],
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
