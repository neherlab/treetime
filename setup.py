"""
Construct time-stamped phylogenies from a precomputed trees
Author: Pavel Sagulenko and Richard Neher
"""
import os
from setuptools import setup

setup(
        name = "treetime",
        version = "0.2.3",
        author = "Pavel Sagulenko and Richard Neher",
        author_email = "pavel.sagulenko@tuebingen.mpg.de",
        description = (""),
        license = "MIT",
        keywords = "Time-stamped phylogenies, evolution",
        url = "https://github.com/neherlab/treetime",
        packages=['treetime'],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Science",
            "License :: OSI Approved :: MIT License",
            ],
        )
