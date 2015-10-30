"""
Optimizing branch lengths of a given tree given the time stamp constraints of 
(some) leaves
Author: Dr. R.Neher's lab
"""
import os
from setuptools import setup

setup(
        name = "tree_time",
        version = "0.1",
        author = "Lab of Dr. R. Neher",
        author_email = "pavel.sagulenko@tuebingen.mpg.de",
        description = (""),
        license = "MIT",
        keywords = "coalescent evolution",
        url = "https://github.com/neherlab/timetree",
        packages=['time_tree'],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Science",
            "License :: OSI Approved :: MIT License",
            ],
        )
