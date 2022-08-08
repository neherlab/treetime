import os
from setuptools import setup

def get_version():
    v = "0.0.0"
    with open('treetime/__init__.py') as ifile:
        for line in ifile:
            if line[:7]=='version':
                v = line.split('=')[-1].strip()[1:-1]
                break
    return v

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name = "phylo-treetime",
        version = get_version(),
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
            'biopython>=1.67,!=1.77,!=1.78',
            'numpy>=1.10.4',
            'pandas>=0.17.1',
            'scipy>=0.16.1'
        ],
        extras_require = {
            ':python_version < "3.6"':['matplotlib>=2.0, ==2.*'],
            ':python_version >= "3.6"':['matplotlib>=2.0'],
        },
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            ],
        entry_points = {
            "console_scripts": [
                "treetime = treetime.__main__:main",
            ]
        }
    )

