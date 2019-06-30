Installation
============

TreeTime is compatible with Python 2.7 upwards and is tested on 2.7, 3.5, and 3.6.  It depends on several Python libraries:

* numpy, scipy, pandas: for all kind of mathematical operations as matrix
  operations, numerical integration, interpolation, minimization, etc.
* BioPython: for parsing multiple sequence alignments and phylogenetic trees
* matplotlib: optional dependency for plotting


Installing from PyPi or Conda
-----------------------------

You can also install TreeTime from PyPi via

.. code:: bash

  pip install phylo-treetime

You might need root privileges for system wide installation.
Similarly, you can install from conda using

.. code:: bash

  conda install -c bioconda treetime


Installing from source
----------------------
Clone or download the source code.

.. code:: bash

  git clone https://github.com/neherlab/treetime.git
  cd treetime
  pip install .

You might need root privileges for system wide installation. Alternatively, you can simply use it TreeTime locally without installation. In this case, just download and unpack it, and then add the TreeTime folder to your $PYTHONPATH.


Building the documentation
--------------------------

The API documentation for the TreeTime package is generated created with Sphinx. The source code for the documentaiton is located in doc folder.

  - sphinx-build to generate static html pages from source. Installed as

  .. code:: bash

	pip install Sphinx

  - basicstrap Html theme for sphinx:

  .. code:: bash

	pip install recommonmark sphinx-argparse

After required packages are installed, navigate to doc directory, and build the docs by typing:

 .. code:: bash

   make html

Instead of html, another target as `latex` or `epub` can be specified to build the docs in the desired format.

