
Ancestral sequence reconstruction using TreeTime
------------------------------------------------

At the core of TreeTime is a class that models how sequences change along the tree.
This class allows to reconstruct likely sequences of internal nodes of the tree.
On the command-line, ancestral reconstruction can be done via the command

.. code-block:: bash

   treetime ancestral --aln data/h3n2_na/h3n2_na_20.fasta --tree data/h3n2_na/h3n2_na_20.nwk --outdir ancestral_results

This command will save a number of files into the directory `ancestral_results` and generate the output

.. code-block::

   Inferred GTR model:
   Substitution rate (mu): 1.0

   Equilibrium frequencies (pi_i):
     A: 0.2983
     C: 0.1986
     G: 0.2353
     T: 0.2579
     -: 0.01

   Symmetrized rates from j->i (W_ij):
     A C G T -
     A 0 0.8273  2.8038  0.4525  1.031
     C 0.8273  0 0.5688  2.8435  1.0561
     G 2.8038  0.5688  0 0.6088  1.0462
     T 0.4525  2.8435  0.6088  0 1.0418
     - 1.031 1.0561  1.0462  1.0418  0

   Actual rates from j->i (Q_ij):
     A C G T -
     A 0 0.2468  0.8363  0.135 0.3075
     C 0.1643  0 0.1129  0.5646  0.2097
     G 0.6597  0.1338  0 0.1432  0.2462
     T 0.1167  0.7332  0.157 0 0.2686
     - 0.0103  0.0106  0.0105  0.0104  0

   --- alignment including ancestral nodes saved as
      ancestral_results/ancestral_sequences.fasta

   --- tree saved in nexus format as
      ancestral_results/annotated_tree.nexus

TreeTime has inferred a GTR model and used it to reconstruct the most likely ancestral sequences.
The reconstructed sequences will be written to a file ending in ``_ancestral.fasta`` and a tree with mutations mapped to branches will be saved in nexus format in a file ending on ``_mutation.nexus``.
Mutations are added as comments to the nexus file like ``[&mutations="G27A,A58G,A745G,G787A,C1155T,G1247A,G1272A"]``.

Amino acid sequences
^^^^^^^^^^^^^^^^^^^^

Ancestral reconstruction of amino acid sequences works analogously to nucleotide sequences.
However, the user has to either explicitly choose an amino acid substitution model (JTT92)

.. code-block:: bash

   treetime ancestral --tree data/h3n2_na/h3n2_na_20.nwk  --aln data/h3n2_na/h3n2_na_20_aa.fasta --gtr JTT92

or specify that this is a protein sequence alignment using the flag ``--aa``\ :

.. code-block:: bash

   treetime ancestral --tree data/h3n2_na/h3n2_na_20.nwk  --aln data/h3n2_na/h3n2_na_20_aa.fasta --aa

VCF files as input
^^^^^^^^^^^^^^^^^^

In addition to standard fasta files, TreeTime can ingest sequence data in form of vcf files which is common for bacterial data sets where short reads are mapped against a reference and only variable sites are reported.
In this case, an additional argument specifying the mapping reference is required.

.. code-block:: bash

   treetime ancestral --aln data/tb/lee_2015.vcf.gz --vcf-reference data/tb/tb_ref.fasta --tree data/tb/lee_2015.nwk

The ancestral reconstruction is saved as a vcf files with the name ``ancestral_sequences.vcf``.

