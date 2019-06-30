
Analyzing homoplasies and recurrent mutations using TreeTime
------------------------------------------------------------

Homoplasies are common if the same mutations are selected independently in many lineages, for example resistance mutations under treatment or cell-culture adaptation.
Homoplasies are also observed if the sequences used to build a tree underwent recombination or if there was contamination during sample prep and sequencing.
For all these reasons, a quick and easy way of gathering statistics on homoplasies and spotting sites and sequences that are disproportionally affected by homoplasies is useful.

After joint reconstruction of ancestral states, it is straight-forward to gather statistics on how often a particular site mutates and how often it hits a particular state.
The ``treetime homoplasy`` does exactly that.
To perform a simple analysis of homoplasies in a set of sequences, call TreeTime as

.. code-block:: bash

   treetime homoplasy --aln data/zika/zika.fasta --tree data/zika/zika.nwk

This command will reconstruct ancestral sequences and count how many time a particular site changed along the tree and how often the exact same mutation was observed.

The basic output looks like this:

.. code-block::

   The TOTAL tree length is 6.617e-02 and 674 mutations were observed.
   Of these 674 mutations,
      - 544 occur 1 times
      - 54 occur 2 times
      - 6 occur 3 times
      - 1 occur 4 times

This block reports how many times specific mutations, e.g. ``G27A``\ , were observed.
The next output block summarizes how often sites in the sequence are hit by mutations.
This can be compared to a null expectation that mutations are uniformly and independently distributed along the sequence according to a Poisson distribution with the same number of mutations.

.. code-block::

   Of the 10807 positions in the genome,
      - 10175 were hit 0 times (expected 10153.59)
      - 480 were hit 1 times (expected 633.25)
      - 69 were hit 2 times (expected 19.75)
      - 12 were hit 3 times (expected 0.41)
      - 5 were hit 4 times (expected 0.01)

This null expectation is almost always going to be rejected since many sites will not tolerate mutations and the evolutionary rate varies considerably across sites.
Nonetheless, this can be a useful comparison.

Lastly, there is a block that highlights mutations that are most homoplasic.
By default, ``treetime homoplasy`` will list the first 10 mutations -- this behavior can be adjusted using the option ``-n``.

.. code-block::

   The ten most homoplasic mutations are:
     mut multiplicity
     G9344A  4
     T428A 3
     T2060C  3
     A2925G  3
     G4319A  3
     T8207C  3
     C9287T  3
     A9C 2
     C11G  2
     A19C  2

Additional output
^^^^^^^^^^^^^^^^^

``treetime homoplasy`` can be run with the flag ``--detailed`` which will print additional statistics, mostly concerning terminal branches.
The first extra outputs are on the mutation statistics on terminal branch in analogy to the statistics for the total tree above.

.. code-block::

   The TERMINAL branch length is 3.962e-02 and 408 mutations were observed.
   Of these 408 mutations,
      - 353 occur 1 times
      - 23 occur 2 times
      - 3 occur 3 times

   The ten most homoplasic mutation on terminal branches are:
     mut multiplicity
     A2925G  3
     G4319A  3
     G9344A  3
     T73G  2
     ...

Lastly, the script outputs a list of sequences that have many mutations on the terminal branches leading up to them that also occur elsewhere in the tree.
This is can often be a sign of contamination or other problematic sequence.

.. code-block::

   Taxons that carry positions that mutated elsewhere in the tree:
     taxon name  #of homoplasic mutations
     Natal_RGN|KU527068|2015-07-01|brazil  10
     Haiti/1/2016|KX051563|2016-02-05|haiti  8
     PuertoRico/ZF8/2016|PuertoRico/ZF8/2016|2016-06-13|puerto_rico  8
     Brazil/PE243/2015|KX197192|2015-07-01|brazil  7
     ....

VCF input and reduced alignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When working with larger genomes, sequences are often represented as variants relative to a reference and stored as a vcf file or only variable positions are collated in a fasta file.
Trees build from such 'variable position only alignments' or SNP matrices have branch length that don't reflect real divergence and need to be rescaled.

If you have a vcf sequence input and a tree build on informative sites only (which will result in much longer branch length), you need to provide a rescaling factor unless branch length in your tree are already rescaled.
In our case, this would be 271/4411532=0.0000614.

.. code-block:: bash

   treetime homoplasy --aln data/tb/lee_2015.vcf.gz --vcf-reference data/tb/tb_ref.fasta --tree data/tb/lee_2015.nwk --rescale 0.0000614

Similarly, if you have an alignment with informative sites only, you should supplement the number of constant sites not included

.. code-block:: bash

   treetime homoplasy --aln data/tb/lee_2015.informative_sites.fasta --tree data/tb/lee_2015.nwk --rescale 0.0000614 --const 4411361

