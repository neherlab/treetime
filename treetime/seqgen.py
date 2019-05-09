from __future__ import division, print_function, absolute_import
from collections import defaultdict
import numpy as np
from treetime import config as ttconf
from .seq_utils import alphabets, profile_maps, alphabet_synonyms, seq2array, seq2prof
from .gtr import GTR
from .treeanc import TreeAnc


class SeqGen(TreeAnc):
    '''
    Evolve sequences along a given tree with a specific GTR model.
    This class inherits from TreeAnc.
    '''

    def __init__(self, *args, **kwargs):
        """Instantiate. Mandatory arguments are a tree and GTR model.
        """
        super(SeqGen, self).__init__(reduce_alignment=False, **kwargs)


    def sample_from_profile(self, p):
        """returns a sequence sampled from a profile (column wise state probabilities)

        Parameters
        ----------
        p : np.array
            sequence profile with dimensions (L,q)

        Returns
        -------
        np.array (character)
            sequence as character array array(['A', 'C', 'G',...])
        """
        cum_p = p.cumsum(axis=1).T

        prand = np.random.random(self.seq_len)
        seq = self.gtr.alphabet[np.argmax(cum_p>prand, axis=0)]
        return seq


    def evolve(self, root_seq=None):
        """Evolve a root sequences along a tree. If no root sequences
        is provided, one will be sampled from the equilibrium
        probabilities of the GTR model

        Parameters
        ----------
        root_seq : numpy character array, optional
            sequence to be used as the root sequence of the tree. if not given,
            will sample a sequence from the equilibrium probabilities of the GTR model.
        """
        self.seq_len = self.gtr.seq_len
        # set root if not given
        if root_seq:
            self.tree.root.sequence = seq2array(root_seq)
        else:
            if len(self.gtr.Pi.shape)==2:
                self.tree.root.sequence = self.sample_from_profile(self.gtr.Pi.T)
            else:
                self.tree.root.sequence = self.sample_from_profile(np.repeat([self.gtr.Pi], self.seq_len, axis=0))

        # generate sequences in preorder
        for n in self.tree.get_nonterminals(order='preorder'):
            profile_p = seq2prof(n.sequence, self.gtr.profile_map)
            for c in n:
                profile = self.gtr.evolve(profile_p, c.branch_length)
                c.sequence = self.sample_from_profile(profile)
        self.make_reduced_alignment()

        # gather mutations
        for n in self.tree.find_clades():
            if n==self.tree.root:
                n.mutations=[]
            else:
                n.mutations = self.get_mutations(n)


    def get_aln(self, internal=False):
        """assemble a multiple sequence alignment from the evolved
        sequences. Optionally in clude internal sequences

        Parameters
        ----------
        internal : bool, optional
            include sequences of internal nodes in the alignment

        Returns
        -------
        Bio.Align.MultipleSeqAlignment
            multiple sequence alignment
        """
        from Bio import SeqRecord, Seq
        from Bio.Align import MultipleSeqAlignment

        tmp = []
        for n in self.tree.get_terminals():
            if n.is_terminal() or internal:
                tmp.append(SeqRecord.SeqRecord(id=n.name, name=n.name, description='', seq=Seq.Seq(''.join(n.sequence))))

        return MultipleSeqAlignment(tmp)


