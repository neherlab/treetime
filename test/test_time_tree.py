"""
Test time tree class
"""

import unittest,os,sys
from Bio import AlignIO

from ..time_tree import tree_anc

resources_dir = os.path.join(path.dirname(__file__), '../data/')

class TestTreeAnc(unittest.TestCase):
    def test_read_newick(self):
        anc_t = tree_anc.TreeAnc.from_file(resources_dir+'PR.B.100.nwk', 'newick')
        assert len(anc_t.tree.get_terminals()) == 100 #  tree read successfully

    def test_set_links_to_parent(self):
        #  up-link installed succesfully
        anc_t = tree_anc.TreeAnc.from_file(resources_dir+'PR.B.100.nwk', 'newick')
        for node in anc_t.tree.get_nonterminals():
            for ch in node.clades:
                assert ch.up == node

    def test_aln_to_leaves(self):
        anc_t = tree_anc.TreeAnc.from_file(resources_dir+'PR.B.100.nwk', 'newick')
        aln = AlignIO.read(resources_dir+'PR.B.100.fasta', 'fasta')

        err = anc_t.set_seqs_to_leaves(aln)
        assert (err==0) # all sequencs were set up successfully




suite = unittest.TestLoader().loadTestsFromTestCase(TestTreeAnc)
unittest.TextTestRunner(verbosity=10).run(suite)

if __name__=='__main__':
    pass


