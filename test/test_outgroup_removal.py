from treetime import TreeTime
from treetime.utils import parse_dates
from Bio import Phylo
import os
from io import StringIO



def test_outgroup_rerooting():
    myTree = TreeTime(gtr='Jukes-Cantor', tree = nwk, use_fft=False,
                    aln = None, verbose = 1, dates = dates, precision=3, debug=True)
    myTree.reroot(root=outgroup, clock_rate=1)
    assert [c.name for c in myTree.tree.root.clades[0].clades] == outgroup

def test_outgroup_removal():
    myTree = TreeTime(gtr='Jukes-Cantor', tree = nwk,
                    aln = fasta, verbose = 1, dates = dates, precision=3, debug=True)
    pre_length = len([n for n in myTree.tree.get_terminals()])
    myTree.run(infer_gtr='JC69', root=outgroup, remove_outgroup=False, **seq_kwargs, **tt_kwargs)
    Phylo.write(myTree.tree, "rerooted.nwk", 'newick')
    myTree = TreeTime(gtr='Jukes-Cantor', tree = nwk,
                    aln = fasta, verbose = 1, dates = dates, precision=3, debug=True)
    pre_length = len([n for n in myTree.tree.get_terminals()])
    myTree.run(infer_gtr='JC69', root=outgroup, remove_outgroup=True, **seq_kwargs, **tt_kwargs)
    assert [n for n in myTree.tree.find_clades() if n.name in outgroup] == []
    assert len([n for n in myTree.tree.get_terminals()]) == pre_length-2
    Phylo.write(myTree.tree, "outgroup_removed.nwk", 'newick')


outgroup = ['A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409', 'A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409']
root_dir = os.path.dirname(os.path.realpath(__file__))
nwk = root_dir + "/treetime_examples/data/h3n2_na/h3n2_na_20.nwk"
dates = parse_dates(root_dir + "/treetime_examples/data/h3n2_na/h3n2_na_20.metadata.csv")  
fasta = root_dir + "/treetime_examples/data/h3n2_na/h3n2_na_20.fasta"

seq_kwargs = {"marginal_sequences":True,
                    "branch_length_mode": 'input',
                    "reconstruct_tip_states":False}
tt_kwargs = {'clock_rate': 0.0028,
                    'time_marginal':'assign'}
                    
test_outgroup_rerooting()
test_outgroup_removal()

