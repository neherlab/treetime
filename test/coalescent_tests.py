import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib import cm

from treetime import TreeTime
from treetime.utils import parse_dates


def get_test_node(tt, pattern):
    return [n for n in tt.tree.find_clades() if pattern in n.name][0]

plt.ion()
ebola=True
if ebola:
    node_pattern = 'Gueckedou-C05'
    base_name = '../treetime_examples/data/ebola/ebola'
else:
    node_pattern = 'Indiana'
    base_name = '../treetime_examples/data/h3n2_na/h3n2_na_20'

dates = parse_dates(base_name+'.metadata.csv')
tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                aln = base_name+'.fasta', verbose = 1, dates = dates, precision=1, debug=True)


fixed_clock_rate = 0.0028
tt.infer_ancestral_sequences(infer_gtr=False, marginal=False)
tt.prune_short_branches()
tt.reroot(root='least-squares', clock_rate=fixed_clock_rate)
#tt.init_date_constraints(clock_rate=fixed_clock_rate)
tt.make_time_tree(clock_rate=fixed_clock_rate, time_marginal='assign')
tree_events_tt = np.array(sorted([(n.time_before_present, len(n.clades)-1)
                            for n in tt.tree.find_clades() if not n.bad_branch],
                            key=lambda x:-x[0]))
Phylo.draw(tt.tree, label_func=lambda x:"")

Tc=0.01
tt.add_coalescent_model(Tc)
tt.make_time_tree(clock_rate=fixed_clock_rate, time_marginal='assign')
tree_events_tt_post_coal = np.array(sorted([(n.time_before_present, len(n.clades)-1)
                            for n in tt.tree.find_clades() if not n.bad_branch],
                            key=lambda x:-x[0]))
Phylo.draw(tt.tree, label_func=lambda x:"")