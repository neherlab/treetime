from __future__ import print_function, division

import numpy as np
from tree_time.tree_time import tree_anc as ta
from tree_time.tree_time import tree_time as tt
import datetime
import os
from Bio import Phylo
import matplotlib.pyplot as plt
plt.ion()


def date_from_seq_name(name):
    date = None
    try:
        date_str = name.split('|')[-1]
        date = datetime.datetime.strptime(date_str, "%Y-%m-%d")
    except:        
        pass
    return date

def binstr(x):
    if x > 100:
        return '>10 or <0'
    else:
        return '<%.2f' % x


if __name__=='__main__':

    
    gtr = ta.GTR.standard()

    root_dir = os.path.dirname(os.path.realpath(__file__)) 
    fasta = os.path.join(root_dir, 'ebola.fasta')
    nwk = os.path.join(root_dir, 'ebola.nwk')

    t = tt.TreeTime.from_files(nwk, fasta)

    #t.reconstruct_anc('fitch')

    
    t.set_dates_to_nodes(date_from_seq_name)

    t.reroot_to_oldest()

    t.reconstruct_anc('ml', model=gtr)
    t.optimize_branch_len(gtr, verbose=10, store_old=False)
    t.prune_short_branches()
    t.reconstruct_anc('ml', model=gtr)
    t.optimize_branch_len(gtr, verbose=10, store_old=False)
    
    t.init_date_constraints(gtr)

    rds = []
    for n in t.tree.find_clades():
        if hasattr(n, 'raw_date') and n.raw_date is not None: 
            rds.append((n.dist2root, n.raw_date ))
    rds = np.array(rds)
    m = rds.max()
    plt.plot(rds[:,1],rds[:, 0], 'o')
    i = t.date2dist.intersect
    s = t.date2dist.slope
    plt.plot([0,m], [i, i + s*m ])


    t.ml_t(gtr)

    # plotting the results
    bins = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 1000.0])
    t._score_branches(bins)
    scores = []
    for n in t.tree.find_clades():
        scores.append(n.score)
    
    k = map(binstr, bins)
    
    plt.figure(1)
    plt.hist(scores)
    plt.xticks(np.linspace(0,0.9,10) + 0.05, list(k))
    plt.ylabel('Branch count')
    plt.xlabel('Branch length deviation from optimum,\nrelated to the optimal branch length')
    plt.title ('The deviation of the branch lengths from their optimal values\nEBOLA')

    #
    t.tree.ladderize()
    Phylo.draw(t.tree, label_func = lambda x:'', show_confidence=False, branch_labels='')
    plt.xlim(0.9999, 1.0025)
    

    plt.show()
    