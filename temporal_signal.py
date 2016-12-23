#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from treetime import TreeTime
from treetime.utils import DateConversion
from Bio import Phylo, AlignIO

if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description="Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. "
                        "It will optionally reroot the tree to maximize the clock-like signal and recalculate branch length.")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")
    parser.add_argument('--dates', required = True, type = str,
                        help ="csv with dates for nodes with 'node_name, date' where date is float (as in 2012.15)")
    parser.add_argument('--aln', required = False, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--infer_gtr', default = False, action='store_true', help='infer substitution model')
    parser.add_argument('--reroot', required = False, action="store_true", default=False,
                        help ="reroot the tree to maximize the correlation of root-to-tip distance with sampling time")
    parser.add_argument('--plot', required = False, action="store_true", default=False,)
    parser.add_argument('--verbose', default = 0, type=int,
                        help='verbosity of output 0-6')
    params = parser.parse_args()

    ###########################################################################
    ### PARSING DATES
    ###########################################################################
    with open(params.dates) as date_file:
        dates = {}
        failed_dates = 0
        for line in date_file:
            try:
                name, date = line.strip().split(',')[:2]
                dates[name] = float(date)
            except:
                failed_dates+=1
                print("couldn't parse date from:",line.strip(),"\n\texpecting float in second column")
        if len(dates)<failed_dates:
            print("\n\nDATE PARSING FAILED, ABORTING...")
            import sys
            sys.exit(1)

    ###########################################################################
    ### FAKING ALIGMENT IF NONE GIVEN
    ###########################################################################
    if params.aln is None:
        from Bio import Seq, SeqRecord, Align
        aln = Align.MultipleSeqAlignment([SeqRecord.SeqRecord(Seq.Seq("AAA"), id=name, name=name)
                                    for node in dates])


    ###########################################################################
    ### ESTIMATE ROOT (if requested) AND DETERMINE TEMPORAL SIGNAL
    ###########################################################################
    base_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])
    myTree = TreeTime(dates=dates, tree=params.tree,
                      aln=params.aln, gtr='JC69', verbose=params.verbose)
    myTree.one_mutation=0.001

    if params.reroot:
        myTree.reroot('best')
        # write rerooted tree to file
        outtree_name = base_name+'_rerooted.newick'
        Phylo.write(myTree.tree, outtree_name, 'newick')


    d2d = DateConversion.from_tree(myTree.tree)
    print('\n',d2d)


    ###########################################################################
    ### PLOT AND SAVE RESULT
    ###########################################################################
    if params.plot:
        import matplotlib.pyplot as plt
        myTree.plot_root_to_tip(label=False)
        t = np.array([np.min(dates.values()), np.max(dates.values())])
        plt.plot(t, t*d2d.slope+d2d.intercept,
                 label='y=%1.4f+%1.5ft, r^2=%1.2f'%(d2d.intercept, d2d.slope, d2d.r_val**2))
        plt.legend(loc=2)
        plt.savefig(base_name+'_root_to_tip_regression.pdf')




