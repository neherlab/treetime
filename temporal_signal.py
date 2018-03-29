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
                        "It will reroot the tree to maximize the clock-like "
                        "signal and recalculate branch length unless run with --keep_root.")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")
    parser.add_argument('--dates', required = True, type = str,
                        help ="csv with dates for nodes with 'node_name, date' where date is float (as in 2012.15)")
    parser.add_argument('--aln', required = False, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--infer_gtr', default = False, action='store_true', help='infer substitution model')
    parser.add_argument('--keep_root', required = False, action="store_true", default=False,
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

        if len(dates)<failed_dates:
            print("\n\nDATE PARSING FAILED, ABORTING...")
            import sys
            sys.exit(1)


    ###########################################################################
    ### FAKING ALIGMENT IF NONE GIVEN
    ###########################################################################
    if params.aln is None:
        from Bio import Seq, SeqRecord, Align
        aln = Align.MultipleSeqAlignment([SeqRecord.SeqRecord(Seq.Seq("AAA"), id=node, name=node)
                                    for node in dates])


    ###########################################################################
    ### ESTIMATE ROOT (if requested) AND DETERMINE TEMPORAL SIGNAL
    ###########################################################################
    base_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])
    myTree = TreeTime(dates=dates, tree=params.tree,
                      aln=aln, gtr='JC69', verbose=params.verbose)

    if not params.keep_root:
        myTree.reroot('best')

    d2d = DateConversion.from_tree(myTree.tree)
    print('\n',d2d)
    print('The R^2 value indicates the fraction of variation in'
          '\nroot-to-tip distance explained by the temporal sampling.'
          '\nHigher values corresponds more clock-like behavior (max 1.0).')

    print('\nThe rate is the slope of the best fit of the date to'
          '\nthe root-to-tip distance and provides an estimate of'
          '\nthe substitution rate. The rate needs to be positive!'
          '\nNegative rates suggest an inappropriate root.\n\n')

    if not params.keep_root:
        # write rerooted tree to file
        outtree_name = base_name+'_rerooted.newick'
        Phylo.write(myTree.tree, outtree_name, 'newick')
        print("--- re-rooted tree written to \n\t %s\n"%outtree_name)

    table_fname = base_name+'_rtt.csv'
    with open(table_fname, 'w') as ofile:
        ofile.write("#name, date, root-to-tip distance\n")
        ofile.write("#Dates of nodes that didn't have a specified date are inferred from the root-to-tip regression.\n")
        for n in myTree.tree.get_terminals():
            if hasattr(n, "numdate_given"):
                ofile.write("%s, %f, %f\n"%(n.name, n.numdate_given, n.dist2root))
            else:
                ofile.write("%s, %f, %f\n"%(n.name, d2d.numdate_from_dist2root(n.dist2root), n.dist2root))
        for n in myTree.tree.get_nonterminals(order='preorder'):
            ofile.write("%s, %f, %f\n"%(n.name, d2d.numdate_from_dist2root(n.dist2root), n.dist2root))
        print("--- wrote dates and root-to-tip distances to \n\t %s\n"%table_fname)


    ###########################################################################
    ### PLOT AND SAVE RESULT
    ###########################################################################
    if params.plot:
        import matplotlib.pyplot as plt
        myTree.plot_root_to_tip(label=False)
        t = np.array([np.min(dates.values()), np.max(dates.values())])
        plt.plot(t, t*d2d.clock_rate+d2d.intercept,
                 label='y=%1.4f+%1.5ft, r^2=%1.2f'%(d2d.intercept, d2d.clock_rate, d2d.r_val**2))
        plt.legend(loc=2)
        plt.savefig(base_name+'_root_to_tip_regression.pdf')
        print("--- root-to-tip plot saved to  \n\t %s_root_to_tip_regression.pdf\n"%base_name)




