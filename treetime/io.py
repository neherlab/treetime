from Bio import Phylo
import numpy as np
import utils


def plot_vs_years(tt, years = 1, ax=None, confidence=None, ticks=True, **kwargs):
    '''
    converts branch length to years and plots the time tree on a time axis.
    Args:
        tt:     treetime object after a time tree is inferred
        years:  width of shaded boxes indicating blocks of years, default 1
        ax:     axis object. will create new axis of none specified
        confidence:     draw confidence intervals. This assumes that marginal
                        time tree inference was run
        **kwargs:   arbitrary kew word arguments that are passed down to Phylo.draw
    '''
    import matplotlib.pyplot as plt
    tt.branch_length_to_years()
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)
    # draw tree
    if "label_func" not in kwargs:
        nleafs = tt.tree.count_terminals()
        kwargs["label_func"] = lambda x:x.name if (x.is_terminal() and nleafs<30) else ""
    Phylo.draw(tt.tree, axes=ax, **kwargs)

    # set axis labels
    offset = tt.tree.root.numdate - tt.tree.root.branch_length
    xticks = ax.get_xticks()
    dtick = xticks[1]-xticks[0]
    shift = offset - dtick*(offset//dtick)
    tick_vals = [x+offset-shift for x in xticks]
    ax.set_xticks(xticks - shift)
    ax.set_xticklabels(map(str, tick_vals))
    ax.set_xlabel('year')
    ax.set_ylabel('')
    ax.set_xlim((0,np.max([n.numdate for n in tt.tree.get_terminals()])+2-offset))

    # put shaded boxes to delineate years
    if years:
        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        if type(years) in [int, float]:
            dyear=years
        from matplotlib.patches import Rectangle
        for yi,year in enumerate(np.arange(tick_vals[0], tick_vals[-1],dyear)):
            pos = year - offset
            r = Rectangle((pos, ylim[1]-5),
                          dyear, ylim[0]-ylim[1]+10,
                          facecolor=[0.7+0.1*(1+yi%2)] * 3,
                          edgecolor=[1,1,1])
            ax.add_patch(r)
            if year in tick_vals and pos>xlim[0] and pos<xlim[1] and ticks:
                ax.text(pos,ylim[0]-0.04*(ylim[1]-ylim[0]),str(int(year)),
                        horizontalalignment='center')
        ax.set_axis_off()

    # add confidence intervals to the tree graph -- grey bars
    if confidence:
        utils.tree_layout(tt.tree)
        if not hasattr(tt.tree.root, "marginal_inverse_cdf"):
            print("marginal time tree reconstruction required for confidence intervals")
        elif len(confidence)==2:
            cfunc = tt.get_confidence_interval
        elif len(confidence)==1:
            cfunc = tt.get_max_posterior_region
        else:
            print("confidence needs to be either a float (for max posterior region) or a two numbers specifying lower and upper bounds")
            return

        for n in tt.tree.find_clades():
            pos = cfunc(n, confidence)
            ax.plot(pos-offset, np.ones(len(pos))*n.ypos, lw=3, c=(0.5,0.5,0.5))


def treetime_to_newick(tt, outf):
    Phylo.write(tt.tree, outf, 'newick')


if __name__=='__main__':
    pass
