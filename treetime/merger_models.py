"""
methods to calculate merger models for a time tree
"""
from __future__ import print_function, division
import numpy as np
from Bio import AlignIO, Phylo
from scipy.interpolate import interp1d
import config as ttconf


class Coalescent(object):
    """docstring for Coalescent"""
    def __init__(self, tree, Tc=0.001):
        super(Coalescent, self).__init__()
        self.tree = tree
        self.Tc = Tc
        self.make_skyline()

    def make_skyline(self):
        # collect all node locations and difference in branch count at that point
        # this is sorted by negative time before present, i.e.the root is first
        tmp = np.array(sorted([(-n.time_before_present, len(n.clades)-1)
                                for n in self.tree.find_clades()], key=lambda x:x[0]))
        # calculate the branch count at each point summing the delta branch counts
        tmp_t, tmp_n = tmp[:,0], np.cumsum(tmp[:,1])+1
        # evaluate the number of branches at unique temporal positions
        tvals = np.unique(tmp_t)
        nbranches = tmp_n[np.searchsorted(tmp_t, tvals)]
        # integrate the piecewise constant branch count function.
        cost = np.concatenate(([0],np.cumsum(np.diff(tvals)*nbranches[1:])))
        # make interpolation objects for the branch count and its integral
        # the latter is scales by 0.5/Tc
        self.nbranches = interp1d(-tvals, nbranches, kind='linear')
        # need to add extra point at very large time before present to prevent 'out of interpolation range' errors
        self.cost_func = interp1d(np.concatenate(([-ttconf.BIG_NUMBER], -tvals,[ttconf.BIG_NUMBER])),
                                  np.concatenate(([cost[-1]], cost,[cost[0]]))*0.5/self.Tc, kind='linear')

        # calculate merger rates
        mergers = np.array(sorted([(-n.time_before_present, len(n.clades)-1)
                                for n in self.tree.get_nonterminals()], key=lambda x:x[0]))

        events_t = -mergers[:,0]
        nlin = self.nbranches(events_t)
        events = 2.0*mergers[:,1]/(nlin*(nlin-1))
        self.normalized_mergers = np.array((events_t, events))

        # smooth this with a Gaussian kernel and add 0.01 of average to fill long gaps
        dt = 0.05*(events_t[0]-events_t[-1])
        windows = np.linspace(events_t[-1], events_t[0]-dt, 100)
        avg = np.sum(events)/np.abs(events_t[0]-events_t[-1])
        smoothing_kernel = lambda x: np.exp(-x**2/2.0/dt**2)/np.sqrt(2.0*np.pi)/dt
        self.Tc_inv = interp1d(windows,
                        [0.01*avg+0.99*np.sum(smoothing_kernel(events_t-w)*events) for w in windows])


    def cost(self, t_node, branch_length):
        # return the cost associated with a branch starting at t_node
        # t_node is time before present, the branch goes back in time
        return self.cost_func(t_node) - self.cost_func(t_node+branch_length)

    def attach_to_tree(self):
        for clade in self.tree.find_clades():
            if clade.up is not None:
                clade.branch_length_interpolator.merger_cost = self.cost


    def skyline(self, gen=1.0, to_numdate=None):
        '''
        return the skyline, i.e., an estimate of the inverse rate of coalesence
        parameters:
            gen -- number of generations per unit of time. Unit of time is branch length,
                   hence this needs to be the inverse substitution rate per generation
            to_numdate -- function to convert time before present to numerical dates
        '''
        if to_numdate is None:
            to_numdate =lambda x:x

        return interp1d(to_numdate(self.Tc_inv.x), gen/self.Tc_inv.y)



def traveling_wave(tree, Tc=None, tau=None):
    '''
    assigns coalescent merger rates to all branches in the tree
    '''
    #
    if tau is None:
        tau = Tc/8.0   # 8 is roughly the factor between total coalescence and expansion rate in realistic populations
    # determine the msg to parents
    for n in tree.find_clades(order='postorder'):
        n._polarizer_to_parent = np.sum([c._polarizer_to_parent for c in n.clades])
        n._polarizer_to_parent*= np.exp(-n.branch_length/tau)
        n._polarizer_to_parent+=(1-np.exp(-n.branch_length/tau))*tau

    # determine the msg to children
    tree.root._polarizer_from_parent = 0.0
    for n in tree.get_nonterminals(order='preorder'):
        tmp = np.sum([c._polarizer_to_parent for c in n.clades]) + n._polarizer_from_parent
        for c in n.clades:
            c._polarizer_from_parent = tmp-c._polarizer_to_parent
            c._polarizer_from_parent*= np.exp(-c.branch_length/tau)
            c._polarizer_from_parent+=(1-np.exp(-c.branch_length/tau))*tau

    # determine the msg to parents
    for n in tree.find_clades(order='postorder'):
        n.lbi = n._polarizer_from_parent + np.sum([c._polarizer_to_parent for c in n.clades])

    # assign those rates to all nodes in the tree
    for n in tree.find_clades():
        n.merger_rate = n.lbi/tau/Tc

