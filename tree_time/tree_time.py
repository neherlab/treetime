"""
Class, which contains methods to optimize branch lengths given the time
constraints set to leaves
"""
# import tree_anc as ta
from __future__ import print_function, division
from .tree_anc import TreeAnc
import numpy as np
from Bio import AlignIO
import datetime
from scipy import stats
import matplotlib.pyplot as plt


class DateConversion(object):

    """
    Small container class to store parameters to convert between branch length
    as it is used in ML computations and the dates of the nodes.
    It is assumed that the conversion formula is 'length = k*date + b'
    """

    def __init__(self):

        self.slope = 0
        self.intersect = 0
        self.r_val = 0
        self.pi_val = 0
        self.sigma = 0

    @classmethod
    def from_tree(cls, t):
        

        dates = []
        for node in self.tree.find_clades():
            if node.date is not None:
                dates.append((node.date, node.dist2root))
        dates = np.array(dates)
        cls.slope,\
            cls.intersect,\
            cls.r_val,\
            cls.pi_val,\
            cls.sigma = stats.linregress(dates[:, 0], dates[:, 1])
        return cls

        # set dates to the internal nodes

        self._ml_t_init(gtr)


    def get_branch_len(self, date1, date2):
        """
        Compute branch length given the dates of the two nodes.

        Args:
         - date1 (int): date of the first node (days before present)
         - date2 (int): date of the second node (days before present)

        Returns:
         - branch length (double): Branch length, assuming that the dependence
         between the node date and the node depth in the the tree is linear.
        """
        return abs(date1 - date2) * self.slope

    def get_date(self, abs_t):
        """
        Get the approximate date of the tree node, assuming that the
        dependence between the node date and the node depth int the tree is
        linear.

        Args:
         - node(Phylo.Tree.Clade): node of the tree. Must be from the TreeAnc
         class (or its derivative), to contain the necessary attributes (
            dist2root).

        """
        year = (self.intersect - abs_t) / self.slope
        if year < 0:
            print ("Warning: got the negative date! Returning the inverse.")
            year = abs(year)
        return year


class TreeTime(TreeAnc, object):

    """
    TreeTime is the main class to perform the optimization of the node
    positions  given the temporal constraints of (some) nodes and leaves.
    To se how to use it, please refer to the examples section.
    """

    def __init__(self, tree):
        super(TreeTime, self).__init__(tree)
        self.date2dist = None  # we do not know anything about the conversion
        self.tree_file = ""

    @classmethod
    def from_files(cls, tree_file, aln_file, **kwargs):
        """
        """
        # parsing kwargs:
        if 'tree_fmt' in kwargs:
            tree_fmt = kwargs['tree_fmt']
        else:
            print ("Trying to read tree using default tree format")
            tree_fmt = 'newick'
        if 'aln_fmt' in kwargs:
            aln_fmt = kwargs['aln_fmt']
        else:
            print ("Trying to read alignment using default format (fasta)")
            aln_fmt = 'fasta'

        tree = cls.from_file(tree_file, tree_fmt)
        tree.tree_file = tree_file

        aln = AlignIO.read(aln_file, aln_fmt)
        tree.set_seqs_to_leaves(aln)

        if 'dates_file' in kwargs:
            dates_file = kwargs['dates_file']
            d_dic = tree._read_dates_file(dates_file)
            tree._set_dates_to_all_nodes(d_dic)
        return tree

    def load_dates(self, dates_file):
        """
        Args:
         - dates_file(str): path to file containing timing info about the tree
         nodes. File should be in csv format: 'nodeName,year'
        """
        d_dic = self._read_dates_file(dates_file)
        if (len(d_dic)) < 1:
            return
        self._set_dates_to_all_nodes(d_dic)

    def _str_to_date(self, instr):
        """
        Convert input string to datetime object.

        Args:
         - instr (str): input string. Accepts one of the formats:
         {YYYY.MM.DD, YYYY.MM, YYYY}.

        Returns:
         - date (datetime.datetime): parsed date object. If the parsing
         failed, None is returned
        """
        # import ipdb; ipdb.set_trace()
        try:
            date = datetime.datetime.strptime(instr, "%Y.%m.%d")
        except ValueError:
            date = None
        if date is not None:
            return date

        try:
            date = datetime.datetime.strptime(instr, "%Y.%m")
        except ValueError:
            date = None

        if date is not None:
            return date

        try:
            date = datetime.datetime.strptime(instr, "%Y")
        except ValueError:
            date = None

        return date

    def _read_dates_file(self, inf, **kwargs):
        """
        Read dates from the file into python dictionary. The input file should
        be in csv format 'node name, date'. The date will be converted to the
        datetime object and added to the dictionary {node name: datetime}
        """

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = 10

        if verbose > 3:
            print ("Reading datetime info for the tree nodes...")
        with open(inf, 'r') as finf:
            all_ss = finf.readlines()
        if verbose > 5:
            print ("Loaded %d lines form dates file" % len(all_ss))
        try:
            dic = {s.split(',')[0]: self._str_to_date(s.split(',')[1].strip())
                   for s in all_ss if not s.startswith("#")}
            if verbose > 3:
                print ("Parsed data in %d lines of %d input, %d corrupted"
                       % (len(dic), len(all_ss), len(all_ss) - len(dic)))
            return dic
        except ValueError:
            # unable to read all dates, the file is corrupted - go one by one
            print ("Unable to perform parsing of the dates file, file is "
                   "corrupted. Return empty dictionary.")
            return {}

    def _set_dates_to_all_nodes(self, dates_dic):
        """
        Set the time information to all nodes.
        Gets the datetime object of the nodes specified, calculate the time
        before present  (in days) for the nodes and sets this parameter (as
            int) to the node.date attribute.
        Args:
         - dates_dic(dic): dictionary with datetime informationfor nodes.
         Format: {node name: datetime object}

        """
        now = datetime.datetime.now()

        for node in self.tree.find_clades(order='preorder'):
            if node.name in dates_dic \
                    and dates_dic[node.name] is not None:
                days_before_present = (now - dates_dic[node.name]).days
                if days_before_present < 0:
                    print ("Cannot set the date ")
                    continue
                node.date = days_before_present
            else:
                node.date = None
        self.reroot_to_oldest()

    def reroot_to_oldest(self):
        """
        Set the root to the oldest node
        """
        # for now, we just reroot the to the most ancient node. Later,
        # it should be done in a cleverer way.

        self.tree.root_with_outgroup(sorted(self.tree.get_terminals(),
                                            key=lambda x: x.date)[-1])
        self.tree.root.date = None
        # fix tree lengths, etc
        self._add_node_params()

    def init_date_constraints(self, gtr):
        """
        Get the conversion coefficients between the dates and the branch
        lengths as they are used in ML computations. The conversion formula is
        assumed to be 'length = k*date + b'. For convenience, these
        coefficients as well as regression parameters are stored in the
        dates2dist object.

        Note: that tree must have dates set to all nodes before calling this
        function. (The latter is accomplished by calling load_dates func).
        """

        self.date2dist = DateConversion()

        dates = []
        for node in self.tree.find_clades():
            if node.date is not None:
                dates.append((node.date, node.dist2root))
        dates = np.array(dates)
        self.date2dist.slope,\
            self.date2dist.intersect,\
            self.date2dist.r_val,\
            self.date2dist.pi_val,\
            self.date2dist.sigma = stats.linregress(dates[:, 0], dates[:, 1])

        # set dates to the internal nodes

        self._ml_t_init(gtr)

    def _ml_t_init(self, gtr):
        """
        Initialize the tree nodes for ML computations with temporal
        constraints.
        Set the absolute positions for the nodes, init grid and constraints,
        set sequence profiles to the nodes.

        Args:
         - gtr(GTR): Evolutionary model, required to compute some node
         parameters.
        """
        if self.date2dist is None:
            print ("error")
            return

        for node in self.tree.find_clades():
            node.ml_t_prefactor = 0.0
            self._set_rotated_profiles(node, gtr)
            if node.date is not None:
                node.abs_t = node.date * self.date2dist.slope + \
                    self.date2dist.intersect
                # constraints mean that the probability distribution is
                # delta-function, so we do not need big grid
                node.grid = np.array([node.abs_t])
                node.prob = np.array([1])
            else:
                node.abs_t = node.dist2root  # not corrected!
                # if there are no constraints - grid will be set on-the-fly
                node.grid = None
                node.prob = None

    def _ml_t_msgup(self, gtr):
        """
        Propagate up- messages for ML computations with temporal constraints.
        To each node, sets the grid and the likelihood distribution on the grid

        Args:
         - gtr(GTR): evolutionary model
        """

        for node in self.tree.find_clades(order='postorder'):  # down->up
            if node.is_terminal():
                continue

            # we already have processed the node
            if hasattr(
                    node,
                    "prob") and node.prob is not None:
                continue

            # children nodes with constraints
            clades = [k for k in node.clades if k.prob is not None]
            if len(clades) < 1:  # we need at least one constrainted
                continue

            max_t = np.min([_c.grid.max()
                            for _c in clades])  # maximal grid node
            min_t = max_t - 1 * np.max([c_.grid.max() - c_.grid.min()
                                        for c_ in clades] + [2 * c_.branch_length for c_ in clades])

            node.grid = np.linspace(min_t, max_t, 100)

            # find probabilites in the grid:

            node.prob = np.ones(node.grid.shape)

            # iterate all children
            for clade in clades:
                # grid shape can changfe when we go from internal nodes
                # to terminal and back. So, we create the array every time.
                grid = np.zeros((clade.grid.shape[0], node.grid.shape[0]))
                grid[:, :] -= node.grid
                grid[:, :] = (grid.T + clade.grid).T

                _prob = self._ml_t_grid_prob(
                    node.prf_r,
                    clade.prf_l,
                    grid,
                    gtr)
                # sum over the child grid to get conditional on terminal
                _prob[:, :] = (_prob.T * clade.prob).T

                node.prob *= _prob.sum(0)
                node.ml_t_prefactor += clade.ml_t_prefactor
            # import ipdb; ipdb.set_trace()
            pre = node.prob.sum()

            node.prob /= pre

            node.ml_t_prefactor += np.log(pre)  # and store log-prefactor

    def _ml_t_msgdwn(self, gtr):
        """
        Propagate down- messages for ML computations with temporal constraints.
        for each node, set the grid and the likelihood distribution of the
        position on the on the grid

        Args:
         - gtr(GTR): evolutionary model
        """
        for node in self.tree.find_clades(order='preorder'):  # up->down
            if node.up is None:  # root node
                node.abs_t = node.grid[node.prob.argmax()]

                node.dist2root = 0.0
                node.date = self.date2dist.intersect

            else:
                grid = (node.grid - node.up.abs_t)
                idx = (grid > 0)

                if idx.sum() == 0:
                    print ("Warning! all desired branch lengths are negative!"
                           " Ignoring the messages from upside node")
                    prob_corr = np.ones(node.prob.shape)

                else:
                    prob_corr = np.array(
                        [-1 * self._neg_prob(t_, node.up.profile, node.profile, gtr) for t_ in grid])

                    prob_corr[~idx] = 0

                node.prob *= prob_corr
                node.abs_t = node.grid[node.prob.argmax()]

                node.branch_length = node.abs_t - node.up.abs_t
                node.dist2root = node.up.dist2root + node.branch_length
                node.date = (
                    node.dist2root - self.date2dist.intersect) / self.date2dist.slope

    def _ml_t_grid_prob(self, p_parent, p_child, grid, gtr):
        """
        Compute probability for 2D grid of times

        Args:

         - p_parent(numpy.array): parent profile (left profile). Shape: axL (a
            - alphabet size, L - sequence length)

         - p_child(numpy.array): child profile (right profile). Shape: axL (a
            - alphabet size, L - sequence length)

         - grid(numpy.array): prepared grid of times (branch lengths) to
         compute the probabilites for double-gridded nodes. The grid must have
         shape of (Lc, Lp, L), where Lc is the length of child grid, Lp is the
         length of the parent grid and L is the length of the sequence.

         - gtr(GTR): model of evolution.

         - res(numpy.array, default None): array to store results of the
         probability computations in order to avoid the construction of big
         arrays. Must have the same dimensions as grid. If set to none, the
         array will be constructed.

         - out_prob(numpy.array, default None): output probability. If
         specified, should have the shape of (Lc, Lp). If None or shape
         mismatch, will be constructed.
        """
        out_prob = np.zeros(grid.shape[:2])
        seq_l = p_parent.shape[0]

        if grid.ndim == 2:
            # grid.shape = (len(child),len(parent), 1)
            grid = grid.reshape(grid.shape + (1,))

        # indexes to exclude overlapping grids:
        idx = (
            grid >= 0).reshape(
            grid.shape[:2])  # idx.shape=(len(child),len(parent))

        # temp result
        # shape=(child_grid,parent_grid,L)
        tmp_res = np.zeros(idx.shape + (seq_l,))

        for state in range(gtr.alphabet.shape[0]):
            egrid = np.tile(np.exp(gtr.eigenmat[state] * grid), (1, 1, seq_l))
            # ma_egrid = np.ma.masked_array(egrid, idx, fill_value=0.0)
            tmp_res += (p_parent[:, state] * p_child[:, state]) * egrid

        out_prob[idx] = tmp_res.prod(-1)[idx]  # multiply along sequence
        # out_prob[~idx] = 0.0

        return out_prob

    def ml_t(self, gtr):
        """
        Perform tree optimization with temporal constarints.
        """
        #  propagate messages up
        self._ml_t_msgup(gtr)

        #  propagate messages down - reconstruct node positions
        self._ml_t_msgdwn(gtr)

    def date2dist_plot(self):
        """
        Plot the dependence between the node depth in the tree and the given
        node date information.
        """
        dates = []
        for node in self.tree.find_clades():
            if node.date is not None:
                dates.append((node.date, node.dist2root))
        dates = np.array(dates)

        if self.date2dist is None:
            self.date2dist = DateConversion()
            self.date2dist.intersect,\
                self.date2dist.slope,\
                self.date2dist.r_val,\
                self.date2dist.pi_val,\
                self.date2dist.sigma = stats.linregress(dates[:, 0],
                        dates[:, 1])

        plt.plot(dates[:, 0], dates[:, 1], 'o', c='r',
                 fillstyle='none', markersize=9, label='Data')

        plt.plot([0, dates[:, 0].max()],
                 [self.date2dist.intersect, self.date2dist.intersect +
                  self.date2dist.slope * dates[:, 0].max()],
                 lw=2, c='r', label='Linear regression')

        plt.grid()
        plt.ylabel("Distance from root (node depth)")
        plt.xlabel("Node date (days before present)")
        plt.title(
            "Dependence between the node depth\nand the given date time"
            " constraints of the node.\n Tree file: %s" %
            self.tree_file)
        plt.legend()

    def _set_rotated_profiles(self, node, gtr):
        """
        Set sequence and its profiles in the eigenspace of the transition
        matrix.
        """
        node.prf_r = node.profile.dot(gtr.v)
        node.prf_l = (gtr.v_inv.dot(node.profile.T)).T
