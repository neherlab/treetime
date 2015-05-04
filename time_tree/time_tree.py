"""
Class, which contains methods to optimize branch lengths given the time
constraints set to leaves
"""
# import tree_anc as ta
from tree_anc import TreeAnc
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
        self.k = 0
        self.b = 0
        self.r = 0
        self.pi = 0
        self.sigma = 0

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
        return abs(date1 - date2) * self.k

    def get_date(self, node):
        """
        Get the approximate date of the tree node, assuming that the
        dependence between the node date and the node depth int the tree is
        linear.

        Args:
         - node(Phylo.Tree.Clade): node of the tree. Must be from the TreeAnc
         class (or its derivative), to contain the necessary attributes (
            dist2root).

        """
        yr = (self.b - node.dist2root) / self.k
        if (yr < 0):
            print ("Warning: got the negative date! Returning the inverse.")
            yr = abs(yr)
        return yr


class TreeTime(TreeAnc, object):

    def __init__(self, tree, dates_file=""):
        super(TreeTime, self).__init__(tree)
        self.date2dist = None  # we do not know anything about the conversion
        self.tree_file = ""

    @classmethod
    def from_files(cls, tree_file, aln_file, dates_file="", tree_fmt='newick', aln_fmt='fasta'):
        t = cls.from_file(tree_file, tree_fmt)
        t.tree_file = tree_file

        aln = AlignIO.read(aln_file, aln_fmt)
        t.set_seqs_to_leaves(aln)
        if dates_file != "":
            d_dic = t._read_dates_file(dates_file)
            t._set_dates_to_all_nodes(d_dic)  # even with empty dictionary we need to set dates as empty objects
        return t

    def load_dates(self, dates_file):
        """
        Args:
         - dates_file(str): path to file containing timing info about the tree nodes. File should be in csv format: 'nodeName,year'
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
         - date (datetime.datetime): parsed date object. If the parsing failed, None is returned
        """
        #import ipdb; ipdb.set_trace()
        try:
            date  = datetime.datetime.strptime(instr,  "%Y.%m.%d")
        except ValueError:
            date = None
        if date is not None:
            return date

        try:
            date  = datetime.datetime.strptime(instr,  "%Y.%m")
        except ValueError:
            date = None

        if date is not None:
            return date

        try:
            date  = datetime.datetime.strptime(instr,  "%Y")
        except ValueError:
            date = None

        return date

    def _read_dates_file(self, inf, **kwargs):
        """
        Read dates from the file into python dictionary. The input file should be in csv format 'node name, date'. The date will be converted to the datetime object and added to the dictionary {node name: datetime}
        """

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = 10

        if (verbose > 3):
            print ("Reading datetime info for the tree nodes...")
        with open(inf, 'r') as finf:
            ss = finf.readlines()
        if verbose > 5:
            print ("Loaded %d lines form dates file" % len(ss))
        try:
            d = {s.split(',')[0]: self._str_to_date(s.split(',')[1].strip())
                    for s in ss if not s.startswith("#")}
            if verbose > 3:
                print ("Parsed data in %d lines of %d input, %d corrupted"
                    % (len(d), len(ss), len(ss) - len(d)))
            return d
        except:
            # unable to read all dates, the file is corrupted - go one by one
            print ("Unable to perform parsing of the dates file, file is corrupted. Return empty dictionary.")
            return {}

    def _set_dates_to_all_nodes(self, dates_dic):
        """
        Set the time information to all nodes.
        Gets the datetime object of the nodes specified, calculate the time before present  (in days) for the nodes and sets this parameter (as int) to the node.date attribute.
        Args:
         - dates_dic(dic): dictionary with datetime informationfor nodes. Format: {node name: datetime object}

        """
        now = datetime.datetime.now()

        for node in self.tree.find_clades(order='preorder'):
            if dates_dic.has_key(node.name) \
                    and dates_dic[node.name] is not None:
                days_before_present = (now - dates_dic[node.name]).days
                if (days_before_present < 0):
                    print ("Cannot set the date ")
                    continue
                node.date = days_before_present
            else:
                node.date = None
        self.reroot_to_oldest()

    def reroot_to_oldest(self):
        # for now, we just reroot the to the most ancient node. Later,
        # it should be done in a cleverer way.

        self.tree.root_with_outgroup(sorted(self.tree.get_terminals(),
            key=lambda x: x.date)[-1])
        self.tree.root.date = None
        # fix tree lengths, etc
        self._add_node_params()

    def init_date_constraints(self, gtr):
        """
        Get the conversion coefficients between the dates and the branch lengths as they are used in ML computations. The conversion formula is assumed to be 'length = k*date + b'. For convenience, these coefficients as well as regression parameters are stored in the dates2dist object.

        Note: that tree must have dates set to all nodes before calling this function. (The latter is accomplished by calling load_dates func).
        """

        self.date2dist = DateConversion()

        d = []
        for n in self.tree.find_clades():
            if n.date is not None:
                d.append((n.date, n.dist2root))
        d = np.array(d)
        self.date2dist.k,\
            self.date2dist.b,\
            self.date2dist.r,\
            self.date2dist.pi,\
            self.date2dist.sigma = stats.linregress(d[:, 0], d[:, 1])

        # set dates to the internal nodes

        self._ml_t_init(gtr)

    def _ml_t_init(self, gtr):
        """
        Initialize the tree nodes for ML computations with temporal constraints.
        Set the absolute positions for the nodes, init grid and constraints, set sequence profiles to the nodes.

        Args:
         - gtr(GTR): Evolutionary model, required to compute some node parameters.
        """
        if self.date2dist is None:
            print ("error")
            return

        for n in self.tree.find_clades():
            n.ml_t_prefactor = 0.0
            self._set_rotated_profiles(n, gtr)
            if n.date is not None:
                n.abs_t = n.date * self.date2dist.k + self.date2dist.b
                # constraints mean that the probability distribution is delta-function, so we do not need big grid
                n.grid = np.array([n.abs_t])
                n.prob = np.array([1])
            else:
                n.abs_t = n.dist2root  # not corrected!
                # if there are no constraints - grid will be set on-the-fly
                n.grid = None
                n.prob = None

    def _ml_t_msgup(self, gtr):
        """
        Propagate up- messages for ML computations with temporal constraints.
        To each node, sets the grid and the likelihood distribution on the grid

        Args:
         - gtr(GTR): evolutionary model
        """
        i = 0
        for n in self.tree.find_clades(order='postorder'): # down->up
            if n.is_terminal():
                continue

            if hasattr(n, "prob") and  n.prob is not None:  # we already have processed the node
                continue

            # children nodes with constraints
            clades = [k for k in n.clades if k.prob is not None]
            if len(clades) < 1:  # we need at least one constrainted
                continue

            max_t = np.min([_c.grid.max() for _c in clades])  # maximal grid node
            min_t = max_t - 1 * np.max([c_.grid.max() - c_.grid.min()
                for c_ in clades] + [2 * c_.branch_length for c_ in clades])



            n.grid = np.linspace(min_t, max_t, 100)

            nn = n.grid.max() - n.grid.min()
            print ("relative grid width: %.4f" %(nn / n.abs_t))

            if nn > 1 :
                return n
            # find probabilites in the grid:

            n.prob = np.ones(n.grid.shape)

            # iterate all children
            for c in clades:
                # grid shape can changfe when we go from internal nodes
                # to terminal and back. So, we create the array every time.
                grid = np.zeros((c.grid.shape[0], n.grid.shape[0]))
                grid[:, :] -= n.grid
                grid [:, :] = (grid.T + c.grid).T

                _prob = self._ml_t_grid_prob(n.prf_r, c.prf_l, grid, gtr)
                # sum over the child grid to get conditional on terminal
                _prob[:, :] = (_prob.T * c.prob).T

                n.prob *= _prob.sum(0)
                n.ml_t_prefactor += c.ml_t_prefactor
            # import ipdb; ipdb.set_trace()
            pre = n.prob.sum()

            n.prob /= pre

            n.ml_t_prefactor += np.log(pre) # and store log-prefactor

    def _ml_t_msgdwn(self, gtr):
        """
        Propagate down- messages for ML computations with temporal constraints.
        for each node, set the grid and the likelihood distribution of the
        position on the on the grid

        Args:
         - gtr(GTR): evolutionary model
        """
        for n in self.tree.find_clades(order='preorder'):  # up->down
            if n.up is None:  # root node
                n.abs_t = n.grid[n.prob.argmax()]

                n.dist2root = 0.0
                n.date = self.date2dist.b

            else:
                grid = (n.grid - n.up.abs_t)
                idx = (grid > 0)

                if idx.sum() == 0:
                    print ("Warning! all desired branch lengths are negative!"
                    " Ignoring the messages from upside node")
                    prob_corr = np.ones(n.prob.shape)

                else:
                    prob_corr = np.array([-1 * self._neg_prob(t_, n.up.profile, n.profile, gtr) for t_ in grid])

                    prob_corr[~idx] = 0

                n.prob *= prob_corr
                n.abs_t = n.grid[n.prob.argmax()]

                n.branch_length = n.abs_t - n.up.abs_t
                n.dist2root = n.up.dist2root + n.branch_length
                n.date = (n.dist2root - self.date2dist.b) / self.date2dist.k

    def _ml_t_grid_prob(self, p_parent, p_child, grid, gtr,
            res=None, out_prob=None):
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

        # FIXME: with masked arrays should be more efficient (?)

        if out_prob is None or out_prob.shape != grid.shape[:2]:
            out_prob = np.zeros(grid.shape[:2])
        else:
            out_prob[:, :] = 0.0

        L = p_parent.shape[0]

        if grid.ndim == 2:
            # grid.shape = (len(child),len(parent), 1)
            grid = grid.reshape(grid.shape + (1,))

        # indexes to exclude overlapping grids:
        idx = (grid >= 0).reshape(grid.shape[:2]) # idx.shape=(len(child),len(parent))

        # temp result
        tmp_res = np.zeros(idx.shape+(L,)) # shape=(child_grid,parent_grid,L)

        for state in xrange(gtr.alphabet.shape[0]):
            # FIXME  time-limiting operation is to copy-tile all values
            egrid = np.tile(np.exp(gtr.eigenmat[state] * grid), (1, 1, L))
            # ma_egrid = np.ma.masked_array(egrid, idx, fill_value=0.0)
            tmp_res += (p_parent[:, state] * p_child[:, state]) * egrid

        out_prob[idx] = tmp_res.prod(-1)[idx]  # multiply along sequence
        # out_prob[~idx] = 0.0

        return out_prob

    def ml_t(self, gtr):
        #  propagate messages up
        self._ml_t_msgup(gtr)

        #  propagate messages down - reconstruct node positions
        self._ml_t_msgdwn(gtr)


    def date2dist_plot(self):
        """
        Plot the dependence between the node depth in the tree and the given node date information.
        """
        d = []
        for n in self.tree.find_clades():
            if n.date is not None:
                d.append((n.date, n.dist2root ))
        d = np.array(d)

        if self.date2dist is None:
            self.date2dist = DateConversion()
            self.date2dist.k,\
                self.date2dist.b,\
                self.date2dist.r,\
                self.date2dist.pi,\
                self.date2dist.sigma = stats.linregress(d[:, 0], d[:, 1])

        plt.plot(d[:, 0], d[:, 1], 'o', c='r',
                 fillstyle='none', markersize=9, label='Data')

        plt.plot([0, d[:, 0].max()],
                 [self.date2dist.b, self.date2dist.b + self.date2dist.k *
                 d[:, 0].max()], lw=2, c='r', label='Linear regression')

        plt.grid()
        plt.ylabel("Distance from root (node depth)")
        plt.xlabel("Node date (days before present)")
        plt.title ("Dependence between the node depth\nand the given date time constraints of the node.\n Tree file: %s" % self.tree_file)
        plt.legend()

    def _set_rotated_profiles(self, node, gtr):
        """
        Set sequence and its profiles in the eigenspace of the transition
        matrix.
        """
        node.prf_r = node.profile.dot (gtr.v)
        node.prf_l = (gtr.v_inv.dot(node.profile.T)).T