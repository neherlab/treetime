"""
Class, which contains methods to optimize branch lengths given the time
constraints set to leaves
"""
from __future__ import print_function, division

from treeanc import TreeAnc
import utils
import config as ttconf

import numpy as np
from Bio import AlignIO, Phylo
import datetime
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import copy
from scipy import optimize as sciopt

class TreeTime(TreeAnc, object):
    """
    TreeTime is the main class to perform the optimization of the node
    positions given the temporal constraints of (some) nodes and leaves.

    The optimization workflow includes the inferrence of the ancestral sequences
    using Fitch's method or maximum-likelihood (ML), followed by the unconstrained optimization
    of the branch lengths with maximum-likelihood method. After the optimization
    is done, the nodes with date-time information are arranged along the time axis,
    the appropriate conversion btween the branch lengths units and the date-time units
    is found. Then, for each internal node, we compute the the probability distribution
    of the node's location conditional on the fixed location of the leaves, which
    have temporal information. In the end, the most probable location of the internal nodes
    is converted to the time of the internal node.
    """

    def __init__(self, gtr):
        super(TreeTime, self).__init__(gtr)
        self.date2dist = None  # we do not know anything about the conversion
        self.tree_file = ""
        self.max_diam = 0.0
        self.debug=False

    @property
    def average_branch_len(self):
        """
        Compute the average branch length of the tree.
        Used to estimate the scale  of the branch-lengths
        """
        return np.mean([n.branch_length for n in self.tree.find_clades()])

    def reroot_to_oldest(self):
        """
        Set the root of the tree to the oldest node.
        """
        def numdate_given(node):
            if not hasattr(node, 'numdate_given') or node.numdate_given is None:
                return 0
            return node.numdate_given

        self.tree.root_with_outgroup(sorted(self.tree.get_terminals(), key=numdate_given)[0])
        self.tree.ladderize()
        og = self.tree.root.clades[0]
        self.tree.root.clades[1].branch_length += og.branch_length
        og.branch_length = self.one_mutation
        self.tree.root.branch_length = self.one_mutation
        self.tree.root.numdate_given = None
        # fix tree lengths, etc
        self.set_additional_tree_params()

    def init_date_constraints(self, slope=None):
        """
        Get the conversion coefficients between the dates and the branch
        lengths as they are used in ML computations. The conversion formula is
        assumed to be 'length = k*numdate_given + b'. For convenience, these
        coefficients as well as regression parameters are stored in the
        dates2dist object.

        Note: that tree must have dates set to all nodes before calling this
        function. (This is accomplished by calling load_dates func).
        """
        self.date2dist = utils.DateConversion.from_tree(self.tree, slope)
        self.max_diam = self.date2dist.intercept

        # set the None  for the date-related attributes in the internal nodes.
        # make interpolation objects for the branches
        self._ml_t_init()

    def _make_branch_len_interpolator(self, node, n=ttconf.BRANCH_GRID_SIZE):
        """
        Makes an interpolation object for probability of branch length. Builds
        an adaptive grid with n points, fine around the optimal branch length,
        and more sparse at the tails. This method does **not**  require the
        branch length to be at their optimal value, as it computes the optimal
        length itself. The method, however **requires** the knowledge of the
        sequences and/or sequence profiles for the node and its parent in order
        to compute distance probability in the scope of the GTR models

        Args:
         - node(Phylo.Clade): tree node. The probability distribution of the branch
           length from the node to its parent will be computed

         - n(int): number of points in the branch length grid.

        Returns:
         - None. The node gets new attribute - the linear interpolation object
           of the branch length probability distribution.

        """
        # no need to optimize the root branch length
        if node.up is None:
            node.branch_neg_log_prob = None
            return None
        if not hasattr(node, 'gamma'):
            node.gamma = 1.0

        parent = node.up
        prof_p = parent.profile
        prof_ch = node.profile

        # optimal branch length
        obl = self.gtr.optimal_t(node.up.profile, node.profile) # not rotated profiles!

        node.opt_branch_length = obl #  need for some computations

        if obl < np.min((1e-5, 0.1*self.one_mutation)): # zero-length

            grid = ttconf.MAX_BRANCH_LENGTH * (np.linspace(0, 1.0 , n/3)**2)

        else: # branch length is not zero

            sigma = obl #np.max([self.average_branch_len, obl])
            # from zero to optimal branch length
            grid_left = obl * (1 - np.linspace(1, 0.0, n/3)**2)
            # from optimal branch length to the right (--> 3*branch lengths),
            grid_right = obl + obl/100 + (3*sigma*(np.linspace(0, 1, n/3)**2))
            # far to the right (3*branch length ---> MAX_LEN), very sparse
            far_grid = grid_right.max() + obl/2 + ttconf.MAX_BRANCH_LENGTH*np.linspace(0, 1, n)**2

            grid = np.concatenate((grid_left,grid_right,far_grid))
            grid.sort() # just for safety

        grid = np.concatenate((
            [ttconf.MIN_T, -ttconf.TINY_NUMBER],
            grid,
            [ttconf.MAX_T]))

        # log-probability of the branch len to be at this value
        logprob = np.concatenate([
            [0., 0.],
            [self.gtr.prob_t(prof_p, prof_ch, node.gamma*t_, return_log=True) for t_ in grid[2:-2]],
            [0., 0.]])

        logprob[((0,1,-2,-1),)] = ttconf.MIN_LOG
        logprob *= -1.0

        dt = np.diff(grid)
        tmp_prob = np.exp(-logprob)
        integral = np.sum(0.5*(tmp_prob[1:]+tmp_prob[:-1])*dt)

        # save raw branch length interpolators without coalescent contribution
        # TODO: better criterion to catch bad branch
        if integral < 1e-200:
            print ("!!WARNING!! Node branch length probability distribution "
                "integral is ZERO. Setting bad_branch flag..."
                "Not accounting for the normalization, coalescence theory.")
            node.bad_branch = True
            node.raw_branch_neg_log_prob = interp1d(grid, logprob, kind='linear')

        else:
            node.raw_branch_neg_log_prob = interp1d(grid, logprob+np.log(integral), kind='linear')


        # add merger rate contribution to the raw branch length
        logprob += node.merger_rate * np.minimum(ttconf.MAX_BRANCH_LENGTH, np.maximum(0,grid))

        # normalize the branch lengths prob distribution
        min_prob = np.min(logprob)
        if np.exp(-1*min_prob) == 0:
            print ("!!Warning!! the branch length probability is zero. "
                   "Are the branches and sequences correct?")
            # setting bad branch flag
            node.bad_branch = True

        logprob -= min_prob
        dt = np.diff(grid)
        tmp_prob = np.exp(-logprob)
        integral = np.sum(0.5*(tmp_prob[1:]+tmp_prob[:-1])*dt)

        if integral < 1e-200:
            print ("!!WARNING!! Node branch length probability distribution "
                "integral is ZERO. Setting bad_branch flag..."
                "Not accounting for the normalization, coalescence theory.")
            node.bad_branch = True
            node.branch_neg_log_prob = interp1d(grid, logprob, kind='linear')

        else:
            node.branch_neg_log_prob = interp1d(grid, logprob+np.log(integral), kind='linear')


        # node gets new attribute
        return None

    def _update_branch_len_interpolators(self):
        """
        reassign interpolator object for branch length after changing the merger_rate
        """
        for node in self.tree.find_clades():
            if node.up is None:
                continue
            # make sure to copy raw_branch_neg_log_prob.y -> otherwise only referenced and modified
            grid,y = node.raw_branch_neg_log_prob.x, np.array(node.raw_branch_neg_log_prob.y)
            y+=node.merger_rate * np.minimum(ttconf.MAX_BRANCH_LENGTH, np.maximum(0,grid))
            dt = np.diff(grid)
            tmp_prob = np.exp(-y)
            integral = np.sum(0.5*(tmp_prob[1:]+tmp_prob[:-1])*dt)
            node.branch_neg_log_prob = interp1d(grid, y + np.log(integral), kind='linear')

    def _ml_t_init(self,ancestral_inference=True, **kwarks):
        """
        Initialize the attributes in all tree nodes that are required
        by the ML algorithm to compute the probablility distribution of the node
        locations. These attributes include the distance from the node postions
        to the present (in branch length units), branch length interpolation
        objects, and the probability distributions for the nodes which have the
        date-time information (these are going to be delta-functions), and
        set the sequence profiles in the eigenspace of the GTR matrix.

        """
        tree = self.tree
        if ancestral_inference:
            self.optimize_seq_and_branch_len(**kwarks)
        print('Initializing branch length interpolation object')
        if self.date2dist is None:
            print ("error - no date to dist conversion set. "
                "Run init_date_constraints and try once more.")
            return

        for node in tree.find_clades():

            if not hasattr(node, 'merger_rate'):
                node.merger_rate=ttconf.BRANCH_LEN_PENALTY

            # make interpolation object for branch lengths
            self._make_branch_len_interpolator(node, n=ttconf.BRANCH_GRID_SIZE)
            # set the profiles in the eigenspace of the GTR matrix
            # in the following, we only use the prf_l and prf_r (left and right
            # profiles in the matrix eigenspace)
            self._set_rotated_profiles(node)

            # node is constrained
            if hasattr(node, 'numdate_given') and node.numdate_given is not None:
                if hasattr(node, 'bad_branch') and node.bad_branch==True:
                    print ("Branch is marked as bad, excluding it from the optimization process"
                        " Will be optimized freely")
                    node.numdate_given = None
                    node.abs_t = None
                    #    if there are no constraints - log_prob will be set on-the-fly
                    node.msg_to_parent = None
                else:

                    # set the absolute time in branch length units
                    # the abs_t zero is today, and the direction is to the past

                    # this is the conversion between the branch-len and the years
                    node.abs_t = (utils.numeric_date() - node.numdate_given) * abs(self.date2dist.slope)
                    node.msg_to_parent = utils.delta_fun(node.abs_t, return_log=True, normalized=False)
            # unconstrained node
            else:
                node.numdate_given = None
                node.abs_t = None
                # if there are no constraints - log_prob will be set on-the-fly
                node.msg_to_parent = None


    def _convolve(self, src_neglogprob, src_branch_neglogprob, inverse_time):
        """
        Compute the convolution of parent (target) and child (source)
        nodes negative log-likelihood distributions.
        Take the source node log-LH distribution, extracts its grid. Based on
        the branch length probability distribution (also neg log-LH), find
        approximate position of the target node. Make the grid for the target
        node, and for each point of this newly generated grid, compute the
        convolution over all possible positions of the source node.

        Args:

        - src_neglogprob (scipy.interpolate.interp1d): neg log-LH
         distribution of the node to be integrated, represented as scipy
         interpolation object

        - src_branch_neglogprob(scipy.interpolate.interp1d): neg log-LH
         distribution of the branch lenghts between the two nodes, represented
         as scipy interpolation object

         - inverse_time (bool): Whether the time should be inversed.
         True if we go from leaves to root (against absolute time scale), and
         the convolution is computed over positions of the child node.
         False if the messages are propagated from root towards leaves (the same
         direction as the absolute time axis), and the convolution is being
         computed over the position of the parent node

        """

        opt_source_pos = utils.min_interp(src_neglogprob)
        opt_branch_len = utils.min_interp(src_branch_neglogprob)
        if inverse_time:
            opt_target_pos = opt_source_pos + opt_branch_len # abs_t
        else:
            opt_target_pos = opt_source_pos - opt_branch_len

        # T
        target_grid = utils.make_node_grid(opt_target_pos)
        target_grid.sort() # redundant
        if hasattr(src_neglogprob, 'delta_pos'): # convolve with delta-fun
            x_axis = target_grid - src_neglogprob.delta_pos
            x_axis[x_axis < ttconf.MIN_T] = ttconf.MIN_T
            x_axis[x_axis > ttconf.MAX_T] = ttconf.MAX_T
            res_y = src_branch_neglogprob(x_axis)
            res = interp1d(target_grid, res_y, kind='linear')
        else: # convolve two different distributions
            pre_b = np.min(src_branch_neglogprob.y)
            pre_n = np.min(src_neglogprob.y)
            src_branch_neglogprob.y -= pre_b
            src_neglogprob.y -= pre_n
            res = utils.convolve(target_grid, src_neglogprob, src_branch_neglogprob)
            src_branch_neglogprob.y += pre_b
            src_neglogprob.y += pre_n
            res.y += pre_b
            res.y += pre_n
        return res

    def _ml_t_leaves_root(self):
        """
        Compute the probability distribution of the internal nodes positions by
        propagating from the tree leaves towards the root. The result of
        this operation are the probability distributions of each internal node,
        conditional on the constraints on leaves in the descendant subtree. The exception
        is the root of the tree, as its subtree includes all the constrained leaves.
        To the final location probability distribution of the internal nodes,
        is calculated via back-propagation in _ml_t_root_to_leaves.

        Args:

         - None: all required parameters are pre-set as the node attributes during
           tree preparation

        Returns:

         - None: Every internal node is assigned the probability distribution in form
           of an interpolation object and sends this distribution further towards the
           root.

        """

        print("Maximum likelihood tree optimization with temporal constraints:"
            " Propagating leaves -> root...")
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents

            if node.is_terminal():
                continue # either have constraints, or will be optimized freely on the way back

            # children nodes with constraints
            msgs_from_clades = [self._convolve(clade.msg_to_parent,
                               clade.branch_neg_log_prob,
                               inverse_time=True)
                               for clade in node.clades if clade.msg_to_parent is not None]
            if len(msgs_from_clades) < 1:  # we need at least one constraint
                continue

            new_neglogprob = utils.multiply_dists(msgs_from_clades)
            node.msg_to_parent = new_neglogprob

    def _ml_t_root_leaves(self):
        """
        Given the location probability distribution, computed by the propagation
        from leaves to root, set the root most-likely location. Estimate the
        tree likelihood. Report the root location probability distribution
        message towards the leaves. For each internal node, compute the final
        location probability distribution based on the pair of messages (from the
        leaves and from the root), and find the most likely position of the
        internal nodes and finally, convert it to the date-time information

        Args:

        - None: all the requires parameters are pre-set in the previous steps.

        Returns:
         - None: all the internal nodes are assigned probability distributions
           of their locations. The branch lengths are updated to reflect the most
           likely node locations.

        """

        print("Maximum likelihood tree optimization with temporal constraints:"
            " Propagating root -> leaves...")
        for node in self.tree.find_clades(order='preorder'):  # ancestors first, msg to children
            if not hasattr(node, "msg_to_parent"):
                print ("ERROR: node has no log-prob interpolation object! "
                    "Aborting.")
            collapse_func = utils.median_interp
            if node.up is None:  # root node
                node.total_prob = utils.delta_fun(collapse_func(node.msg_to_parent),
                                                  return_log=True,normalized=False)
                #node.total_prefactor = node.msg_to_parent_prefactor
                #node.msg_from_parent_prefactor = 0
                self._set_final_date(node)
                continue

            if node.msg_to_parent is not None: # constrained terminal
                                              # and all internal nodes


                if not hasattr(node.up.total_prob ,'delta_pos'):
                    print ("Cannot infer the position of the node: the position "
                           "of the parent is not delta function")
                    continue

                node_grid = node.up.total_prob.delta_pos - node.branch_neg_log_prob.x
                node_grid[node_grid < ttconf.MIN_T/2] = ttconf.MIN_T
                node_grid[node_grid > ttconf.MAX_T/2] = ttconf.MAX_T
                node.msg_from_parent = interp1d(node_grid, node.branch_neg_log_prob.y, kind='linear')

                final_prob = utils.multiply_dists((node.msg_from_parent, node.msg_to_parent))

                if collapse_func(final_prob) > node.up.abs_t + 1e-9:
                    # must never happen, just for security
                    # I think this can sometimes happen when using median collapsing
                    if self.debug: import ipdb; ipdb.set_trace()
                    node.total_prob = utils.delta_fun(node.up.abs_t, return_log=True, normalized=False)
                    print ("Warn: the child node wants to be {0} earlier than "
                        "the parent node. Setting the child location to the parent's "
                        "one.".format((collapse_func(final_prob) - node.up.abs_t)))

                else:
                    node.total_prob = utils.delta_fun(collapse_func(final_prob),
                        return_log=True, normalized=False)

            else: # unconstrained terminal nodes
                node_grid = node.up.total_prob.delta_pos - node.branch_neg_log_prob.x
                node_grid[node_grid < ttconf.MIN_T/2] = ttconf.MIN_T
                node_grid[node_grid > ttconf.MAX_T/2] = ttconf.MAX_T

                node.msg_from_parent = interp1d(node_grid, node.branch_neg_log_prob.y, kind='linear')
                #final_prob = utils.multiply_dists((node.msg_from_parent, node.msg_to_parent))
                #node.msg_from_parent = msg_from_parent

                node.total_prob = utils.delta_fun(collapse_func(node.msg_from_parent),

                        return_log=True, normalized=False)

            self._set_final_date(node)

    def _set_final_date(self, node):
        """
        Given the location of the node in branch length units, convert it to the
        date-time information.

        Args:
         - node(Phylo.Clade): tree node. NOTE the node should have the abs_t attribute
         to have a valid value. This is automatically taken care of in the
         procedure to get the node location probability distribution.

        """
        node.abs_t = utils.min_interp(node.total_prob)
        if node.up is not None:
            node.branch_length = node.up.abs_t - node.abs_t
            node.dist2root = node.up.dist2root + node.branch_length
        else:
            node.branch_length = self.one_mutation
            node.dist2root = 0.0



        node.years_bp = self.date2dist.get_date(node.abs_t)
        if node.years_bp < 0:
            if not hasattr(node, "bad_branch") or node.bad_branch==False:
                raise ArithmeticError("The node is later than today, but it is not"
                    "marked as \"BAD\", which indicates the error in the "
                    "likelihood optimization.")
            else:
                print ("Warning! node, which is marked as \"BAD\" optimized "
                    "later than present day")

        now = utils.numeric_date()
        node.numdate = now - node.years_bp

        # set the human-readable date
        days = 365.25 * (node.numdate - int(node.numdate))
        year = int(node.numdate)
        try:
            n_date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=days)
            node.date = datetime.datetime.strftime(n_date, "%Y-%m-%d")
        except:
            # this is the approximation
            node.date = str(year) + "-" + str( int(days / 30)) + "-" + str(int(days % 30))


    def coalescent_model(self, Tc=None, optimize_Tc = False,**kwarks):
        """
        This is a wrapper function doing the full inference of node placing and
        ancestral sequences. In addition to standard branch length probabilie,
        a coalescent model prior is used. The branch length probabilities are adjusted to
        reflect the number of concurrent branches with which the branch could
        potentially merge.
        Args:
            - Tc(float): coalescent time scale, if None, half the depth of the
              tree is used. This is expected to be accurate when the tree is not
              ladder like
            - optimize_Tc(bool):  adjust the coalescence time scale to optimize the
              likelihood.
        Returns:
            - None. attributes are added to the nodes
        """
        from merger_models import coalescent
        #  propagate messages up and down and reconstruct ancestral states
        self.ml_t(**kwarks)

        # if no coalescence time scale is provided, use half the root time
        if Tc is None:
            Tc = 0.5*self.tree.root.abs_t

        # resolve polytomies if there are any
        coalescent(self.tree, Tc=Tc)
        self._update_branch_len_interpolators()
        self.resolve_polytomies()

        # if desired, optimize the coalescence time scale
        if optimize_Tc:
            def tmpTotalLH(Tc):
                coalescent(self.tree, Tc=Tc)
                self._update_branch_len_interpolators()
                self.ml_t()
                print("Tc:",Tc)
                self.print_lh()
                return -self.total_LH()
            sol = sciopt.minimize_scalar(tmpTotalLH, bounds=[0, Tc*5], method='bounded')
            if sol['success']:
                self.Tc_opt = sol['x']
                print('coalescent time scale optimization successful, Tc_opt=',self.Tc_opt)
                # final run with optimal Tc
                Tc = self.Tc_opt
            else:
                print('coalescent time scale optimization failed')


    def ml_t(self, max_iter = 3):
        """
        Perform the maximum-likelihood -- based optimization of the tree with temporal
        constraints of (some) internal nodes.

        Args:
         - None
        Returns:
         - None: Updates the tree, its branch lengths and information about the
         internal nodes.
        """
        #  propagate messages up
        self._ml_t_leaves_root()
        #  propagate messages down - reconstruct node positions
        self._ml_t_root_leaves()
        Ndiff = self.reconstruct_anc(method='ml')

        niter=0
        while Ndiff>0 and niter<max_iter:
            print('rerunning treetime inference iteration', niter+1, 'number of state changes observed:',Ndiff)
            self._ml_t_init(ancestral_inference=False)
            self._ml_t_leaves_root()
            self._ml_t_root_leaves()
            Ndiff = self.reconstruct_anc(method='ml')
            niter+=1
        print ("Done tree optimization after",niter+1,"iterations, final state changes:",Ndiff)

    def _set_rotated_profiles(self, node):
        """
        Set sequence and its profiles in the eigenspace of the transition
        matrix.
        """
        node.prf_r = node.profile.dot(self.gtr.v)
        node.prf_l = (self.gtr.v_inv.dot(node.profile.T)).T

    def _score_branch(self, node):
        """
        Auxilliary function to see how well is the particular branch optimized
        (how far it is from its optimal value)
        """
        cmap = mpl.cm.get_cmap ()
        def dev(n):
            sign = np.sign(node.branch_length - utils.opt_branch_len(node))
            opt_bl = sign * abs(node.branch_neg_log_prob(utils.opt_branch_len(n))
                                - node.branch_neg_log_prob(node.branch_length))
            return opt_bl

        node._score = dev(node)
        return None

    def _score_branches(self):
        """
        Set score to the branch. The score is how far is the branch length from
        its optimal value
        """
        for n in self.tree.find_clades():
            if n.up is not None:
                self._score_branch(n)

    def log_lh(self, node):
        """
        Get log-likelihood of the tree given the constrained leaves.
        """
        if hasattr(node, 'lh_prefactor') and hasattr(node, 'msg_to_parent_prefactor'):
            return -node.root.msg_to_parent_prefactor + node.lh_prefactor.sum()
        else:
            return ttconf.MIN_LOG


    def resolve_polytomies(self, merge_compressed=False):
        """
        Resolve the polytomies on the tree.
        The function scans the tree, resolves polytomies in case there are any,
        and re-optimizes the tree with new topology.
        Args:
            - merge_compressed(bool): whether to keep compressed branches as
              polytomies or return a strictly binary tree.
        """
        print('resolving polytomies')
        for n in self.tree.find_clades():
            if len(n.clades) > 3: self._poly(n, merge_compressed)

        # reoptimize branch length and sequences after topology changes
        self.optimize_branch_len()
        self.optimize_seq_and_branch_len(prune_short=False)
        self._ml_t_init()
        self.ml_t()
        self.tree.ladderize()


    def _poly(self, clade, merge_compressed, verbose=1):
        """
        Function to resolve polytomies for a given parent node. If the number of the
        direct decendants is less than three (not a polytomy), does nothing.
        Otherwise, for each pair of nodes, assess the possible LH increase which could be
        gained by merging the two nodes. The increase in the LH is basically the
        tradeoff between the gain of the LH due to the changing the branch lenghts towardsthe optimal
        values and the decrease due to the introduction of the new branch with zero
        optimal length. After the cost gains been determined,
        """

        # TODO coefficient from the gtr
        zero_branch_slope = self.one_mutation / 0.8

        def _c_gain(t, n1, n2, parent):
            """
            cost gain if nodes n1, n2 are joined and their parent is placed at time t

            cost gain = (LH loss now) - (LH loss when placed at time t)
                      = [LH(opt) - LH(now)] - [LH(opt) - LH(t)] approx.=
                      approx.= LH(branch_len(now)) - LH (branch_len(t))

            """
            cg2 = n2.branch_neg_log_prob(parent.abs_t - n2.abs_t ) - n2.branch_neg_log_prob (t - n2.abs_t)
            cg1 = n1.branch_neg_log_prob(parent.abs_t - n1.abs_t ) - n1.branch_neg_log_prob (t - n1.abs_t)
            cg_new = - zero_branch_slope * (parent.abs_t - t) # loss in LH due to the new branch
            return -(cg2+cg1+cg_new)

        def cost_gain(n1, n2, parent):
            """
            cost gained if the two nodes would have been connected.
            """
            cg = sciopt.minimize_scalar(_c_gain,
                    bounds=[np.max(n1.abs_t,n2.abs_t), parent.abs_t],
                    method='Bounded',args=(n1,n2, parent))
            return cg['x'], - cg['fun']

        def merge_nodes(source_arr):
            mergers = np.array([[cost_gain(n1,n2, clade) for n1 in source_arr]for n2 in source_arr])
            while len(source_arr) > 1:
                #print (len(source_arr))

                LH = 0

                # max possible gains of the cost when connecting the nodes:
                # this is only a rough approximation because it assumes the new node positions
                # to be optimal
                new_positions = mergers[:,:,0]
                cost_gains = mergers[:,:,1]
                np.fill_diagonal(cost_gains, -1e9)

                idxs = np.unravel_index(cost_gains.argmax(),cost_gains.shape)
                try:
                    assert (idxs[0] != idxs[1])
                except:
                    if self.debug:
                        import ipdb; ipdb.set_trace()
                    else:
                        print("problem merging nodes")
                n1, n2 = source_arr[idxs[0]], source_arr[idxs[1]]
                if self.debug:
                    print (n1,n2)
                    print ("Delta-LH = " + str(cost_gains[idxs].round(3)))
                LH += cost_gains[idxs]

                new_node = Phylo.BaseTree.Clade()

                # fix positions and branch lengths
                new_node.abs_t = new_positions[idxs] # (n1.abs_t + tree.opt_branch_len(n1) + n2.abs_t + tree.opt_branch_len(n2))/2
                new_node.branch_length = clade.abs_t - new_node.abs_t
                new_node.clades = [n1,n2]
                n1.branch_length = new_node.abs_t - n1.abs_t
                n2.branch_length = new_node.abs_t - n2.abs_t

                # set parameters for the new node
                new_node.up = clade
                n1.up = new_node
                n2.up = new_node
                new_node.sequence = clade.sequence
                new_node.profile = clade.profile
                new_node.mutations = []
                new_node.merger_rate = clade.merger_rate
                self._make_branch_len_interpolator(new_node, n=ttconf.BRANCH_GRID_SIZE)
                clade.clades.remove(n1)
                clade.clades.remove(n2)
                clade.clades.append(new_node)

                # and modify source_arr array for the next loop
                if len(source_arr)>3: # if more than 3 nodes in polytomy, replace row/column
                    for ii in np.sort(idxs)[::-1]:
                        tmp_ind = np.arange(mergers.shape[0])!=ii
                        mergers = mergers[tmp_ind].swapaxes(0,1)
                        mergers = mergers[tmp_ind].swapaxes(0,1)

                    source_arr.remove(n1)
                    source_arr.remove(n2)
                    new_gains = np.array([[cost_gain(n1,new_node, clade) for n1 in source_arr]])
                    mergers = np.vstack((mergers, new_gains)).swapaxes(0,1)

                    source_arr.append(new_node)
                    new_gains = np.array([[cost_gain(n1,new_node, clade) for n1 in source_arr]])
                    mergers = np.vstack((mergers, new_gains)).swapaxes(0,1)
                else: # otherwise just recalculate matrix
                    source_arr.remove(n1)
                    source_arr.remove(n2)
                    source_arr.append(new_node)
                    mergers = np.array([[cost_gain(n1,n2, clade) for n1 in source_arr]
                                       for n2 in source_arr])

            return LH

        stretched = [c for c  in clade.clades if utils.opt_branch_len(c) < c.branch_length]
        compressed = [c for c in clade.clades if c not in stretched]

        if verbose>5:
            print (stretched)
        LH = 0.0

        if len(stretched)==1:
            return LH

        merge_nodes(stretched)
        if merge_compressed:
            merge_nodes(compressed)

        return LH

    def print_lh(self):
        """
        Print the total likelihood of the tree given the constrained leaves
        """
        s_lh = -self.tree.sequence_LH
        t_lh = self.tree.root.msg_to_parent.y.min()

        print ("###  Tree Likelihood  ###\n"
                " Seq log-LH:      {0}\n"
                " Temp.Pos log-LH: {1}\n"
                " Total log-LH:    {2}\n"
               "#########################".format(s_lh, t_lh, s_lh+t_lh))

    def total_LH(self):
        s_lh = self.tree.sequence_LH
        t_lh = -self.tree.root.msg_to_parent.y.min()
        return s_lh+t_lh


    def relaxed_clock(self, slack=None, coupling=None):
        """
        Allow the mutation rate to vary on the tree (relaxed molecular clock).
        Changes of the mutation rates from one branch to another are penalized.
        In addition, deviations of the mutation rate from the mean rate are
        penalized.
        """
        c=1.0
        if slack is None: slack=ttconf.MU_ALPHA
        if coupling is None: coupling=ttconf.MU_BETA
        stiffness = (self.tree.count_terminals()/self.tree.total_branch_length())
        for node in self.tree.find_clades(order='postorder'):
            if node.up is None:
                opt_len = node.branch_length
            else:
                opt_len = self.gtr.optimal_t(node.profile, node.up.profile)
                #opt_len = 1.0*len(node.mutations)/node.profile.shape[0]
            # contact term: stiffness*(g*bl - bl_opt)^2 + slack(g-1)^2 =
            #               (slack+bl^2) g^2 - 2 (bl*bl_opt+1) g + C= k2 g^2 + k1 g + C
            node._k2 = slack + stiffness*node.branch_length**2/(c*opt_len+self.one_mutation)
            node._k1 = -2*(stiffness*node.branch_length*opt_len/(c*opt_len+self.one_mutation) + slack)
            # coupling term: \sum_c coupling*(g-g_c)^2 + Cost_c(g_c|g)
            # given g, g_c needs to be optimal-> 2*coupling*(g-g_c) = 2*child.k2 g_c  + child.k1
            # hence g_c = (coupling*g - 0.5*child.k1)/(coupling+child.k2)
            # substituting yields
            for child in node.clades:
                denom = coupling+child._k2
                node._k2 += coupling*(1.0-coupling/denom)**2 + child._k2*coupling**2/denom**2
                node._k1 += (coupling*(1.0-coupling/denom)*child._k1/denom \
                            - coupling*child._k1*child._k2/denom**2 \
                            + coupling*child._k1/denom)


        all_gammas = []
        for node in self.tree.find_clades(order='preorder'):
            if node.up is None:
                node.gamma = - 0.5*node._k1/node._k2
            else:
                node.gamma = (coupling*node.up.gamma - 0.5*node._k1)/(coupling+node._k2)
            all_gammas.append(node.gamma)
        # normalize avg gamma values to avoid drift in overall mutation rate.
        avg_gamma = np.mean(all_gammas)
        for node in self.tree.find_clades(order='preorder'):
            node.gamma/=avg_gamma

        print('reevaluating branch length interpolators')
        self.init_date_constraints()

    def autocorr_molecular_clock(self, slack=None, coupling=None):
        """
        Allow the mutation rate to vary on the tree (relaxed molecular clock).
        Changes of the mutation rates from one branch to another are penalized.
        In addition, deviations of the mutation rate from the mean rate are
        penalized.
        """
        if slack is None: slack=ttconf.MU_ALPHA
        if coupling is None: coupling=ttconf.MU_BETA
        def opt_mu(node):
            if node.up is None: return mu_0
            mu = (node.up.sequence!=node.sequence).mean()/(node.numdate-node.up.numdate)
            #print (mu)
            return mu

        def get_mu_avg():
            muts = 0.0
            years = 0.0
            L = self.tree.get_terminals()[0].sequence.shape[0]
            for node in self.tree.find_clades(order="preorder"):
                if node.up is None: continue
                muts +=    (node.up.sequence!=node.sequence).sum()
                years += node.numdate-node.up.numdate

            return muts/years/L

        mu_0 = get_mu_avg()
        MAX_N = 1000
        D_MU = 1e-10

        def init_iterative():
            for node in self.tree.find_clades(order="preorder"):
                denom = 1 + slack + coupling * (1 + len(node.clades))
                node._Cn = (opt_mu(node) + slack * mu_0) / denom
                node._Bn = coupling / denom
                node._mu_n = mu_0
                node._mu_n1 = 0.0

        init_iterative()
        converged = False
        N = 0
        while not converged:
            delta_mu = 0.0

            # first pass, we set the mu values at N+1 step
            for node in self.tree.find_clades(order="preorder"):
                if node.up is None: continue
                node._mu_n1 = node._Cn + node._Bn * node.up._mu_n + node._Bn * np.sum([0.0] + [k._mu_n for k in node.clades])
                delta_mu += (node._mu_n1 - node._mu_n)**2

            # update the Nth mu value
            for node in self.tree.find_clades(order="preorder"):
                node._mu_n = node._mu_n1

            N += 1

            if N > MAX_N:
                print ("The autocorrelated molecular clock failed to converge.")
                break

            converged = delta_mu < D_MU

        if converged:
            print ("Autocorrelated molecular clock was computed in " + str(N+1)
                + " steps")

            for node in self.tree.find_clades(order="preorder"):
                denom = 1 + slack + coupling * (1 + len(node.clades))
                node._mu_n1 /= mu_0

        else:
            print ("Autocorrelated molecular clock computation has not converged "
                "after " + str(N) + "steps. Computation failed. The mutation rates will be purged now...")

            for node in self.tree.find_clades(order="preorder"):
                denom = 1 + slack + coupling * (1 + len(node.clades))
                del(node._Cn)
                del(node._Bn)
                del(node._mu_n  )
                del(node._mu_n1 )

    def find_best_root_and_regression(self):
        """
        Find the best root for the tree in linear time, given the timestamps of
        the leaves. The branch lengths should be optimized prior to the run;
        the terminal nodes should have the timestamps assigned as numdate_given
        attribute.
        """

        sum_ti = np.sum([node.numdate_given for node in self.tree.get_terminals() if node.numdate_given is not None])
        sum_ti2 = np.sum([node.numdate_given ** 2 for node in self.tree.get_terminals() if node.numdate_given is not None])
        N = len(self.tree.get_terminals())

        #  fill regression terms for one of the two subtrees
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents
            if node.is_terminal():  # skip leaves
                #  will not rely on the standard func - count terminals directly
                node._st_n_leaves = 1
                node._st_sum_di = 0.0
                node._st_sum_diti = 0.0
                node._st_sum_di2 = 0.0

                if node.numdate_given is not None:
                    node._st_sum_ti = node.numdate_given

                node._sum_ti = sum_ti
                node._sum_ti2 = sum_ti2

                continue

            #  theese all account for subtree only (except for the root, which collects whole tree)
            node._st_sum_ti = np.sum([k._st_sum_ti for k in node.clades])
            node._st_n_leaves = np.sum([k._st_n_leaves for k in node.clades])
            node._st_sum_di   = np.sum([k._st_sum_di + k._st_n_leaves * k.branch_length for k in node.clades])
            node._st_sum_diti = np.sum([k._st_sum_diti + k.branch_length * k._st_sum_ti for k in node.clades])
            node._st_sum_di2  = np.sum([k._st_sum_di2 + k._st_sum_di * 2 * k.branch_length + k._st_n_leaves * k.branch_length ** 2 for k in node.clades])

            node._sum_ti = sum_ti
            node._sum_ti2 = sum_ti2


        best_root = self.tree.root

        for node in self.tree.find_clades(order='preorder'):  # root first

            if node.up is None:
                # assign the values for the root node
                node._sum_di   = node._st_sum_di
                node._sum_diti = node._st_sum_diti
                node._sum_di2  = node._st_sum_di2

            else: # basing on the parent, compute the values for regression
                #  NOTE order of the values computation matters
                node._sum_di = node.up._sum_di + (1.0*N-2.0*node._st_n_leaves) * node.branch_length
                node._sum_di2 = node.up._sum_di2 - 4.0*node.branch_length*node._st_sum_di + 2.0 * node.branch_length * node._sum_di + (1.0*N - 2.0 * node._st_n_leaves) * node.branch_length**2
                node._sum_diti = node.up._sum_diti + node.branch_length * (sum_ti - 2.0 * node._st_sum_ti)

            node._R2 = ((1.0*N * node._sum_diti - sum_ti * node._sum_di) / (np.sqrt((1.0*N * sum_ti2 - node._sum_ti**2)*(1.0*N * node._sum_di2 - node._sum_di**2))))**2
            node._beta = ( 1.0 * N * node._sum_diti - sum_ti * node._sum_di ) / (1.0*N*sum_ti2 - sum_ti**2)
            node._alpha = 1.0 / N *node._sum_di - 1.0 / N * node._beta * sum_ti

            if (node._R2 > best_root._R2 and node._beta>0) or best_root._beta<0:
                best_root = node
                print("Better root found: R2:", best_root._R2, " slope:", best_root._beta)

        return best_root, best_root._alpha, best_root._beta

    def reroot_to_best_root(self,infer_gtr = False, **kwarks):
        '''
        determine the node that, when the tree is rooted on this node, results
        in the best regression of temporal constraints and root to tip distances
        '''
        best_root, a, b = self.find_best_root_and_regression()
        # first, re-root the tree
        self.tree.root_with_outgroup(best_root)
        # set the date2dist params
        self.date2dist = utils.DateConversion()
        self.date2dist.slope = best_root._beta
        self.date2dist.intercept = best_root._alpha
        self.date2dist.r_val = best_root._R2
        # finally, re-compute the basic tree params as we changed the root
        self.tree.ladderize()
        self.tree.root.branch_length = self.one_mutation
        self.tree.root.numdate_given = None
        # fix tree lengths, etc
        self.set_additional_tree_params()
        if infer_gtr:
            self.infer_gtr()
        self.init_date_constraints()



