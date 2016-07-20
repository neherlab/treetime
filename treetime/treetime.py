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
import json
import copy
from scipy import optimize as sciopt
from scipy.ndimage import binary_dilation
from weakref import WeakKeyDictionary

class _Descriptor_Distribution(object):
    """
    Descriptor to manage the settings, common for the LH distributions of different types
    (Branch lengths, messages, LH distributions for the positions, etc)
    It controls the values of the distributions, and, when setting the distribution
    as the attribute to the node, it automatically assigns some additional parameters
    to the same node (incl the sigma)
    """
    def __init__(self, sigma_name, default=None):
        self.default = default
        self._sigma_name = sigma_name
        self.data = WeakKeyDictionary()

    def __get__(self, instance, owner):
        # we get here when someone calls x.d, and d is a NonNegative instance
        # instance = x
        # owner = type(x)
        return self.data.get(instance, self.default)

    def __set__(self, instance, value):
        # we get here when someone calls x.d = val, and d is a NonNegative instance
        # instance = x
        # value = val
        if value is None:
            self.data[instance] = self.default
            return

        if not isinstance(value, interp1d):
            raise TypeError("Cannot set the branch length LH distribution property."
                "The interpolation object is expected")

        self.data[instance] = value
        # compute the connected parameters:
        instance.__setattr__(self._sigma_name, _Descriptor_Distribution._logprob_sigma(value))

    @staticmethod
    def _logprob_sigma(logprob):
        """
        Assess the width of the probability distribution.
        """
        if logprob is None: return 0.0

        if not isinstance(logprob, interp1d):
            raise TypeError("Error in computing the sigma for the branch length LH distribution."
                "Interpolation object expected")
        ymin = logprob.y.min()
        real_prob = np.exp(-(logprob.y-ymin))
        xs = logprob.x[real_prob > (real_prob.max() - real_prob.min()) / 2]
        return xs.max() - xs.min()

##  set the necessary descriptors to the Phylo.Clade objects
##  (NOTE all descriptors assigned at object level)
Phylo.BaseTree.Clade.branch_neg_log_prob = _Descriptor_Distribution("branch_len_sigma")
Phylo.BaseTree.Clade.msg_to_parent = _Descriptor_Distribution("msg_to_parent_sigma")

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

    def init_date_constraints(self, slope=None, **kwarks):
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
        self._ml_t_init(**kwarks)

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


        if not hasattr(node, 'gamma'):
            node.gamma = 1.0

        if not hasattr(node, 'merger_rate') or node.merger_rate is None:
            node.merger_rate = ttconf.BRANCH_LEN_PENALTY

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

        min_prob = np.min(logprob)
        logprob -= min_prob

        tmp_prob = np.exp(-logprob)
        integral = self._integral(grid, tmp_prob)

        # save raw branch length interpolators without coalescent contribution
        # TODO: better criterion to catch bad branch
        if integral < 1e-200:
            print ("!!WARNING!!", node.name, " branch length probability distribution",
                "integral is ZERO. Setting bad_branch flag...")
            node.bad_branch = True
            node.raw_branch_neg_log_prob = interp1d(grid, logprob, kind='linear')

        else:
            node.raw_branch_neg_log_prob = interp1d(grid, logprob+np.log(integral), kind='linear')


        # add merger rate contribution to the raw branch length
        logprob += node.merger_rate * np.minimum(ttconf.MAX_BRANCH_LENGTH, np.maximum(0,grid))

        tmp_prob = np.exp(-logprob)
        integral = self._integral(grid, tmp_prob)

        if integral < 1e-200:
            print ("!!WARNING!! Node branch length probability distribution "
                "integral is ZERO. Setting bad_branch flag...")
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

        if ttconf.BRANCH_LEN_PENALTY is None:
            ttconf.BRANCH_LEN_PENALTY = 0.0

        if ancestral_inference:
            self.optimize_seq_and_branch_len(**kwarks)

        print('Initializing branch length interpolation objects')
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

    def _integral(self, x, y):
        """
        compute the integral of y(x) value over all range of x
        """
        dx = np.diff(x)
        res = (0.5*(y[1:]+y[:-1])*dx).sum()
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
        def _send_message(node, **kwargs):
            """
            Calc the desired LH distribution of the parent
            """
            if hasattr(node.msg_to_parent, 'delta_pos'): # convolve distribution  with delta-fun
                target_grid = node.branch_neg_log_prob.x + node.msg_to_parent.delta_pos
                target_grid[target_grid < ttconf.MIN_T/2] = ttconf.MIN_T
                target_grid[target_grid > ttconf.MAX_T/2] = ttconf.MAX_T
                res = interp1d(target_grid, node.branch_neg_log_prob.y, kind='linear')
            else: # convolve two distributions
                target_grid =  self._conv_grid(node.msg_to_parent, node.branch_neg_log_prob,
                    base_s=node.msg_to_parent_sigma, prop_s=node.branch_len_sigma, inverse_time=True)
                    #self._conv_grid(node)
                res = self._convolve(target_grid, node.msg_to_parent, node.branch_neg_log_prob, inverse_time=True)
            return res

        print("Maximum likelihood tree optimization with temporal constraints:"
            " Propagating leaves -> root...")


        # go through the nodes from leaves towards the root:
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents

            if node.is_terminal():
                node.msgs_from_leaves = {}
                continue # either have constraints, or will be optimized freely on the way back

            # save all messages from the children nodes with constraints
            # store as dictionary to exclude nodes from the set when necessary
            # (see below)
            node.msgs_from_leaves = {clade: _send_message(clade) for clade in node.clades
                                            if clade.msg_to_parent is not None}


            if len(node.msgs_from_leaves) < 1:  # we need at least one constraint
                continue

            # this is what the node sends to the parent
            node_grid = self._make_node_grid(node)
            node.msg_to_parent = self._multiply_dists(node.msgs_from_leaves.values(), node_grid)


    def _convolve(self, time, f_func, g_func, inverse_time=None, cutoff=1e7, n_integral=100):

        """
        Compute convolution of the functions f, g:= backwards and forward in time:
        if inverse_time is False, then the function computes  the proper convolution:
            F(t) = (f*g) = integral {f(t-tau) g(tau) d_tau}.
        Otherwise, the result  is the "adapted" convolution to account the backwards message passing:
            F(t) = integral {f(t+tau) g(tau) d_tau}.
        """
        if inverse_time is None:
            raise RunTimeError("temporal direction of convolution not specified")

        # get the support ranges for the raw functions
        #include first points below the cutoff to have at least three points in the integration region
        frange = binary_dilation((f_func.y - f_func.y.min()) < cutoff)
        grange = binary_dilation((g_func.y - g_func.y.min()) < cutoff)
        # nothing to convolve, one of the distribution is flat zero
        if (frange.sum() == 0 or grange.sum() == 0):
            # we lost the probability
            # NOTE binary_dilation does not extend the False array
            print ("Function F values: \n" + "\t".join(map(str,f.y)) + "\n")
            print ("Function G values: \n" + "\t".join(map(str,g.y)) + "\n")
            raise ArithmeticError("Cannot convolve functions. At least one of the "
                                  "probability distributions has been lost!")
        fx_min = f_func.x[frange].min()
        fx_max = f_func.x[frange].max()
        gx_min = g_func.x[grange].min()
        gx_max = g_func.x[grange].max()
        # resulting convolution
        res = np.ones(time.shape[0]) * 1e8
        for idx, ti in enumerate(time):
            if inverse_time:
                tau_min = np.max((ti-fx_max, gx_min))
                tau_max = np.min((ti-fx_min, gx_max))
            else:
                tau_min = np.max((fx_min - ti, gx_min))
                tau_max = np.min((fx_max - ti, gx_max))
            if (tau_min > tau_max):
                # functions not overlap
                continue

            # get the step for the grid
            #dtau = np.min((
            #    (f_func.x[frange][-1] - f_func.x[frange][0]) / 10, # if f sharp - at least 10 points cover f range
            #    (g_func.x[grange][-1] - g_func.x[grange][0]) / 10, # if g sharp - at least 10 points cover g range
            #    (tau_max - tau_min) / 100.0)) # normal situation, regular grid of 100 points

            tau = np.linspace(tau_min, tau_max, n_integral)
            # include the distributions extremum positions to the grid to avoid round error:
            #tau = np.concatenate(((f_func.x[f_func.y.argmin()], g_func.x[g_func.y.argmin()]), tau))
            # make sure the values are unique (NOTE unique method also sorts in-place)
            #tau = np.unique(tau)

            if len(tau) < 2:
                #print "Cannot convolve the distributions: functions do not overlap!"
                continue
            if inverse_time:
                fg = f_func(ti-tau) + g_func(tau)
            else:
                fg = f_func(ti+tau) + g_func(tau)

            min_fg = fg.min() # exponent pre-factor
            expfg = np.exp(-1*(fg-min_fg))
            dtau = np.diff(tau)
            integral = (0.5*(expfg[1:]+expfg[:-1])*dtau).sum()
            res[idx] = -np.log(integral) + min_fg

        res[-1] = 1e8
        res[-2] = 1e8
        res = interp1d(time, res, kind='linear')
        return res


    def _conv_grid(self, base_dist, propagator, base_s=None, prop_s=None,
                    inverse_time=None, sigma_factor=4, n_points=600, **kwargs):
        """
        Make the grid for the two convolving functions
        # NOTE! n_points largely defines the smoothness of the final distribution
        """
        if inverse_time is None:
            raise RunTimeError("temporal direction of convolution not specified")

        ##  compute the width of the message from parent
        if base_s is None:
            base_s = _Descriptor_Distribution._logprob_sigma(base_dist)
        if prop_s is None:
            prop_s = _Descriptor_Distribution._logprob_sigma(propagator)

        # FIXME not min of the dist, but average (especially for two zero-max dists)
        T0_base  = utils.min_interp(base_dist)
        T0_prop = utils.min_interp(propagator)

        # rough estimate for the convolution maximum (max overlap of the distributions)
        if inverse_time: # we move towards the root, but the time values (years before present) increase
            T0 = T0_base + T0_prop
        else: # we move towards leaves (years before present decrease)
            T0 = T0_base - T0_prop

        # rough estimate for the convolution width (sum of dists widths)
        T_sigma = sigma_factor * (base_s + prop_s)

        # TODO sometimes, the distribution is "cut" from one side,
        # we do not need the symmetric grid in this case
        return self._make_grid(T0, T_sigma, n_points)


    def _make_node_grid(self, node, sigma_factor=6, n_points=300, cutoff=1e7, **kwargs):

        if hasattr(node, 'msgs_from_leaves'):

            xmin = np.max([scx.x[binary_dilation(scx.y < cutoff).argmax()] for scx in node.msgs_from_leaves.values()])

            pos = np.array([utils.min_interp(msg) for msg in node.msgs_from_leaves.values()])
            sigmas = np.array([_Descriptor_Distribution._logprob_sigma(msg) for msg in node.msgs_from_leaves.values()])
            steps = np.array([sigmas[i] / ((msg.x > pos[i] -sigmas[i]) & (msg.x < pos[i] + sigmas[i])).sum()
                        for i,msg in enumerate(node.msgs_from_leaves.values())])

        else:

            xmin = None

            pos = np.array([])
            sigmas = np.array([])
            densities = np.array([])

        # account for the node position, which has bee inferred from the molecular clock
        _pos, _sigma = self._opt_node_pos(node)
        pos = np.concatenate((pos, [_pos]))
        sigmas = np.concatenate((sigmas,[_sigma]))
        steps = np.concatenate((steps, [_sigma * sigma_factor / n_points]))
        steps = steps[steps>0]


        # choose the finest grid in the selected region
        extreme_pos = np.concatenate((pos-sigmas, pos+sigmas))

        Npoints = (extreme_pos.max() - extreme_pos.min())/steps.min()

        return self._make_grid((extreme_pos.max() + extreme_pos.min()) / 2,
           (extreme_pos.max() - extreme_pos.min()) , Npoints, xmin)

    def _opt_node_pos(self, node):
        """
        Make a guess about the optimal position of a node.
        Returns:
         - opt_pos: guess for the optimal position
         - sigma: guess for the opt_pos RMSE
        """

        numdate = (node.dist2root - self.date2dist.intercept) / self.date2dist.slope
        opt_pos = self.date2dist.get_abs_t(numdate)
        sigma = self.date2dist.rms
        return opt_pos, sigma

    def _make_grid(self, center, sigma, N, xmin=None, xmax=None):

        alpha=1.0
        grid_center = center + sigma * np.sign(np.linspace(-1, 1, N/2)) * np.abs(np.linspace(-1, 1, N/2)**alpha)

        # points for the "wings" grid
        pts = np.arange(0, int(N/2))
        N_pts = pts.shape[0]


        # derivatives and values of node grid position as a function of grid index
        center_deriv_r = grid_center[-1]-grid_center[-2]
        start_point_r = grid_center[-1]
        end_point_r = ttconf.MAX_BRANCH_LENGTH

        center_deriv_l = grid_center[1]-grid_center[0]
        start_point_l = grid_center[0]
        end_point_l =  - ttconf.MAX_BRANCH_LENGTH

        grid_wings_r = pts**2 * (end_point_r - center_deriv_r*N_pts - start_point_r) / N_pts**2 + center_deriv_r * pts + start_point_r
        grid_wings_l = pts**2 * (end_point_l - center_deriv_l*N_pts - start_point_l) / N_pts**2 + center_deriv_l * pts + start_point_l

        if xmin is not None:
            grid_center = grid_center[grid_center > xmin]
            grid_wings_r = grid_wings_r[grid_wings_r > xmin]
            grid_wings_l = grid_wings_l[grid_wings_l > xmin]

        if xmax is not None:
            grid_center = grid_center[grid_center < xmax]
            grid_wings_r = grid_wings_r[grid_wings_r < xmax]
            grid_wings_l = grid_wings_l[grid_wings_l < xmax]

        #grid_wings_r = center + sigma + ttconf.MAX_BRANCH_LENGTH * (np.linspace(0, 1, N/2)**2) [1:]
        #grid_wings_l = center - sigma - ttconf.MAX_BRANCH_LENGTH * (np.linspace(0, 1, N/2)**2) [1:]

        grid = np.unique(np.concatenate( # unique will also sort the array
            ([ttconf.MIN_T],
             grid_center,
             grid_wings_r,
             grid_wings_l,
             [ttconf.MAX_T])))

        return grid

    def _multiply_dists(self, dists, grid):

        res = np.sum((k(grid) for k in dists))
        res[((0,-1),)] = -1 * ttconf.MIN_LOG
        res[((1,-2),)] = -1 * ttconf.MIN_LOG / 2
        interp = interp1d(grid, res, kind='linear')
        return interp

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

        collapse_func = utils.median_interp

        def _send_message(msg_parent_to_node, branch_lh, msg_parent_s=None, branch_lh_s=None, **kwargs):

            """
            Compute the 'back-convolution' when propagating from the root to the leaves.
            The message from the parent to the child node is the
            $$
            Pr(T_c | comp_tips_pos) = \int_{d\tau} Pr(L = \tau) Pr(T_p = T_c + \tau)
            $$

            """
            if hasattr(msg_parent_to_node, 'delta_pos'): # convolve distribution  with delta-fun

                #grid will be same as for the branch len
                target_grid = msg_parent_to_node.delta_pos - branch_lh.x
                target_grid[target_grid < ttconf.MIN_T/2] = ttconf.MIN_T
                target_grid[target_grid > ttconf.MAX_T/2] = ttconf.MAX_T
                res = interp1d(target_grid, branch_lh.y, kind='linear')

            else: # all other cases
                # make the grid for the node
                target_grid = self._conv_grid(msg_parent_to_node, branch_lh,
                                msg_parent_s, branch_lh_s, inverse_time=False)
                res = self._convolve(target_grid, msg_parent_to_node, branch_lh, inverse_time = False)

            return res

        def _set_joint_lh_pos(node):
            """
            Compute the joint LH distribution and collapse it to the delta-function,
            which is later will be converted to the node's position
            """
            if not hasattr(node.up, 'joint_lh_pos') or \
               not hasattr(node.up.joint_lh_pos ,'delta_pos'):
               #  the joint lh is missing in the parent node, or it is not the delta function
                print ("Cannot infer the joint distribution of the node: {0}. "
                       "The position "
                       "of the parent is not delta function").format(node.name)
                node.joint_lh_pos =  utils.delta_fun(node.up.abs_t, return_log=True, normalized=False)
                return;

            # compute the jont LH pos
            joint_msg_from_parent = _send_message(node.up.joint_lh_pos, node.branch_neg_log_prob)

            if node.msg_to_parent is not None:
                node.joint_lh = self._multiply_dists((joint_msg_from_parent, node.msg_to_parent),
                                                    node.msg_to_parent.x)
            else: #unconstrained terminal node
                node.joint_lh = joint_msg_from_parent


            # median/mean of LH dist is the node's date
            node_date = collapse_func(node.joint_lh)

            if node_date > node.up.abs_t + 1e-9:
                if self.debug: import ipdb; ipdb.set_trace()
                # collapse to the parent's position
                node.joint_lh_pos = utils.delta_fun(node.up.abs_t, return_log=True, normalized=False)
                print ("Warn: the child node wants to be {0} earlier than "
                    "the parent node. Setting the child location to the parent's "
                    "one. Node name: {1}".format(node_date - node.up.abs_t, node.name))
            else:
                node.joint_lh_pos = utils.delta_fun(node_date, return_log=True, normalized=False)

            return;

        def _set_marginal_lh_dist(node):
            # Compute the marginal distribution for the internal node
            # set the node's distribution
            parent = node.up
            assert (parent is not None)

            # get the set of the messages from the complementary subtree of the node
            # these messages are ready-to-use (already convolved with the branch lengths)
            complementary_msgs = [parent.msgs_from_leaves[k]
                for k in parent.msgs_from_leaves
                if k != node]

            # prepare the message that will be sent from the root to the node
            if parent.msg_from_parent is not None: # the parent is not root => got something from the parent
                complementary_msgs += [parent.msg_from_parent]

            # prepare the message, which will be propagated to the child.
            # NOTE the message is created on the parental grid
            # we reuse the msg_to_parent grid to save some computations
            msg_parent_to_node = self._multiply_dists(complementary_msgs, parent.msg_to_parent.x)

            # propagate the message from the parent to the node:
            node.msg_from_parent = _send_message(msg_parent_to_node, node.branch_neg_log_prob, None, node.branch_len_sigma)

            # finally, set the node marginal LH distribution:

            if node.msg_to_parent is not None:
                # we again reuse the msg_to_parent grid
                node.marginal_lh =  self._multiply_dists((node.msg_from_parent, node.msg_to_parent), node.msg_to_parent.x)
            else: # terminal node without constraints
                # the smoothness of the dist is defined by the grid of the message,
                # no need to create new grid
                node.marginal_lh = node.msg_from_parent


        # Main method - propagate from root to the leaves and set the LH distributions
        # to each node
        for node in self.tree.find_clades(order='preorder'):  # ancestors first, msg to children
            if not hasattr(node, "msg_to_parent"):
                print ("ERROR: node has no log-prob interpolation object! "
                    "Aborting.")

            ## This is the root node
            if node.up is None:
                node.marginal_lh = node.msg_to_parent
                node.msg_from_parent = None # nothing beyond the root
                node.joint_lh_pos = utils.delta_fun(collapse_func(node.marginal_lh),
                                                  return_log=True,normalized=False)
            else:
                _set_joint_lh_pos(node)
                _set_marginal_lh_dist(node)
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
        node.abs_t = utils.min_interp(node.joint_lh_pos)
        if node.up is not None:
            node.branch_length = node.up.abs_t - node.abs_t
            node.dist2root = node.up.dist2root + node.branch_length
        else:
            node.branch_length = self.one_mutation
            node.dist2root = 0.0

        node.years_bp = self.date2dist.get_date(node.abs_t)
        if node.years_bp < 0:
            if not hasattr(node, "bad_branch") or node.bad_branch==False:
                #import ipdb; ipdb.set_trace()
                print("ERROR: The node is later than today, but it is not"
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
        self.ml_t(max_iter=1, **kwarks)

        # if no coalescence time scale is provided, use half the root time
        if Tc is None:
            Tc = 0.5*self.tree.root.abs_t

        # resolve polytomies if there are any
        coalescent(self.tree, Tc=Tc)
        self._update_branch_len_interpolators()
        self.resolve_polytomies(rerun=False)
        self.init_date_constraints(ancestral_inference=True, prune_short=False)
        self.ml_t(max_iter=1)

        # if desired, optimize the coalescence time scale
        if optimize_Tc:
            def tmpTotalLH(Tc):
                coalescent(self.tree, Tc=Tc)
                self._update_branch_len_interpolators()
                self.ml_t(max_iter=1)
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

    def ml_t(self, max_iter = 3,**kwarks):
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

        niter=1
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
        from matplotlib import cm
        cmap = cm.get_cmap ()
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

    def resolve_polytomies(self, merge_compressed=False, rerun=True):
        """
        Resolve the polytomies on the tree.
        The function scans the tree, resolves polytomies in case there are any,
        and re-optimizes the tree with new topology.
        Args:
            - merge_compressed(bool): whether to keep compressed branches as
              polytomies or return a strictly binary tree.
        """
        print('resolving polytomies')
        poly_found=False
        for n in self.tree.find_clades():
            if len(n.clades) > 2:
                self._poly(n, merge_compressed)
                poly_found=True
                #import ipdb; ipdb.set_trace()


        print('Checking for obsolete nodes')
        obsolete_nodes = [n for n in self.tree.find_clades() if len(n.clades)==1 and n.up is not None]
        for node in obsolete_nodes:
            print('remove obsolete node',node.name)
            if node.up is not None:
                self.tree.collapse(node)
        # reoptimize branch length and sequences after topology changes
        if rerun and poly_found:
            print("topology of the tree has changed, will rerun inference...")
            self.optimize_branch_len()
            self.optimize_seq_and_branch_len(prune_short=False)
            self.init_date_constraints(ancestral_inference=False)
            self.ml_t()

        else:
            self._set_each_node_params() # set node info to the new nodes

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

        def merge_nodes(source_arr, isall=False):
            mergers = np.array([[cost_gain(n1,n2, clade) for n1 in source_arr]for n2 in source_arr])
            while len(source_arr) > 1 + int(isall):
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
                if len(source_arr)>2: # if more than 3 nodes in polytomy, replace row/column
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

        if len(stretched)==1 and merge_compressed==False:
            return LH

        merge_nodes(stretched, isall=len(stretched)==len(clade.clades))
        if merge_compressed and len(compressed)>1:
            merge_nodes(compressed, isall=len(compressed)==len(clade.clades))

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
        self.init_date_constraints(ancestral_inference=False)

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
        sum_ti2 = np.sum([node.numdate_given**2 for node in self.tree.get_terminals() if node.numdate_given is not None])
        N = 1.0*len([x for x in self.tree.get_terminals() if x.numdate_given is not None])
        Ninv = 1.0/N
        time_variance = (N*sum_ti2 - sum_ti**2)*Ninv**2

        #  fill regression terms for one of the two subtrees
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents
            if node.is_terminal():  # inititalize the leaves
                #  will not rely on the standard func - count terminals directly
                node._st_n_leaves = 0 if node.numdate_given is None else 1
                node._st_di = 0.0
                node._st_diti = 0.0
                node._st_di2 = 0.0

                if node.numdate_given is not None:
                    node._st_ti = node.numdate_given
                else:
                    node._st_ti = 0

                node._ti = sum_ti
            else:
                #  for non-terminal nodes,
                node._st_ti = np.sum([k._st_ti for k in node.clades])
                node._st_n_leaves = np.sum([k._st_n_leaves for k in node.clades])
                node._st_di   = np.sum([k._st_di + k._st_n_leaves*k.branch_length for k in node.clades])
                node._st_diti = np.sum([k._st_diti + k.branch_length*k._st_ti for k in node.clades])
                node._st_di2  = np.sum([k._st_di2 + 2*k._st_di*k.branch_length + k._st_n_leaves*k.branch_length**2 for k in node.clades])
                node._ti = sum_ti

        best_root = self.tree.root
        for node in self.tree.find_clades(order='preorder'):  # root first

            if node.up is None:
                # assign the values for the root node
                node._di   = node._st_di
                node._diti = node._st_diti
                node._di2  = node._st_di2

                # TODO
                dist_variance = (N*node._di2 - node._di**2)*(Ninv**2)
                disttime_cov = (N*node._diti - sum_ti*node._di)*(Ninv**2)
                time_variance = time_variance

                node._beta = disttime_cov/time_variance
                node._alpha = (node._di - node._beta*sum_ti)/N
                node._R2 = disttime_cov**2/(time_variance*dist_variance)
                node._R2_delta_x = 0 # there is no branch to move the root

            else: # based on the parent, compute the values for regression
                #  NOTE order of the values computation matters
                n_up = N - node._st_n_leaves
                n_down = node._st_n_leaves
                node._di = node.up._di + (n_up-n_down)*node.branch_length
                node._di2 = (node.up._di2 + 2*node.branch_length*node.up._di
                            - 4*(node.branch_length*(node._st_di + n_down*node.branch_length))
                            + N*node.branch_length**2)
                node._diti = node.up._diti + node.branch_length*(sum_ti - 2*node._st_ti)

                L = node.branch_length

                ## Express Node's sum_Di as the function of parent's sum_Di
                # and **displacement from parent's node x** :
                # sum_Di = A1 + A2 * x
                A1 = node.up._di
                A2 = n_up - node._st_n_leaves

                ## Express Node's sum_Di**2 as the function of parent's params
                # and **displacement from parent's node x** :
                # sum_Di = B1 + B2 * x + B3 * x**2
                B1 = node.up._di2
                B2 = 2 * (node.up._di - 2 * node._st_di - 2 * node.branch_length * node._st_n_leaves )
                B3 = N

                ## Express Node's sum_DiTi as the function of parent's params
                # and **displacement from parent's node x** :
                # sum_DiTi = C1 + C2 * x
                C1 = node.up._diti
                C2 = sum_ti - 2 * node._st_ti

                ## Substituting Ai, Bi, Ci to the expression for R2, and
                ## making all the algebra, we get the R2 as the function of the
                ## displacement from the parent's node x:
                # R2(x) = CONST * (alpha * x**2 + beta * x+ gamma) / (mu * x**2 + nu * x + delta)
                # Expressions for alpha, beta, etc through Ai, Bi, Ci:
                alpha = (N * C2 - sum_ti * A2)**2
                beta = 2 * (N*C2 - sum_ti*A2) * (N*C1 - sum_ti*A1)
                gamma = (N*C1 - sum_ti*A1)**2
                mu = N * B3 - A2**2
                nu = N * B2 - 2 * A1 * A2
                delta = N * B1 - A1**2

                # Search for the maximum of R2 in the middle of the branch.
                # Eq: dR2/dx = 0 -> square equation:
                # x**2*(alpha*nu - beta *  mu) + 2x*(alpha*delta-mu*gamma) + (beta*delta - nu*gamma) = 0
                # look for the root(s):
                # Determinant is
                D2 =  (alpha * delta - mu * gamma) ** 2 - (alpha * nu - beta * mu) * (beta * delta - nu * gamma)

                if D2 < 0:
                    # somehow there is no extremum for the R2(x) function
                    x1 = -1 # any arbitrary value out of range [0, L], see below
                    x2 = -1
                else:
                    # actual roots - the extrema for the R2(x) function
                    x1 = (-1 * (alpha * delta - mu * gamma) + D2 **0.5) / (alpha * nu - beta * mu)
                    x2 = (-1 * (alpha * delta - mu * gamma) - D2 **0.5) / (alpha * nu - beta * mu)

                # possible positions, where the new root can possibly be located
                # (restrict to the branch length)
                max_points = [k for k in (x1,x2,L) if k >= 0 and k <= L]
                # values of the R2 at these positions
                R2s = [(alpha * x**2 + beta * x + gamma) / (mu * x**2 + nu * x + delta) / time_variance / N**2 for x in max_points]
                # choose the best R2
                node._R2 = np.max(R2s)
                # and set the position for the best R2 value
                node._R2_delta_x = L - max_points[np.argmax(R2s)]

                # for this position, define the slope and intercept:
                node._beta = ((L - node._R2_delta_x) * (N * C2 - sum_ti * A2) + (N*C1-sum_ti*A1)) / time_variance / N**2
                node._alpha = (L - node._R2_delta_x) * A2 / N  + (A1 - node._beta * sum_ti) / N

            if node.up is None:
                print("Initial root: R2:", best_root._R2, " slope:", best_root._beta)

            if (node._R2 > best_root._R2 and node._beta>0) or best_root._beta<0:
            #if (node._beta>best_root._beta and node._beta>0) or best_root._beta<0:
                best_root = node
                print("Better root found: R2:", best_root._R2,
                    " slope:", best_root._beta,
                    " branch_displacement: ", (best_root._R2_delta_x) / ( node.branch_length + self.one_mutation))

        return best_root, best_root._alpha, best_root._beta

    def reroot_to_best_root(self,infer_gtr = False, n_iqd = None, **kwarks):
        '''
        determine the node that, when the tree is rooted on this node, results
        in the best regression of temporal constraints and root to tip distances
        '''
        best_root, a, b = self.find_best_root_and_regression()
        # first, re-root the tree

        if hasattr(best_root, "_R2_delta_x") and  best_root._R2_delta_x > 0 and best_root.up is not None:
            # create new node in the branch and root the tree to it
            new_node = copy.copy(best_root)
            new_node.clades = [best_root]
            new_node.branch_length = best_root.branch_length - best_root._R2_delta_x
            best_root.branch_length = best_root._R2_delta_x
            best_root.up.clades = [k if k != best_root else new_node for k in best_root.up.clades]
            self.tree.root_with_outgroup(new_node)

        else:
            # simply use the existing node as the new root
            self.tree.root_with_outgroup(best_root)

        if n_iqd is not None:
            root_to_tip = self.tree.depths()
            self.tree.root.up=None
            self.tree.root.branch_length = self.one_mutation
            res = {}
            for node in self.tree.get_terminals():
                if hasattr(node, 'numdate_given') and node.numdate_given is not None:
                    res[node] = root_to_tip[node] - best_root._beta*node.numdate_given - best_root._alpha
            residuals = np.array(res.values())
            iqd = np.percentile(residuals,75) - np.percentile(residuals,25)
            for node,r in res.iteritems():
                if r>n_iqd*iqd:
                    print('marking',node.name,'as outlier, residual',r/iqd, 'interquartile distances')
                    node.bad_branch=True
                    node.numdate_given = None
            # redo root estimation after outlier removal
            best_root, a, b = self.find_best_root_and_regression()
            # first, re-root the tree
            self.tree.root_with_outgroup(best_root)
        print('Checking for obsolete nodes')
        obsolete_nodes = [n for n in self.tree.find_clades() if len(n.clades)==1]
        for node in obsolete_nodes:
            print('remove obsolete node',node.name)
            self.tree.collapse(node)

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
        self.init_date_constraints(ancestral_inference=True, **kwarks)



