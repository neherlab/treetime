import utils
import numpy as np
import config as ttconf
from treeanc import TreeAnc
from distribution import Distribution
from branch_len_interpolator import BranchLenInterpolator
from node_interpolator import NodeInterpolator


class ClockTree(TreeAnc):
    """
    Class to produce molecular clock trees.
    """

    def __init__(self,  dates=None,*args, **kwargs):
        super(ClockTree, self).__init__(*args, **kwargs)
        if dates is None:
            raise("ClockTree requires date contraints!")
        self.date_dict = dates
        self.date2dist = None  # we do not know anything about the conversion
        self.max_diam = 0.0
        self.debug=False

        for node in self.tree.find_clades():
            if node.name in self.date_dict:
                node.numdate_given = self.date_dict[node.name]
            else:
                node.numdate_given = None


    @property
    def date2dist(self):
        return self._date2dist

    @date2dist.setter
    def date2dist(self, val):
        if val is None:
            self._date2dist = None
            return
        else:
            self.logger("TreeTime.date2dist: Setting new date to branchlength conversion. slope=%f, R2=%.4f"%(val.slope, val.r_val), 2)
            self._date2dist = val

    def init_date_constraints(self, ancestral_inference=True, slope=None, **kwarks):
        """
        Get the conversion coefficients between the dates and the branch
        lengths as they are used in ML computations. The conversion formula is
        assumed to be 'length = k*numdate_given + b'. For convenience, these
        coefficients as well as regression parameters are stored in the
        dates2dist object.

        Note: that tree must have dates set to all nodes before calling this
        function. (This is accomplished by calling load_dates func).
        """
        self.logger("ClockTree.init_date_constraints...",2)

        if ancestral_inference or (not hasattr(self.tree.root, 'sequence')):
            self.optimize_seq_and_branch_len(**kwarks)

        # set the None  for the date-related attributes in the internal nodes.
        # make interpolation objects for the branches
        print('\n----- Initializing branch length interpolation objects...\n')
        for node in self.tree.find_clades():
            if node.up is not None:
                node.branch_length_interpolator = BranchLenInterpolator(node, self.gtr, one_mutation=self.one_mutation)
            else:
                node.branch_length_interpolator = None
        self.date2dist = utils.DateConversion.from_tree(self.tree, slope)
        self.max_diam = self.date2dist.intercept

        # make node distribution objects
        for node in self.tree.find_clades():
            # node is constrained
            if hasattr(node, 'numdate_given') and node.numdate_given is not None:
                if hasattr(node, 'bad_branch') and node.bad_branch==True:
                    print ("Branch is marked as bad, excluding it from the optimization process"
                        " Will be optimized freely")
                    node.numdate_given = None
                    node.abs_t = None
                    # if there are no constraints - log_prob will be set on-the-fly
                    node.msg_to_parent = None
                else:
                    # set the absolute time before present in branch length units
                    node.abs_t = (utils.numeric_date() - node.numdate_given) * abs(self.date2dist.slope)
                    node.msg_to_parent = NodeInterpolator.delta_function(node.abs_t, weight=1)

            else: # node without sampling date set
                node.numdate_given = None
                node.abs_t = None
                # if there are no constraints - log_prob will be set on-the-fly
                node.msg_to_parent = None


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
            if node.msg_to_parent.is_delta:
                res = Distribution.shifted_x(node.branch_length_interpolator, node.msg_to_parent.peak_pos)
            else: # convolve two distributions
                res =  NodeInterpolator.convolve(node.msg_to_parent, node.branch_length_interpolator)
                # TODO deal with grid size explosion
            return res

        self.logger("ClockTree: Maximum likelihood tree optimization with temporal constraints:",1)
        self.logger("ClockTree: Propagating leaves -> root...", 2)
        # go through the nodes from leaves towards the root:
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents
            if node.is_terminal():
                node.msgs_from_leaves = {}
            else:
                # save all messages from the children nodes with constraints
                # store as dictionary to exclude nodes from the set when necessary
                # (see below)
                node.msgs_from_leaves = {clade: _send_message(clade) for clade in node.clades
                                                if clade.msg_to_parent is not None}

                if len(node.msgs_from_leaves) < 1:  # we need at least one constraint
                    continue
                # this is what the node sends to the parent
                node.msg_to_parent = NodeInterpolator.multiply(node.msgs_from_leaves.values())



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

        Npoints = int(min(600,(extreme_pos.max() - extreme_pos.min())/steps.min()))

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

            # compute the joint LH pos
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
            # FIXME
            grid = np.unique(np.concatenate((np.concatenate([cm.x for cm in complementary_msgs]), parent.msg_to_parent.x)))
            #msg_parent_to_node = self._multiply_dists(complementary_msgs, parent.msg_to_parent.x)
            msg_parent_to_node = self._multiply_dists(complementary_msgs, parent.msg_to_parent.x)

            # propagate the message from the parent to the node:
            node.msg_from_parent = _send_message(msg_parent_to_node, node.branch_neg_log_prob, None, node.branch_len_sigma)

            # finally, set the node marginal LH distribution:

            if node.msg_to_parent is not None:

                grid = np.unique(np.concatenate((node.msg_to_parent.x, node.msg_from_parent.x)))
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


if __name__=="__main__":
    import matplotlib.pyplot as plt
    plt.ion()

    with open('data/H3N2_NA_allyears_NA.20.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    from Bio import Phylo
    tree = Phylo.read("data/H3N2_NA_allyears_NA.20.nwk", 'newick')
    tree.root_with_outgroup([n for n in tree.get_terminals()
                              if n.name=='A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409'][0])
    myTree = ClockTree(gtr='Jukes-Cantor', tree = tree,
                        aln = 'data/H3N2_NA_allyears_NA.20.fasta', verbose = 6, dates = dates)

    myTree.init_date_constraints()
    myTree._ml_t_leaves_root()

    plt.figure()
    x = np.linspace(0,0.05,100)
    for node in myTree.tree.find_clades():
        if node.up is not None:
            print(node.branch_length_interpolator.peak_val, node.mutations)
            plt.plot(x, node.branch_length_interpolator.prob(x))
    plt.yscale('log')

    plt.figure()
    x = np.linspace(0,0.2,100)
    for node in myTree.tree.find_clades():
        if node.up is not None:
            print(node.branch_length_interpolator.peak_val, node.mutations)
            plt.plot(x, node.msg_to_parent.prob_relative(x))
    #plt.yscale('log')

