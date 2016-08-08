
class ClockTree(TreeAnc):
    """
    Class to produce molecular clock trees.
    """

    def __init__(self, *args,**kwargs):
        super(TreeTime, self).__init__(*args, **kwargs)
        self.date2dist = None  # we do not know anything about the conversion
        self.max_diam = 0.0
        self.debug=False

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
        self.logger("TreeTime.init_date_constraints...",2)
        self.date2dist = utils.DateConversion.from_tree(self.tree, slope)
        self.max_diam = self.date2dist.intercept

        # set the None  for the date-related attributes in the internal nodes.
        # make interpolation objects for the branches
        self._ml_t_init(**kwarks)

    def _make_branch_len_interpolator(self, node, n=ttconf.BRANCH_GRID_SIZE, ignore_gaps=False):
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

        if not hasattr(node, 'merger_rate') or node.merger_rate is None:
            node.merger_rate = ttconf.BRANCH_LEN_PENALTY

        # optimal branch length
        mutation_length = node.mutation_length

        if mutation_length < np.min((1e-5, 0.1*self.one_mutation)): # zero-length
            grid = ttconf.MAX_BRANCH_LENGTH * (np.linspace(0, 1.0 , n)**2)

        else: # branch length is not zero

            sigma = mutation_length #np.max([self.average_branch_len, mutation_length])
            # from zero to optimal branch length
            grid_left = mutation_length * (1 - np.linspace(1, 0.0, n/3)**2.0)
            # from optimal branch length to the right (--> 3*branch lengths),
            grid_right = mutation_length + 1e-3 * self.one_mutation + (3*sigma*(np.linspace(0, 1, n/3)**2))
            # far to the right (3*branch length ---> MAX_LEN), very sparse
            far_grid = grid_right.max() + mutation_length/2 + ttconf.MAX_BRANCH_LENGTH*np.linspace(0, 1, n/3)**2

            grid = np.concatenate((grid_left,grid_right,far_grid))
            grid.sort() # just for safety


        if hasattr(node, 'compressed_sequence'):
            log_prob = np.array([-self.gtr.prob_t_compressed(node.compressed_sequence['pair'],
                                                node.compressed_sequence['multiplicity'],
                                                k,
                                                return_log=True)
                    for k in grid])
        else:
            log_prob = np.array([-self.gtr.prob_t(node.sequence,
                                     node.up.sequence,
                                     k,
                                     return_log=True,
                                     ignore_gaps=self.ignore_gaps)
                    for k in grid])

        log_prob [np.isnan(log_prob)] = ttconf.BIG_NUMBER
        tmp_log_prob = np.exp(-log_prob)
        integral = self._integral(grid, tmp_log_prob)

        # save raw branch length interpolators without coalescent contribution
        # TODO: better criterion to catch bad branch
        if integral < 1e-200:
            print ("!!WARNING!!", node.name, " branch length probability distribution",
                "integral is ZERO. Setting bad_branch flag...")
            node.bad_branch = True
            node.raw_branch_neg_log_prob = Distribution(grid, log_prob, is_log=True, kind='linear')

        else:
            node.raw_branch_neg_log_prob = Distribution(grid, log_prob+np.log(integral), is_log=True, kind='linear')

        # add merger rate contribution to the raw branch length
        log_prob += node.merger_rate * np.minimum(ttconf.MAX_BRANCH_LENGTH, np.maximum(0,grid))

        tmp_log_prob = np.exp(-log_prob)
        integral = self._integral(grid, tmp_log_prob)

        if integral < 1e-200:
            print ("!!WARNING!! Node branch length probability distribution "
                "integral is ZERO. Setting bad_branch flag...")
            node.bad_branch = True
            node.branch_neg_log_prob = Distribution(grid, log_prob, is_log=True, kind='linear')

        else:
            node.branch_neg_log_prob = Distribution(grid, log_prob+np.log(integral), is_log=True, kind='linear')

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
           node.branch_neg_log_prob = Distribution(grid, y + np.log(integral), is_log=True, kind='linear')
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

        print('\n----- Initializing branch length interpolation objects...\n')
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
                    node.msg_to_parent = Distribution.delta_function(node.abs_t, normalized=False)

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

            if node.msg_to_parent.is_delta:
                res = Distribution.shifted_x(node.branch_neg_log_prob, node.msg_to_parent.peak_pos)

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

        tau = g_func.x
        for t_idx, t_val in enumerate(time):
            F = Distribution(t_val - tau, f_func.y)
            FG = F * g_func
            res[t_idx] = FG.integrate(tau.min(), tau.max(), n_integral)

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

            dtau = np.min((self.one_mutation, (tau_max - tau_min) / n_integral))

            tau = np.arange(tau_min, tau_max, dtau)

            if inverse_time:
                fg = f_func(ti-tau) + g_func(tau) # f,g - log probabilities
            else:
                fg = f_func(ti+tau) + g_func(tau)

            min_fg = fg.min() # exponent pre-factor
            expfg = np.exp(-1*(fg-min_fg))

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
