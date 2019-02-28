from __future__ import print_function, division, absolute_import
import numpy as np
from treetime import config as ttconf
from .treeanc import TreeAnc
from .utils import numeric_date, DateConversion
from .distribution import Distribution
from .branch_len_interpolator import BranchLenInterpolator
from .node_interpolator import NodeInterpolator

class ClockTree(TreeAnc):
    """
    ClockTree is the main class to perform the optimization of the node
    positions given the temporal constraints of (some) leaves.

    The optimization workflow includes the inference of the ancestral sequences
    and branch length optimization using TreeAnc. After the optimization
    is done, the nodes with date-time information are arranged along the time axis,
    the conversion between the branch lengths units and the date-time units
    is determined. Then, for each internal node, we compute the the probability distribution
    of the node's location conditional on the fixed location of the leaves, which
    have temporal information. In the end, the most probable location of the internal nodes
    is converted to the most likely time of the internal nodes.
    """

    def __init__(self,  dates=None, debug=False, real_dates=True, precision='auto',
                 branch_length_mode='joint', use_covariation=False, *args, **kwargs):

        """
        ClockTree constructor

        Parameters
        ----------

         dates : dict
            :code:`{leaf_name:leaf_date}` dictionary

         debug : bool
            If True, the debug mode is ON, which means no or less clean-up of
            obsolete parameters to control program execution in intermediate
            states. In debug mode, the python debugger is also allowed to interrupt
            program execution with intercative shell if an error occurs.

         real_dates : bool
            If True, some additional checks for the input dates sanity will be
            performed.

         precision : int
            Precision can be 0 (rough), 1 (default), 2 (fine), or 3 (ultra fine).
            This parameter determines the number of grid points that are used
            for the evaluation of the branch length interpolation objects.
            When not specified, this will default to 1 for short sequences and 2
            for long sequences with L>1e4

         branch_length_mode : str
            determines whether branch length are calculated using the 'joint' ML,
            'marginal' ML, or branch length of the input tree ('input').

         use_covariation : bool
            determines whether root-to-tip regression accounts for covariance
            introduced by shared ancestry.

         **kwargs:
            Key word argments needed to construct parent class (TreeAnc)

        """
        super(ClockTree, self).__init__(*args, **kwargs)
        if dates is None:
            raise ValueError("ClockTree requires date constraints!")

        self.debug=debug
        self.real_dates = real_dates
        self.date_dict = dates
        self._date2dist = None  # we do not know anything about the conversion
        self.tip_slack = ttconf.OVER_DISPERSION  # extra number of mutations added
                                                 # to terminal branches in covariance calculation
        self.rel_tol_prune = ttconf.REL_TOL_PRUNE
        self.rel_tol_refine = ttconf.REL_TOL_REFINE
        self.branch_length_mode = branch_length_mode
        self.clock_model=None
        self.use_covariation=use_covariation # if false, covariation will be ignored in rate estimates.
        self._set_precision(precision)
        if self._assign_dates()==ttconf.ERROR:
            raise ValueError("ClockTree requires date constraints!")


    def _assign_dates(self):
        """assign dates to nodes

        Returns
        -------
        str
            success/error code
        """
        if self.tree is None:
            self.logger("ClockTree._assign_dates: tree is not set, can't assign dates", 0)
            return ttconf.ERROR

        bad_branch_counter = 0
        for node in self.tree.find_clades(order='postorder'):
            if node.name in self.date_dict:
                tmp_date = self.date_dict[node.name]
                if np.isscalar(tmp_date) and np.isnan(tmp_date):
                    self.logger("WARNING: ClockTree.init: node %s has a bad date: %s"%(node.name, str(tmp_date)), 2, warn=True)
                    node.raw_date_constraint = None
                    node.bad_branch = True
                else:
                    try:
                        tmp = np.mean(tmp_date)
                        node.raw_date_constraint = tmp_date
                        node.bad_branch = False
                    except:
                        self.logger("WARNING: ClockTree.init: node %s has a bad date: %s"%(node.name, str(tmp_date)), 2, warn=True)
                        node.raw_date_constraint = None
                        node.bad_branch = True
            else: # nodes without date contraints

                node.raw_date_constraint = None

                if node.is_terminal():
                    # Terminal branches without date constraints marked as 'bad'
                    node.bad_branch = True
                else:
                    # If all branches dowstream are 'bad', and there is no date constraint for
                    # this node, the branch is marked as 'bad'
                    node.bad_branch = np.all([x.bad_branch for x in node])

            if node.is_terminal() and node.bad_branch:
                bad_branch_counter += 1

        if bad_branch_counter>self.tree.count_terminals()-3:
            self.logger("ERROR: ALMOST NO VALID DATE CONSTRAINTS, EXITING", 1, warn=True)
            return ttconf.ERROR

        return ttconf.SUCCESS


    def _set_precision(self, precision):
        '''
        function that sets precision to an (hopfully) reasonable guess based
        on the length of the sequence if not explicitly set
        '''
        # if precision is explicitly specified, use it.

        if self.one_mutation:
            self.min_width = 10*self.one_mutation
        else:
            self.min_width = 0.001
        if precision in [0,1,2,3]:
            self.precision=precision
            if self.one_mutation and self.one_mutation<1e-4 and precision<2:
                self.logger("ClockTree._set_precision: FOR LONG SEQUENCES (>1e4) precision>=2 IS RECOMMENDED."
                            " \n\t **** precision %d was specified by the user"%precision, level=0)
        else:
            # otherwise adjust it depending on the minimal sensible branch length
            if self.one_mutation:
                if self.one_mutation>1e-4:
                    self.precision=1
                else:
                    self.precision=2
            else:
                self.precision=1
            self.logger("ClockTree: Setting precision to level %s"%self.precision, 2)

        if self.precision==0:
            self.node_grid_points = ttconf.NODE_GRID_SIZE_ROUGH
            self.branch_grid_points = ttconf.BRANCH_GRID_SIZE_ROUGH
            self.n_integral = ttconf.N_INTEGRAL_ROUGH
        elif self.precision==2:
            self.node_grid_points = ttconf.NODE_GRID_SIZE_FINE
            self.branch_grid_points = ttconf.BRANCH_GRID_SIZE_FINE
            self.n_integral = ttconf.N_INTEGRAL_FINE
        elif self.precision==3:
            self.node_grid_points = ttconf.NODE_GRID_SIZE_ULTRA
            self.branch_grid_points = ttconf.BRANCH_GRID_SIZE_ULTRA
            self.n_integral = ttconf.N_INTEGRAL_ULTRA
        else:
            self.node_grid_points = ttconf.NODE_GRID_SIZE
            self.branch_grid_points = ttconf.BRANCH_GRID_SIZE
            self.n_integral = ttconf.N_INTEGRAL

    @property
    def date2dist(self):
        return self._date2dist

    @date2dist.setter
    def date2dist(self, val):
        if val is None:
            self._date2dist = None
        else:
            self.logger("ClockTree.date2dist: Setting new molecular clock."
                        " rate=%.3e, R^2=%.4f"%(val.clock_rate, val.r_val**2), 2)
            self._date2dist = val


    def setup_TreeRegression(self, covariation=True):
        """instantiate a TreeRegression object and set its tip_value and branch_value function
        to defaults that are sensible for treetime instances.

        Parameters
        ----------
        covariation : bool, optional
            account for phylogenetic covariation
        Returns
        -------
        TreeRegression
            a TreeRegression instance with self.tree attached as tree.
        """
        from .treeregression import TreeRegression
        tip_value = lambda x:np.mean(x.raw_date_constraint) if (x.is_terminal() and (x.bad_branch is False)) else None
        branch_value = lambda x:x.mutation_length
        if covariation:
            om = self.one_mutation
            branch_variance = lambda x:((x.clock_length if hasattr(x,'clock_length') else x.mutation_length)
                                        +(self.tip_slack**2*om if x.is_terminal() else 0.0))*om
        else:
            branch_variance = lambda x:1.0 if x.is_terminal() else 0.0

        Treg = TreeRegression(self.tree, tip_value=tip_value,
                             branch_value=branch_value, branch_variance=branch_variance)
        Treg.valid_confidence = covariation
        return Treg


    def get_clock_model(self, covariation=True, slope=None):
        Treg = self.setup_TreeRegression(covariation=covariation)
        self.clock_model = Treg.regression(slope=slope)
        if not np.isfinite(self.clock_model['slope']):
            raise ValueError("Clock rate estimation failed. If your data lacks temporal signal, please specify the rate explicitly!")

        if not Treg.valid_confidence or (slope is not None):
            if 'cov' in self.clock_model:
                self.clock_model.pop('cov')
            self.clock_model['valid_confidence']=False
        else:
            self.clock_model['valid_confidence']=True
        self.clock_model['r_val'] = Treg.explained_variance()
        self.date2dist = DateConversion.from_regression(self.clock_model)


    def init_date_constraints(self, ancestral_inference=False, clock_rate=None, **kwarks):
        """
        Get the conversion coefficients between the dates and the branch
        lengths as they are used in ML computations. The conversion formula is
        assumed to be 'length = k*numdate + b'. For convenience, these
        coefficients as well as regression parameters are stored in the
        'dates2dist' object.

        .. Note::
            The tree must have dates set to all nodes before calling this
            function.

        Parameters
        ----------

         ancestral_inference: bool
            If True, reinfer ancestral sequences
            when ancestral sequences are missing

         clock_rate: float
            If specified, timetree optimization will be done assuming a
            fixed clock rate as specified

        """
        self.logger("ClockTree.init_date_constraints...",2)
        self.tree.coalescent_joint_LH = 0
        if self.aln and (ancestral_inference or (not hasattr(self.tree.root, 'sequence'))):
            self.infer_ancestral_sequences('probabilistic', marginal=self.branch_length_mode=='marginal',
                                            sample_from_profile='root',**kwarks)

        # set the None  for the date-related attributes in the internal nodes.
        # make interpolation objects for the branches
        self.logger('ClockTree.init_date_constraints: Initializing branch length interpolation objects...',3)
        has_clock_length = []
        for node in self.tree.find_clades(order='postorder'):
            if node.up is None:
                node.branch_length_interpolator = None
            else:
                has_clock_length.append(hasattr(node, 'clock_length'))
                # copy the merger rate and gamma if they are set
                if hasattr(node,'branch_length_interpolator') and node.branch_length_interpolator is not None:
                    gamma = node.branch_length_interpolator.gamma
                    merger_cost = node.branch_length_interpolator.merger_cost
                else:
                    gamma = 1.0
                    merger_cost = None

                if self.branch_length_mode=='marginal':
                    node.profile_pair = self.marginal_branch_profile(node)

                node.branch_length_interpolator = BranchLenInterpolator(node, self.gtr,
                            pattern_multiplicity = self.multiplicity, min_width=self.min_width,
                            one_mutation=self.one_mutation, branch_length_mode=self.branch_length_mode)

                node.branch_length_interpolator.merger_cost = merger_cost
                node.branch_length_interpolator.gamma = gamma

        # use covariance in clock model only after initial timetree estimation is done
        use_cov = (np.sum(has_clock_length) > len(has_clock_length)*0.7) and self.use_covariation
        self.get_clock_model(covariation=use_cov, slope=clock_rate)

        # make node distribution objects
        for node in self.tree.find_clades(order="postorder"):
            # node is constrained
            if hasattr(node, 'raw_date_constraint') and node.raw_date_constraint is not None:
                # set the absolute time before present in branch length units
                if np.isscalar(node.raw_date_constraint):
                    tbp = self.date2dist.get_time_before_present(node.raw_date_constraint)
                    node.date_constraint = Distribution.delta_function(tbp, weight=1.0, min_width=self.min_width)
                else:
                    tbp = self.date2dist.get_time_before_present(np.array(node.raw_date_constraint))
                    node.date_constraint = Distribution(tbp, np.ones_like(tbp), is_log=False, min_width=self.min_width)

                if hasattr(node, 'bad_branch') and node.bad_branch is True:
                    self.logger("ClockTree.init_date_constraints -- WARNING: Branch is marked as bad"
                                ", excluding it from the optimization process.\n"
                                "\t\tDate constraint will be ignored!", 4, warn=True)
            else: # node without sampling date set
                node.raw_date_constraint = None
                node.date_constraint = None


    def make_time_tree(self, time_marginal=False, clock_rate=None, **kwargs):
        '''
        Use the date constraints to calculate the most likely positions of
        unconstrained nodes.

        Parameters
        ----------

         time_marginal : bool
            If true, use marginal reconstruction for node positions

         **kwargs
            Key word arguments to initialize dates constraints

        '''
        self.logger("ClockTree: Maximum likelihood tree optimization with temporal constraints",1)

        self.init_date_constraints(clock_rate=clock_rate, **kwargs)

        if time_marginal:
            self._ml_t_marginal(assign_dates = time_marginal=="assign")
        else:
            self._ml_t_joint()

        self.convert_dates()


    def _ml_t_joint(self):
        """
        Compute the joint maximum likelihood assignment of the internal nodes positions by
        propagating from the tree leaves towards the root. Given the assignment of parent nodes,
        reconstruct the maximum-likelihood positions of the child nodes by propagating
        from the root to the leaves. The result of this operation is the time_before_present
        value, which is the position of the node, expressed in the units of the
        branch length, and scaled from the present-day. The value is assigned to the
        corresponding attribute of each node of the tree.

        Returns
        -------

         None
            Every internal node is assigned the probability distribution in form
            of an interpolation object and sends this distribution further towards the
            root.

        """

        def _cleanup():
            for node in self.tree.find_clades():
                del node.joint_pos_Lx
                del node.joint_pos_Cx


        self.logger("ClockTree - Joint reconstruction:  Propagating leaves -> root...", 2)
        # go through the nodes from leaves towards the root:
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents
            # Lx is the maximal likelihood of a subtree given the parent position
            # Cx is the branch length corresponding to the maximally likely subtree
            if node.bad_branch:
                # no information at the node
                node.joint_pos_Lx = None
                node.joint_pos_Cx = None
            else: # all other nodes
                if node.date_constraint is not None and node.date_constraint.is_delta: # there is a time constraint
                    # subtree probability given the position of the parent node
                    # Lx.x is the position of the parent node
                    # Lx.y is the probablity of the subtree (consisting of one terminal node in this case)
                    # Cx.y is the branch length corresponding the optimal subtree
                    bl = node.branch_length_interpolator.x
                    x = bl + node.date_constraint.peak_pos
                    node.joint_pos_Lx = Distribution(x, node.branch_length_interpolator(bl),
                                                     min_width=self.min_width, is_log=True)
                    node.joint_pos_Cx = Distribution(x, bl, min_width=self.min_width) # map back to the branch length
                else: # all nodes without precise constraint but positional information
                    msgs_to_multiply = [node.date_constraint] if node.date_constraint is not None else []
                    msgs_to_multiply.extend([child.joint_pos_Lx for child in node.clades
                                             if child.joint_pos_Lx is not None])

                    # subtree likelihood given the node's constraint and child messages
                    if len(msgs_to_multiply) == 0: # there are no constraints
                        node.joint_pos_Lx = None
                        node.joint_pos_Cx = None
                        continue
                    elif len(msgs_to_multiply)>1: # combine the different msgs and constraints
                        subtree_distribution = Distribution.multiply(msgs_to_multiply)
                    else: # there is exactly one constraint.
                        subtree_distribution = msgs_to_multiply[0]

                    if node.up is None: # this is the root, set dates
                        subtree_distribution._adjust_grid(rel_tol=self.rel_tol_prune)
                        # set root position and joint likelihood of the tree
                        node.time_before_present = subtree_distribution.peak_pos
                        node.joint_pos_Lx = subtree_distribution
                        node.joint_pos_Cx = None
                        node.clock_length = node.branch_length
                    else: # otherwise propagate to parent
                        res, res_t = NodeInterpolator.convolve(subtree_distribution,
                                        node.branch_length_interpolator,
                                        max_or_integral='max',
                                        inverse_time=True,
                                        n_grid_points = self.node_grid_points,
                                        n_integral=self.n_integral,
                                        rel_tol=self.rel_tol_refine)

                        res._adjust_grid(rel_tol=self.rel_tol_prune)

                        node.joint_pos_Lx = res
                        node.joint_pos_Cx = res_t


        # go through the nodes from root towards the leaves and assign joint ML positions:
        self.logger("ClockTree - Joint reconstruction:  Propagating root -> leaves...", 2)
        for node in self.tree.find_clades(order='preorder'):  # root first, msgs to children

            if node.up is None: # root node
                continue # the position was already set on the previous step

            if node.joint_pos_Cx is None: # no constraints or branch is bad - reconstruct from the branch len interpolator
                node.branch_length = node.branch_length_interpolator.peak_pos

            elif isinstance(node.joint_pos_Cx, Distribution):
                # NOTE the Lx distribution is the likelihood, given the position of the parent
                # (Lx.x = parent position, Lx.y = LH of the node_pos given Lx.x,
                # the length of the branch corresponding to the most likely
                # subtree is node.Cx(node.time_before_present))
                subtree_LH = node.joint_pos_Lx(node.up.time_before_present)
                node.branch_length = node.joint_pos_Cx(max(node.joint_pos_Cx.xmin,
                                            node.up.time_before_present)+ttconf.TINY_NUMBER)

            node.time_before_present = node.up.time_before_present - node.branch_length
            node.clock_length = node.branch_length

            # just sanity check, should never happen:
            if node.branch_length < 0 or node.time_before_present < 0:
                if node.branch_length<0 and node.branch_length>-ttconf.TINY_NUMBER:
                    self.logger("ClockTree - Joint reconstruction: correcting rounding error of %s"%node.name, 4)
                    node.branch_length = 0

        self.tree.positional_joint_LH = self.timetree_likelihood()
        # cleanup, if required
        if not self.debug:
            _cleanup()


    def timetree_likelihood(self):
        '''
        Return the likelihood of the data given the current branch length in the tree
        '''
        LH = 0
        for node in self.tree.find_clades(order='preorder'):  # sum the likelihood contributions of all branches
            if node.up is None: # root node
                continue
            LH -= node.branch_length_interpolator(node.branch_length)

        # add the root sequence LH and return
        if self.aln:
            LH += self.gtr.sequence_logLH(self.tree.root.cseq, pattern_multiplicity=self.multiplicity)
        return LH


    def _ml_t_marginal(self, assign_dates=False):
        """
        Compute the marginal probability distribution of the internal nodes positions by
        propagating from the tree leaves towards the root. The result of
        this operation are the probability distributions of each internal node,
        conditional on the constraints on all leaves of the tree, which have sampling dates.
        The probability distributions are set as marginal_pos_LH attributes to the nodes.

        Parameters
        ----------

         assign_dates : bool, default False
            If True, the inferred dates will be assigned to the nodes as
            :code:`time_before_present' attributes, and their branch lengths
            will be corrected accordingly.
            .. Note::
                Normally, the dates are assigned by running joint reconstruction.

        Returns
        -------

         None
            Every internal node is assigned the probability distribution in form
            of an interpolation object and sends this distribution further towards the
            root.

        """

        def _cleanup():
            for node in self.tree.find_clades():
                try:
                    del node.marginal_pos_Lx
                    del node.subtree_distribution
                    del node.msg_from_parent
                    #del node.marginal_pos_LH
                except:
                    pass


        self.logger("ClockTree - Marginal reconstruction:  Propagating leaves -> root...", 2)
        # go through the nodes from leaves towards the root:
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents
            if node.bad_branch:
                # no information
                node.marginal_pos_Lx = None
            else: # all other nodes
                if node.date_constraint is not None and node.date_constraint.is_delta: # there is a time constraint
                    # initialize the Lx for nodes with precise date constraint:
                    # subtree probability given the position of the parent node
                    # position of the parent node is given by the branch length
                    # distribution attached to the child node position
                    node.subtree_distribution = node.date_constraint
                    bl = node.branch_length_interpolator.x
                    x = bl + node.date_constraint.peak_pos
                    node.marginal_pos_Lx = Distribution(x, node.branch_length_interpolator(bl),
                                                        min_width=self.min_width, is_log=True)

                else: # all nodes without precise constraint but positional information
                      # subtree likelihood given the node's constraint and child msg:
                    msgs_to_multiply = [node.date_constraint] if node.date_constraint is not None else []
                    msgs_to_multiply.extend([child.marginal_pos_Lx for child in node.clades
                                             if child.marginal_pos_Lx is not None])

                    # combine the different msgs and constraints
                    if len(msgs_to_multiply)==0:
                        # no information
                        node.marginal_pos_Lx = None
                        continue
                    elif len(msgs_to_multiply)==1:
                        node.subtree_distribution = msgs_to_multiply[0]
                    else: # combine the different msgs and constraints
                        node.subtree_distribution = Distribution.multiply(msgs_to_multiply)

                    if node.up is None: # this is the root, set dates
                        node.subtree_distribution._adjust_grid(rel_tol=self.rel_tol_prune)
                        node.marginal_pos_Lx = node.subtree_distribution
                        node.marginal_pos_LH = node.subtree_distribution
                        self.tree.positional_marginal_LH = -node.subtree_distribution.peak_val
                    else: # otherwise propagate to parent
                        res, res_t = NodeInterpolator.convolve(node.subtree_distribution,
                                        node.branch_length_interpolator,
                                        max_or_integral='integral',
                                        n_grid_points = self.node_grid_points,
                                        n_integral=self.n_integral,
                                        rel_tol=self.rel_tol_refine)
                        res._adjust_grid(rel_tol=self.rel_tol_prune)
                        node.marginal_pos_Lx = res

        self.logger("ClockTree - Marginal reconstruction:  Propagating root -> leaves...", 2)
        from scipy.interpolate import interp1d
        for node in self.tree.find_clades(order='preorder'):

            ## The root node
            if node.up is None:
                node.msg_from_parent = None # nothing beyond the root
            # all other cases (All internal nodes + unconstrained terminals)
            else:
                parent = node.up
                # messages from the complementary subtree (iterate over all sister nodes)
                complementary_msgs = [sister.marginal_pos_Lx for sister in parent.clades
                                            if (sister != node) and (sister.marginal_pos_Lx is not None)]

                # if parent itself got smth from the root node, include it
                if parent.msg_from_parent is not None:
                    complementary_msgs.append(parent.msg_from_parent)
                elif parent.marginal_pos_Lx is not None:
                    complementary_msgs.append(parent.marginal_pos_LH)

                if len(complementary_msgs):
                    msg_parent_to_node = NodeInterpolator.multiply(complementary_msgs)
                    msg_parent_to_node._adjust_grid(rel_tol=self.rel_tol_prune)
                else:
                    x = [parent.numdate, numeric_date()]
                    msg_parent_to_node = NodeInterpolator(x, [1.0, 1.0],min_width=self.min_width)

                # integral message, which delivers to the node the positional information
                # from the complementary subtree
                res, res_t = NodeInterpolator.convolve(msg_parent_to_node, node.branch_length_interpolator,
                                                    max_or_integral='integral',
                                                    inverse_time=False,
                                                    n_grid_points = self.node_grid_points,
                                                    n_integral=self.n_integral,
                                                    rel_tol=self.rel_tol_refine)

                node.msg_from_parent = res
                if node.marginal_pos_Lx is None:
                    node.marginal_pos_LH = node.msg_from_parent
                else:
                    node.marginal_pos_LH = NodeInterpolator.multiply((node.msg_from_parent, node.subtree_distribution))

                self.logger('ClockTree._ml_t_root_to_leaves: computed convolution'
                                ' with %d points at node %s'%(len(res.x),node.name),4)

                if self.debug:
                    tmp = np.diff(res.y-res.peak_val)
                    nsign_changed = np.sum((tmp[1:]*tmp[:-1]<0)&(res.y[1:-1]-res.peak_val<500))
                    if nsign_changed>1:
                        import matplotlib.pyplot as plt
                        plt.ion()
                        plt.plot(res.x, res.y-res.peak_val, '-o')
                        plt.plot(res.peak_pos - node.branch_length_interpolator.x,
                                 node.branch_length_interpolator(node.branch_length_interpolator.x)-node.branch_length_interpolator.peak_val, '-o')
                        plt.plot(msg_parent_to_node.x,msg_parent_to_node.y-msg_parent_to_node.peak_val, '-o')
                        plt.ylim(0,100)
                        plt.xlim(-0.05, 0.05)
                        import ipdb; ipdb.set_trace()

            # assign positions of nodes and branch length only when desired
            # since marginal reconstruction can result in negative branch length
            if assign_dates:
                node.time_before_present = node.marginal_pos_LH.peak_pos
                if node.up:
                    node.clock_length = node.up.time_before_present - node.time_before_present
                    node.branch_length = node.clock_length

            # construct the inverse cumulant distribution to evaluate confidence intervals
            if node.marginal_pos_LH.is_delta:
                node.marginal_inverse_cdf=interp1d([0,1], node.marginal_pos_LH.peak_pos*np.ones(2), kind="linear")
            else:
                dt = np.diff(node.marginal_pos_LH.x)
                y = node.marginal_pos_LH.prob_relative(node.marginal_pos_LH.x)
                int_y = np.concatenate(([0], np.cumsum(dt*(y[1:]+y[:-1])/2.0)))
                int_y/=int_y[-1]
                node.marginal_inverse_cdf = interp1d(int_y, node.marginal_pos_LH.x, kind="linear")
                node.marginal_cdf = interp1d(node.marginal_pos_LH.x, int_y, kind="linear")

        if not self.debug:
            _cleanup()

        return


    def convert_dates(self):
        '''
        This function converts the estimated "time_before_present" properties of all nodes
        to numerical dates stored in the "numdate" attribute. This date is further converted
        into a human readable date string in format %Y-%m-%d assuming the usual calendar.

        Returns
        -------
         None
            All manipulations are done in place on the tree

        '''
        from datetime import datetime, timedelta
        now = numeric_date()
        for node in self.tree.find_clades():
            years_bp = self.date2dist.to_years(node.time_before_present)
            if years_bp < 0 and self.real_dates:
                if not hasattr(node, "bad_branch") or node.bad_branch is False:
                    self.logger("ClockTree.convert_dates -- WARNING: The node is later than today, but it is not "
                        "marked as \"BAD\", which indicates the error in the "
                        "likelihood optimization.",4 , warn=True)
                else:
                    self.logger("ClockTree.convert_dates -- WARNING: node which is marked as \"BAD\" optimized "
                        "later than present day",4 , warn=True)

            node.numdate = now - years_bp

            # set the human-readable date
            year = np.floor(node.numdate)
            days = max(0,365.25 * (node.numdate - year)-1)
            try:  # datetime will only operate on dates after 1900
                n_date = datetime(year, 1, 1) + timedelta(days=days)
                node.date = datetime.strftime(n_date, "%Y-%m-%d")
            except:
                # this is the approximation not accounting for gap years etc
                n_date = datetime(1900, 1, 1) + timedelta(days=days)
                node.date = "%04d-%02d-%02d"%(year, n_date.month, n_date.day)


    def branch_length_to_years(self):
        '''
        This function sets branch length to reflect the date differences between parent and child
        nodes measured in years. Should only be called after :py:meth:`timetree.ClockTree.convert_dates` has been called.

        Returns
        -------
         None
            All manipulations are done in place on the tree

        '''
        self.logger('ClockTree.branch_length_to_years: setting node positions in units of years', 2)
        if not hasattr(self.tree.root, 'numdate'):
            self.logger('ClockTree.branch_length_to_years: infer ClockTree first', 2,warn=True)
        self.tree.root.branch_length = 0.1
        for n in self.tree.find_clades(order='preorder'):
            if n.up is not None:
                n.branch_length = n.numdate - n.up.numdate


    def calc_rate_susceptibility(self, rate_std=None, params=None):
        """return the time tree estimation of evolutionary rates +/- one
        standard deviation form the ML estimate.

        Returns
        -------
        TreeTime.return_code : str
            success or failure
        """
        params = params or {}
        if rate_std is None:
            if not (self.clock_model['valid_confidence'] and 'cov' in self.clock_model):
                self.logger("ClockTree.calc_rate_susceptibility: need valid standard deviation of the clock rate to estimate dating error.", 1, warn=True)
                return ttconf.ERROR
            rate_std = np.sqrt(self.clock_model['cov'][0,0])

        current_rate = np.abs(self.clock_model['slope'])
        upper_rate = self.clock_model['slope'] + rate_std
        lower_rate = max(0.1*current_rate, self.clock_model['slope'] - rate_std)
        for n in self.tree.find_clades():
            if n.up:
                n._orig_gamma = n.branch_length_interpolator.gamma
                n.branch_length_interpolator.gamma = n._orig_gamma*upper_rate/current_rate

        self.logger("###ClockTree.calc_rate_susceptibility: run with upper bound of rate estimate", 1)
        self.make_time_tree(**params)
        self.logger("###ClockTree.calc_rate_susceptibility: rate: %f, LH:%f"%(upper_rate, self.tree.positional_joint_LH), 2)
        for n in self.tree.find_clades():
            n.numdate_rate_variation = [(upper_rate, n.numdate)]
            if n.up:
                n.branch_length_interpolator.gamma = n._orig_gamma*lower_rate/current_rate

        self.logger("###ClockTree.calc_rate_susceptibility: run with lower bound of rate estimate", 1)
        self.make_time_tree(**params)
        self.logger("###ClockTree.calc_rate_susceptibility: rate: %f, LH:%f"%(lower_rate, self.tree.positional_joint_LH), 2)
        for n in self.tree.find_clades():
            n.numdate_rate_variation.append((lower_rate, n.numdate))
            if n.up:
                n.branch_length_interpolator.gamma  = n._orig_gamma

        self.logger("###ClockTree.calc_rate_susceptibility: run with central rate estimate", 1)
        self.make_time_tree(**params)
        self.logger("###ClockTree.calc_rate_susceptibility: rate: %f, LH:%f"%(current_rate, self.tree.positional_joint_LH), 2)
        for n in self.tree.find_clades():
            n.numdate_rate_variation.append((current_rate, n.numdate))
            n.numdate_rate_variation.sort(key=lambda x:x[1]) # sort estimates for different rates by numdate

        return ttconf.SUCCESS


    def date_uncertainty_due_to_rate(self, node, interval=(0.05, 0.095)):
        """use previously calculated variation of the rate to estimate
        the uncertainty in a particular numdate due to rate variation.

        Parameters
        ----------
        node : PhyloTree.Clade
            node for which the confidence interval is to be calculated
        interval : tuple, optional
            Array of length two, or tuple, defining the bounds of the confidence interval

        """
        if hasattr(node, "numdate_rate_variation"):
            from scipy.special import erfinv
            nsig = [np.sqrt(2.0)*erfinv(-1.0 + 2.0*x) if x*(1.0-x) else 0
                    for x in interval]
            l,c,u = [x[1] for x in node.numdate_rate_variation]
            return np.array([c + x*np.abs(y-c) for x,y in zip(nsig, (l,u))])

        else:
            return None

    def combine_confidence(self, center, limits, c1=None, c2=None):
        if c1 is None and c2 is None:
            return np.array(limits)
        elif c1 is None:
            min_val,max_val = c2
        elif c2 is None:
            min_val,max_val = c1
        else:
            min_val = center - np.sqrt((c1[0]-center)**2 + (c2[0]-center)**2)
            max_val = center + np.sqrt((c1[1]-center)**2 + (c2[1]-center)**2)

        return np.array([max(limits[0], min_val),
                         min(limits[1], max_val)])



    def get_confidence_interval(self, node, interval = (0.05, 0.95)):
        '''
        If temporal reconstruction was done using the marginal ML mode, the entire distribution of
        times is available. This function determines the 90% (or other) confidence interval, defined as the
        range where 5% of probability is below and above. Note that this does not necessarily contain
        the highest probability position.
        In absense of marginal reconstruction, it will return uncertainty based on rate
        variation. If both are present, the wider interval will be returned.

        Parameters
        ----------

         node : PhyloTree.Clade
            The node for which the confidence interval is to be calculated

         interval : tuple, list
            Array of length two, or tuple, defining the bounds of the confidence interval

        Returns
        -------

         confidence_interval : numpy array
            Array with two numerical dates delineating the confidence interval

        '''
        rate_contribution = self.date_uncertainty_due_to_rate(node, interval)

        if hasattr(node, "marginal_inverse_cdf"):
            min_date, max_date = [self.date2dist.to_numdate(x) for x in
                                  (node.marginal_pos_LH.xmax, node.marginal_pos_LH.xmin)]
            if node.marginal_inverse_cdf=="delta":
                return np.array([node.numdate, node.numdate])
            else:
                mutation_contribution = self.date2dist.to_numdate(node.marginal_inverse_cdf(np.array(interval))[::-1])
        else:
            min_date, max_date = [-np.inf, np.inf]

        return self.combine_confidence(node.numdate, (min_date, max_date),
                                  c1=rate_contribution, c2=mutation_contribution)

    def get_max_posterior_region(self, node, fraction = 0.9):
        '''
        If temporal reconstruction was done using the marginal ML mode, the entire distribution of
        times is available. This function determines the interval around the highest
        posterior probability region that contains the specified fraction of the probability mass.
        In absense of marginal reconstruction, it will return uncertainty based on rate
        variation. If both are present, the wider interval will be returned.

        Parameters
        ----------

         node : PhyloTree.Clade
            The node for which the posterior region is to be calculated

         interval : float
            Float specifying who much of the posterior probability is
            to be contained in the region

        Returns
        -------
         max_posterior_region : numpy array
            Array with two numerical dates delineating the high posterior region

        '''
        if node.marginal_inverse_cdf=="delta":
            return np.array([node.numdate, node.numdate])


        min_max = (node.marginal_pos_LH.xmin, node.marginal_pos_LH.xmax)
        min_date, max_date = [self.date2dist.to_numdate(x) for x in min_max][::-1]
        if node.marginal_pos_LH.peak_pos == min_max[0]: #peak on the left
            return self.get_confidence_interval(node, (0, fraction))
        elif node.marginal_pos_LH.peak_pos == min_max[1]: #peak on the right
            return self.get_confidence_interval(node, (1.0-fraction ,1.0))
        else: # peak in the center of the distribution
            rate_contribution = self.date_uncertainty_due_to_rate(node, ((1-fraction)*0.5, 1.0-(1.0-fraction)*0.5))

            # construct height to position interpolators left and right of the peak
            # this assumes there is only one peak --- might fail in odd cases
            from scipy.interpolate import interp1d
            from scipy.optimize import minimize_scalar as minimize
            pidx = np.argmin(node.marginal_pos_LH.y)
            pval = np.min(node.marginal_pos_LH.y)
            left =  interp1d(node.marginal_pos_LH.y[:(pidx+1)]-pval, node.marginal_pos_LH.x[:(pidx+1)],
                            kind='linear', fill_value=min_max[0], bounds_error=False)
            right = interp1d(node.marginal_pos_LH.y[pidx:]-pval, node.marginal_pos_LH.x[pidx:],
                            kind='linear', fill_value=min_max[1], bounds_error=False)

            # function to minimize -- squared difference between prob mass and desired fracion
            def func(x, thres):
                interval = np.array([left(x), right(x)]).squeeze()
                return (thres - np.diff(node.marginal_cdf(np.array(interval))))**2

            # minimza and determine success
            sol = minimize(func, bracket=[0,10], args=(fraction,))
            if sol['success']:
                mutation_contribution = self.date2dist.to_numdate(np.array([right(sol['x']), left(sol['x'])]).squeeze())
            else: # on failure, return standard confidence interval
                mutation_contribution = None

            return self.combine_confidence(node.numdate, (min_date, max_date),
                                      c1=rate_contribution, c2=mutation_contribution)


if __name__=="__main__":

    pass
