import utils
import numpy as np
import config as ttconf
from treeanc import TreeAnc
from distribution import Distribution
from branch_len_interpolator import BranchLenInterpolator
from node_interpolator import NodeInterpolator
import collections

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

    def __init__(self,  dates=None, debug=False, real_dates=True, *args, **kwargs):
        """
        ClockTree constructor

        Parameters
        ----------

         dates : dic, None
            {leaf_name:leaf_date} dictionary

         debug : bool
            If True, the debug mode is ON, which means no or less clean-up of
            obsolete parameters to control program  execution in intermediate
            states. In debug mode, the python debugger is also allowed to interrupt
            program execution with intercative shell if an error occurs.

         real_dates : bool
            If True, some additional checks for the input dates sanity will be
            performed.

        Keyword Args
        ------------
            Kwargs needed to construct parent class (TreeAnc)

        """
        super(ClockTree, self).__init__(*args, **kwargs)
        if dates is None:
            raise ValueError("ClockTree requires date constraints!")

        self.debug=debug
        self.real_dates = real_dates
        self.date_dict = dates
        self.date2dist = None  # we do not know anything about the conversion
        self.n_integral = ttconf.NINTEGRAL
        self.rel_tol_prune = ttconf.REL_TOL_PRUNE
        self.rel_tol_refine = ttconf.REL_TOL_REFINE

        for node in self.tree.find_clades(order='postorder'):
            if node.name in self.date_dict:
                try:
                    tmp = np.mean(self.date_dict[node.name])
                    node.numdate_given = self.date_dict[node.name]
                    node.bad_branch = False
                except:
                    self.logger("WARNING: ClockTree.init: node %s has a bad date: %s"%(node.name, str(self.date_dict[node.name])), 2, warn=True)
                    node.numdate_given = None
                    node.bad_branch = True
            else: # nodes without date contraints

                node.numdate_given = None

                if node.is_terminal():
                    # Terminal branches without date constraints marked as 'bad'
                    node.bad_branch = True
                else:
                    # If all branches dowstream are 'bad', and there is no date constraint for
                    # this node, the branch is marked as 'bad'
                    node.bad_branch = np.all([x.bad_branch for x in node])


    @property
    def date2dist(self):
        return self._date2dist

    @date2dist.setter
    def date2dist(self, val):
        if val is None:
            self._date2dist = None
        else:
            self.logger("ClockTime.date2dist: Setting new molecular clock. rate=%.3e, R^2=%.4f"%(val.clock_rate, val.r_val**2), 2)
            self._date2dist = val


    def init_date_constraints(self, ancestral_inference=False, clock_rate=None, **kwarks):
        """
        Get the conversion coefficients between the dates and the branch
        lengths as they are used in ML computations. The conversion formula is
        assumed to be 'length = k*numdate + b'. For convenience, these
        coefficients as well as regression parameters are stored in the
        dates2dist object.

        ..Note:: that tree must have dates set to all nodes before calling this
        function. (This is accomplished by calling load_dates func).

        Parameters
        ----------

         ancestral_inference: bool
            Whether or not to reinfer ancestral sequences
            done by default when ancestral sequences are missing

        """
        self.logger("ClockTree.init_date_constraints...",2)
        self.tree.coalescent_joint_LH = 0
        if ancestral_inference or (not hasattr(self.tree.root, 'sequence')):
            self.infer_ancestral_sequences('ml', marginal=False, sample_from_profile='root',**kwarks)

        # set the None  for the date-related attributes in the internal nodes.
        # make interpolation objects for the branches
        self.logger('ClockTree.init_date_constraints: Initializing branch length interpolation objects...',3)
        for node in self.tree.find_clades(order='postorder'):
            if node.up is None:
                node.branch_length_interpolator = None
            else:
                # copy the merger rate and gamma if they are set
                if hasattr(node,'branch_length_interpolator') and node.branch_length_interpolator is not None:
                    gamma = node.branch_length_interpolator.gamma
                    merger_cost = node.branch_length_interpolator.merger_cost
                else:
                    gamma = 1.0
                    merger_cost = None
                node.branch_length_interpolator = BranchLenInterpolator(node, self.gtr, one_mutation=self.one_mutation)
                node.branch_length_interpolator.merger_cost = merger_cost
                node.branch_length_interpolator.gamma = gamma
        self.date2dist = utils.DateConversion.from_tree(self.tree, clock_rate)

        # make node distribution objects
        for node in self.tree.find_clades(order="postorder"):
            # node is constrained
            if hasattr(node, 'numdate_given') and node.numdate_given is not None:
                # set the absolute time before present in branch length units
                if np.isscalar(node.numdate_given):
                    tbp = self.date2dist.get_time_before_present(node.numdate_given)
                    node.date_constraint = Distribution.delta_function(tbp, weight=1.0)
                else:
                    tbp = self.date2dist.get_time_before_present(np.array(node.numdate_given))
                    node.date_constraint = Distribution(tbp, np.ones_like(tbp), is_log=False)

                if hasattr(node, 'bad_branch') and node.bad_branch==True:
                    self.logger("ClockTree.init_date_constraints -- WARNING: Branch is marked as bad"
                                ", excluding it from the optimization process"
                                " Will be optimized freely", 4, warn=True)
            else: # node without sampling date set
                node.numdate_given = None
                node.date_constraint = None


    def make_time_tree(self, time_marginal=False, **kwargs):
        '''
        Use the date constraints to calculate the most likely positions of
        unconstrained nodes.

        Parameters
        ----------

         time_marginal : bool
            Whether use marginal reconstruction for node positions or not

        Keyword Args
        ------------

         Kwargs needed to initialize dates constraints

        '''
        self.logger("ClockTree: Maximum likelihood tree optimization with temporal constraints:",1)
        self.init_date_constraints(**kwargs)

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
                    node.joint_pos_Lx = Distribution(x, node.branch_length_interpolator(bl), is_log=True)
                    node.joint_pos_Cx = Distribution(x, bl) # map back to the branch length
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
        return LH + self.gtr.sequence_logLH(self.tree.root.cseq, pattern_multiplicity=self.multiplicity)


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
                    node.marginal_pos_Lx = Distribution(x, node.branch_length_interpolator(bl), is_log=True)

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
                    from utils import numeric_date
                    x = [parent.numdate, numeric_date()]
                    msg_parent_to_node = NodeInterpolator(x, [1.0, 1.0])

                # integral message, which delivers to the node the positional information
                # from the complementary subtree
                res, res_t = NodeInterpolator.convolve(msg_parent_to_node, node.branch_length_interpolator,
                                                    max_or_integral='integral',
                                                    inverse_time=False,
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
        this fucntion converts the estimated "time_before_present" properties of all nodes
        to numerical dates stored in the "numdate" attribute. This date is further converted
        into a human readable date string in format %Y-%m-%d assuming the usual calendar


        Returns
        -------
         None
            All manipulations are done in place on the tree

        '''
        from datetime import datetime, timedelta
        now = utils.numeric_date()
        for node in self.tree.find_clades():
            years_bp = self.date2dist.to_years(node.time_before_present)
            if years_bp < 0 and self.real_dates:
                if not hasattr(node, "bad_branch") or node.bad_branch==False:
                    self.logger("ClockTree.convert_dates -- WARNING: The node is later than today, but it is not "
                        "marked as \"BAD\", which indicates the error in the "
                        "likelihood optimization.",4 , warn=True)
                else:
                    self.logger("ClockTree.convert_dates -- WARNING: node which is marked as \"BAD\" optimized "
                        "later than present day",4 , warn=True)

            node.numdate = now - years_bp

            # set the human-readable date
            days = 365.25 * (node.numdate - int(node.numdate))
            year = int(node.numdate)
            try:  # datetime will only operate on dates after 1900
                n_date = datetime(year, 1, 1) + timedelta(days=days)
                node.date = datetime.strftime(n_date, "%Y-%m-%d")
            except:
                # this is the approximation not accounting for gap years etc
                n_date = datetime(1900, 1, 1) + timedelta(days=days)
                node.date = str(year) + "-" + str(n_date.month) + "-" + str(n_date.day)


    def branch_length_to_years(self):
        '''
        This function sets branch length to reflect the date differences between parent and child
        nodes measured in years. Should only be called after convert_dates has been called

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


    def get_confidence_interval(self, node, interval = (0.05, 0.95)):
        '''
        If temporal reconstruction was done using the marginal ML mode, the entire distribution of
        times is available. this function here determines the 90%( or other) confidence interval defines as the
        range where 5% of probability are below and above. Note that this does not necessarily contain
        the highest probability position.

        Parameters
        ----------

         node : PhyloTree.Clade
            The node for which the confidence interval is to be calculated

         interval : tuple, list default=(0.05, 0.95)
            Array like of length two defining the bounds of the confidence interval

        Returns
        -------

         confidence_interval : numpy array
            Array with two numerical dates delineating the confidence interval

        '''
        if hasattr(node, "marginal_inverse_cdf"):
            if node.marginal_inverse_cdf=="delta":
                return np.array([node.numdate, node.numdate])
            else:
                return self.date2dist.to_numdate(node.marginal_inverse_cdf(np.array(interval))[::-1])
        else:
            return np.array([np.nan, np.nan])


    def get_max_posterior_region(self, node, fraction = 0.9):
        '''
        If temporal reconstruction was done using the marginal ML mode, the entire distribution of
        times is available. this function here determines the 95% confidence interval defines as the
        range where 5% of probability are below and above. Note that this does not necessarily contain
        the highest probability position.

        Parameters
        ----------

         node : PhyloTree.Clade
            The node for which the confidence region is to be calculated

         interval : float
            Float specifying who much of the posterior probability is
            to be contained in the region

        Returns
        -------
         max_posterior_region : numpy array
            Array with two numerical dates delineating the high posterior region

        '''
        if not hasattr(node, "marginal_inverse_cdf"):
            return np.array([np.nan, np.nan])

        if node.marginal_inverse_cdf=="delta":
            return np.array([node.numdate, node.numdate])

        min_max = (node.marginal_pos_LH.xmin, node.marginal_pos_LH.xmax)
        if node.marginal_pos_LH.peak_pos == min_max[0]: #peak on the left
            return self.get_confidence_interval(node, (0, fraction))
        elif node.marginal_pos_LH.peak_pos == min_max[1]: #peak on the right
            return self.get_confidence_interval(node, (1.0-fraction ,1.0))
        else: # peak in the center of the distribution

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
                return self.date2dist.to_numdate(np.array([right(sol['x']), left(sol['x'])]).squeeze())
            else: # on failure, return standard confidence interval
                return self.get_confidence_interval(node, (0.5*(1-fraction), 1-0.5*(1-fraction)))


if __name__=="__main__":

    pass
