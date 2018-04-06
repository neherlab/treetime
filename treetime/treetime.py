from __future__ import print_function, division
from clock_tree import ClockTree
import utils as ttutils
import config as ttconf
import numpy as np
from scipy import optimize as sciopt
from Bio import Phylo
from version import tt_version as __version__


class TreeTime(ClockTree):
    """
    TreeTime is a wrapper class to ClockTree that adds additional functionality
    such as reroot, detection and exclusion of outliers, resolution of polytomies
    using temporal information, and relaxed molecular clock models
    """

    def __init__(self, *args,**kwargs):
        """
        TreeTime constructor

        Parameters
        -----------
         Arguments to construct ClockTree

        Keyword Args
        ------------
         Kwargs to construct ClockTree

        """
        super(TreeTime, self).__init__(*args, **kwargs)
        self.n_iqd = ttconf.NIQD

    def run(self, root=None, infer_gtr=True, relaxed_clock=None, n_iqd = None,
            resolve_polytomies=True, max_iter=0, Tc=None, fixed_clock_rate=None,
            time_marginal=False, use_input_branch_length = False, **kwargs):

        """
        Run TreeTime reconstruction. Based on the input parameters, it divides
        the analysis into semi-independent jobs and conquers them one-by one
        gradually optimizing the tree given the temporal constarints and leaf
        nodes sequences.

        Parameters
        ----------

         root : str, None
            Try to find better root position on a given tree. If string is passed,
            the root will be searched according to the specified method. Available
            reroot methods are: 'best', 'oldest', '<leaf_name>'

            If None, use tree as-is.

         infer_gtr : bool default True
            Should infer GTR model?

         relaxed_clock : dic, None
            If not None, use autocorrelated molecular clock model. Specify the
            clock parameters as {slack:<slack>, coupling:<coupling>} dictionary.

         n_iqd : int, None
            If not None, filter tree nodes, which do not obey molecular clock
            for the particular tree. The nodes, which deviate more than
            :code:`n_iqd` interquantile intervals from the molecular clock
            regression will be marked as 'BAD' and not account in the TreeTime
            analysis

         resolve_polytomies : bool
            Should attempt to resolve multiple mergers?

         max_iter : int
            Maximum number of iterations to optimize the tree

         Tc : float, str, None
            If not None, use coalescent model to correct the branch lengths by
            introducing merger costs.

            If Tc is float, it is interpreted as the coalescence time scale

            If Tc is str, it should be one of (:code:`opt`, :code:`skyline`)

         fixed_clock_rate : float, None
            If None, infer clock rate from the molecular clock

            If float, use this rate

         time_marginal : bool default False
            Should perform marginal reconstruction of the node's positions?

         use_input_branch_length : bool
            If True, rely on the branch lengths in the imput tree and skip directly
            to the maximum-likelihood ancestral sequence reconstruction.
            Otherwise, perform preliminary sequence reconstruction using parsimony
            algorithm and do branch length optimization

        Keyword Args
        ------------

         Additional arguments needed by the dowstream funcitons


        """
        # determine how to reconstruct and sample sequences
        seq_kwargs = {"marginal":False, "sample_from_profile":"root"}
        if "fixed_pi" in kwargs:
            seq_kwargs["fixed_pi"] = kwargs["fixed_pi"]
        if "do_marginal" in kwargs:
            time_marginal=kwargs["do_marginal"]

        # initially, infer ancestral sequences and infer gtr model if desired
        if use_input_branch_length:
            self.infer_ancestral_sequences(infer_gtr=infer_gtr, **seq_kwargs)
            self.prune_short_branches()
        else:
            self.optimize_sequences_and_branch_length(infer_gtr=infer_gtr,
                                                  max_iter=2, prune_short=True, **seq_kwargs)
        avg_root_to_tip = np.mean([x.dist2root for x in self.tree.get_terminals()])

        # optionally reroot the tree either by oldest, best regression or with a specific leaf
        if n_iqd or root=='clock_filter':
            if "plot_rtt" in kwargs and kwargs["plot_rtt"]:
                plot_rtt=True
            else:
                plot_rtt=False
            self.clock_filter(reroot='best' if root=='clock_filter' else root,
                              n_iqd=n_iqd, plot=plot_rtt)
        elif root is not None:
            self.reroot(root=root)

        if use_input_branch_length:
            self.infer_ancestral_sequences(**seq_kwargs)
        else:
            self.optimize_sequences_and_branch_length(max_iter=2,**seq_kwargs)

        # infer time tree and optionally resolve polytomies
        self.logger("###TreeTime.run: INITIAL ROUND",0)
        self.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False, **kwargs)

        self.LH = [[self.tree.sequence_marginal_LH if seq_kwargs['marginal'] else self.tree.sequence_joint_LH,
                    self.tree.positional_joint_LH, 0.0]]

        # iteratively reconstruct ancestral sequences and re-infer
        # time tree to ensure convergence.
        niter = 0
        while niter < max_iter:

            self.logger("###TreeTime.run: ITERATION %d out of %d iterations"%(niter+1,max_iter),0)
            # add coalescent prior
            if Tc and (Tc is not None):
                from merger_models import Coalescent
                self.logger('TreeTime.run: adding coalescent prior with Tc='+str(Tc),1)
                self.merger_model = Coalescent(self.tree, Tc=avg_root_to_tip,
                                               date2dist=self.date2dist, logger=self.logger)

                if Tc=='skyline' and niter==max_iter-1: # restrict skyline model optimization to last iteration
                    self.merger_model.optimize_skyline(**kwargs)
                    self.logger("optimized a skyline ", 2)
                else:
                    if Tc in ['opt', 'skyline']:
                        self.merger_model.optimize_Tc()
                        self.logger("optimized Tc to %f"%self.merger_model.Tc.y[0], 2)
                    else:
                        try:
                            self.merger_model.set_Tc(Tc)
                        except:
                            self.logger("setting of coalescent time scale failed", 1, warn=True)

                self.merger_model.attach_to_tree()

            # estimate a relaxed molecular clock
            if relaxed_clock:
                self.relaxed_clock(**relaxed_clock)

            n_resolved=0
            if resolve_polytomies:
                # if polytomies are found, rerun the entire procedure
                n_resolved = self.resolve_polytomies()
                if n_resolved:
                    self.prepare_tree()
                    # when using the input branch length, only infer ancestral sequences
                    if use_input_branch_length:
                        self.infer_ancestral_sequences(**seq_kwargs)
                    else: # otherwise reoptimize branch length while preserving branches without mutations
                        self.optimize_sequences_and_branch_length(prune_short=False,
                                                                  max_iter=0, **seq_kwargs)

                    self.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False, **kwargs)
                    ndiff = self.infer_ancestral_sequences('ml',**seq_kwargs)
                else:
                    ndiff = self.infer_ancestral_sequences('ml',**seq_kwargs)
                    self.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False, **kwargs)
            elif (Tc and (Tc is not None)) or relaxed_clock: # need new timetree first
                self.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False, **kwargs)
                ndiff = self.infer_ancestral_sequences('ml',**seq_kwargs)
            else: # no refinements, just iterate
                ndiff = self.infer_ancestral_sequences('ml',**seq_kwargs)
                self.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False, **kwargs)

            self.tree.coalescent_joint_LH = self.merger_model.total_LH() if Tc else 0.0

            self.LH.append([self.tree.sequence_marginal_LH if seq_kwargs['marginal'] else self.tree.sequence_joint_LH,
                            self.tree.positional_joint_LH, self.tree.coalescent_joint_LH])
            niter+=1

            if ndiff==0 & n_resolved==0:
                self.logger("###TreeTime.run: CONVERGED",0)
                break

        # if marginal reconstruction requested, make one more round with marginal=True
        # this will set marginal_pos_LH, which to be used as error bar estimations
        if time_marginal:
            self.logger("###TreeTime.run: FINAL ROUND - confidence estimation via marginal reconstruction", 0)
            self.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=time_marginal, **kwargs)




    def clock_filter(self, reroot='best', n_iqd=None, plot=False):
        '''
        Labels outlier branches that don't seem to follow a molecular clock
        and excludes them from subsequent the molecular clock estimate and
        the timetree propagation

        Parameters
        ----------
         reroot : str, None
            Method to find the best root in the tree.

         n_iqd : int, None
            Number of iqd intervals. The outlier nodes are those which do not fall
            into :math:`IQD\cdot n_iqd` interval (:math:`IQD` is the interval between
            75 and 25 percentiles)

            if None, the default (3) assumed

         plot : bool
            Should plot the reults?

        '''
        if n_iqd is None:
            n_iqd = self.n_iqd

        terminals = self.tree.get_terminals()
        if reroot:
            self.reroot(root=reroot)
            icpt, clock_rate = self.tree.root._alpha, self.tree.root._beta
        else:
            tmp_date2dist = ttutils.DateConversion.from_tree(self.tree)
            icpt, clock_rate = tmp_date2dist.intercept, tmp_date2dist.clock_rate

        res = {}
        for node in terminals:
            if hasattr(node, 'numdate_given') and  (node.numdate_given is not None):
                res[node] = node.dist2root - clock_rate*np.mean(node.numdate_given) - icpt
        residuals = np.array(res.values())
        iqd = np.percentile(residuals,75) - np.percentile(residuals,25)
        for node,r in res.iteritems():
            if abs(r)>n_iqd*iqd and node.up.up is not None:
                self.logger('TreeTime.ClockFilter: marking %s as outlier, residual %f interquartile distances'%(node.name,r/iqd), 3)
                node.bad_branch=True
            else:
                node.bad_branch=False

        if plot:
            self.plot_root_to_tip()

        # redo root estimation after outlier removal
        if reroot:
            self.reroot(root=reroot)


    def plot_root_to_tip(self, add_internal=False, label=True, ax=None, **kwargs):
        """
        Plot root-to-tip regression

        Parameters
        ----------

         add_internal : bool
            Should plot internal node positoins?

         label : bool
            Should label the plots?

         ax: matplotlib axes, None
            If not None, use the provided matplotlib axes to plot the results

        Keyword Args
        ------------
         Additional arguments for matplotlib.pyplot.scatter function

        """
        import matplotlib.pyplot as plt
        tips = self.tree.get_terminals()
        internal = self.tree.get_nonterminals()
        if ax is None:
            plt.figure()
            ax=plt.subplot(111)
        dates = np.array([np.mean(n.numdate_given) for n in tips if n.numdate_given is not None])
        dist = np.array([n.dist2root for n in tips if n.numdate_given is not None])
        ind = np.array([n.bad_branch for n in tips if n.numdate_given is not None])
        # plot tips
        ax.scatter(dates[ind], dist[ind]  , c='r', label="bad tips" if label else "" , **kwargs)
        ax.scatter(dates[~ind], dist[~ind], c='g', label="good tips" if label else "", **kwargs)
        if add_internal and hasattr(self.tree.root, "numdate"):
            dates = np.array([n.numdate for n in internal])
            dist = np.array([n.dist2root for n in internal])
            ind = np.array([n.bad_branch for n in internal])
            # plot internal
            ax.scatter(dates[~ind], dist[~ind], c='b', marker='<', label="internal" if label else "", **kwargs)

        if label:
            ax.legend(loc=2)
        ax.set_ylabel('root-to-tip distance')
        ax.set_xlabel('date')
        ax.ticklabel_format(useOffset=False)
        plt.tight_layout()


    def reroot(self,root='best'):
        """
        Find best root and re-root the tree to the new root

        Parameters
        ----------

         root : str
            Which method should be used to find the best root. Available methods are:

            :code:`best` - maximize root-to-tip regression coefficient

            :code:`oldest` - choose the oldest node

            :code:`<node_name>` - reroot to the node with name :code:`<node_name>`

            :code:`[<node_name1>, <node_name2>, ...]` - reroot to the MRCA of these nodes
        """
        self.logger("TreeTime.reroot: with method or node: %s"%root,1)
        for n in self.tree.find_clades():
            n.branch_length=n.mutation_length
        from Bio import Phylo
        if isinstance(root,Phylo.BaseTree.Clade):
            new_root = root
        elif isinstance(root, list):
            new_root = self.tree.common_ancestor(*root)
        elif root in self._leaves_lookup:
            new_root = self._leaves_lookup[root]
        elif root=='oldest':
            new_root = sorted([n for n in self.tree.get_terminals()
                               if n.numdate_given is not None],
                               key=lambda x:np.mean(x.numdate_given))[0]
        elif root=='best':
            new_root = self.reroot_to_best_root(criterium='residual')
        elif root=='rsq':
            new_root = self.reroot_to_best_root(criterium='rsq')
        elif root=='residual':
            new_root = self.reroot_to_best_root(criterium='residual')
        elif root=='min_dev':
            new_root = self.reroot_to_best_root(criterium='min_dev')
        else:
            self.logger('TreeTime.reroot -- WARNING: unsupported rooting mechanisms or root not found',2,warn=True)
            return

        self.logger("TreeTime.reroot: Tree is being re-rooted to node "
                    +('new_node' if new_root.name is None else new_root.name), 2)
        if isinstance(root, list):
            #this forces a bifurcating root, as we want. Branch lengths will be reoptimized anyway.
            #(Without outgroup_branch_length, gives a trifurcating root, but this will mean
            #mutations may have to occur multiple times.)
            self.tree.root_with_outgroup(new_root, outgroup_branch_length=new_root.branch_length/2)
        else:
            self.tree.root_with_outgroup(new_root)
       # new nodes are produced when rooting with a terminal node, copy this clock info

        if new_root.is_terminal():
            if hasattr(new_root, "_alpha"):
                self.tree.root._alpha = new_root._alpha
            if hasattr(new_root, "_beta"):
                self.tree.root._beta = new_root._beta
            if hasattr(new_root, "_R2"):
                self.tree.root._R2 = new_root._R2
            if hasattr(new_root, "_residual"):
                self.tree.root._residual = new_root._residual

        self.tree.root.branch_length = self.one_mutation
        for n in self.tree.find_clades():
            n.mutation_length=n.branch_length
        self.tree.root.numdate_given = None
        # set root.gamma bc root doesn't have a branch_length_interpolator but gamma is needed
        if not hasattr(self.tree.root, 'gamma'):
            self.tree.root.gamma = 1.0
        self.prepare_tree()


    def resolve_polytomies(self, merge_compressed=False, rerun=True):
        """
        Resolve the polytomies on the tree.
        The function scans the tree, resolves polytomies in case there are any,
        and re-optimizes the tree with new topology. Note that polytomies are only
        resolved if that would result in higher likelihood. Sometimes, stretching
        two or more branches that carry several mutations are less costly than
        an additional branch with zero mutations (long branches are not stiff,
        short branches are.)

        Parameters
        ----------
         merge_compressed : bool
            Whether to keep compressed branches as polytomies or
            return a strictly binary tree.

        """
        self.logger("TreeTime.resolve_polytomies: resolving multiple mergers...",1)

        poly_found=0
        for n in self.tree.find_clades():
            if len(n.clades) > 2:
                prior_n_clades = len(n.clades)
                self._poly(n, merge_compressed)
                poly_found+=prior_n_clades - len(n.clades)

        obsolete_nodes = [n for n in self.tree.find_clades() if len(n.clades)==1 and n.up is not None]
        for node in obsolete_nodes:
            self.logger('TreeTime.resolve_polytomies: remove obsolete node '+node.name,4)
            if node.up is not None:
                self.tree.collapse(node)

        if poly_found:
            self.logger('TreeTime.resolve_polytomies: introduces %d new nodes'%poly_found,3)
        else:
            self.logger('TreeTime.resolve_polytomies: No more polytomies to resolve',3)
        return poly_found


    def _poly(self, clade, merge_compressed, verbose=1):

        """
        Function to resolve polytomies for a given parent node. If the
        number of the direct decendants is less than three (not a polytomy), does
        nothing. Otherwise, for each pair of nodes, assess the possible LH increase
        which could be gained by merging the two nodes. The increase in the LH is
        basically the tradeoff between the gain of the LH due to the changing the
        branch lenghts towards the optimal values and the decrease due to the
        introduction of the new branch with zero optimal length.
        """

        from branch_len_interpolator import BranchLenInterpolator
        from Bio import Phylo

        zero_branch_slope = self.gtr.mu*self.seq_len

        def _c_gain(t, n1, n2, parent):
            """
            cost gain if nodes n1, n2 are joined and their parent is placed at time t
            cost gain = (LH loss now) - (LH loss when placed at time t)
            """
            cg2 = n2.branch_length_interpolator(parent.time_before_present - n2.time_before_present) - n2.branch_length_interpolator(t - n2.time_before_present)
            cg1 = n1.branch_length_interpolator(parent.time_before_present - n1.time_before_present) - n1.branch_length_interpolator(t - n1.time_before_present)
            cg_new = - zero_branch_slope * (parent.time_before_present - t) # loss in LH due to the new branch
            return -(cg2+cg1+cg_new)

        def cost_gain(n1, n2, parent):
            """
            cost gained if the two nodes would have been connected.
            """
            try:
                cg = sciopt.minimize_scalar(_c_gain,
                    bounds=[max(n1.time_before_present,n2.time_before_present), parent.time_before_present],
                    method='Bounded',args=(n1,n2, parent))
                return cg['x'], - cg['fun']
            except:
                self.logger("TreeTime._poly.cost_gain: optimization of gain failed", 3, warn=True)
                return parent.time_before_present, 0.0


        def merge_nodes(source_arr, isall=False):
            mergers = np.array([[cost_gain(n1,n2, clade) if i1<i2 else (0.0,-1.0)
                                    for i1,n1 in enumerate(source_arr)]
                                for i2, n2 in enumerate(source_arr)])
            LH = 0
            while len(source_arr) > 1 + int(isall):
                # max possible gains of the cost when connecting the nodes:
                # this is only a rough approximation because it assumes the new node positions
                # to be optimal
                new_positions = mergers[:,:,0]
                cost_gains = mergers[:,:,1]
                # set zero to large negative value and find optimal pair
                np.fill_diagonal(cost_gains, -1e11)
                idxs = np.unravel_index(cost_gains.argmax(),cost_gains.shape)
                if (idxs[0] == idxs[1]) or cost_gains.max()<0:
                    self.logger("TreeTime._poly.merge_nodes: node is not fully resolved "+clade.name,4)
                    return LH

                n1, n2 = source_arr[idxs[0]], source_arr[idxs[1]]
                LH += cost_gains[idxs]

                new_node = Phylo.BaseTree.Clade()

                # fix positions and branch lengths
                new_node.time_before_present = new_positions[idxs]
                new_node.branch_length = clade.time_before_present - new_node.time_before_present
                new_node.clades = [n1,n2]
                n1.branch_length = new_node.time_before_present - n1.time_before_present
                n2.branch_length = new_node.time_before_present - n2.time_before_present

                # set parameters for the new node
                new_node.up = clade
                n1.up = new_node
                n2.up = new_node
                new_node.cseq = clade.cseq
                self._store_compressed_sequence_to_node(new_node)

                new_node.mutations = []
                new_node.mutation_length = 0.0
                new_node.branch_length_interpolator = BranchLenInterpolator(new_node, self.gtr, one_mutation=self.one_mutation)
                clade.clades.remove(n1)
                clade.clades.remove(n2)
                clade.clades.append(new_node)
                self.logger('TreeTime._poly.merge_nodes: creating new node as child of '+clade.name,3)
                self.logger("TreeTime._poly.merge_nodes: Delta-LH = " + str(cost_gains[idxs].round(3)), 3)

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

        stretched = [c for c  in clade.clades if c.mutation_length < c.clock_length]
        compressed = [c for c in clade.clades if c not in stretched]

        if len(stretched)==1 and merge_compressed==False:
            return 0.0

        LH = merge_nodes(stretched, isall=len(stretched)==len(clade.clades))
        if merge_compressed and len(compressed)>1:
            LH += merge_nodes(compressed, isall=len(compressed)==len(clade.clades))

        return LH


    def print_lh(self, joint=True):
        """
        Print the total likelihood of the tree given the constrained leaves

        Parameters
        ----------

         joint : bool
            Whether joint or marginal LH should be printed

        """
        try:
            u_lh = self.tree.unconstrained_sequence_LH
            if joint:
                s_lh = self.tree.sequence_joint_LH
                t_lh = self.tree.positional_joint_LH
                c_lh = self.tree.coalescent_joint_LH
            else:
                s_lh = self.tree.sequence_marginal_LH
                t_lh = self.tree.positional_marginal_LH
                c_ls = 0

            print ("###  Tree Log-Likelihood  ###\n"
                " Sequence log-LH without constraints: \t%1.3f\n"
                " Sequence log-LH with constraints:    \t%1.3f\n"
                " TreeTime sequence log-LH:            \t%1.3f\n"
                " Coalescent log-LH:                   \t%1.3f\n"
               "#########################"%(u_lh, s_lh,t_lh, c_lh))
        except:
            print("ERROR. Did you run the corresponding inference (joint/marginal)?")


    def relaxed_clock(self, slack=None, coupling=None, **kwargs):
        """
        Allow the mutation rate to vary on the tree (relaxed molecular clock).
        Changes of the mutation rates from one branch to another are penalized.
        In addition, deviation of the mutation rate from the mean rate are
        penalized.

        Parameters
        ----------
         slack : float
            Maximum change in substitution rate between parent and child nodes

         coupling : float
            Maximum difference in substitution rates in sibling nodes

        """
        if slack is None: slack=ttconf.MU_ALPHA
        if coupling is None: coupling=ttconf.MU_BETA
        self.logger("TreeTime.relaxed_clock: slack=%f, coupling=%f"%(slack, coupling),2)

        c=1.0/self.one_mutation
        for node in self.tree.find_clades(order='postorder'):
            opt_len = node.mutation_length

            # opt_len \approx 1.0*len(node.mutations)/node.profile.shape[0] but calculated via gtr model
            # contact term: stiffness*(g*bl - bl_opt)^2 + slack(g-1)^2 =
            #               (slack+bl^2) g^2 - 2 (bl*bl_opt+1) g + C= k2 g^2 + k1 g + C
            node._k2 = slack + c*node.branch_length**2/(opt_len+self.one_mutation)
            node._k1 = -2*(c*node.branch_length*opt_len/(opt_len+self.one_mutation) + slack)
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

        for node in self.tree.find_clades(order='preorder'):
            if node.up is None:
                node.gamma =- 0.5*node._k1/node._k2
            else:
                if node.up.up is None:
                    g_up = node.up.gamma
                else:
                    g_up = node.up.branch_length_interpolator.gamma
                node.branch_length_interpolator.gamma = (coupling*g_up - 0.5*node._k1)/(coupling+node._k2)

###############################################################################
### rerooting
###############################################################################
    def find_best_root_and_regression(self, criterium='rsq'):
        """
        Find the best root for the tree in linear time, given the timestamps of
        the leaves. The branch lengths should be optimized prior to the run;
        the terminal nodes should have the timestamps assigned as numdate_given
        attribute.
        """
        sum_ti =  np.sum([np.mean(node.numdate_given)*node.count for node in self.tree.get_terminals() if (not node.bad_branch)])
        sum_ti2 = np.sum([np.mean(node.numdate_given)**2*node.count for node in self.tree.get_terminals() if (not node.bad_branch)])
        N = 1.0*np.sum([node.count for node in self.tree.get_terminals() if not node.bad_branch])
        tip_count = 1.0*np.sum([1.0 for node in self.tree.get_terminals() if not node.bad_branch])
        if tip_count<2:
            self.logger("****ERROR: TreeTime.find_best_root_and_regression: need at least two dates to reroot!", 0, warn=True)
            self.logger("****ERROR: only %d tips have valid dates!"%N, 0, warn=True)
            return selt.tree.root, np.nan, np.nan

        Ninv = 1.0/N
        time_variance = (N*sum_ti2 - sum_ti**2)*Ninv**2

        #  fill regression terms for one of the two subtrees
        for node in self.tree.find_clades(order='postorder'):  # children first, msg to parents
            if node.is_terminal():  # inititalize the leaves
                #  will not rely on the standard func - count terminals directly
                node._st_n_leaves = 0 if node.bad_branch else node.count
                node._st_di = 0.0
                node._st_diti = 0.0
                node._st_di2 = 0.0

                if node.bad_branch:
                    node._st_ti = 0
                else:
                    node._st_ti = np.mean(node.numdate_given)*node.count

                node._ti = sum_ti
            else:
                #  for non-terminal nodes,
                node._st_ti = np.sum([k._st_ti for k in node.clades])
                node._st_n_leaves = np.sum([k._st_n_leaves for k in node.clades])
                node._st_di   = np.sum([k._st_di + k._st_n_leaves*k.branch_length for k in node.clades])
                node._st_diti = np.sum([k._st_diti + k.branch_length*k._st_ti for k in node.clades])
                node._st_di2  = np.sum([k._st_di2 + 2*k._st_di*k.branch_length + k._st_n_leaves*k.branch_length**2
                                       for k in node.clades])
                node._ti = sum_ti
                node.bad_branch = np.all([x.bad_branch for x in node])

        best_root = self.tree.root
        best_root_any = self.tree.root
        for node in self.tree.find_clades(order='preorder'):  # root first

            if node.up is None:
                # assign the values for the root node
                node._di   = node._st_di
                node._diti = node._st_diti
                node._di2  = node._st_di2

                dist_variance = (N*node._di2 - node._di**2)*(Ninv**2)
                disttime_cov = (N*node._diti - sum_ti*node._di)*(Ninv**2)
                time_variance = time_variance

                node._beta = disttime_cov/time_variance
                node._alpha = (node._di - node._beta*sum_ti)/N
                node._residual = (node._di2 - 2*node._beta*node._diti - 2*node._alpha*node._di
                                   + node._beta**2*sum_ti2 + 2*node._alpha*node._beta*sum_ti + node._alpha**2*N)

                node._R2 = disttime_cov**2/(time_variance*dist_variance)
                node._R2_delta_x = 0.0 # there is no branch to move the root
            elif node.bad_branch:
                node._beta = np.nan
                node._alpha = np.nan
                node._R2 = 0.0
                node._R2_delta_x = 0.0
                node._residual = np.inf

            elif criterium=='rsq': # calculate the r^2 of the root to tip regression and pick the best intermediate position on the branch
                #  NOTE order of these computation matters
                n_up = N - node._st_n_leaves
                n_down = node._st_n_leaves
                L = node.branch_length
                node._di = node.up._di + (n_up-n_down)*L
                node._di2 = (node.up._di2 + 2*L*node.up._di
                            - 4*(L*(node._st_di + n_down*L))
                            + N*L**2)
                node._diti = node.up._diti + L*(sum_ti - 2*node._st_ti)


                ## Express Node's sum_Di as the function of parent's sum_Di
                # and **displacement from parent's node x** :
                # sum_Di = A1 + A2 * x
                A1 = node.up._di
                A2 = n_up - n_down

                ## Express Node's sum_Di**2 as the function of parent's params
                # and **displacement from parent's node x** :
                # sum_Di2 = B1 + B2 * x + B3 * x**2
                B1 = node.up._di2
                B2 = 2 * (node.up._di - 2 * node._st_di - 2 * L * n_down )
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
                    x1 = -1.0 # any arbitrary value out of range [0, L], see below
                    x2 = -1.0
                else:
                    # actual roots - the extrema for the R2(x) function
                    if np.abs(alpha * nu - beta * mu)>0:
                        x1 = (-1 * (alpha * delta - mu * gamma) + D2**0.5) / (alpha * nu - beta * mu)
                        x2 = (-1 * (alpha * delta - mu * gamma) - D2**0.5) / (alpha * nu - beta * mu)
                    else:
                        x1 = -(beta*delta - nu*gamma)/(alpha * delta - mu * gamma)
                        x2 = x1

                # possible positions, where the new root can possibly be located
                # (restrict to the branch length)
                max_points = [k for k in (x1,x2,L) if k >= 0 and k <= L]
                # values of the R2 at these positions
                R2s = [(alpha * x**2 + beta * x + gamma) / (mu * x**2 + nu * x + delta) / time_variance / N**2 for x in max_points]
                # choose the best R2
                node._R2 = np.max(R2s)
                # and set the position for the best R2 value
                node._R2_delta_x = L - max_points[np.argmax(R2s)]

                # for this position, define the clock_rate and intercept:
                node._beta = ((L - node._R2_delta_x) * (N * C2 - sum_ti * A2) + (N*C1-sum_ti*A1)) / time_variance / N**2
                node._alpha = (L - node._R2_delta_x) * A2 / N  + (A1 - node._beta * sum_ti) / N
            elif criterium in ['residual', 'min_dev']: # calculate the squared residuals and minimize as rooting criterium
                L = node.branch_length
                # number of nodes descendent and outgrouping this node
                n_up = N - node._st_n_leaves
                n_down = node._st_n_leaves
                # sum of branch length of the descendent tree and the outgrouping one
                node._di = node.up._di + (n_up-n_down)*L
                nd_down = node._st_di
                nd_up = node._di - node._st_di
                # sum of times of descendent tips and outgrouping ones
                nt_down = node._st_ti
                nt_up = sum_ti - node._st_ti

                # sum of squared branch length of the descendent tree and the outgrouping one
                node._di2 = (node.up._di2 + 2*L*node.up._di
                            - 4*(L*(node._st_di + n_down*L))
                            + N*L**2)
                nd2_down = node._st_di2
                nd2_up = node._di2 - node._st_di2

                # sum of timexbranch length of the descendent tree and the outgrouping one
                node._diti = node.up._diti + L*(sum_ti - 2*node._st_ti)
                ndt_down = node._st_diti
                ndt_up = node._diti - node._st_diti

                disttime_cov = (N*node._diti - sum_ti*node._di)*(Ninv**2)

                # decompose expression for alpha and beta into parts independent and linear in epsilon
                # where epsilon is the shift in root position along the branch
                beta_0 = disttime_cov/time_variance
                b = L*(nt_down - nt_up - (n_down-n_up)*sum_ti*Ninv)*Ninv/time_variance
                alpha_0 = (node._di - beta_0*sum_ti)*Ninv
                a = (-b*sum_ti + L*(n_down-n_up))*Ninv
                eps = -(((nd_down - beta_0*nt_down - n_down*alpha_0) - (nd_up - beta_0*nt_up - n_up*alpha_0))/
                      ((n_down*L - b*nt_down  - a*n_down) - (-n_up*L - b*nt_up - a*n_up)))

                # only shifts between 0 and 1 are admissible (where 0 is the node itself and 1 is the parent)
                eps = min(1,max(0,eps))
                beta = beta_0 + eps*b
                alpha = alpha_0 + eps*a

                # calculate the residual and assign the regression coefficients
                node._residual = (node._di2 - 2*beta*node._diti - 2*alpha*node._di
                                   + beta**2*sum_ti2 + 2*alpha*beta*sum_ti + alpha**2*N)
                node._residual += N*(L*eps)**2
                node._residual += 2*eps*L*(nd_down - beta*nt_down - n_down*alpha) - 2*eps*L*(nd_up - beta*nt_up - n_up*alpha)
                node._alpha = alpha
                node._beta = beta
                node._R2_delta_x = eps*L
            else:
                self.logger("TreeTime.find_best_root_and_regression: unknown criterium",0)

            if criterium=='rsq':
                if node.up is None:
                    self.logger("TreeTime.find_best_root_and_regression: Initial root: R2:%f\tclock_rate:%.3e"%(best_root._R2, best_root._beta),3)
                if  node._R2 > best_root_any._R2:
                    best_root_any = node
                if  (node._R2 > best_root._R2 and node._beta>0) or best_root._beta<0:
                    best_root = node
                    self.logger("TreeTime.find_best_root_and_regression: Better root found: R2:%f\tclock_rate:%.3e\tbranch_displacement:%f"
                            %(best_root._R2, best_root._beta, (best_root._R2_delta_x) / ( best_root.branch_length + self.one_mutation)),4)
            elif criterium in ['residual', 'min_dev']:
                if node.up is None:
                    self.logger("TreeTime.find_best_root_and_regression: Initial root: residual:%.3e\tclock_rate:%.3e"%(best_root._residual, best_root._beta),3)
                if  node._residual < best_root_any._residual:
                    best_root_any = node
                if (node._residual < best_root._residual and node._beta>0) or best_root._beta<0:
                    best_root = node
                    self.logger("TreeTime.find_best_root_and_regression: Better root found: residual:%.3e\tclock_rate:%.3e\tbranch_displacement:%f"
                            %(best_root._residual, best_root._beta, (best_root._R2_delta_x) / ( best_root.branch_length + self.one_mutation)),4)


        if criterium=='rsq':
            if (best_root_any._R2 > best_root._R2):
                self.logger("WARNING: TreeTime.find_best_root_and_regression: optimal regression has negative rate: R2:%f\tclock_rate:%.3e"
                        %(best_root_any._R2, best_root_any._beta), 1)
            self.logger("TreeTime.find_best_root_and_regression: Best root: R2:%f\tclock_rate:%.3e\tbranch_displacement:%f"
                        %(best_root._R2, best_root._beta, (best_root._R2_delta_x) / ( best_root.branch_length + 0.001*self.one_mutation)),3)
        elif criterium=='residual':
            if (best_root_any._residual < best_root._residual):
                self.logger("WARNING: TreeTime.find_best_root_and_regression: optimal regression has negative rate: residual:%.3e\tclock_rate:%.3e"
                        %(best_root_any._residual, best_root_any._beta), 1)
            self.logger("TreeTime.find_best_root_an_R2_delta_xd_regression: Best root: residual:%.3e\tclock_rate:%.3e\tbranch_displacement:%f"
                        %(best_root._residual, best_root._beta, (best_root._R2_delta_x) / ( best_root.branch_length + 0.001*self.one_mutation)),3)
        elif criterium=='min_dev':
            if (best_root_any._residual < best_root._residual):
                self.logger("WARNING: TreeTime.find_best_root_and_regression: optimal regression has negative rate",1)
                best_root = best_root_any
            self.logger("TreeTime.find_best_root_an_R2_delta_xd_regression: Best root: residual:%.3e\tclock_rate:%.3e\tbranch_displacement:%f"
                        %(best_root._residual, best_root._beta, (best_root._R2_delta_x) / ( best_root.branch_length + 0.001*self.one_mutation)),3)

        return best_root, best_root._alpha, best_root._beta


    def reroot_to_best_root(self,infer_gtr = False, criterium='rsq', **kwarks):
        '''
        determine the node that, when the tree is rooted on this node, results
        in the best regression of temporal constraints and root to tip distances

        Parameters
        ----------

         infer_gtr : bool
            Should infer new GTR model after re-root?

        '''
        from Bio import Phylo
        self.logger("TreeTime.reroot_to_best_root: searching for the best root position...",2)
        best_root, a, b = self.find_best_root_and_regression(criterium=criterium)
        # first, re-root the tree

        if hasattr(best_root, "_R2_delta_x") and  best_root._R2_delta_x > 0 and best_root.up is not None:

            # create new node in the branch and root the tree to it
            new_node = Phylo.BaseTree.Clade()

            # insert the new node in the middle of the branch
            # by simple re-wiring the links on the both sides of the branch
            # and fix the branch lengths
            new_node.branch_length = best_root.branch_length - best_root._R2_delta_x
            new_node.up = best_root.up
            new_node._alpha = a
            new_node._beta = b
            new_node.clades = [best_root]
            if hasattr(best_root, "_R2"):
                new_node._R2 = best_root._R2
            if hasattr(best_root, "_residual"):
                new_node._residual = best_root._residual
            new_node.up.clades = [k if k != best_root else new_node
                                  for k in best_root.up.clades]

            best_root.branch_length = best_root._R2_delta_x
            best_root.up = new_node
            self.logger("TreeTime.reroot_to_best_root:"
                        " branch length of children of new root:"
                        " %.3e, %.3e"%(best_root.branch_length,
                                       new_node.branch_length),2)
            return new_node
        else:
            # simply use the existing node as the new root
            return best_root


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
        ttutils.tree_layout(tt.tree)
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


if __name__=="__main__":
    pass


