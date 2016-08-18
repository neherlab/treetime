"""
Class, which contains methods to optimize branch lengths given the time
constraints set to leaves
"""
from __future__ import print_function, division
from clock_tree import ClockTree
import config as ttconf
import numpy as np
from scipy import optimize as sciopt


class TreeTime(ClockTree):
    """
    TreeTime is a wrapper class to ClockTree that adds additional functionality
    such as reroot, detection and exclusion of outliers, resolution of polytomies
    using temporal information, and relaxed molecular clock models
    """
    def __init__(self, *args,**kwargs):
        super(TreeTime, self).__init__(*args, **kwargs)


    def run(self, root=None, infer_gtr=True, relaxed_clock=False, resolve_polytomies=True, max_iter=0, Tc=None):
        # initially, infer ancestral sequences and infer gtr model if desired
        self.optimize_seq_and_branch_len(infer_gtr=infer_gtr, prune_short=True)

        # optionally reroot the tree either by oldest, best regression or with a specific leaf
        if root is not None:
            self.reroot(root=root)

        # infer time tree and optionally resolve polytomies
        self.make_time_tree()
        if resolve_polytomies:
            # if polytomies are found, rerun the entire procedure
            if self.resolve_polytomies():
                self.prepare_tree()
                self.optimize_seq_and_branch_len(prune_short=False)
                if root=='best':
                    self.reroot(root=root)
                self.make_time_tree()

        # set root.gamma bc root doesn't have a branch_length_interpolator but gamma is needed
        if not hasattr(self.tree.root, 'gamma'):
            self.tree.root.gamma = 1.0
        # add coalescent prior
        if Tc is not None:
            from merger_models import coalescent
            self.logger('TreeTime.run: adding coalescent prior',2)
            coalescent(self.tree, Tc=Tc)
            self.make_time_tree()

        if relaxed_clock  and len(relaxed_clock)==2:
            slack, coupling = relaxed_clock
        # iteratively reconstruct ancestral sequences and re-infer
        # time tree to ensure convergence.
        niter = 0
        while niter<max_iter:
            if relaxed_clock:
                # estimate a relaxed molecular clock
                self.relaxed_clock(slack=slack, coupling=coupling)

            ndiff = self.reconstruct_anc('ml')
            if ndiff==0:
                break
            self.make_time_tree()
            niter+=1


    def reroot(self,root='best'):
        from Bio import Phylo
        if isinstance(root,Phylo.BaseTree.Clade):
            new_root = root
        elif root in self._leaves_lookup:
            new_root = self._leaves_lookup[root]
        elif root=='oldest':
            new_root = sorted([n for n in self.tree.get_terminals()
                               if n.numdate_given is not None],
                               key=lambda x:x.numdate_given)[0]
        elif root=='best':
            new_root = self.reroot_to_best_root()
        else:
            self.logger('TreeTime.reroot -- WARNING: unsupported rooting mechanisms or root not found',2,warn=True)
            return

        self.logger("TreeTime.reroot: Tree is being re-rooted to node "
                    +('new_node' if new_root.name is None else new_root.name), 2)
        self.tree.root_with_outgroup(new_root)
        self.tree.root.branch_length = self.one_mutation
        self.tree.root.numdate_given = None
        self.prepare_tree()


    def resolve_polytomies(self, merge_compressed=False, rerun=True):
        """
        Resolve the polytomies on the tree.
        The function scans the tree, resolves polytomies in case there are any,
        and re-optimizes the tree with new topology.
        Args:
            - merge_compressed(bool): whether to keep compressed branches as
              polytomies or return a strictly binary tree.
        """
        self.logger("TreeTime.resolve_polytomies: resolving multiple mergers...",2)

        poly_found=False
        for n in self.tree.find_clades():
            if len(n.clades) > 2:
                self._poly(n, merge_compressed)
                poly_found=True

        obsolete_nodes = [n for n in self.tree.find_clades() if len(n.clades)==1 and n.up is not None]
        for node in obsolete_nodes:
            self.logger('TreeTime.resolve_polytomies: remove obsolete node '+node.name,4)
            if node.up is not None:
                self.tree.collapse(node)

        return poly_found


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
        from branch_len_interpolator import BranchLenInterpolator
        from Bio import Phylo
        # TODO coefficient from the gtr
        zero_branch_slope = self.one_mutation / 0.8

        def _c_gain(t, n1, n2, parent):
            """
            cost gain if nodes n1, n2 are joined and their parent is placed at time t

            cost gain = (LH loss now) - (LH loss when placed at time t)
                      = [LH(opt) - LH(now)] - [LH(opt) - LH(t)] approx.=
                      approx.= LH(branch_len(now)) - LH (branch_len(t))

            """
            cg2 = n2.branch_length_interpolator(parent.time_before_present - n2.time_before_present) - n2.branch_length_interpolator(t - n2.time_before_present)
            cg1 = n1.branch_length_interpolator(parent.time_before_present - n1.time_before_present) - n1.branch_length_interpolator(t - n1.time_before_present)
            cg_new = - zero_branch_slope * (parent.time_before_present - t) # loss in LH due to the new branch
            return -(cg2+cg1+cg_new)

        def cost_gain(n1, n2, parent):
            """
            cost gained if the two nodes would have been connected.
            """
            cg = sciopt.minimize_scalar(_c_gain,
                    bounds=[max(n1.time_before_present,n2.time_before_present), parent.time_before_present],
                    method='Bounded',args=(n1,n2, parent))
            return cg['x'], - cg['fun']

        def merge_nodes(source_arr, isall=False):
            mergers = np.array([[cost_gain(n1,n2, clade) for n1 in source_arr]for n2 in source_arr])
            LH = 0
            while len(source_arr) > 1 + int(isall):
                # max possible gains of the cost when connecting the nodes:
                # this is only a rough approximation because it assumes the new node positions
                # to be optimal
                new_positions = mergers[:,:,0]
                cost_gains = mergers[:,:,1]
                np.fill_diagonal(cost_gains, -1e11)

                idxs = np.unravel_index(cost_gains.argmax(),cost_gains.shape)
                if (idxs[0] == idxs[1]) or cost_gains.max()<0:
                    if self.debug:
                        import ipdb; ipdb.set_trace()
                    else:
                        self.logger("TreeTime._poly.merge_nodes: problem merging child nodes of "+clade.name,4,warn=True)
                        continue

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
                new_node.sequence = clade.sequence
                new_node.profile = clade.profile
                self.store_compressed_sequence_to_node(new_node)

                new_node.mutations = []
                new_node.mutation_length = new_node.branch_length
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
        t_lh = self.tree.root.msg_to_parent.peak_val

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
        if slack is None: slack=ttconf.MU_ALPHA
        if coupling is None: coupling=ttconf.MU_BETA
        self.logger("TreeTime.relaxed_clock: slack=%f, coupling=%f"%(slack, coupling),2)

        c=1.0/self.one_mutation
        for node in self.tree.find_clades(order='postorder'):
            opt_len = node.mutation_length

            #opt_len = 1.0*len(node.mutations)/node.profile.shape[0]
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
                node._st_di   = np.sum([k._st_di + k._st_n_leaves*k.mutation_length for k in node.clades])
                node._st_diti = np.sum([k._st_diti + k.mutation_length*k._st_ti for k in node.clades])
                node._st_di2  = np.sum([k._st_di2 + 2*k._st_di*k.mutation_length + k._st_n_leaves*k.mutation_length**2 for k in node.clades])
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
                node._di = node.up._di + (n_up-n_down)*node.mutation_length
                node._di2 = (node.up._di2 + 2*node.mutation_length*node.up._di
                            - 4*(node.mutation_length*(node._st_di + n_down*node.mutation_length))
                            + N*node.mutation_length**2)
                node._diti = node.up._diti + node.mutation_length*(sum_ti - 2*node._st_ti)

                L = node.mutation_length

                ## Express Node's sum_Di as the function of parent's sum_Di
                # and **displacement from parent's node x** :
                # sum_Di = A1 + A2 * x
                A1 = node.up._di
                A2 = n_up - node._st_n_leaves

                ## Express Node's sum_Di**2 as the function of parent's params
                # and **displacement from parent's node x** :
                # sum_Di = B1 + B2 * x + B3 * x**2
                B1 = node.up._di2
                B2 = 2 * (node.up._di - 2 * node._st_di - 2 * node.mutation_length * node._st_n_leaves )
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
                self.logger("TreeTime.find_best_root_and_regression: Initial root: R2:%f\tslope:%f"%(best_root._R2, best_root._beta),3)
            elif (node._R2 > best_root._R2 and node._beta>0) or best_root._beta<0:
                best_root = node
                self.logger("TreeTime.find_best_root_and_regression: Better root found: R2:%f\tslope:%f\tbranch_displacement:%f"
                            %(best_root._R2, best_root._beta, (best_root._R2_delta_x) / ( best_root.mutation_length + self.one_mutation)),4)

        self.logger("TreeTime.find_best_root_and_regression: Best root: R2:%f\tslope:%f\tbranch_displacement:%f"
                    %(best_root._R2, best_root._beta, (best_root._R2_delta_x) / ( best_root.mutation_length + self.one_mutation)),3)
        return best_root, best_root._alpha, best_root._beta

    def reroot_to_best_root(self,infer_gtr = False, n_iqd = None, **kwarks):
        '''
        determine the node that, when the tree is rooted on this node, results
        in the best regression of temporal constraints and root to tip distances
        '''
        from Bio import Phylo
        self.logger("TreeTime.reroot_to_best_root: searching for the best root position...",2)

        best_root, a, b = self.find_best_root_and_regression()
        # first, re-root the tree

        if hasattr(best_root, "_R2_delta_x") and  best_root._R2_delta_x > 0 and best_root.up is not None:

            # create new node in the branch and root the tree to it
            new_node = Phylo.BaseTree.Clade()

            # insert the new node in the middle of the branch
            # by simple re-wiring the links on the both sides of the branch
            # and fix the branch lengths
            new_node.clades = [best_root]
            new_node.mutation_length = best_root.mutation_length - best_root._R2_delta_x
            best_root.mutation_length = best_root._R2_delta_x
            best_root.up.clades = [k if k != best_root else new_node for k in best_root.up.clades]
            new_node.branch_length = new_node.mutation_length
            return new_node
        else:
            # simply use the existing node as the new root
            return best_root

if __name__=="__main__":
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('whitegrid')
    from Bio import Phylo
    plt.ion()
    base_name = 'data/H3N2_NA_allyears_NA.20'
    with open(base_name+'.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    myTree = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 4, dates = dates)

    myTree.run(root='best', relaxed_clock=False, max_iter=1, resolve_polytomies=True, Tc=0.01) #(1.0,1.0), max_iter=1)

    plt.figure()
    x = np.linspace(0,0.05,100)
    leaf_count=0
    for node in myTree.tree.find_clades(order='postorder'):
        if node.up is not None:
            plt.plot(x, node.branch_length_interpolator.prob_relative(x))
        if node.is_terminal():
            leaf_count+=1
            node.ypos = leaf_count
        else:
            node.ypos = np.mean([c.ypos for c in node.clades])
    plt.yscale('log')
    plt.ylim([0.01,1.2])

    from matplotlib import cm
    fig, axs = plt.subplots(2,1, sharex=True, figsize=(8,12))
    x = np.linspace(-0.25,0.05,10000)+ myTree.tree.root.time_before_present
    for n in myTree.tree.find_clades():
        if n.up is None:
            g = n.gamma
        else:
            g = n.branch_length_interpolator.gamma
        n.color = [int(_x*255) for _x in cm.coolwarm(max(0,min(1,(g-1)+0.5)))[:3]]

    Phylo.draw(myTree.tree, axes=axs[0], show_confidence=False)
    offset = myTree.tree.root.time_before_present + myTree.tree.root.branch_length
    cols = sns.color_palette()
    depth = myTree.tree.depths()
    for ni,node in enumerate(myTree.tree.find_clades()):
        if (not node.is_terminal()):
            axs[1].plot(offset-x, node.marginal_lh.prob_relative(x), '-', c=cols[ni%len(cols)])
            axs[1].plot(offset-x, node.joint_lh.prob_relative(x), '--', c=cols[ni%len(cols)])
        if node.up is not None:
            x_branch = np.linspace(depth[node]-2*node.branch_length-0.005,depth[node],100)
            axs[0].plot(x_branch, node.ypos - 0.7*node.branch_length_interpolator.prob_relative(depth[node]-x_branch), '-', c=cols[ni%len(cols)])
    axs[1].set_yscale('log')
    axs[1].set_ylim([0.01,1.2])
    axs[0].set_xlabel('')
    plt.tight_layout()

    r2tip = np.array([[n.numdate, n.dist2root] for n in myTree.tree.get_terminals()])
    r2int = np.array([[n.numdate, n.dist2root] for n in myTree.tree.get_nonterminals()])
    plt.figure()
    plt.scatter(r2tip[:,0], r2tip[:,1], c='g', label='terminal nodes', s=30)
    plt.scatter(r2int[:,0], r2int[:,1], c='b', label='internal nodes', s=30)
    plt.ylabel('root to node distance')
    plt.xlabel('date')
    plt.legend(loc=2)

