import utils
import numpy as np
import config as ttconf
from treeanc import TreeAnc
from distribution import Distribution
from branch_len_interpolator import BranchLenInterpolator
from node_interpolator import NodeInterpolator


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

    def __init__(self,  dates=None,*args, **kwargs):
        super(ClockTree, self).__init__(*args, **kwargs)
        if dates is None:
            raise("ClockTree requires date contraints!")
        self.date_dict = dates
        self.date2dist = None  # we do not know anything about the conversion
        self.debug=False
        self.n_integral = ttconf.NINTEGRAL
        self.rel_tol_prune = ttconf.REL_TOL_PRUNE
        self.rel_tol_refine = ttconf.REL_TOL_REFINE
        self.merger_rate_default = ttconf.BRANCH_LEN_PENALTY

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
        else:
            self.logger("ClockTime.date2dist: Setting new date to branchlength conversion. slope=%f, R^2=%.4f"%(val.slope, val.r_val**2), 2)
            self._date2dist = val


    def init_date_constraints(self, ancestral_inference=False, slope=None, **kwarks):
        """
        Get the conversion coefficients between the dates and the branch
        lengths as they are used in ML computations. The conversion formula is
        assumed to be 'length = k*numdate_given + b'. For convenience, these
        coefficients as well as regression parameters are stored in the
        dates2dist object.

        Note: that tree must have dates set to all nodes before calling this
        function. (This is accomplished by calling load_dates func).

        Params:
            ancestral_inference: bool -- whether or not to reinfer ancestral sequences
                                 done by default when ancestral sequences are missing

        """
        self.logger("ClockTree.init_date_constraints...",2)

        if ancestral_inference or (not hasattr(self.tree.root, 'sequence')):
            self.infer_ancestral_sequences('ml',sample_from_profile='root',**kwarks)

        # set the None  for the date-related attributes in the internal nodes.
        # make interpolation objects for the branches
        self.logger('ClockTree.init_date_constraints: Initializing branch length interpolation objects...',3)
        for node in self.tree.find_clades():
            if node.up is None:
                node.branch_length_interpolator = None
            else:
                # copy the merger rate and gamma if they are set
                if hasattr(node,'branch_length_interpolator'):
                    gamma = node.branch_length_interpolator.gamma
                    merger_rate = node.branch_length_interpolator.merger_rate
                else:
                    gamma = 1.0
                    merger_rate = self.merger_rate_default
                node.branch_length_interpolator = BranchLenInterpolator(node, self.gtr, one_mutation=self.one_mutation)
                node.branch_length_interpolator.merger_rate = merger_rate
                node.branch_length_interpolator.gamma = gamma
        self.date2dist = utils.DateConversion.from_tree(self.tree, slope)

        # make node distribution objects
        for node in self.tree.find_clades():
            # node is constrained
            if hasattr(node, 'numdate_given') and node.numdate_given is not None:
                if hasattr(node, 'bad_branch') and node.bad_branch==True:
                    self.logger("ClockTree.init_date_constraints -- WARNING: Branch is marked as bad"
                                ", excluding it from the optimization process"
                                " Will be optimized freely", 4, warn=True)
                    node.numdate_given = None
                    node.time_before_present = None
                    # if there are no constraints - log_prob will be set on-the-fly
                    node.msg_to_parent = None
                else:
                    # set the absolute time before present in branch length units
                    node.time_before_present = self.date2dist.get_time_before_present(node.numdate_given)
                    node.msg_to_parent = NodeInterpolator.delta_function(node.time_before_present, weight=1)

            else: # node without sampling date set
                node.numdate_given = None
                node.time_before_present = None
                # if there are no constraints - log_prob will be set on-the-fly
                node.msg_to_parent = None


    def make_time_tree(self):
        '''
        use the date constraints to calculate the most likely positions of
        unconstraint nodes.
        '''
        self.logger("ClockTree: Maximum likelihood tree optimization with temporal constraints:",1)
        self.init_date_constraints()
        self._ml_t_leaves_root()
        self._ml_t_root_leaves()
        self._set_final_dates()
        self.convert_dates()


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
                res =  NodeInterpolator.convolve(node.msg_to_parent,
                            node.branch_length_interpolator, n_integral=self.n_integral,
                            rel_tol=self.rel_tol_refine)
            self.logger("ClockTree._ml_t_leaves_root._send_message: "
                        "computed convolution with %d points at node %s"%(len(res.x),node.name),4)
            return res

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
                node.msg_to_parent._adjust_grid(rel_tol=self.rel_tol_prune)


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
        self.logger("ClockTree: Propagating root -> leaves...", 2)
        # Main method - propagate from root to the leaves and set the LH distributions
        # to each node
        for node in self.tree.find_clades(order='preorder'):  # ancestors first, msg to children
            ## This is the root node
            if node.up is None:
                node.msg_from_parent = None # nothing beyond the root
            elif node.msg_to_parent.is_delta:
                node.msg_from_parent = None
            else:
                parent = node.up
                complementary_msgs = [parent.msgs_from_leaves[k]
                                      for k in parent.msgs_from_leaves
                                      if k != node]

                if parent.msg_from_parent is not None: # the parent is not root => got something from the parent
                    complementary_msgs.append(parent.msg_from_parent)

                msg_parent_to_node = NodeInterpolator.multiply(complementary_msgs)
                msg_parent_to_node._adjust_grid(rel_tol=self.rel_tol_prune)
                res = NodeInterpolator.convolve(msg_parent_to_node, node.branch_length_interpolator,
                                                inverse_time=False, n_integral=self.n_integral,
                                                rel_tol=self.rel_tol_refine)
                node.msg_from_parent = res
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
                                 node.branch_length_interpolator.y-node.branch_length_interpolator.peak_val, '-o')
                        plt.plot(msg_parent_to_node.x,msg_parent_to_node.y-msg_parent_to_node.peak_val, '-o')
                        plt.ylim(0,100)
                        plt.xlim(-0.01, 0.01)
                        import ipdb; ipdb.set_trace()


    def _set_final_dates(self):
        """
        Given the location of the node in branch length units, convert it to the
        date-time information.

        Args:
         - node(Phylo.Clade): tree node. NOTE the node should have the abs_t attribute
         to have a valid value. This is automatically taken care of in the
         procedure to get the node location probability distribution.

        """
        self.logger("ClockTree: Setting dates and node distributions...", 2)
        def collapse_func(dist):
            if dist.is_delta:
                return dist.peak_pos
            else:
                return dist.peak_pos


        for node in self.tree.find_clades(order='preorder'):  # ancestors first, msg to children
            # set marginal distribution
            ## This is the root node
            if node.up is None:
                node.marginal_lh = node.msg_to_parent
            elif node.msg_to_parent.is_delta:
                node.marginal_lh = node.msg_to_parent
            else:
                node.marginal_lh = NodeInterpolator.multiply((node.msg_from_parent, node.msg_to_parent))

            if node.up is None:
                node.joint_lh = node.msg_to_parent
                node.time_before_present = collapse_func(node.joint_lh)
                node.branch_length = self.one_mutation
            else:
                # shift position of parent node (time_before_present) by the branch length
                # towards the present. To do so, add branch length to negative time_before_present
                # and rescale the resulting distribution by -1.0
                res = Distribution.shifted_x(node.branch_length_interpolator, -node.up.time_before_present)
                res.x_rescale(-1.0)
                # multiply distribution from parent with those from children and determine peak
                if node.msg_to_parent is not None:
                    node.joint_lh = NodeInterpolator.multiply((node.msg_to_parent, res))
                else:
                    node.joint_lh = res
                node.time_before_present = collapse_func(node.joint_lh)

                node.branch_length = node.up.time_before_present - node.time_before_present
            node.clock_length = node.branch_length


    def convert_dates(self):
        from datetime import datetime, timedelta
        now = utils.numeric_date()
        for node in self.tree.find_clades():
            years_bp = self.date2dist.get_date(node.time_before_present)
            if years_bp < 0:
                if not hasattr(node, "bad_branch") or node.bad_branch==False:
                    self.logger("ClockTree.convert_dates -- WARNING: The node is later than today, but it is not"
                        "marked as \"BAD\", which indicates the error in the "
                        "likelihood optimization.",4 , warn=True)
                else:
                    self.logger("ClockTree.convert_dates -- WARNING: node which is marked as \"BAD\" optimized "
                        "later than present day",4 , warn=True)

            node.numdate = now - years_bp

            # set the human-readable date
            days = 365.25 * (node.numdate - int(node.numdate))
            year = int(node.numdate)
            try:
                n_date = datetime(year, 1, 1) + timedelta(days=days)
                node.date = datetime.strftime(n_date, "%Y-%m-%d")
            except:
                # this is the approximation
                n_date = datetime(1900, 1, 1) + timedelta(days=days)
                node.date = str(year) + "-" + str(n_date.month) + "-" + str(n_date.day)


    def branch_length_to_years(self):
        self.logger('ClockTree.branch_length_to_years: setting node positions in units of years', 2)
        if not hasattr(self.tree.root, 'numdate'):
            self.logger('ClockTree.branch_length_to_years: infer ClockTree first', 2,warn=True)
        self.tree.root.branch_length = 1.0
        for n in self.tree.find_clades(order='preorder'):
            if n.up is not None:
                n.branch_length = n.numdate - n.up.numdate


if __name__=="__main__":
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.ion()
    base_name = 'data/H3N2_NA_allyears_NA.20'
    #root_name = 'A/Hong_Kong/JY2/1968|CY147440|1968|Hong_Kong||H3N2/8-1416'
    root_name = 'A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409'
    with open(base_name+'.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    from Bio import Phylo
    tree = Phylo.read(base_name + ".nwk", 'newick')
    tree.root_with_outgroup([n for n in tree.get_terminals() if n.name==root_name][0])
    myTree = ClockTree(gtr='Jukes-Cantor', tree = tree,
                        aln = base_name+'.fasta', verbose = 6, dates = dates)

    myTree.optimize_seq_and_branch_len(prune_short=True)
    myTree.init_date_constraints()
    myTree.make_time_tree()

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

    fig, axs = plt.subplots(2,1, sharex=True, figsize=(8,12))
    x = np.linspace(-0.1,0.05,1000)+ myTree.tree.root.time_before_present
    Phylo.draw(tree, axes=axs[0], show_confidence=False)
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

    myTree.branch_length_to_years()
    Phylo.draw(myTree.tree)
    plt.xlim(myTree.tree.root.numdate-1,
             np.max([x.numdate for x in myTree.tree.get_terminals()])+1)
