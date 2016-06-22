from __future__ import print_function, division
from Bio import Phylo
from Bio import AlignIO
import   numpy as np
from gtr import GTR
import seq_utils
try:
    from itertools import izip
except ImportError:  #python3.x
    izip = zip

import json
from weakref import WeakKeyDictionary
from scipy.interpolate import interp1d


class _Descriptor_PhyloTree(object):
    """
    Descriptor to manage the Phylo.Tree object of the TreeTime class.
    When reading the new tree, it prepares the general properties of the tree for
    the further use. In particular, it sets the descriptors to the Clade object
    to manipulate some sensible properties, which depend on each other.
    """
    def __init__(self, default=Phylo.BaseTree.Tree()):
        self.default = default
        self.data = WeakKeyDictionary()

    def __get__(self, instance, owner):
        # we get here when someone calls treetime.tree and tree is a Descriptor instance
        # instance = treetime (real object)
        # owner = type(treetime) (should be TreeTime)
        return self.data.get(instance, self.default) # if the data not present - default will be returned

    def __set__(self, instance, value):
        # we get here when someone calls x.d = val, and d is a Descriptor instance
        # instance = x
        # value = val
        if value is None:
            self.data[instance] = self.default
            return

        if not isinstance(value, Phylo.BaseTree.Tree):
            raise TypeError("The TreeTime object can only assign the BioPython.BaseTree objects.")

        self.data[instance] = value


class TreeAnc(object):
    """
    Class defines simple tree object with basic interface methods: reading and
    saving from/to files, initializing leaves with sequences from the
    alignment, making ancestral state inferrence
    """

    tree = _Descriptor_PhyloTree() # Always set the descriptors at the class level

    class DisplayAttr(object):

        def __init__(self, name, attr):
            self._name = name
            self._attr = attr

        @property
        def name(self):
            return self._name

        def attr(self, node):

            if callable(self._attr):
                try:
                    return self._attr(node)
                except:
                    return ""
            elif isinstance(self._attr, str):
                if hasattr(node, self._attr):
                    return node.__dict__[self._attr]
                else:
                    return ""
            else:
                return  ""


    @property
    def leaves_lookup(self):
        return self._leaves_lookup

    @property
    def gtr(self):
        return self._gtr

    @gtr.setter
    def gtr_setter(self, value):
        pass



    def __init__(self, gtr):
        assert(isinstance(gtr, GTR))
        self.one_mutation = 1.0
        self._max_node_num = 0
        self._gtr = gtr
        self.tree = None
        self._leaves_lookup = {}
        self._internal_metadata_names = [
                    self.DisplayAttr("numdate", "numdate"),
                    self.DisplayAttr("mutation_rate/avg", "gamma"),
                    self.DisplayAttr("branch_len/opt", lambda n: (abs(n.branch_length - n.branch_neg_log_prob.x[(n.branch_neg_log_prob.y.argmin())]) / self.one_mutation)),
                    self.DisplayAttr("time_since_MRCA (yr)", "tvalue")
                ]
        self._terminal_metadata_names = self._internal_metadata_names

        # self.set_additional_tree_params()

    def infer_gtr(self, print_raw=False, **kwargs):
        self._ml_anc(**kwargs)
        alpha = list(self.gtr.alphabet)
        n=len(alpha)
        nij = np.zeros((n,n))
        Ti = np.zeros(n)
        for node in self.tree.find_clades():
            if hasattr(node,'mutations'):
                for a,pos, d in node.mutations:
                    i,j = alpha.index(a), alpha.index(d)
                    nij[i,j]+=1
                    Ti[i] += 0.5*node.branch_length
                    Ti[j] -= 0.5*node.branch_length
                for nuc in node.sequence:
                    i = alpha.index(nuc)
                    Ti[i]+=node.branch_length
        if print_raw:
            print('alphabet:',alpha)
            print('n_ij:', nij)
            print('T_i:', Ti)
        root_state = np.array([np.sum(self.tree.root.sequence==nuc) for nuc in alpha])
        self._gtr = GTR.infer(nij, Ti, root_state, pc=5.0, alphabet=self.gtr.alphabet)
        return self._gtr

    def set_additional_tree_params(self):
        """
        Set link to parent and net distance to root for all tree nodes.
        Should be run once the tree is read and after every tree topology or branch
        lengths optimizations.
        """
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        self._set_each_node_params()
        self._leaves_lookup = {node.name:node for node in self.tree.get_terminals()}


    def set_metadata_to_node(self, node, **metadata):
        """
        Set the metadata to the given tree node from the given dictionary.

        Args:
         - node(Phylo.Clade): node the metadata should be assigned to

        KWargs:
         - metadata: dictionary for the values to be set as attributes.

        Returns:
         - None
        """

        if isinstance(node, Phylo.BaseTree.Clade):

            for key in metadata:
                if key != "name": #  filter name node if any  (must be already set)
                    setattr(node, key, metadata[key])



        elif isinstance(node, str):

            if node not in  self._leaves_lookup:
                print ("Cannot set metadata to the node: node not found")
                return

            node = self._leaves_lookup[node]
            for key in metadata:
                if key != "name": #  filter name node if any  (must be already set)
                    setattr(node, key, metadata[key])




        else:
            print ("Cannot set metadata to node. Input node must be "
                "either tree node instance, or name of a node.")

    def set_metadata(self, **all_metadata):
        """
        Set metadata from dictionary to all nodes
        """

        metadata_list_set = False

        for node_key in all_metadata:
            if node_key not in self._leaves_lookup:
                print ("Cannot set metadata to the tree node: node name not found")
                print (node_key)
                continue

            self.set_metadata_to_node(node_key, **all_metadata[node_key])
            if not metadata_list_set:
                self._terminal_metadata_names = [
                        self.DisplayAttr(k, k) for k in all_metadata[node_key]
                    ]
                metadata_list_set = True



    def _set_each_node_params(self):

        """
        Set auxilliary parameters to every node of the tree.
        """
        self.tree.root.dist2root_0 = 0.0

        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            for c in clade.clades:
                c.up = clade
                if c.up.name is None:
                    c.up.name = "NODE_" + format(self._max_node_num, '07d')
                    self._max_node_num += 1
                c.dist2root = c.up.dist2root + c.branch_length
                c.dist2root_0 = c.dist2root #  store the values used later for date-branchLen conversion
        return

    def reconstruct_anc(self, method, infer_gtr=False, **kwargs):
        """
        Reconstruct ancestral states
        Args:
         - method(str): method to use. Supported values are "fitch" and "ml"
        KWargs:
         - model(TMat): model to use. required for maximum-likelihood ("ml")

        Returns:
         - N_diff(int): number of nucleotides different from the previous
         reconstruction. If there were no pre-set sequences, returns N*L

        """
        if infer_gtr:
            self.infer_gtr(**kwargs)
            N_diff = self._ml_anc(**kwargs)
        else:
            if method == 'fitch':
                N_diff = self._fitch_anc(**kwargs)
            elif method == 'ml':
                N_diff = self._ml_anc(**kwargs)
            else:
                raise NotImplementedError("The reconstruction method %s is not supported. " % method)
        return N_diff


    def _fitch_anc(self, **kwargs):
        """
        Reconstruct ancestral states using Fitch's algorithm. The method requires
        sequences to be assigned to leaves. It implements the iteration from
        leaves to the root constructing the Fitch profiles for each character of
        the sequence, and then by propagating from the root to the leaves,
        reconstructs the sequences of the internal nodes.

        KWargs:
         -

        Returns:
         - Ndiff (int): number of the characters that changed since the previous
         reconstruction. These changes are determined from the pre-set sequence attributes
         of the nodes. If there are no sequences available (i.e., no reconstruction
         has been made before), returns the total number of characters in the tree.

        """
        # set fitch profiiles to each terminal node
        for l in self.tree.get_terminals():
            l.state = [[k] for k in l.sequence]

        print ("Walking up the tree, creating the Fitch profiles")
        for node in self.tree.get_nonterminals(order='postorder'):
            node.state = [self._fitch_state(node, k) for k in range(self.L)]

        ambs = [i for i in range(self.L) if len(self.tree.root.state[i])>1]
        if len(ambs) > 0:
            for amb in ambs:
                print ("Ambiguous state of the root sequence "
                                    "in the position %d: %s, "
                                    "choosing %s" % (amb, str(self.tree.root.state[amb]),
                                                     self.tree.root.state[amb][0]))
        self.tree.root.sequence = np.array([k[np.random.randint(len(k)) if len(k)>1 else 0]
                                           for k in self.tree.root.state])



        print ("Walking down the self.tree, generating sequences from the "
                         "Fitch profiles.")
        N_diff = 0
        for node in self.tree.get_nonterminals(order='preorder'):
            if node.up != None: # not root
                sequence =  np.array([node.up.sequence[i]
                        if node.up.sequence[i] in node.state[i]
                        else node.state[i][0] for i in range(self.L)])
                if hasattr(node, 'sequence'):
                    N_diff += (sequence!=node.sequence).sum()
                else:
                    N_diff += self.L
                node.sequence = sequence

            node.profile = seq_utils.seq2prof(node.sequence, self.gtr.profile_map)
            del node.state # no need to store Fitch states
        print ("Done ancestral state reconstruction")
        for node in self.tree.get_terminals():
            node.profile = seq_utils.seq2prof(node.sequence, self.gtr.profile_map)
        return N_diff

    def _fitch_state(self, node, pos):
        """
        Determine the Fitch profile for a single character of the node's sequence.
        The profile is essentially the intersection between the children's
        profiles or, if the former is empty, the union of the profiles.

        Args:
         - node (Phylo.Node) internal node which the profiles are to be
         determined

         - pos (int): position in the node's sequence which the profiles should
         be determinedf for.

        Return:
         - state(numpy array): Fitch profile for the character at position pos
         of the given node.
        """
        state = self._fitch_intersect([k.state[pos] for k in node.clades])
        if len(state) == 0:
            state = np.concatenate([k.state[pos] for k in node.clades])
        return state

    def _fitch_intersect(self, arrays, assume_unique=False):
        """
        Find the intersection of any number of 1D arrays.
        Return the sorted, unique values that are in all of the input arrays.
        Adapted from numpy.lib.arraysetops.intersect1d
        """
        N = len(arrays)
        arrays = list(arrays) # allow assignment
        if not assume_unique:
            for i, arr in enumerate(arrays):
                arrays[i] = np.unique(arr)
        aux = np.concatenate(arrays) # one long 1D array
        aux.sort() # sorted
        shift = N-1
        # if an element is in all N arrays, is shows up N consecutive times in the sorted
        # concatenation. those elements can be found by comparing the array shifted by N-1
        # since the initital arrays are unique, only the correct elements are found this way.
        return aux[aux[shift:] == aux[:-shift]]

    def _ml_anc(self, marginal=False, verbose=0, **kwargs):
        """
        Perform ML reconstruction of the ancestral states
        KWargs:
         - store_lh (bool): if True, all likelihoods will be stored for all nodes.
           Useful for testing, diagnostics and if special post-processing is required.
         - verbose (int): how verbose the output should be
        """
        tree = self.tree
        # number of nucleotides changed from prev reconstruction
        N_diff = 0
        if 'store_lh' in kwargs:
            store_lh = kwargs['store_lh'] == True

        L = tree.get_terminals()[0].sequence.shape[0]
        n_states = self.gtr.alphabet.shape[0]
        if verbose > 2:
            print ("Walking up the tree, computing likelihoods... type of reconstruction:", 'marginal' if marginal else "joint")
        for leaf in tree.get_terminals():
            # in any case, set the profile
            leaf.profile = seq_utils.seq2prof(leaf.sequence, self.gtr.profile_map)
            leaf.lh_prefactor = np.zeros(L)
        for node in tree.get_nonterminals(order='postorder'): #leaves -> root
            # regardless of what was before, set the profile to ones
            node.lh_prefactor = np.zeros(L)
            node.profile = np.ones((L, n_states)) # we will multiply it
            for ch in node.clades:
                ch.seq_msg_to_parent = self.gtr.propagate_profile(ch.profile,
                    ch.branch_length,
                    rotated=False, # use unrotated
                    return_log=False) # raw prob to transfer prob up
                node.profile *= ch.seq_msg_to_parent
                node.lh_prefactor += ch.lh_prefactor
            pre = node.profile.sum(axis=1) #sum over nucleotide states

            node.profile = (node.profile.T/pre).T # normalize so that the sum is 1
            node.lh_prefactor += np.log(pre) # and store log-prefactor
        if (verbose > 2):
            print ("Walking down the tree, computing maximum likelihood sequences...")

        # extract the likelihood from the profile
        tree.root.lh_prefactor += np.log(tree.root.profile.max(axis=1))
        tree.anc_LH = tree.root.lh_prefactor.sum()
        tree.sequence_LH = 0
        # reset profile to 0-1 and set the sequence
        tree.root.profile *= np.diag(self.gtr.Pi) # Msg to the root from the distant part (equ frequencies)
        tree.root.sequence, tree.root.profile = \
            seq_utils.prof2seq(tree.root.profile, self.gtr, correct_prof=not marginal)
        tree.root.seq_msg_from_parent = np.repeat([self.gtr.Pi.diagonal()], len(tree.root.sequence), axis=0)

        for node in tree.find_clades(order='preorder'):
            if node.up is None: # skip if node is root
                continue
            # integrate the information coming from parents with the information
            # of all children my multiplying it to the prev computed profile
            if marginal:
                tmp_msg = np.copy(node.up.seq_msg_from_parent)
                for c in node.up.clades:
                    if c != node:
                        tmp_msg*=c.seq_msg_to_parent
                node.seq_msg_from_parent = self.gtr.propagate_profile(tmp_msg,
                            node.branch_length,
                            rotated=False, # use unrotated
                            return_log=False)
                node.profile *= node.seq_msg_from_parent
            else:
                node.seq_msg_from_parent = self.gtr.propagate_profile(node.up.profile,
                            node.branch_length,
                            rotated=False, # use unrotated
                            return_log=False)
                node.profile *= node.seq_msg_from_parent

            # reset the profile to 0-1 and  set the sequence
            sequence, profile = seq_utils.prof2seq(node.profile, self.gtr, correct_prof=not marginal)
            node.mutations = [(anc, pos, der) for pos, (anc, der) in
                            enumerate(izip(node.up.sequence, sequence)) if anc!=der]

            # this needs fixing for marginal reconstruction
            if not marginal:
                tree.sequence_LH += np.sum(np.log(node.seq_msg_from_parent[profile>0.9]))
            if hasattr(node, 'sequence') and node.sequence is not None:
                try:
                    N_diff += (sequence!=node.sequence).sum()
                except:
                    import ipdb; ipdb.set_trace()
            else:
                N_diff += L
            node.sequence = sequence
            node.profile = profile
        return N_diff


    def calc_branch_twopoint_functions(self):
        '''
        attaches L x n x n (L sequence length, n alphabet size) to each branch containing
        the two-node profile, i.e. the probability that the parent node is in state
        s1 and the child node is in state s2. This allows probabilistic counting of mutations.
        '''
        def get_two_point_func(p1,p2,T):
            tmp = np.outer(p1, p2)*T
            tmp/=tmp.sum()
            return tmp

        for node in self.tree.get_nonterminals():
            for child in node:
                transition_matrix = self.gtr.v.dot(np.dot(self.gtr._exp_lt(child.branch_length), self.gtr.v_inv))
                if child.is_terminal():
                    from_children=child.profile
                else:
                    from_children = np.prod([c.seq_msg_to_parent for c in child], axis=0)
                to_parent = np.prod([node.seq_msg_from_parent] +
                                [c.seq_msg_to_parent for c in node if c!=child], axis=0)

                child.mutation_matrix=np.array([get_two_point_func(upmsg, downmsg, transition_matrix)
                                          for upmsg, downmsg in zip(to_parent,from_children)])

    def optimize_branch_len(self, **kwargs):
        """
        Perform ML optimization for the branch lengths of the whole tree or any
        subtree. **Note** this method assumes that each node stores information
        about its sequence as numpy.array object (node.sequence attribute).
        Therefore, before calling this method, sequence reconstruction with
        either of the available models must be performed.

        KWargs:
         - verbose (int): output detalization
         - store_old (bool): if True, the old lenths will be saved in
         node._old_dist attribute. Useful for testing, and special post-processing.
        Returns:
         - None, the phylogenetic tree is modified in-place.
        """

        verbose = 0
        store_old_dist = False

        if 'verbose' in kwargs:
            verbose = int(kwargs['verbose'])
        if 'store_old' in kwargs:
            store_old_dist = kwargs['store_old'] == True

        if verbose > 3:
            print ("Walking up the tree, computing likelihood distributions")

        for node in self.tree.find_clades(order='postorder'):
            parent = node.up
            if parent is None: continue # this is the root
            prof_p = parent.profile
            prof_ch = node.profile

            if store_old_dist:
                node._old_length = node.branch_length

            # optimization method
            #import ipdb; ipdb.set_trace()
            new_len = self.gtr.optimal_t(prof_p, prof_ch) # not rotated profiles!
            if new_len < 0:
                continue

            if verbose > 5:
                print ("Optimization results: old_len=%.4f, new_len=%.4f "
                        " Updating branch length..."
                        %(node.branch_length, new_len))

            node.branch_length = new_len

        # as branch lengths changed, the params must be fixed
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        self._set_each_node_params()
        return

    def prune_short_branches(self):
        """
        If the branch length is less than the minimal value, remove the branch
        from the tree. **Requires** the ancestral sequence reconstruction
        """
        print("pruning short branches (max prob at zero)")
        for node in self.tree.find_clades():
            if node.up is None:
                continue
            # probability of the two seqs separated by zero time is not zero
            if self.gtr.prob_t(node.up.profile, node.profile, 0.0) > 0.1:
                if node.is_terminal(): # leaf stays as is
                    continue
                # re-assign the node children directly to its parent
                node.up.clades = [k for k in node.up.clades if k != node] + node.clades
                if hasattr(node, "lh_prefactor"):
                    node.up.lh_prefactor += node.lh_prefactor
                for clade in node.clades:
                    clade.up = node.up

    def optimize_seq_and_branch_len(self,reuse_branch_len=True,prune_short=True, **kwargs):
        """
        Iteratively set branch lengths and reconstruct ancestral sequences until
        the values of either former or latter do not change. The algorithm assumes
        knowing only the topology of the tree, and requires that sequences are assigned
        to all leaves of the tree. The first step is to pre-reconstruct ancestral
        states using Fitch reconstruction algorithm or ML using existing branch length
        estimates. Then, optimize branch lengths and re-do reconstruction until
        convergence using ML method.

        Args:
         - reuse_branch_len(bool, default True): if True, rely on the initial
         branch lenghts, and start with the Maximum-likelihood ancestral sequence
         inference using existing branch lengths.
         Otherwise, initial reconstruction of ancestral states with Fitch algorithm,
         which uses only the tree topology.

         - prune_short (bool, default True): If True, the branches with zero
         optimal length will be pruned from the tree hence creating polytomies.
         The polytomies could be further processde using resolve_polytomies from
         the TreeTime class.
        """

        if reuse_branch_len:
            N_diff = self.reconstruct_anc('ml', **kwargs)
        else:
            N_diff = self.reconstruct_anc(method='fitch', **kwargs)
        n = 0
        while True: # at least one cycle must be done

            n += 1

            self.optimize_branch_len(verbose=0, store_old=False)
            if prune_short:
                self.prune_short_branches()
            N_diff = self.reconstruct_anc('ml')

            print ("Optimizing ancestral states and branch lengths. Round %d."
                   " #Nuc changed since prev reconstructions: %d" %(n, N_diff))

            if N_diff < 1:
                break

            if n > 100:
                print ("sequences and branch lengths optimization did not"
                       "converge in 100 cycles, aborting.")
                break

        self._set_each_node_params() # fix dist2root and up-links after reconstruction
        print("Unconstrained sequence LH:",self.tree.sequence_LH)
        return
