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


class TreeAnc(object):
    """
    Class defines simple tree object with basic interface methdos: reading and
    saving from/to files, initializing leaves with sequences from the
    alignment, making ancestral state inferrence
    """

    def __init__(self):
        self.tree = None
        # self.set_additional_tree_params()

    def set_additional_tree_params(self):
        """
        Set link to parent and net distance to root for all tree nodes.
        Should be run once the tree is read and after every tree topology or branch
        lengths optimizations.
        """
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        self._set_each_node_params()

    def _set_each_node_params(self):
        """
        Set auxilliary parameters to every node of the tree.
        """
        for clade in self.tree.get_nonterminals(order='preorder'): # up->down
            for c in clade.clades:
                c.up = clade
                c.dist2root = c.up.dist2root + c.branch_length
        return

    def reconstruct_anc(self, method, **kwargs):
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
        if method == 'fitch':
            N_diff = self._fitch_anc(**kwargs)
        elif method == 'ml':

            if ('model' not in kwargs):
                print("Warning: You chose Maximum-likelihood reconstruction,"
                    " but did not specified any model. Jukes-Cantor will be used as default.")
                gtr = GTR.standard(model='Jukes-Cantor')
            else:
                gtr = kwargs.pop('model')

            N_diff = self._ml_anc(gtr, **kwargs)
        else:
            raise NotImplementedError("The reconstruction method %s is not supported. " % method)

        return N_diff

    def _fitch_anc(self, **kwargs):
        """
        Reconstruct ancestral states using Fitch algorithm. The method reequires
        the leaves sequences to be assigned. It implements the iteration from
        leaves to the root constructing the Fitch profiles for each character of
        the sequence, and then by propagating from the root to the leaves,
        reconstructs the sequences of the internal nodes.

        KWargs:
         -

        Returns:
         - Ndiff (int): number of the nodes changed since the previous
         reconstruction. These changes are deermined from the pre-set sequence attributes
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

            node.profile = seq_utils.seq2prof(node.sequence)
            #if np.sum([k not in alphabet for k in node.sequence]) > 0:
            #    import ipdb; ipdb.set_trace()
            del node.state # no need to store Fitch states
        print ("Done ancestral state reconstruction")
        for node in self.tree.get_terminals():
            node.profile = seq_utils.seq2prof(node.sequence)
        return N_diff

    def _fitch_state(self, node, pos):
        """
        Determine the Fitch porfile for a single checracter of the node's sequence.
        The profile is essentially the  intersection between  the children's
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
        state = self._fitch_intersept([k.state[pos] for k in node.clades])
        if len(state) == 0:
            state = np.concatenate([k.state[pos] for k in node.clades])
        return state

    def _fitch_intersept(self, arrays, assume_unique=False):
        """
        Find the interseption of any number of 1D arrays.
        Return the sorted, unique values that are in all of the input arrays.
        Adapted from numpy.lib.arraysetops.intersept1d
        """
        N = len(arrays)
        arrays = list(arrays) # allow assignment
        if not assume_unique:
            for i, arr in enumerate(arrays):
                arrays[i] = np.unique(arr)
        aux = np.concatenate(arrays) # one long 1D array
        aux.sort() # sorted
        shift = N-1
        return aux[aux[shift:] == aux[:-shift]]

    def _ml_anc(self, gtr, **kwargs):
        """
        Perform ML reconstruction for the ancestral states
        Args:
         - model (GTR): General time-reversible model of evolution.
        KWargs:
         - store_lh (bool): if True, all likelihoods will be stored for all nodes. Useful for testing, diagnostics and if special post-processing is required.
         - verbose (int): how verbose the output should be
        """
        tree = self.tree
        # number of nucleotides changed from prev reconstruction
        N_diff = 0
        verbose = 0 # how verbose to be at the output
        if 'store_lh' in kwargs:
            store_lh = kwargs['store_lh'] == True
        if 'verbose' in kwargs:
            try:
                verbose = int(kwargs['verbose'])
            except:
                print ("ML ERROR in input: verbose param must be int")
        L = tree.get_terminals()[0].sequence.shape[0]
        a = gtr.alphabet.shape[0]
        if verbose > 2:
            print ("Walking up the tree, computing joint likelihoods...")
        for leaf in tree.get_terminals():
            # in any case, set the profile
            leaf.profile = seq_utils.seq2prof(leaf.sequence, gtr.alphabet_type)
            leaf.lh_prefactor = np.zeros(L)
        for node in tree.get_nonterminals(order='postorder'): #leaves -> root
            # regardless of what was before, set the profile to ones
            node.lh_prefactor = np.zeros(L)
            node.profile = np.ones((L, a)) # we will multiply it
            for ch in node.clades:
                ch.seq_msg_to_parent = gtr.propagate_profile(ch.profile,
                    ch.branch_length,
                    rotated=False, # use unrotated
                    return_log=False) # raw prob to transfer prob up
                node.profile *= ch.seq_msg_to_parent
                node.lh_prefactor += ch.lh_prefactor
            pre = node.profile.sum(axis=1) #sum over nucleotide states

            node.profile = (node.profile.T/pre).T # normalize so that the sum is 1
            node.lh_prefactor += np.log(pre) # and store log-prefactor
        if (verbose > 2):
            print ("Walking down the tree, computing maximum likelihood     sequences...")
        tree.root.profile *= np.diag(gtr.Pi) # Msg to the root from the distant part (equ frequcies)

        # extract the likelihood from the profile
        tree.root.lh_prefactor += np.log(tree.root.profile.max(axis=1))
        tree.anc_LH = tree.root.lh_prefactor.sum()
        tree.sequence_LH = 0
        # reset profile to 0-1 and set the sequence
        tree.root.sequence, tree.root.profile = \
            seq_utils.prof2seq(tree.root.profile, gtr, True)


        for node in tree.find_clades(order='preorder'):
            if node.up is None: # skip if node is root
                continue
            # integrate the information coming from parents with the information
            # of all children my multiplying it to the prev computed profile
            node.seq_msg_from_parent = gtr.propagate_profile(node.up.profile,
                            node.branch_length,
                            rotated=False, # use unrotated
                            return_log=False)
            node.profile *= node.seq_msg_from_parent

            # reset the profile to 0-1 and  set the sequence
            sequence, profile = seq_utils.prof2seq(node.profile, gtr, True)
            node.mutations = [(anc, pos, der) for pos, (anc, der) in
                            enumerate(izip(node.up.sequence, sequence)) if anc!=der]

            tree.sequence_LH += np.sum(np.log(node.seq_msg_from_parent[profile>0.9]))

            if hasattr(node, 'sequence'):
                N_diff += (sequence!=node.sequence).sum()
            else:
                N_diff += self.L
            node.sequence = sequence
            node.profile = profile
        return N_diff

    def optimize_branch_len(self, model, tree=None, **kwargs):
        """
        Perform ML optimization for the branch lengths of the whole tree or any
        subtree. **Note** this method assumes that each node stores information
        about its sequence as numpy.array object (node.sequence attribute).
        Therefore, before calling this method, sequence reconstruction with
        either of the available models must be performed.

        Args:
         - model(GTR): evolutionary model to be used for ML optimization of the
         branch lengths.
         - tree(None or Phylo.Clade): the root of the subtree to be optimized.
         if None, the optimization is being performed for the whole tree.

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
            print ("Walking up the tree, computing likelihood distrubutions")

        if tree is None:
            tree = self.tree

        for node in tree.find_clades(order='postorder'):
            parent = node.up
            if parent is None: continue # this is the root
            prof_p = parent.profile
            prof_ch = node.profile

            if store_old_dist:
                node._old_length = node.branch_length

            # optimization method
            #import ipdb; ipdb.set_trace()
            new_len = model.optimal_t(prof_p, prof_ch) # not rotated profiles!
            if new_len < 0:
                continue

            if verbose > 5:
                print ("Optimization results: old_len=%.4f, new_len=%.4f "
                        " Updating branch length..."
                        %(node.branch_length, new_len))

            node.branch_length = new_len

        # as branch lengths changed, the params must be fixed
        tree.root.up = None
        tree.root.dist2root = 0.0
        self._set_each_node_params()
        return

    def prune_short_branches(self,gtr):
        """
        If the branch length is less than the minimal value, remove the branch
        from the tree. **Requires** the ancestral seequence reconstruction
        """
        for node in self.tree.find_clades():
            if node.up is None:
                continue
            # probability of the two seqs separated by zero time is not zero
            if gtr.prob_t(node.up.profile, node.profile, 0.0) > 0.1:
                #FIXME: Why don't allow merging with the root?
                if node.is_terminal(): # or (node.up == self.tree.root): # leaf stays as is
                    continue
                # re-assign the node children directly to its parent
                node.up.clades = [k for k in node.up.clades if k != node] + node.clades
                if hasattr(node, "lh_prefactor"):
                    node.up.lh_prefactor += node.lh_prefactor
                for clade in node.clades:
                    clade.up = node.up

    def optimize_seq_and_branch_len(self,gtr,tree=None,reuse_branch_len=True,prune_short=True):
        """
        Iteratively set branch lengths and reconstruct ancestral sequences until
        the values of either former or latter do not change. The algorithm assumes
        knowing only the topology of the tree, and requires the  sequences assigned
        to all leaves of the tree. The first step is to pre-reconstruct ancestral
        states using Fitch reconstruction algorithm. Then, optimize branch lengths
        and re-do reconstruction until convergence using ML method.
        Args:

         - gtr(GTR): general time-reversible model to be used by every ML algorithm
        """

        if reuse_branch_len:
            N_diff = self.reconstruct_anc('ml', model=gtr)
        else:
            N_diff = self.reconstruct_anc(method='fitch')
        n = 0
        while True: # at least one cycle must be done

            n += 1

            self.optimize_branch_len(gtr, tree=tree, verbose=0, store_old=False)
            if prune_short:
                self.prune_short_branches(gtr)
            N_diff = self.reconstruct_anc('ml', model=gtr)

            print ("Optimizing ancestral states and branch lengths. Round %d."
                   " #Nuc changed since prev reconstructions: %d" %(n, N_diff))

            if N_diff < 1:
                break

            if n > 100:
                print ("sequences and branch lengths optimization did not"
                       "converge in 100 cycles, aborting.")

        self._set_each_node_params() # fix dist2root and up-links after reconstruction
        print("Unconstrained sequence LH:",self.tree.sequence_LH)
        return
