from __future__ import print_function, division
import time
import config as ttconf
from Bio import Phylo
from Bio import AlignIO
import numpy as np
from gtr import GTR
import seq_utils
try:
    from itertools import izip
except ImportError:  #python3.x
    izip = zip

min_branch_length = 1e-3


class TreeAnc(object):
    """
    Class defines simple tree object with basic interface methods: reading and
    saving from/to files, initializing leaves with sequences from the
    alignment, making ancestral state inferrence
    """

    def __init__(self, tree=None, aln=None, gtr=None, verbose = ttconf.VERBOSE):
        if tree is None:
            raise("TreeAnc requires a tree!")
        self.t_start = time.time()
        self.verbose = verbose
        self.logger("TreeAnc: set-up",1)
        self._internal_node_count = 0
        self.one_mutation = None
        # TODO: set explicitly
        self.ignore_gaps = True
        if gtr is not None:
            self.gtr = gtr
        if aln is not None:
            self.aln = aln
        self.tree = tree
        if self.tree is None:
            self.logger("TreeAnc: tree loading failed! exiting",0)
            return
        if aln is not None:
            self.attach_sequences_to_nodes()

    def logger(self, msg, level, warn=False):
        if level<self.verbose or warn:
            dt = time.time() - self.t_start
            outstr = '\n' if level<2 else ''
            outstr+=format(dt, '4.2f')+'\t'
            outstr+= level*'-'
            outstr+=msg
            print(outstr)


####################################################################
## SET-UP
####################################################################
    @property
    def leaves_lookup(self):
        return self._leaves_lookup

    @property
    def gtr(self):
        return self._gtr
    @gtr.setter
    def gtr(self, in_gtr):
        if type(in_gtr)==str:
            self._gtr = GTR.standard(model=in_gtr, logger=self.logger)
        elif isinstance(in_gtr, GTR):
            self._gtr = in_gtr
            self._gtr.logger=self.logger
        else:
            self.logger("TreeAnc.gtr_setter: can't interpret GTR model", 1, warn=True)

    @property
    def tree(self):
        return self._tree
    @tree.setter
    def tree(self, in_tree):
        from os.path import isfile
        if isinstance(in_tree, Phylo.BaseTree.Tree):
            self._tree = in_tree
        elif type(in_tree)==str and isfile(in_tree):
            self._tree=Phylo.read(in_tree, 'newick')
        else:
            self.logger('TreeAnc: could not load tree! input was '+in_tree,1)
            self._tree = None
            return

        for node in self._tree.find_clades():
            node.original_length = node.branch_length
            node.mutation_length = node.branch_length
        self.prepare_tree()

    @property
    def aln(self):
        return self._aln
    @aln.setter
    def aln(self,in_aln):
        # load alignment from file if necessary
        from os.path import isfile
        from Bio.Align import MultipleSeqAlignment
        if isinstance(in_aln, MultipleSeqAlignment):
            self._aln = in_aln
        elif type(in_aln)==str and isfile(in_aln):
            self._aln=AlignIO.read(in_aln, 'fasta')

        if hasattr(self, '_tree'):
            self.attach_sequences_to_nodes()
        else:
            self.logger("TreeAnc.aln: sequences not yet attached to tree",3,warn=True)

    def attach_sequences_to_nodes(self):
        # loop over tree,
        failed_leaves= 0
        dic_aln = {k.name: seq_utils.seq2array(k.seq) for k in self.aln} #
        for l in self.tree.get_terminals():
            if l.name in dic_aln:
                l.state_seq = dic_aln[l.name]
                l.sequence=l.state_seq
            else:
                self.logger("TreeAnc.attach_sequences_to_nodes: Cannot find sequence for leaf: %s" % l.name, 4, warn=True)
                failed_leaves += 1
                if failed_leaves == 100:
                    self.logger("Error: cannot set sequences to the terminal nodes.\n", 2, warn=True)
                    self.logger("Are you sure the alignment belongs to the tree?", 2, warn=True)
                    break
        self.seq_len = self.aln.get_alignment_length()
        self.one_mutation = 1.0/self.seq_len


    def prepare_tree(self):
        """
        Set link to parent and net distance to root for all tree nodes.
        Should be run once the tree is read and after every tree topology or branch
        lengths optimizations.
        """
        if self.one_mutation is None:
            self.tree.root.branch_length = 0.001
        else:
            self.tree.root.branch_length = self.one_mutation
        self.tree.root.mutation_length = self.tree.root.branch_length
        self.tree.root.mutations = []
        self.tree.ladderize()
        self._prepare_nodes()
        self._leaves_lookup = {node.name:node for node in self.tree.get_terminals()}


    def _prepare_nodes(self):
        """
        Set auxilliary parameters to every node of the tree.
        """
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            if clade.name is None:
                clade.name = "NODE_" + format(self._internal_node_count, '07d')
                self._internal_node_count += 1
            for c in clade.clades:
                c.up = clade
                if not hasattr(c, 'mutation_length'):
                    c.mutation_length=c.branch_length
                c.dist2root = c.up.dist2root + c.mutation_length

####################################################################
## END SET-UP
####################################################################

    def infer_gtr(self, print_raw=False, **kwargs):

        self.logger("TreeAnc inferring the GTR model from the tree...", 1)
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
                    Ti[i] += 0.5*node.mutation_length
                    Ti[j] -= 0.5*node.mutation_length
                for nuc in node.sequence:
                    i = alpha.index(nuc)
                    Ti[i] += node.mutation_length
        if print_raw:
            print('alphabet:',alpha)
            print('n_ij:', nij)
            print('T_i:', Ti)
        root_state = np.array([np.sum(self.tree.root.sequence==nuc) for nuc in alpha])
        self._gtr = GTR.infer(nij, Ti, root_state, pc=5.0, alphabet=self.gtr.alphabet, logger=self.logger)
        return self._gtr


###################################################################
### ancestral reconstruction
###################################################################
    def reconstruct_anc(self, method='ml', infer_gtr=False, **kwargs):
        """
        Reconstruct ancestral states
        Args:
         - method(str): method to use. Supported values are "fitch" and "ml"

        Returns:
         - N_diff(int): number of nucleotides different from the previous
         reconstruction. If there were no pre-set sequences, returns N*L

        """

        self.logger("TreeAnc reconstructing ancestral states with method: "+method,1)

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


###################################################################
### FITCH
###################################################################
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

        self.logger("TreeAnc._fitch_anc: Walking up the tree, creating the Fitch profiles",2)
        for node in self.tree.get_nonterminals(order='postorder'):
            node.state = [self._fitch_state(node, k) for k in range(self.L)]

        ambs = [i for i in range(self.L) if len(self.tree.root.state[i])>1]
        if len(ambs) > 0:
            for amb in ambs:
                self.logger("Ambiguous state of the root sequence "
                                    "in the position %d: %s, "
                                    "choosing %s" % (amb, str(self.tree.root.state[amb]),
                                                     self.tree.root.state[amb][0]), 4)
        self.tree.root.sequence = np.array([k[np.random.randint(len(k)) if len(k)>1 else 0]
                                           for k in self.tree.root.state])



        self.logger("TreeAnc._fitch_anc: Walking down the self.tree, generating sequences from the "
                         "Fitch profiles.", 2)
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
        self.logger("Done ancestral state reconstruction",3)
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


###################################################################
### Maximum Likelihood
###################################################################
    def branch_length_to_gtr(self, node):
        return max(min_branch_length*self.one_mutation, node.mutation_length)

    def _ml_anc(self, marginal=False, verbose=0, store_compressed=True, **kwargs):
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
        self.logger("TreeAnc._ml_anc: type of reconstruction:"+ ('marginal' if marginal else "joint"), 2)
        self.logger("Walking up the tree, computing likelihoods... ", 3)
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
                    self.branch_length_to_gtr(ch), return_log=False) # raw prob to transfer prob up
                node.profile *= ch.seq_msg_to_parent
                node.lh_prefactor += ch.lh_prefactor

            pre = node.profile.max(axis=1) #sum over nucleotide states
            node.profile = (node.profile.T/pre).T # normalize so that the sum is 1
            node.lh_prefactor += np.log(pre) # and store log-prefactor


        self.logger("Walking down the tree, computing maximum likelihood sequences...",3)

        # extract the likelihood from the profile
        tree.root.profile *= np.diag(self.gtr.Pi) # Msg to the root from the distant part (equ frequencies)
        pre=tree.root.profile.sum(axis=1)
        tree.root.profile = (tree.root.profile.T/pre).T
        tree.root.lh_prefactor += np.log(pre)

        tree.anc_LH = tree.root.lh_prefactor.sum()
        tree.sequence_LH = 0
        # reset profile to 0-1 and set the sequence
        tree.root.sequence, tree.root.profile = \
            seq_utils.prof2seq(tree.root.profile, self.gtr, sample_from_prof=True, collapse_prof=not marginal)
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
                                                self.branch_length_to_gtr(node), return_log=False)
                node.profile *= node.seq_msg_from_parent
            else:
                node.seq_msg_from_parent = self.gtr.propagate_profile(node.up.profile,
                                                self.branch_length_to_gtr(node), return_log=False)
                node.profile *= node.seq_msg_from_parent

            # reset the profile to 0-1 and  set the sequence
            sequence, profile = seq_utils.prof2seq(node.profile, self.gtr, sample_from_prof=False, collapse_prof=not marginal)
            node.mutations = [(anc, pos, der) for pos, (anc, der) in
                            enumerate(izip(node.up.sequence, sequence)) if anc!=der]

            # this needs fixing for marginal reconstruction
            if not marginal:
                tree.sequence_LH += np.sum(np.log(node.seq_msg_from_parent[profile>0.9]))

            if hasattr(node, 'sequence') and node.sequence is not None:
                N_diff += (sequence!=node.sequence).sum()
            else:
                N_diff += L

            node.sequence = sequence
            node.profile = profile

        # note that the root doesn't contribute to N_diff (intended, since root sequence is often ambiguous)
        self.logger("TreeAnc._ml_anc: ...done",3)
        if store_compressed:
            self.store_compressed_sequence_pairs()
        return N_diff

    def store_compressed_sequence_to_node(self, node):
            seq_pairs, multiplicity = self.gtr.compress_sequence_pair(node.up.sequence,
                                                                      node.sequence,
                                                                      ignore_gaps = self.ignore_gaps)
            node.compressed_sequence = {'pair':seq_pairs, 'multiplicity':multiplicity}

    def store_compressed_sequence_pairs(self):
        self.logger("TreeAnc.store_compressed_sequence_pairs...",2)
        for node in self.tree.find_clades():
            if node.up is None:
                continue
            self.store_compressed_sequence_to_node(node)


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
                transition_matrix = self.gtr.expQt(self.branch_length_to_gtr(child))
                if child.is_terminal():
                    from_children=child.profile
                else:
                    from_children = np.prod([c.seq_msg_to_parent for c in child], axis=0)
                to_parent = np.prod([node.seq_msg_from_parent] +
                                [c.seq_msg_to_parent for c in node if c!=child], axis=0)

                child.mutation_matrix=np.array([get_two_point_func(upmsg, downmsg, transition_matrix)
                                          for upmsg, downmsg in zip(to_parent,from_children)])


###################################################################
### Branch length
###################################################################
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


        self.logger("TreeAnc.optimize_branch_length: running branch lengths optimization...",1)

        verbose = 0
        store_old_dist = False

        if 'verbose' in kwargs:
            verbose = int(kwargs['verbose'])
        if 'store_old' in kwargs:
            store_old_dist = kwargs['store_old'] == True

        for node in self.tree.find_clades(order='postorder'):
            if node.up is None: continue # this is the root
            if store_old_dist:
                node._old_length = node.branch_length

            new_len = self.optimal_branch_length(node)

            if new_len < 0:
                continue

            self.logger("Optimization results: old_len=%.4f, new_len=%.4f "
                   " Updating branch length..."%(node.branch_length, new_len), 5)

            node.branch_length = new_len
            node.mutation_length=new_len

        # as branch lengths changed, the params must be fixed
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        self._prepare_nodes()


    def optimal_branch_length(self, node):
        '''
        calculate optimal branch length given the sequences of node and parent
        IMPORTANTLY: this needs to use sequences and 0-1 profiles!
        '''
        if node.up is None:
            return self.one_mutation

        parent = node.up
        if hasattr(node, 'compressed_sequence'):
            new_len = self.gtr.optimal_t_compressed(node.compressed_sequence['pair'],
                                                    node.compressed_sequence['multiplicity'])
        else:
            new_len = self.gtr.optimal_t(parent.sequence, node.sequence,
                                         ignore_gaps=self.ignore_gaps)
        return new_len


    def prune_short_branches(self):
        """
        If the branch length is less than the minimal value, remove the branch
        from the tree. **Requires** the ancestral sequence reconstruction
        """
        self.logger("TreeAnc.prune_short_branches: pruning short branches (max prob at zero)...", 1)
        for node in self.tree.find_clades():
            if node.up is None or node.is_terminal():
                continue

            # probability of the two seqs separated by zero time is not zero
            if self.gtr.prob_t(node.up.sequence, node.sequence, 0.0) > 0.1:
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
        self.logger("TreeAnc.optimize_seq_and_branch_len: ...", 1)
        if reuse_branch_len:
            N_diff = self.reconstruct_anc(method='ml', **kwargs)
        else:
            N_diff = self.reconstruct_anc(method='fitch', **kwargs)
        n = 0
        while True: # at least one cycle must be done
            n += 1

            self.optimize_branch_len(verbose=0, store_old=False)
            if prune_short:
                self.prune_short_branches()
            N_diff = self.reconstruct_anc(method='ml')

            self.logger("TreeAnc.optimize_seq_branch_length: Iteration %d."
                   " #Nuc changed since prev reconstructions: %d" %(n, N_diff), 2)

            if N_diff < 1:
                break
            elif n > 10:
                self.logger("sequences and branch lengths optimization did not"
                       "converge in 10 cycles, aborting.", 4, warn=True)
                break

        self._prepare_nodes() # fix dist2root and up-links after reconstruction
        self.logger("TreeAnc.optimize_seq_and_branch_len: Unconstrained sequence LH:%f"%self.tree.sequence_LH, 2)
        return

###############################################################################
### Utility functions
###############################################################################
    def get_reconstructed_alignment(self):
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        self.logger("TreeAnc.get_reconstructed_alignment ...",2)
        if not hasattr(self.tree.root, 'sequence'):
            self.logger("TreeAnc.reconstructed_alignment... reconstruction not yet done",3)
            self.reconstruct_anc('ml')

        new_aln = MultipleSeqAlignment([SeqRecord(id=n.name, seq=Seq("".join(n.sequence)), description="")
                                        for n in self.tree.find_clades()])

        return new_aln

if __name__=="__main__":
    pass
