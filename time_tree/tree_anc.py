from Bio import Phylo
from Bio import AlignIO
import numpy as np

# FIXME when reading the tree, scale all teh branches so that the maximal branch length is not more than some arbitrary value (needed for branch length optimization, not to make the brent algorithm crazy)

def _seq2idx(x, alph):
    """
    Get alphabet index of a character.
    Args:
     - x(char): character in the alphabet
     - alph(numpy.array): alphabet
    Returns:
     - idx(int): position of the x in the alph
    Throws:
     - ValueError: if the character was not found in the alphabet
    """
    if x not in alph:
        raise ValueError("character %s not found in the specified alphabet %s"
            % (x, alph))
    return (x == alph).argmax()

seq2idx = np.vectorize(_seq2idx)
seq2idx.excluded.add(1)

def seq2prof(x,alph):
    """
    Convert the given character into the profile.
    Args:
     - x(char): character in the alphabet
     - alph(numpy.array): alphabet
    Returns:
     - idx(numpy.array): profile for the character, zero array if the character not found
    """
    prof = np.zeros((alph.shape[0], x.shape[0]))
    for pos,char in enumerate(alph):
        prof[pos, :] = x==char
    return prof

class GTR(object):
    """
    Defines General tme reversible model of character evolution.
    """
    def __init__(self, alphabet):
        """
        Initialize empty evolutionary model.
        Args:
         - alphabet (numpy.array): alphabet of the sequence.
        """
        self.alphabet = alphabet
        # general rate matrix
        self.W = np.zeros((alphabet.shape[0], alphabet.shape[0]))
        # stationary states of the characters
        self.Pi = np.zeros((alphabet.shape[0], alphabet.shape[0]))
        # mutation rate, scaling factor
        self.mu = 1.0
        # eigendecomposition of the GTR matrix
        # Pi.dot(W) = v.dot(eigenmat).dot(v_inv)
        self.v = np.zeros((alphabet.shape[0], alphabet.shape[0]))
        self.v_inv = np.zeros((alphabet.shape[0], alphabet.shape[0]))
        self.eigenmat = np.zeros((alphabet.shape[0], alphabet.shape[0]))

    @classmethod
    def standard(cls, model='Jukes-Cantor', **kwargs):
        if model=='Jukes-Cantor':
            # read kwargs
            if 'alphabet' in kwargs:
                 alphabet = kwargs['alphabet']
            else:
                alphabet = np.array(['A', 'C', 'G', 'T'])
            if 'mu' in kwargs:
                mu = kwargs['mu']
            else:
                mu = 1.0

            gtr = cls(alphabet)
            gtr.W = np.ones((alphabet.shape[0], alphabet.shape[0]))
            np.fill_diagonal(gtr.W, - ((gtr.W).sum(0) - 1))
            gtr.Pi = np.zeros(gtr.W.shape)
            np.fill_diagonal(gtr.Pi, 0.25)
            sqrtPi = np.sqrt(gtr.Pi)
            sqrtPi_inv = np.linalg.inv(sqrtPi)
            W = (sqrtPi.dot(((gtr.Pi).dot(gtr.W)))).dot(sqrtPi_inv)
            eigvals, eigvecs = np.linalg.eig(W)
            gtr.v = sqrtPi.dot(eigvecs)
            gtr.v_inv = np.linalg.inv(gtr.v)
            gtr.eigenmat = np.diagflat(eigvals)
            gtr.mu = mu
            return gtr
        else:
            raise NotImplementedError("The specified evolutionary model is unsupported!")

class TreeAnc(object):
    """
    Class defines simple tree object with basic interface methdos: reading and
    saving from/to files, initializing leaves with sequences from the
    alignment, making ancestral state inferrence
    """

    def __init__(self, tree):
        self.tree = tree
        self._add_node_params()

    # input stuff
    @classmethod
    def from_file(cls, inf, fmt='newick'):
        if fmt == 'newick':
            return cls._from_newick(inf)
        elif fmt == 'json':
            return cls._from_json(inf)
        else:
            raise NotImplementedError("The specified format is unsupported")

    @classmethod
    def _from_newick(cls, inf):
        tree = Phylo.read(inf, 'newick')
        tanc = cls(tree)
        return tanc

    # FIXME
    @classmethod
    def _from_json(cls, inf):
        raise NotImplementedError("This functionality is under development")
        pass

    def _add_node_params(self):
        """
        Set link to parent and net distance to root for all tree nodes
        """
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            for c in clade.clades:
                c.up = clade
                c.dist2root = c.up.dist2root + c.branch_length
        return

    # ancestral state reconstruction
    def set_seqs_to_leaves(self, aln):
        """
        Set sequences from the alignment to the leaves of the tree.

        Args:
         - aln(Bio.MultipleSequenceAlignment): alignment ogbject

        Returns:
         - failed_leaves(int): number of leaves which could not be assigned with sequences.

        Note:
         - If there are more than 100 leaves failed to get sequences, the function breaks, returning 100.
        """
        failed_leaves= 0
        dic_aln = {k.name: np.array(k.seq) for k in aln} #
        for l in self.tree.get_terminals():
            if dic_aln.has_key(l.name):
                l.sequence = dic_aln[l.name]
            else:
                print ("Cannot find sequence for leaf: %s" % l.name)
                failed_leaves += 1
                if failed_leaves == 100:
                    print ("Error: cannot set sequences to the terminal nodes.\n"
                        "Are you sure the alignment belongs to the tree?")
                    break
        return failed_leaves

    def reconstruct_anc(self, method, **kwargs):
        """
        Reconstruct ancestral states
        Args:
         - method(str): method to use. Supported values are "fitch" and "ml"
        KWargs:
         - model(TMat): model to use. required for maximum-likelihood ("ml")
        """
        if method == 'fitch':
            self._fitch_anc(**kwargs)
        elif method == 'ml':
            if ('model' not in kwargs):
                print("Warning: You chose Maximum-likelihood reconstruction,"
                    " but did not specified any model. Jukes-Cantor will be used as default.")
                kwargs['model'] = GTR.standard(model='Jukes-Cantor')

            self._ml_anc(**kwargs)
        else:
            raise NotImplementedError("The reconstruction method %s is not supported. " % method)

    def _fitch_anc(self, **kwargs):
        """
        Reconstruct ancestral states using Fitch algorithm
        """
        # set fitch profiiles to each terminal node
        L = len(self.tree.get_terminals()[0].sequence)
        for l in self.tree.get_terminals():
            l.state = [[k] for k in l.sequence]

        print ("Walking up the tree, creating the Fitch profiles")
        for node in self.tree.get_nonterminals(order='postorder'):
            node.state = [self._fitch_state(node, k) for k in xrange(L)]

        ambs = [i for i in xrange(L) if len(self.tree.root.state[i])>1]
        if len(ambs) > 0:
            for amb in ambs:
                print ("Ambiguous state of the root sequence "
                                    "in the position %d: %s, "
                                    "choosing %s" % (amb, str(self.tree.root.state[amb]),
                                                     self.tree.root.state[amb][0]))
        self.tree.root.sequence = np.array([k[0] for k in self.tree.root.state])


        print ("Walking down the self.tree, generating sequences from the "
                         "Fitch profiles and filling the transition matrix.")
        n0 = 0
        for node in self.tree.get_nonterminals(order='preorder'):
            if node.up != None: # not root
               node.sequence =  np.array([node.up.sequence[i]
                       if node.up.sequence[i] in node.state[i]
                       else node.state[i][0] for i in xrange(L)])
            #if np.sum([k not in alphabet for k in node.sequence]) > 0:
            #    import ipdb; ipdb.set_trace()
            del node.state # no need to store Fitch states
        print ("Done ancestral state reconstruction")
        return

    def _fitch_state(self, node, pos):
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
        return aux[aux[shift:] == aux[:-shift]]

    def _ml_anc(self, model, **kwargs):
        """
        Perform ML reconstruction for the ancestral states
        Args:
         - model (GTR): General time-reversible model of evolution.
        KWargs:
         - store_lh (bool): if True, all likelihoods will be stored for all nodes. Useful for testing, diagnostics and if special post-processing is required.
         - verbose (int): how verbose the output should be
        """
        store_lh = False # store intermediate computations in the tree
        verbose = 1 # how verbose to be at the output
        if 'store_lh' in kwargs:
            store_lh = kwargs['store_lh'] == True
        if 'verbose' in kwargs:
            try:
                verbose = int(kwargs['verbose'])
            except:
                print ("ML ERROR in input: verbose param must be int")
                verbose = 5

        L = self.tree.get_terminals()[0].sequence.shape[0]
        alphabet = model.alphabet

        if verbose > 2:
            print ("Walking up the tree, computing joint likelihoods...")

        for leaf in self.tree.get_terminals():

            leaf.lh = np.ones((L, alphabet.shape[0]))
            leaf.pre = np.zeros(L)
            leaf.c = -1*np.ones((L, alphabet.shape[0]),dtype=int) # state of the node
            for pos,char in enumerate(alphabet):
                # FIXME should be index of the character everywhere
                leaf.c[leaf.sequence == alphabet[pos], :] = pos
            leaf.c[leaf.c.sum(1)==-4] = np.arange(alphabet.shape[0])

            eQT = np.diagflat(np.diag(np.exp(model.mu * leaf.branch_length * model.eigenmat)))
            P_all = (model.v).dot(eQT).dot(model.v_inv)
            leaf.lh[:] = P_all[leaf.c[:, 0], :]


        prob_profile = np.zeros((L, alphabet.shape[0], alphabet.shape[0]))
        for node in self.tree.get_nonterminals(order='postorder'): #leaves -> root

            node.lh = np.ones((L, alphabet.shape[0]))
            node.pre = np.zeros(L)
            node.c = np.zeros((L, alphabet.shape[0]),dtype=int) # state of the node
            if node.up is None: # we are at the root
                node.lh[:, :] = 1
                for ch in node.clades:
                    node.lh[:, :] *= ch.lh

            else: # internal node

                eQT = np.diagflat(np.diag(np.exp(model.mu * node.branch_length * model.eigenmat)))
                # probability to be in state i conditional on parent
                P_all = (model.v).dot(eQT).dot(model.v_inv)
                prob_profile[:] = P_all.T
                for ch in node.clades:
                    prob_profile *= ch.lh.reshape((L, len(alphabet), 1))
                node.lh[:, :] = prob_profile.max(1)
                node.c[:, :] = prob_profile.argmax(1)


        if (verbose > 2):
            print ("Walking down the tree, computing maximum likelihood     sequences...")

        self.tree.root.state = self.tree.root.lh.argmax(-1)
        self.tree.root.sequence = alphabet[self.tree.root.state]
        l_idx = np.arange(L)

        for node in self.tree.find_clades(order='preorder'):
            if node.up != None: # not root
                node.state = node.c[np.arange(L), node.up.state]
                node.sequence = alphabet[node.state]

    # TODO testing
    def optimize_branch_len(self, model, **kwargs):
        """
        Perform ML optimization for the tree branch length. **Note** this method assumes that each node stores information about its sequence as numpy.array object (variable node.sequence). Therefore, before calling this method, sequence reconstruction with either of the available models must be performed.

        Args:
         - model(GTR): evolutionary model

        KWargs:
         - verbose (int): output detalization
         - store_old (bool): if True, the old lenths will be saved in node.old_dist parameter. Useful for testing, and special post-processing.

        Returns:
         - None
        """
        from  scipy import optimize
        verbose = 1
        store_old_dist = False

        if 'verbose' in kwargs:
            verbose = int(kwargs['verbose'])
        if 'store_old' in kwargs:
            store_old_dist = kwargs['store_old'] == True

        if verbose > 3:
            print ("Walking up the tree, computing likelihood distrubutions")

        for node in self.tree.get_nonterminals(order='postorder'):
            parent = node.up
            if parent is None: continue # this is the root
            seq_p = parent.sequence
            seq_ch = node.sequence

            if (seq_p!=seq_ch).sum() == 0:
                if store_old_dist:
                    node.old_length = node.branch_length
                node.branch_length = 0
                if verbose > 5:
                    print ("Parent and child sequences are equal, setting branch len = 0")
            else:
                try:
                    opt = optimize.minimize_scalar(self._neg_prob,
                        bounds=[0,2],
                        method='Bounded',
                        args=(seq_p, seq_ch, model))
                except:
                    # FIXME error is caused by the unknown characters
                    # Introduce characters in the alphabet as profiles
                    continue

                new_len = opt["x"]

                if verbose > 5:
                    print ("Optimization results: old_len=%.4f, new_len=%.4f. "
                    " Updating branch length..." %(node.branch_length, new_len))

                if store_old_dist:
                    node.old_length = node.branch_length

                if new_len > 1.8 or opt["message"] != "Solution found.":
                    if verbose > 0:
                        print ("Cannot optimize branch length, minimization failed. Skipping")
                    continue
                else:
                    node.branch_length = new_len
        # as branch lengths changed, the params must be fixed
        self._add_node_params()
        return

    def _neg_prob(self, t, parent, child, tm):
        """
        Probability to observe child given the the parent state, transition matrix and the time of evolution (branch length).

        Args:
         - t(double): branch length (time between sequences)
         - parent (numpy.array): parent sequence
         - child(numpy.array): child sequence
         - tm (GTR): model of evolution

        Returns:
         - prob(double): negative probability of the two given sequences to be separated by the time t.
        """

        L = len(parent)
        if len(parent) != len(child):
            raise ValueError("Sequence lengths do not match!")

        eQT = np.diagflat(np.diag(np.exp(tm.mu * t * tm.eigenmat)))
        P_all = (tm.v).dot(eQT).dot(tm.v_inv)

        irow = seq2idx(child, tm.alphabet) # child states
        icol = seq2idx(parent, tm.alphabet) # parent states
        prob = P_all[irow, icol].prod()
        return -prob




