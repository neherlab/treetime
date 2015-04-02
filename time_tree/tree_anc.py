"""
Class defines simple tree object with basic interface methdos: reading and
saving from/to files, initializing leaves with sequences from the alignment,
making ancestral state inferrence
"""

from Bio import Phylo
from Bio import AlignIO
import numpy as np

class TreeAnc(object):
    def __init__(self, tree):
        self.tree = tree
        self._set_links_to_parent()

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

    @classmethod
    def _from_json(cls, inf):
        raise NotImplementedError("This functionality is under development")
        pass

    def _set_links_to_parent(self):
        """
        set link to parent for all tree nodes
        """
        self.tree.root.up = None
        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            for c in clade.clades:
                c.up = clade
        return

    # ancestral state reconstruction
    def set_seqs_to_leaves(self, aln):
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
        pass

    def _fitch_anc(self):
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

    def _ml_anc(self, model):
        pass

    # FIXME
    def optimize_branch_len(self, model):
        """
        Perform ML optimization for the tree branch length
        Args:
         - model(TMat): evolutionary model
        KWargs:
         - verbose (int): output detalization
         - store_old (bool): if True, the old lenths will be saved in node.old_dist parameter.
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
            if parent is None: continue
            seq_p = parent.sequence
            seq_ch = node.sequence
            old_len = node.dist

            if (seq_p!=seq_ch).sum() == 0:
                if store_old_dist:
                    node.old_dist = node.dist
                node.dist = 0
                if verbose > 5:
                    print ("Parent and child sequences are equal, setting branch len = 0")
            else:
                opt = optimize.minimize_scalar(_neg_prob, bounds=[0,2],
                    method='Bounded',
                    args=(seq_p, seq_ch, tm))

                new_len = opt["x"]

                if verbose > 5:
                    print ("Optimization results: old_len=%.4f, new_len=%.4f. "
                    " Updating branch length..." %(node.dist, new_len))

                if store_old_dist:
                    node.old_dist = node.dist

                if new_len > 1.8 or opt["message"] != "Solution found.":
                    if verbose > 0:
                        print ("Cannot optiimize tree branch, minimization failed. Skipping")
                    continue
                else:
                    node.dist = new_len
        return

    # FIXME
    def _neg_prob(self, _time, seq_p, seq_ch, tm):
        """
        Probability to observe child given the the parent state, transition matrix
        and the time of evolution
        """
        L = len(seq_p)
        eQT = np.array([np.diagflat(
                            np.exp(tm._mu[pos]*_time*tm._EigVals[pos, :]))
                        for pos in xrange(tm._EigVals.shape[0])])

        P_all = np.einsum('ijk,ikl,ilm->ijm', tm._V, eQT, tm._Vinv)

        irow = seq2idx(seq_ch, tm.alphabet)
        icol = seq2idx(seq_p, tm.alphabet)
        ipos = np.arange(tm.L)
        prob = np.prod(P_all[ipos, irow, icol])
        return -prob




