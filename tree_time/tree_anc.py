from Bio import Phylo
from Bio import AlignIO
import numpy as np
from scipy import optimize as sciopt

alphabets = {
            "nuc": np.array(['A', 'C', 'G', 'T', '-']),
            "aa": np.array(['-'])}

_full_nc_profile = {
    'A': np.array([1, 0, 0, 0, 0], dtype='float'),
    'C': np.array([0, 1, 0, 0, 0], dtype='float'),
    'G': np.array([0, 0, 1, 0, 0], dtype='float'),
    'T': np.array([0, 0, 0, 1, 0], dtype='float'),
    '-': np.array([0, 0, 0, 0, 1], dtype='float'),
    'N': np.array([1, 1, 1, 1, 1], dtype='float'),
    'R': np.array([1, 0, 1, 0, 0], dtype='float'),
    'Y': np.array([0, 1, 0, 1, 0], dtype='float'),
    'S': np.array([0, 1, 1, 0, 0], dtype='float'),
    'W': np.array([1, 0, 0, 0, 1], dtype='float'),
    'K': np.array([0, 0, 0, 1, 1], dtype='float'),
    'M': np.array([1, 1, 0, 0, 0], dtype='float'),
    'D': np.array([1, 0, 1, 1, 0], dtype='float'),
    'H': np.array([1, 1, 0, 1, 0], dtype='float'),
    'B': np.array([0, 1, 1, 1, 0], dtype='float'),
    'V': np.array([1, 1, 1, 0, 0], dtype='float')}


def prepare_seq(seq):
    """
    Args:
     - seq:  sequence as an object of SeqRecord, string or iterable
    Returns:
     - sequence as numpy array
    """
    try:
        sequence = ''.join(seq)
    except TypeError:
        sequence = seq
    sequence = sequence.upper()

    return np.array(list(sequence))


def seq2prof(x, aln_type='nuc'):
    """
    Convert the given character into the profile.
    Args:
     - x(numpy.array): sequence to be converted to the profile
     - alph(str): alphabet type. Can be either 'nuc' for nucleotide alphabet, or 'aa' for amino acid alphabet
    Returns:
     - idx(numpy.array): profile for the character, zero array if the character not found
    """
    if aln_type=='nuc':
        prof = np.array([_full_nc_profile[k]
            if k in _full_nc_profile else _full_nc_profile['N'] for k in x])
        err = ((prof == 0.2).sum(1) != 0).sum()
        if err>0:
            print ("Seq2Profile: %d of %d characters were not identified or"
                    " not sequenced." % (err, prof.shape[0]))
    elif aln_type=='aa':
        raise NotImplementedError("Amino-acid alphabet is under development.")
    else:
        raise TypeError("Alignment type cane be either 'nuc' or 'aa'")
    return prof


class GTR(object):
    """
    Defines General tme reversible model of character evolution.
    """
    def __init__(self, alphabet_type):
        """
        Initialize empty evolutionary model.
        Args:
         - alphabet (numpy.array): alphabet of the sequence.
        """
        if not alphabet_type in alphabets:
            raise AttributeError("Unknown alphabet type specified")

        self.alphabet_type = alphabet_type
        alphabet = alphabets[alphabet_type]
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
        self.eigenmat = np.zeros(alphabet.shape[0])

    @classmethod
    def standard(cls, model='Jukes-Cantor', **kwargs):
        if 'alphabet' in kwargs and alphabet in alphabets.keys():
            alphabet = kwargs['alphabet']
        else:
            print ("No alphabet specified. Using default nucleotide.")
            alphabet = 'nuc'
        if 'mu' in kwargs:
            mu = kwargs['mu']
        else:
            mu = 1.0

        if model=='Jukes-Cantor':

            gtr = cls('nuc')
            gtr.mu = mu
            a = gtr.alphabet.shape[0]

            # flow matrix
            gtr.W = np.ones((a,a))
            np.fill_diagonal(gtr.W, - ((gtr.W).sum(0) - 1))

            # equilibrium concentrations matrix
            gtr.Pi = np.zeros(gtr.W.shape)
            np.fill_diagonal(gtr.Pi, 1.0/a)

            gtr._check_fix_Q() # make sure the main diagonal is correct
            gtr._eig() # eigendecompose the rate matrix
            return gtr

        elif model=='random':
            gtr = cls(alphabet)
            a = gtr.alphabet.shape[0]

            gtr.mu = mu

            Pi = 1.0*np.random.randint(0,100,size=(a))
            Pi /= Pi.sum()
            gtr.Pi = np.diagflat(Pi)

            W = 1.0*np.random.randint(0,100,size=(a,a)) # with gaps
            gtr.W = W+W.T

            gtr._check_fix_Q()
            gtr._eig()
            return gtr
        else:
            raise NotImplementedError("The specified evolutionary model is unsupported!")

    def _check_fix_Q(self):
        """
        Check the main diagonal of Q and fix it in case it does not corresond the definition of Q.
        """
        Q = self.Pi.dot(self.W)
        if (Q.sum(0) < 1e-10).sum() < self.alphabet.shape[0]: # at least one rate is wrong
            # fix Q
            self.Pi /= self.Pi.sum() # correct the Pi manually
            # fix W
            np.fill_diagonal(self.W, 0)
            Wdiag = -((self.W.T*np.diagonal(self.Pi)).T).sum(0)/ \
                    np.diagonal(self.Pi)
            np.fill_diagonal(self.W, Wdiag)
            Q1 = self.Pi.dot(self.W)
            if (Q1.sum(0) < 1e-10).sum() <  self.alphabet.shape[0]: # fix failed
                raise ArithmeticError("Cannot fix the diagonal of the GTR rate matrix.")
        return

    def _eig(self):
        """
        Perform eigendecompositon of the rate matrix
        """
        # eigendecomposition of the rate matrix
        eigvals, eigvecs = np.linalg.eig(self.Pi.dot(self.W))
        self.v = eigvecs
        self.v_inv = np.linalg.inv(self.v)
        self.eigenmat = eigvals
        return

    def prob_t(self, profile_p, profile_ch, t, rotated=False, return_log=False):
        """
        Compute the probability of the two profiles to be separated by the time t.
        Args:
         - profile_p(np.array): parent profile of shape (L, a), where L - length of the sequence, a - alpphabet size.

         - profile_ch(np.array): child profile of shape (L, a), where L - length of the sequence, a - alpphabet size.

         - t (double): time (branch len), separating the profiles.

         - rotated (bool, default False): if True, assume that the supplied profiles are already rotated.

         - return_log(bool, default False): whether return log-probability.

        Returns:
         - prob(np.array): resulting probability.
        """
        
        if t < 0:
            if return_log:
                return -500
            else:
                return 0.0

        L = profile_p.shape[0]
        if L != profile_ch.shape[0]:
            raise ValueError("Sequence lengths do not match!")
        eLambdaT = self._exp_lt(t)
        if not rotated: # we need to rotate first
            p1 = profile_p.dot(self.v) # (L x a).dot(a x a) = (L x a) - prof in eigenspace
            p2 = (self.v_inv.dot(profile_ch.T)).T # (L x a).dot(a x a) = (L x a) - prof in eigenspace
        else:
            p1 = profile_p
            p2 = profile_ch
            #prob = (profile_p*eLambdaT*profile_ch).sum(1) # sum over the alphabet

        prob = (p1*eLambdaT*p2).sum(1) # sum_i (p1_i * exp(l_i*t) * p_2_i) result = vector lenght L
        prob[prob<0] = 0.0 # avoid rounding instability

        if return_log:
            prob = (np.log(prob + 1e-50)).sum() # sum all sites
        else:
            prob = prob.prod() # prod of all sites
        return prob

    def propagate_profile(self, profile, t, rotated=False, return_log=False):
        """
        Compute the probability of the sequence state (profile) at time (t+t0), given the sequence state (profile) at time t0.
        Args:
         - profile(numpy.array): sequence profile. Shape = (L, a), where L - sequence length, a - alphabet size.

         - t(doble): time to propagate

         - rotated(bool default False): whether the supplied profile is in the GTR matrix eigenspace

         - return log (bool, default False): whether to return log-probability

        Returns:
         - res(np.array): profile of the sequence after time t. Shape = (L, a), where L - sequence length, a - alphabet size.
        """
        eLambdaT = self._exp_lt(t) # vector lenght = a

        if not rotated:
            # rotate
            p = self.v_inv.dot(profile.T).T
        else:
            p = profile

        res = (self.v.dot((eLambdaT * p).T)).T

        if not return_log:
            return res

        else:
            return np.log(res)

    def _exp_lt(self, t):
        """
        Returns:
         - exp_lt(numpy.array): array of values exp(lambda(i) * t), where (i) - alphabet index (the eigenvalue number).
        """
        return np.exp(self.mu * t * self.eigenmat)


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

    # def has_attr(self, node, attr):
    #    return node.__dict__.has_key(attr)

    def _add_node_params(self):
        """
        Set link to parent and net distance to root for all tree nodes
        """
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        for clade in self.tree.get_nonterminals(order='preorder'): # up->down
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
            if l.name in dic_aln:
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
                gtr = GTR.standard(model='Jukes-Cantor')
            else:
                gtr = kwargs.pop('model')

            self._ml_anc(gtr, **kwargs)
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
            node.state = [self._fitch_state(node, k) for k in range(L)]

        ambs = [i for i in range(L) if len(self.tree.root.state[i])>1]
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
                       else node.state[i][0] for i in range(L)])
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

    def _ml_anc(self, gtr, **kwargs):
        """
        Perform ML reconstruction for the ancestral states
        Args:
         - model (GTR): General time-reversible model of evolution.
        KWargs:
         - store_lh (bool): if True, all likelihoods will be stored for all nodes. Useful for testing, diagnostics and if special post-processing is required.
         - verbose (int): how verbose the output should be
        """

        verbose = 0 # how verbose to be at the output
        if 'store_lh' in kwargs:
            store_lh = kwargs['store_lh'] == True
        if 'verbose' in kwargs:
            try:
                verbose = int(kwargs['verbose'])
            except:
                print ("ML ERROR in input: verbose param must be int")
                verbose = 5

        L = self.tree.get_terminals()[0].sequence.shape[0]
        a = gtr.alphabet.shape[0]

        if verbose > 2:
            print ("Walking up the tree, computing joint likelihoods...")

        for leaf in self.tree.get_terminals():

            # in any case, set the profile
            leaf.profile = seq2prof(leaf.sequence, gtr.alphabet_type)
            leaf.lh_prefactor = np.zeros(L)

        for node in self.tree.get_nonterminals(order='postorder'): #leaves -> root
            # regardless of what was before, set the profile to zeros
            node.profile = np.ones((L, a)) # we will multiply it
            node.lh_prefactor = np.zeros(L)

            for ch in node.clades:

                node.profile *= gtr.propagate_profile(ch.profile,
                    ch.branch_length,
                    rotated=False, # use unrotated
                    return_log=False) # raw prob to transfer prob up

                node.lh_prefactor += ch.lh_prefactor

            pre = node.profile.sum(1)

            node.profile /= pre.reshape((pre.shape[0], 1)) # normalize so that the sum is 1
            node.lh_prefactor += np.log(pre) # and store log-prefactor

        if (verbose > 2):
            print ("Walking down the tree, computing maximum likelihood     sequences...")

        self.tree.root.profile *= np.diag(gtr.Pi) # enable time-reversibility
        self.tree.root.sequence, self.tree.root.profile = \
            self._prof_to_seq(self.tree.root.profile, gtr, True)

        for node in self.tree.find_clades(order='preorder'):
            if node.up != None: # not root
                node.profile *= gtr.propagate_profile(node.up.profile,
                                node.branch_length,
                                rotated=False, # use unrotated
                                return_log=False)
                # actually, infer sequence
                node.sequence,node.profile=self._prof_to_seq(node.profile,gtr)

    def _prof_to_seq(self, profile, gtr, correct_prof=True):
        seq = gtr.alphabet[profile.argmax(1)]  # max LH over the alphabet
        if correct_prof:  # max profile value to one, others - zeros
            am = profile.argmax(1)
            profile[:, :] = 0.0
            profile[range(profile.shape[0]), am] = 1.0
        return seq, profile

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

        verbose = 1
        store_old_dist = False

        if 'verbose' in kwargs:
            verbose = int(kwargs['verbose'])
        if 'store_old' in kwargs:
            store_old_dist = kwargs['store_old'] == True

        if verbose > 3:
            print ("Walking up the tree, computing likelihood distrubutions")

        for node in self.tree.find_clades(order='postorder'):
            parent = node.up
            if parent is None: continue # this is the root
            prof_p = parent.profile
            prof_ch = node.profile

            if store_old_dist:
                node.old_length = node.branch_length

            # optimization method
            #import ipdb; ipdb.set_trace()
            new_len = self._opt_len(prof_p, prof_ch, model)
            if new_len < 0:
                continue

            if verbose > 5:
                print ("Optimization results: old_len=%.4f, new_len=%.4f "
                        " Updating branch length..."
                        %(node.branch_length, new_len))

            node.branch_length = new_len

        # as branch lengths changed, the params must be fixed
        self._add_node_params()
        return

    def _opt_len(self, seq_p, seq_ch, gtr, verbose=10):
        opt = sciopt.minimize_scalar(self._neg_prob,
                bounds=[0,.2],
                method='Bounded',
                args=(seq_p, seq_ch, gtr))

        new_len = opt["x"]

        if new_len > .18 or opt["success"] != True:
            if verbose > 0:
                print ("Cannot optimize branch length, minimization failed.")
            
            return -1.0
        else:
            return  new_len

    def _neg_prob(self, t, parent, child, gtr):
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
        return - gtr.prob_t (parent, child, t, rotated=False, return_log=False)


