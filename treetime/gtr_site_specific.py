from __future__ import division, print_function, absolute_import
from collections import defaultdict
import numpy as np
from treetime import config as ttconf
from .seq_utils import alphabets, profile_maps, alphabet_synonyms


class GTR_site_specific(object):
    """
    Defines General-Time-Reversible model of character evolution.
    """
    def __init__(self, seq_len=1, alphabet='nuc', prof_map=None, logger=None):
        """
        Initialize empty evolutionary model.

        Parameters
        ----------

         alphabet : str, numpy.array
            Alphabet of the sequence. If a string is passed, it is understood as
            an alphabet name. In this case, the alphabet and its profile map are pulled
            from :py:obj:`treetime.seq_utils`.

            If a numpy array of characters is passed, a new alphabet is constructed,
            and the default profile map is atached to it.

         prof_map : dict
            Dictionary linking characters in the sequence to the likelihood
            of observing characters in the alphabet. This is used to
            implement ambiguous characters like 'N'=[1,1,1,1] which are
            equally likely to be any of the 4 nucleotides. Standard profile_maps
            are defined in file seq_utils.py. If None is provided, no ambigous
            characters are supported.

         logger : callable
            Custom logging function that should take arguments (msg, level, warn=False),
            where msg is a string and level an integer to be compared against verbose.

        """
        self.debug=False
        if isinstance(alphabet, str):
            if alphabet not in alphabet_synonyms:
                raise AttributeError("Unknown alphabet type specified")
            else:
                tmp_alphabet = alphabet_synonyms[alphabet]
                self.alphabet = alphabets[tmp_alphabet]
                self.profile_map = profile_maps[tmp_alphabet]
        else:
            # not a predefined alphabet
            self.alphabet = alphabet
            if prof_map is None: # generate trivial unambiguous profile map is none is given
                self.profile_map = {s:x for s,x in zip(self.alphabet, np.eye(len(self.alphabet)))}
            else:
                self.profile_map = prof_map

        self.seq_len = seq_len
        if logger is None:
            def logger_default(*args,**kwargs):
                """standard logging function if none provided"""
                if self.debug:
                    print(*args)
            self.logger = logger_default
        else:
            self.logger = logger
        n_states = len(self.alphabet)

        self.logger("GTR: with alphabet: "+str(self.alphabet),1)
        # determine if a character exists that corresponds to no info, i.e. all one profile
        if any([x.sum()==n_states for x in self.profile_map.values()]):
            amb_states = [c for c,x in self.profile_map.items() if x.sum()==n_states]
            self.ambiguous = 'N' if 'N' in amb_states else amb_states[0]
            self.logger("GTR: ambiguous character: "+self.ambiguous,2)
        else:
            self.ambiguous=None

        # check for a gap symbol
        try:
            self.gap_index = list(self.alphabet).index('-')
        except:
            self.logger("GTR: no gap symbol!", 4, warn=True)
            self.gap_index=-1


        # NEEDED TO BREAK RATE MATRIX DEGENERACY AND FORCE NP TO RETURN REAL ORTHONORMAL EIGENVECTORS
        # ugly hack, but works and shouldn't affect results
        tmp_rng_state = np.random.get_state()
        np.random.seed(12345)
        self.break_degen = np.random.random(size=(n_states, n_states))*1e-6
        np.random.set_state(tmp_rng_state)

        # init all matrices with dummy values
        self.logger("GTR: init with dummy values!", 3)
        self.v = None # right eigenvectors
        self.v_inv = None # left eigenvectors
        self.eigenvals = None # eigenvalues
        self.assign_rates()


    def assign_rates(self, mu=1.0, pi=None, W=None):
        """
        Overwrite the GTR model given the provided data

        Parameters
        ----------

         mu : float
            Substitution rate

         W : nxn matrix
            Substitution matrix

         pi : n vector
            Equilibrium frequencies

        """
        n = len(self.alphabet)
        self.mu = mu

        if pi is not None and pi.shape[0]==n:
            Pi = np.array(pi)
        else:
            if pi is not None and len(pi)!=n:
                self.logger("length of equilibrium frequency vector does not match alphabet length", 4, warn=True)
                self.logger("Ignoring input equilibrium frequencies", 4, warn=True)
            Pi = np.ones(shape=(n,seq_len))

        self.Pi = Pi/np.sum(Pi, axis=0)

        if W is None or W.shape!=(n,n):
            if (W is not None) and W.shape!=(n,n):
                self.logger("Substitution matrix size does not match alphabet size", 4, warn=True)
                self.logger("Ignoring input substitution matrix", 4, warn=True)
            # flow matrix
            W = np.ones((n,n))
            np.fill_diagonal(W, 0.0)
            np.fill_diagonal(W, -W.sum(axis=0))
        else:
            W=np.array(W)

        self.W = 0.5*(W+W.T)
        self._eig()


    @classmethod
    def random(cls, mu=1.0, alphabet='nuc'):
        """
        Creates a random GTR model

        Parameters
        ----------

         mu : float
            Substitution rate

         alphabet : str
            Alphabet name (should be standard: 'nuc', 'nuc_gap', 'aa', 'aa_gap')


        """

        alphabet=alphabets[alphabet]
        gtr = cls(alphabet)
        n = gtr.alphabet.shape[0]
        pi = 1.0*np.random.random(size=(n,self.seq_len))
        pi /= pi.sum(axis=0)
        tmp = 1.0*np.random.random(size=(n,n)) # with gaps
        W = tmp + tmp.T
        mu = 1.0 + 0.05*np.random.normal(size=(self.seq_len)) # with gaps

        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr

    @classmethod
    def custom(cls, mu=1.0, pi=None, W=None, **kwargs):
        """
        Create a GTR model by specifying the matrix explicitly

        Parameters
        ----------

         mu : float
            Substitution rate

         W : nxn matrix
            Substitution matrix

         pi : n vector
            Equilibrium frequencies

         **kwargs:
            Key word arguments to be passed

        Keyword Args
        ------------

         alphabet : str
            Specify alphabet when applicable. If the alphabet specification is
            required, but no alphabet is specified, the nucleotide alphabet will be used as
            default.

        """
        gtr = cls(**kwargs)
        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr


    def _eig_single_site(self, W, p):
        tmpp = np.sqrt(p)
        symQ = W*np.outer(tmpp, tmpp)
        np.fill_diagonal(symQ, -np.sum(W*p, axis=1))
        eigvals, eigvecs = np.linalg.eigh(symQ)
        tmp_v = eigvecs.T*tmpp
        one_norm = np.sum(np.abs(tmp_v), axis=1)
        return eigvals, tmp_v.T/one_norm, (eigvecs*one_norm).T/tmpp


    def _eig(self):
        eigvals, vec, vec_inv = [], [], []
        for pi in range(self.seq_len):
            if len(self.W.shape)>2:
                W = np.copy(self.W[:,:,pi])
                np.fill_diagonal(W, 0)
            elif pi==1:
                np.fill_diagonal(self.W, 0)
                W=self.W

            ev, evec, evec_inv = self._eig_single_site(W,self.Pi[:,pi])
            eigvals.append(ev)
            vec.append(evec)
            vec_inv.append(evec_inv)

        self.eigenvals = np.array(eigvals).T
        self.v = np.swapaxes(vec,0,-1)
        self.v_inv = np.swapaxes(vec_inv, 0,-1)


    def expQt(self, t):
        eLambdaT = np.exp(self.mu*t*self.eigenvals)
        return np.einsum('jia,ja,kja->ika',self.v, eLambdaT, self.v_inv)

    def prop_t_compressed(self, seq_pair, multiplicity, t, return_log=False):
        print("NOT IMPEMENTED")


    def prob_t_profiles(self, profile_pair, multiplicity, t, return_log=False, ignore_gaps=True):
        '''
        Calculate the probability of observing a node pair at a distance t

        Parameters
        ----------

          profile_pair: numpy arrays
            Probability distributions of the nucleotides at either
            end of the branch. pp[0] = parent, pp[1] = child

          multiplicity : numpy array
            The number of times an alignment pattern is observed

          t : float
            Length of the branch separating parent and child

          ignore_gaps: bool
            If True, ignore mutations to and from gaps in distance calculations

          return_log : bool
            Whether or not to exponentiate the result

        '''
        if t<0:
            logP = -ttconf.BIG_NUMBER
        else:
            Qt = self.expQt(t)
            res = np.einsum('ai,ija,aj->a', profile_pair[1], Qt, profile_pair[0])
            if ignore_gaps: # calculate the probability that neither outgroup/node has a gap
                non_gap_frac = (1-profile_pair[0][:,self.gap_index])*(1-profile_pair[1][:,self.gap_index])
                # weigh log LH by the non-gap probability
                logP = np.sum(multiplicity*np.log(res)*non_gap_frac)
            else:
                logP = np.sum(multiplicity*np.log(res))

        return logP if return_log else np.exp(logP)


    def propagate_profile(self, profile, t, return_log=False):
        """
        Compute the probability of the sequence state of the parent
        at time (t+t0, backwards), given the sequence state of the
        child (profile) at time t0.

        Parameters
        ----------

         profile : numpy.array
            Sequence profile. Shape = (L, a),
            where L - sequence length, a - alphabet size.

         t : double
            Time to propagate

         return_log: bool
            If True, return log-probability

        Returns
        -------

         res : np.array
            Profile of the sequence after time t in the past.
            Shape = (L, a), where L - sequence length, a - alphabet size.

        """
        Qt = self.expQt(t)
        res = np.einsum('ai,ija->aj', profile, Qt)

        return np.log(res) if return_log else res


    def prob_t(self, seq_p, seq_ch, t, pattern_multiplicity = None, return_log=False, ignore_gaps=True):
        """
        Compute the probability to observe seq_ch (child sequence) after time t starting from seq_p
        (parent sequence).

        Parameters
        ----------

         seq_p : character array
            Parent sequence

         seq_c : character array
            Child sequence

         t : double
            Time (branch len) separating the profiles.

         pattern_multiplicity : numpy array
            If sequences are reduced by combining identical alignment patterns,
            these multplicities need to be accounted for when counting the number
            of mutations across a branch. If None, all pattern are assumed to
            occur exactly once.

         return_log : bool
            It True, return log-probability.

        Returns
        -------
         prob : np.array
            Resulting probability

        """
        if t<0:
            logP = -ttconf.BIG_NUMBER
        else:
            tmp_eQT = self.expQt(t)
            bad_indices=(tmp_eQT==0)
            logQt = np.log(tmp_eQT + ttconf.TINY_NUMBER*(bad_indices))
            logQt[np.isnan(logQt) | np.isinf(logQt) | bad_indices] = -ttconf.BIG_NUMBER
            seq_indices_c = np.zeros(len(seq_c), dtype=int)
            seq_indices_p = np.zeros(len(seq_p), dtype=int)
            for ai, a in self.alphabet:
                seq_indices_p[seq_p==a] = ai
                seq_indices_c[seq_c==a] = ai

            if len(logQt.shape)==2:
                logP = np.sum(logQt[seq_indices_p, seq_indices_c]*pattern_multiplicity)
            else:
                logP = np.sum(logQt[seq_indices_p, seq_indices_c, np.arange(len(seq_c))]*pattern_multiplicity)

        return logP if return_log else np.exp(logP)
