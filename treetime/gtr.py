#!/usr/local/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
import numpy as np
import config as ttconf
from seq_utils import alphabets, profile_maps
from aa_models  import JTT92
from nuc_models import JC69, K80, F81, HKY85, T92, TN93

class GTR(object):
    """
    Defines General-Time-Reversible model of character evolution.
    """
    def __init__(self, alphabet='nuc', prof_map=None, logger=None):
        """
        Initialize empty evolutionary model.
        Args:
         - alphabet (numpy.array): alphabet of the sequence.
        """
        self.debug=False
        if type(alphabet)==str:
            if (alphabet not in alphabets):
                raise AttributeError("Unknown alphabet type specified")
            else:
                self.alphabet = alphabets[alphabet]
                self.profile_map = profile_maps[alphabet]
        else:
            self.alphabet = alphabet
            if prof_map is None:
                self.profile_map = {s:x for s,x in zip(self.alphabet, np.eye(len(self.alphabet)))}
            else:
                self.profile_map = prof_map

        if logger is None:
            def logger(*args,**kwargs):
                if self.debug:
                    print(*args)
            self.logger = logger
        else:
            self.logger = logger
        n_states = len(self.alphabet)

        self.logger("GTR: with alphabet: "+str(self.alphabet),1)
        if any([x.sum()==n_states for x in self.profile_map.values()]):
            self.ambiguous = [c for c,x in self.profile_map.iteritems() if x.sum()==n_states][0]
            self.logger("GTR: ambiguous character: "+self.ambiguous,2)
        else:
            self.ambiguous=None
        try:
            self.gap_index = list(self.alphabet).index('-')
        except:
            self.logger("GTR: no gap symbol!", 4, warn=True)
            self.gap_index=-1
        # general rate matrix
        self.W = np.zeros((n_states, n_states))
        # stationary states of the characters
        self.Pi = np.zeros(n_states)

        # mutation rate, scaling factor
        self.mu = 1.0
        # eigendecomposition of the GTR matrix
        # Pi.dot(W) = v.dot(eigenvals).dot(v_inv)
        self.v = np.zeros((n_states, n_states))
        self.v_inv = np.zeros((n_states, n_states))
        self.eigenvals = np.zeros(n_states)
        # NEEDED TO BREAK RATE MATRIX DEGENERACY AND FORCE NP TO RETURN REAL ORTHONORMAL EIGENVECTORS
        tmp_rng_state = np.random.get_state()
        np.random.seed(12345)
        self.break_degen = np.random.random(size=self.W.shape)*1e-6
        np.random.set_state(tmp_rng_state)

        # distance matrix (needed for topology optimization and for NJ)
        self.dm = None

    @property
    def Q(self):
        return (self.W*self.Pi).T


######################################################################
## constructor methods
######################################################################

    def __str__(self):
        '''
        string representation of the GTR model for pretty printing
        '''
        eq_freq_str = "Mutation rate (mu): "+str(np.round(self.mu,6))+'\n'

        eq_freq_str += "\nEquilibrium frequencies (pi_i):\n"
        for a,p in zip(self.alphabet, self.Pi):
            eq_freq_str+=str(a)+': '+str(np.round(p,4))+'\n'

        W_str = "\nSymmetrized rates from j->i (W_ij):\n"
        W_str+='\t'+'\t'.join(map(str, self.alphabet))+'\n'
        for a,Wi in zip(self.alphabet, self.W):
            W_str+=str(a)+'\t'+'\t'.join([str(np.round(max(0,p),4)) for p in Wi])+'\n'

        Q_str = "\nActual rates from j->i (Q_ij):\n"
        Q_str+='\t'+'\t'.join(map(str, self.alphabet))+'\n'
        for a,Qi in zip(self.alphabet, self.Q):
            Q_str+=str(a)+'\t'+'\t'.join([str(np.round(max(0,p),4)) for p in Qi])+'\n'

        return eq_freq_str + W_str + Q_str

    def assign_rates(self, mu=1.0, pi=None, W=None):
        n = len(self.alphabet)
        self.mu = mu

        if pi is not None and len(pi)==n:
            Pi = np.array(pi)
        else:
            if pi is not None and len(pi)!=n:
                self.logger("length of equilibrium frequency vector does not match alphabet length", 4, warn=True)
                self.logger("Ignoring input equilibrium frequencies", 4, warn=True)
            Pi = np.ones(shape=(n,))

        self.Pi = Pi/np.sum(Pi)

        if W is None or W.shape!=(n,n):
            if (W is not None) and W.shape!=(n,n):
                self.logger("Mutation matrix size does not match alphabet size", 4, warn=True)
                self.logger("Ignoring input mutation matrix", 4, warn=True)
            # flow matrix
            W = np.ones((n,n))
            np.fill_diagonal(self.W, - ((self.W).sum(axis=0) - 1))
        else:
            W=np.array(W)

        self.W = 0.5*(W+W.T)
        self._check_fix_Q(fixed_mu=True)
        self._eig()


    @classmethod
    def custom(cls, mu=1.0, pi=None, W=None, **kwargs):
        """
        Create a GTR model by specifying the matrix explicitly

        Args:
         - mu (float): mutation rate
         - W (nxn matrix): mutation matrix
         - pi (n vector): equilibrium frequencies

        KWargs:
         - alphabet(str): specify alphabet when applicable. If the alphabet specification
         is required, but no alphabet specified, the nucleotide will be used as default.
        """
        gtr = cls(**kwargs)
        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr

    @staticmethod
    def standard(model, **kwargs):
        """
        Create standard model of molecular evolution.

        Args:
         - model(str): model to create. List of the available models is
         (see description below): ['JC69',

        Kwargs:
         - model arguments

        Available models:

         --- JC69:
            Jukes-Cantor 1969 model. This model assumes equal concentrations
            of the nucleotides and equal transition rates between nucleotide states.
            For more info, see: Jukes and Cantor (1969). Evolution of Protein Molecules. New York: Academic Press. pp. 21–132.

            To create this model, use:

            Args: model='JC69' - keyword to specify the model

            Kwargs:

             - mu(float)  - specify the mutation rate

             - alphabet(str) - alphabet. By default 'nuc' is used (all gaps are ignored).
             specify 'nuc_gap' to treat gaps as characters


         --- K80:
            Kimura 1980 model. Assumes equal concentrations across nucleotides, but
            allows different rates between transitions and transversions. The ratio
            of the transversion/transition rates is given by kappa parameter.
            For more info, see
            Kimura (1980),  J. Mol. Evol. 16 (2): 111–120. doi:10.1007/BF01731581.

            Current implementation of the model does not account for the gaps.

            Kwargs:

             - mu(float): overall mutation rate

             - kappa(float): ratio of transversion/transition rates


         --- F81:
            Felsenstein 1981 model. Assumes non-equal concentrations across nucleotides,
            but the transition rate between all states is assumed to be equal. See
            Felsenstein (1981), J. Mol. Evol. 17  (6): 368–376. doi:10.1007/BF01734359
            for details.

            Current implementation of the model does not account for the gaps (treatment of
            gaps as characters is possible if specify alphabet='nuc_gap').

            Args:

             - mu(float): mutation rate

             - pi(numpy array): nucleotide concentrations

             - alphabet(str): alphabet to use. Default 'nuc', which discounts al gaps.
             'nuc-gap' alphabet enables treatmen of gaps as characters.



         --- HKY85:
            Hasegawa, Kishino and Yano 1985 model. Allows different concentrations of the
            nucleotides (as in F81) + distinguishes between transition/transversionmutations
            (similar to K80). Link:
            Hasegawa, Kishino, Yano (1985), J. Mol. Evol. 22 (2): 160–174. doi:10.1007/BF02101694

            Current implementation of the model does not account for the gaps

            Args:

             - mu(float): mutation rate

             - pi(numpy array): nucleotide concentrations

             - kappa(float): ratio of transversion/transition mutation rates



         --- T92:
            Tamura 1992 model. Extending Kimura  (1980) model for the case where a
            G+C-content bias exists. Link:
            Tamura K (1992),  Mol.  Biol. Evol. 9 (4): 678–687.  DOI: 10.1093/oxfordjournals.molbev.a040752

            Current implementation of the model does not account for the gaps

            Args:

             - mu(float): mutation rate

             - pi_GC(float): relative GC content

             - kappa(float): relative transversion/transition rate



         --- TN93:
            Tamura and Nei 1993. The model distinguishes between the two different types of
            transition: (A <-> G) is allowed to have a different rate to (C<->T).
            Transversions have the same rate. The frequencies of the nucleotides are allowed
            to be different. Link:
            Tamura, Nei (1993), MolBiol Evol. 10 (3): 512–526. DOI:10.1093/oxfordjournals.molbev.a040023

            Args:

             - mu(float): mutaion rate

             - kappa1(float): relative A<-->C, A<-->T, T<-->G and G<-->C rates

             - kappa2(float): relative C<-->T rate

            Note:

             - Rate of A<-->G mutation is set to one. All other rates (kappa1, kappa2)
            are specified relative to this rate
        """

        if model.lower() in ['jc', 'jc69', 'jukes-cantor', 'jukes-cantor69', 'jukescantor', 'jukescantor69']:
            return JC69(**kwargs)
        elif model.lower() in ['k80', 'kimura80', 'kimura1980']:
            return K80(**kwargs)
        elif model.lower() in ['f81', 'felsenstein81', 'felsenstein1981']:
            return F81(**kwargs)
        elif model.lower() in ['hky', 'hky85', 'hky1985']:
            return HKY85(**kwargs)
        elif model.lower() in ['t92', 'tamura92', 'tamura1992']:
            return T92(**kwargs)
        elif model.lower() in ['tn93', 'tamura_nei_93', 'tamuranei93']:
            return TN93(**kwargs)
        elif model.lower() in ['jtt', 'jtt92']:
            return JTT92(**kwargs)
        else:
            raise KeyError("The GTR model '{}' is not in the list of available models."
                "".format(model))

    @classmethod
    def random(cls, mu=1.0, alphabet='nuc'):

        alphabet=alphabets[alphabet]
        gtr = cls(alphabet)
        n = gtr.alphabet.shape[0]
        pi = 1.0*np.random.randint(0,100,size=(n))
        W = 1.0*np.random.randint(0,100,size=(n,n)) # with gaps

        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr

    # @classmethod
    # def standard(cls, model='Jukes-Cantor', **kwargs):
    #     """
    #     Create one of the standard GTR models.

    #     Args:

    #      - model (str): type of the model. Currently supported models are:
    #      Jukes-Cantor, random.

    #     KWargs:

    #      - alphabet(str): specify alphabet when applicable. If the alphabet specification
    #      is requred, but no alphabet specified, the nucleotide will be used as default.

    #      - mu(double): general mutation rate. **NOTE** that the mutation rate is the
    #      only object which sets the time-scale to the GTR model.
    #      In other words, the unit of branch length and the unit of time are
    #      connected through this variable. By default set to 1.
    #     """
    #     if 'alphabet' in kwargs and kwargs['alphabet'] in alphabets.keys():
    #         alphabet = kwargs['alphabet']
    #     elif model in ['aa', 'prot']:
    #         alphabet='aa'
    #     else:
    #         alphabet = 'nuc'

    #     if 'mu' in kwargs:
    #         mu = kwargs['mu']
    #     else:
    #         mu = 1.0

    #     gtr = cls(alphabet)
    #     n = gtr.alphabet.shape[0]

    #     if model in ['Jukes-Cantor', 'JC69', 'aa', 'prot']:
    #         n = gtr.alphabet.shape[0]
    #         W, pi = np.ones((n,n)), np.ones(n)
    #     elif model=='random':
    #         pi = 1.0*np.random.randint(0,100,size=(n))
    #         W = 1.0*np.random.randint(0,100,size=(n,n)) # with gaps
    #     else:
    #         raise NotImplementedError("The specified evolutionary model is unsupported!")
    #     gtr.assign_rates(mu=mu, pi=pi, W=W)
    #     return gtr

    @classmethod
    def infer(cls, nij, Ti, root_state, fixed_pi=None, pc=5.0, **kwargs):
        """
        Infer a GTR model by specifying the number of transitions and time spent in each
        character. The basic equation that is being solved is
            n_ij = pi_i W_ij T_j
        where n_ij are the transitions, pi_i are the equilibrium state frequencies,
        W_ij is the "mutation attempt matrix", while T_i is the time on the tree
        spent in character state i. To regularize the process, we add pseudocounts and
        also need to account for the fact that the root of the tree is in a particular
        state. the modified equation is
            n_ij + pc = pi_i W_ij (T_j+pc+root_state)

        Args:
         - nij (nxn matrix): the number of times a change in character state is observed
            between state i and j
         - Ti (n vector): the time spent in each character state
         - root_state( n vector): the number of characters in state i in the sequence
            of the root node.
         - pc (float): pseudocounts, this determines the lower cutoff on the rate when
            no mutation are observed
        KWargs:
         - alphabet(str): specify alphabet when applicable. If the alphabet specification
           is required, but no alphabet specified, the nucleotide will be used as default.
        """
        from scipy import linalg as LA
        gtr = cls(**kwargs)
        gtr.logger("GTR: model inference ",1)
        dp = 1e-5
        Nit = 40
        pc_mat = pc*np.ones_like(nij)
        np.fill_diagonal(pc_mat, 0.0)
        count = 0
        pi_old = np.zeros_like(Ti)
        if fixed_pi is None:
            pi = np.ones_like(Ti)
        else:
            pi = np.copy(fixed_pi)
        pi/=pi.sum()
        W_ij = np.ones_like(nij)
        mu = nij.sum()/Ti.sum()
        # if pi is fixed, this will immediately converge
        while LA.norm(pi_old-pi) > dp and count < Nit:
            gtr.logger(' '.join(map(str, ['GTR inference iteration',count,'change:',LA.norm(pi_old-pi)])), 3)
            count += 1
            pi_old = np.copy(pi)
            W_ij = (nij+nij.T+2*pc_mat)/mu/(np.outer(pi,Ti) + np.outer(Ti,pi)
                                                    + ttconf.TINY_NUMBER + 2*pc_mat)

            np.fill_diagonal(W_ij, 0)
            Wdiag = (((W_ij.T*pi).T).sum(axis=0)+ttconf.TINY_NUMBER)/(pi+ttconf.TINY_NUMBER)
            np.fill_diagonal(W_ij, Wdiag)
            Q1 = np.diag(pi).dot(W_ij)
            scale_factor = np.sum(np.diagonal(Q1*np.diag(pi)))
            np.fill_diagonal(W_ij, 0)

            W_ij = W_ij/scale_factor
            if fixed_pi is None:
                pi = (np.sum(nij+pc_mat,axis=1)+root_state)/(ttconf.TINY_NUMBER + mu*np.dot(W_ij,Ti)+root_state.sum()+np.sum(pc_mat, axis=1))
                pi /= pi.sum()
            mu = nij.sum()/(ttconf.TINY_NUMBER + np.sum(pi * (W_ij.dot(Ti))))
        if count >= Nit:
            gtr.logger('WARNING: maximum number of iterations has been reached in GTR inference',3, warn=True)
            np.min(pi.sum(axis=0)), np.max(pi.sum(axis=0))
            if LA.norm(pi_old-pi) > dp:
                gtr.logger('the iterative scheme has not converged',3,warn=True)
            elif np.abs(1-np.max(pi.sum(axis=0))) > dp:
                gtr.logger('the iterative scheme has converged, but proper normalization was not reached',3,warn=True)

        gtr.assign_rates(mu=mu, W=W_ij, pi=pi)
        return gtr

########################################################################
### prepare model
########################################################################
    def _check_fix_Q(self, fixed_mu=False):
        """
        Check the main diagonal of Q and fix it in case it does not corresond
        the definition of the rate matrix. Should be run every time when creating
        custom GTR model.
        """
        # fix Q
        self.Pi /= self.Pi.sum() # correct the Pi manually
        # NEEDED TO BREAK RATE MATRIX DEGENERACY AND FORCE NP TO RETURN REAL ORTHONORMAL EIGENVECTORS
        self.W += self.break_degen + self.break_degen.T
        # fix W
        np.fill_diagonal(self.W, 0)
        Wdiag = -(self.Q).sum(axis=0)/self.Pi
        np.fill_diagonal(self.W, Wdiag)
        scale_factor = -np.sum(np.diagonal(self.Q)*self.Pi)
        self.W /= scale_factor
        if not fixed_mu:
            self.mu *= scale_factor
        if (self.Q.sum(axis=0) < 1e-10).sum() <  self.alphabet.shape[0]: # fix failed
            print ("Cannot fix the diagonal of the GTR rate matrix. Should be all zero", self.Q.sum(axis=0))
            import ipdb; ipdb.set_trace()
            raise ArithmeticError("Cannot fix the diagonal of the GTR rate matrix.")


    def _eig(self):
        """
        Perform eigendecompositon of the rate matrix and stores the left- and right-
        matrices to convert the sequence profiles to the GTR matrix eigenspace
        and hence to speed-up the computations.
        """
        # eigendecomposition of the rate matrix
        eigvals, eigvecs = np.linalg.eig(self.Q)
        self.v = np.real(eigvecs)
        self.v_inv = np.linalg.inv(self.v)
        self.eigenvals = np.real(eigvals)
        return


    def compress_sequence_pair(self, seq_p, seq_ch, pattern_multiplicity=None, ignore_gaps=False):
        '''
        make a compressed representation of a pair of sequences only counting
        the number of times a particular pair of states (e.g. (A,T)) is observed
        the the aligned sequences of parent and child
        Args:
          - seq_p:  parent sequence as numpy array
          - seq_ch: child sequence as numpy array
          - ignore_gap: whether or not to include gapped positions of the alignment
                        in the multiplicity count
        Returns:
          - seq_pair: [(0,1), (2,2), (3,4)] list of parent_child state pairs
                      as indices in the alphabet
          - multiplicity: number of times a particular pair is observed
        '''
        if pattern_multiplicity is None:
            pattern_multiplicity = np.ones_like(seq_p, dtype=float)

        from collections import Counter
        if seq_ch.shape != seq_ch.shape:
            raise ValueError("GTR.compress_sequence_pair: Sequence lengths do not match!")

        if len(self.alphabet)<10: # for small alphabet, repeatedly check array for all state pairs
            pair_count = []
            bool_seqs_p = []
            bool_seqs_ch = []
            for seq, bs in [(seq_p,bool_seqs_p), (seq_ch, bool_seqs_ch)]:
                for ni,nuc in enumerate(self.alphabet):
                    bs.append(seq==nuc)

            for n1,nuc1 in enumerate(self.alphabet):
                if (n1!=self.gap_index or (not ignore_gaps)):
                    for n2,nuc2 in enumerate(self.alphabet):
                        if (n2!=self.gap_index or (not ignore_gaps)):
                            count = ((bool_seqs_p[n1]&bool_seqs_ch[n2])*pattern_multiplicity).sum()
                            if count: pair_count.append(((n1,n2), count))
        else: # enumerate state pairs of the sequence for large alphabets
            num_seqs = []
            for seq in [seq_p, seq_ch]: # for each sequence (parent and child) construct a numerical sequence [0,5,3,1,2,3...]
                tmp = np.ones_like(seq, dtype=int)
                for ni,nuc in enumerate(self.alphabet):
                    tmp[seq==nuc] = ni  # set each position corresponding to a state to the corresponding index
                num_seqs.append(tmp)
            if ignore_gaps:  # if gaps are ingnored skip positions where one or the other sequence is gapped
                pair_count = Counter([x for x in zip(num_seqs[0], num_seqs[1])
                                      if (self.gap_index not in x)])
            else: # otherwise, just count
                pair_count = Counter(zip(num_seqs[0], num_seqs[1]))
            pair_count = pair_count.items()

        return (np.array([x[0] for x in pair_count], dtype=int),    # [(child_nuc, parent_nuc),()...]
                np.array([x[1] for x in pair_count], dtype=int))    # multiplicity of each parent/child nuc pair


########################################################################
### evolution functions
########################################################################
    def prob_t_compressed(self, seq_pair, multiplicity, t, return_log=False):
        '''
        calculate the probability of observing a sequence pair at a distance t
        Args:
          - seq_pair:     np.array([(0,1), (2,2), ()..]) as indicies of
                pairs of aligned positions. (e.g. 'A'==0, 'C'==1 etc)
                this only lists all occuring parent-child state pairs, order is irrelevant
          - multiplicity: the number of times a parent-child state pair is observed
                this allows to compress the sequence representation
          - t:            length of the branch separating parent and child
          - return_log:   whether or not to exponentiate the result
        '''
        if (t<0):
            logP = -ttconf.BIG_NUMBER
        else:
            tmp_eQT = self.expQt(t)
            bad_indices=(tmp_eQT==0)
            logQt = np.log(tmp_eQT + ttconf.TINY_NUMBER*(bad_indices))
            logQt[np.isnan(logQt) | np.isinf(logQt) | bad_indices] = -ttconf.BIG_NUMBER
            logP = np.sum(logQt[seq_pair[:,1], seq_pair[:,0]]*multiplicity)
            if return_log:
                return logP
            else:
                return np.exp(logP)

    def prob_t(self, seq_p, seq_ch, t, pattern_multiplicity = None, return_log=False, ignore_gaps=True):
        """
        Compute the probability to observe seq_ch after time t starting from seq_p.
        Args:
         - profile_p(np.array): parent profile of shape (L, a), where
         L - length of the sequence, a - alphabet size.

         - profile_ch(np.array): child profile of shape (L, a), where
         L - length of the sequence, a - alphabet size.

         - t (double): time (branch len), separating the profiles.

         - return_log(bool, default False): whether return log-probability.

        Returns:
         - prob(np.array): resulting probability.
        """
        seq_pair, multiplicity = self.compress_sequence_pair(seq_p, seq_ch,
                                        pattern_multiplicity=pattern_multiplicity, ignore_gaps=ignore_gaps)
        return self.prob_t_compressed(seq_pair, multiplicity, t, return_log=return_log)


    def optimal_t(self, seq_p, seq_ch, pattern_multiplicity=None, ignore_gaps=False):
        '''
        Find the optimal distance between the two sequences
        '''
        seq_pair, multiplicity = self.compress_sequence_pair(seq_p, seq_ch,
                                                            multiplicity = pattern_multiplicity,
                                                            ignore_gaps=ignore_gaps)
        return self.optimal_t_compressed(seq_pair, multiplicity)


    def optimal_t_compressed(self, seq_pair, multiplicity):
        """
        Find the optimal distance between the two sequences
        """

        def _neg_prob(t, seq_pair, multiplicity):
            """
            Probability to observe child given the the parent state, transition
            matrix and the time of evolution (branch length).

            Args:
             - t(double): branch length (time between sequences)
             - parent (numpy.array): parent sequence
             - child(numpy.array): child sequence
             - tm (GTR): model of evolution

            Returns:
             - prob(double): negative probability of the two given sequences
               to be separated by the time t.
            """
            return -1.0*self.prob_t_compressed(seq_pair, multiplicity,t, return_log=True)

        try:
            from scipy.optimize import minimize_scalar
            opt = minimize_scalar(_neg_prob,
                    bounds=[0,ttconf.MAX_BRANCH_LENGTH],
                    method='bounded',
                    args=(seq_pair, multiplicity), options={'xatol':1e-8})
            new_len = opt["x"]
        except:
            import scipy
            print('legacy scipy', scipy.__version__)
            from scipy.optimize import fminbound
            new_len = fminbound(_neg_prob,
                    0,ttconf.MAX_BRANCH_LENGTH,
                    args=(seq_pair, multiplicity))
            opt={'success':True}

        if new_len > .9 * ttconf.MAX_BRANCH_LENGTH:
            self.logger("WARNING: GTR.optimal_t_compressed -- The branch length seems to be very long!", 4, warn=True)

        if opt["success"] != True:
            # return hamming distance: number of state pairs where state differs/all pairs
            new_len =  np.sum(multiplicity[seq_pair[:,1]!=seq_pair[:,0]])/np.sum(multiplicity)

        return new_len


    def propagate_profile(self, profile, t, return_log=False):
        """
        Compute the probability of the sequence state (profile) at time (t+t0),
        given the sequence state (profile) at time t0.
        Args:
         - profile(numpy.array): sequence profile. Shape = (L, a),
         where L - sequence length, a - alphabet size.

         - t(double): time to propagate

         - return log (bool, default False): whether to return log-probability

        Returns:
         - res(np.array): profile of the sequence after time t.
         Shape = (L, a), where L - sequence length, a - alphabet size.
        """
        Qt = self.expQt(t)
        res = profile.dot(Qt)

        if return_log:
            return np.log(res)
        else:
            return res


    def _exp_lt(self, t, mu_prefactor=1.0):
        """
        Returns:
         - exp_lt(numpy.array): array of values exp(lambda(i) * t),
         where (i) - alphabet index (the eigenvalue number).
        """
        return np.exp(self.mu * t * self.eigenvals)

    def expQt(self, t):
        eLambdaT = np.diag(self._exp_lt(t)) # vector length = a
        Qt = self.v.dot(eLambdaT.dot(self.v_inv))   # This is P(nuc1 | given nuc_2)
        return np.maximum(0,Qt)

    def sequence_logLH(self,seq, pattern_multiplicity=None):
        """
        returns the loglikelihood of sampling a sequence from equilibrium frequency
        expects a sequence as numpy array
        """
        if pattern_multiplicity is None:
            pattern_multiplicity = np.ones_like(seq, dtype=float)
        return np.sum([np.sum((seq==state)*pattern_multiplicity*np.log(self.Pi[si]))
                      for si,state in enumerate(self.alphabet)])


    def save_to_npz(self, outfile):
        full_gtr = self.mu * np.dot(self.Pi, self.W)
        desc=np.array(["GTR matrix description\n", "Mutation rate: " + str(self.mu)])
        np.savez(outfile,   description=desc,
                            full_gtr=full_gtr,
                            char_dist=self.Pi,
                            flow_matrix=self.W)

    def save_to_json(self, zip):
        d = {
        "full_gtr": self.mu * np.dot(self.Pi, self.W),
        "Mutation rate" : mu,
        "Equilibrium character composition": self.Pi,
        "Flow rate matrix": self.W
        }


if __name__ == "__main__":
     pass
