from __future__ import division, print_function
import numpy as np
import config as ttconf
from seq_utils import alphabets, profile_maps

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
        self.break_degen = np.random.random(size=self.W.shape)*0.0001
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
            Pi = np.ones(size=(n))
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
    def custom(cls,mu=1.0, pi=None, W=None, **kwargs):
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


    @classmethod
    def standard(cls, model='Jukes-Cantor', **kwargs):
        """
        Create one of the standard GTR models.

        Args:

         - model (str): type of the model. Currently supported models are:
         Jukes-Cantor, random.

        KWargs:

         - alphabet(str): specify alphabet when applicable. If the alphabet specification
         is requred, but no alphabet specified, the nucleotide will be used as default.

         - mu(double): general mutation rate. **NOTE** that the mutation rate is the
         only object which sets the time-scale to the GTR model.
         In other words, the unit of branch length and the unit of time are
         connected through this variable. By default set to 1.
        """
        if 'alphabet' in kwargs and kwargs['alphabet'] in alphabets.keys():
            alphabet = kwargs['alphabet']
        elif model in ['aa', 'prot']:
            alphabet='aa'
        else:
            alphabet = 'nuc'

        if 'mu' in kwargs:
            mu = kwargs['mu']
        else:
            mu = 1.0

        gtr = cls(alphabet)
        n = gtr.alphabet.shape[0]

        if model in ['Jukes-Cantor', 'JC69', 'aa', 'prot']:
            n = gtr.alphabet.shape[0]
            W, pi = np.ones((n,n)), np.ones(n)
        elif model=='random':
            pi = 1.0*np.random.randint(0,100,size=(n))
            W = 1.0*np.random.randint(0,100,size=(n,n)) # with gaps
        else:
            raise NotImplementedError("The specified evolutionary model is unsupported!")
        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr

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
        gtr = cls(**kwargs)
        gtr.logger("GTR: model inference ",1)
        dp = 1e-5
        Nit = 40
        pc_mat = pc*np.ones_like(nij)
        pc_mat -= np.diag(np.diag(pc_mat))
        from scipy import linalg as LA
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
            Wdiag = ((W_ij.T*pi).T).sum(axis=0)/pi
            np.fill_diagonal(W_ij, Wdiag)
            Q1 = np.diag(pi).dot(W_ij)
            scale_factor = np.sum(np.diagonal(Q1*np.diag(pi)))
            np.fill_diagonal(W_ij, 0)

            W_ij = W_ij/scale_factor
            if fixed_pi is None:
                pi = (np.sum(nij+pc_mat,axis=1)+root_state)/(mu*np.dot(W_ij,Ti)+root_state.sum()+np.sum(pc_mat, axis=1))
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


    def compress_sequence_pair(self, seq_p, seq_ch, ignore_gaps=False):
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
                            count = (bool_seqs_p[n1]&bool_seqs_ch[n2]).sum()
                            if count: pair_count.append(((n1,n2), count))
        else: # enumerate state pairs of the sequence for large alphabets
            num_seqs = []
            for seq in [seq_p, seq_ch]:
                tmp = np.ones_like(seq, dtype=int)
                for ni,nuc in enumerate(self.alphabet):
                    tmp[seq==nuc] = ni
                    num_seqs.append(tmp)
                    if ignore_gaps:
                        pair_count = Counter([x for x in zip(num_seqs[0], num_seqs[1])
                                              if (self.gap_index not in x)])
                    else:
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

    def prob_t(self, seq_p, seq_ch, t, return_log=False, ignore_gaps=True):
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
        seq_pair, multiplicity = self.compress_sequence_pair(seq_p, seq_ch, ignore_gaps=ignore_gaps)
        return self.prob_t_compressed(seq_pair, multiplicity, t, return_log=return_log)


    def optimal_t(self, seq_p, seq_ch, ignore_gaps=False):
        '''
        Find the optimal distance between the two sequences
        '''
        seq_pair, multiplicity = self.compress_sequence_pair(seq_p, seq_ch, ignore_gaps=ignore_gaps)
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

    def sequence_logLH(self,seq):
        """
        returns the loglikelihood of sampling a sequence from equilibrium frequency
        expects a sequence as numpy array
        """
        return np.sum([np.sum((seq==state)*np.log(self.Pi[si]))
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
