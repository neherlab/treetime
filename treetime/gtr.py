from __future__ import division, print_function
import numpy as np
import config as ttconf
from seq_utils import alphabets, profile_maps

class GTR(object):
    """
    Defines General-Time-Reversible model of character evolution.
    """
    def __init__(self, alphabet='nuc', prof_map=None):
        """
        Initialize empty evolutionary model.
        Args:
         - alphabet (numpy.array): alphabet of the sequence.
        """
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
        n_states = len(self.alphabet)

        # general rate matrix
        self.W = np.zeros((n_states, n_states))
        # stationary states of the characters
        self.Pi = np.zeros((n_states, n_states))

        # mutation rate, scaling factor
        self.mu = 1.0
        # eigendecomposition of the GTR matrix
        # Pi.dot(W) = v.dot(eigenmat).dot(v_inv)
        self.v = np.zeros((n_states, n_states))
        self.v_inv = np.zeros((n_states, n_states))
        self.eigenmat = np.zeros(n_states)
        # NEEDED TO BREAK RATE MATRIX DEGENERACY AND FORCE NP TO RETURN REAL ORTHONORMAL EIGENVECTORS
        self.break_degen = np.random.random(size=self.W.shape)*0.0001

        # distance matrix (needed for topology optimization and for NJ)
        self.dm = None

    def __str__(self):
        eq_freq_str = "Equilibrium frequencies (pi_i):\n"
        for a,p in zip(self.alphabet, np.diagonal(self.Pi)):
            eq_freq_str+=str(a)+': '+str(np.round(p,4))+'\n'

        W_str = "\nSymmetrized rates from j->i (W_ij):\n"
        W_str+='\t'+'\t'.join(map(str, self.alphabet))+'\n'
        for a,Wi in zip(self.alphabet, self.W):
            W_str+=str(a)+'\t'+'\t'.join([str(np.round(max(0,p),4)) for p in Wi])+'\n'

        Q_str = "\nActual rates from j->i (Q_ij):\n"
        Q_str+='\t'+'\t'.join(map(str, self.alphabet))+'\n'
        for a,Qi in zip(self.alphabet, self.Pi.dot(self.W)):
            Q_str+=str(a)+'\t'+'\t'.join([str(np.round(max(0,p),4)) for p in Qi])+'\n'

        return eq_freq_str + W_str + Q_str

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
         is requred, but no alphabet specified, the nucleotide will be used as default.
        """
        gtr = cls(**kwargs)
        n = gtr.alphabet.shape[0]

        gtr.mu = mu
        if pi is not None and len(pi)==n:
            Pi = pi
        else:
            if pi is not None and len(pi)!=n:
                print("length of equilibrium frequency vector does not match alphabet length"
                      "Ignoring input equilibrium frequencies")
            Pi = np.ones(size=(n))
        Pi /= Pi.sum()
        gtr.Pi = np.diagflat(Pi)

        if W is None or W.shape!=(n,n):
            if (W is not None) and W.shape!=(n,n):
                print("Mutation matrix size does not match alphabet size"
                      "Ignoring input mutation matrix")
            # flow matrix
            gtr.W = np.ones((a,a))
            np.fill_diagonal(gtr.W, - ((gtr.W).sum(axis=0) - 1))

        gtr.W = 0.5*(W+W.T)
        gtr._check_fix_Q()
        gtr._eig()
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
            np.fill_diagonal(gtr.W, - ((gtr.W).sum(axis=0) - 1))

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

    @classmethod
    def infer(cls, nij, Ti, root_state, pc=5.0, **kwargs):
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
         is requred, but no alphabet specified, the nucleotide will be used as default.
        """
        gtr = cls(**kwargs)
        dp = 1e-5
        Nit = 40
        pc_mat = pc*np.ones_like(nij)
        pc_mat -= np.diag(np.diag(pc_mat))
        from scipy import linalg as LA
        count = 0
        pi_old = np.zeros_like(Ti)
        pi = np.ones_like(Ti)
        pi/=pi.sum()
        W_ij = np.ones_like(nij)
        mu = nij.sum()/Ti.sum()
        while LA.norm(pi_old-pi) > dp and count < Nit:
            print('GTR inference iteration ',count,'change:',LA.norm(pi_old-pi))
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
            pi = (np.sum(nij+pc_mat,axis=1)+root_state)/(mu*np.dot(W_ij,Ti)+root_state.sum()+np.sum(pc_mat, axis=1))
            pi /= pi.sum()
            mu = nij.sum()/(ttconf.TINY_NUMBER + np.sum(pi * (W_ij.dot(Ti))))
        if count >= Nit:
            print ('WARNING: maximum number of iterations has been reached in GTR inference')
            np.min(pi.sum(axis=0)), np.max(pi.sum(axis=0))
            if LA.norm(pi_old-pi) > dp:
                print ('    the iterative scheme has not converged')
            elif np.abs(1-np.max(pi.sum(axis=0))) > dp:
                print ('    the iterative scheme has converged, but proper normalization was not reached')
        gtr.W = W_ij
        gtr.Pi = np.diag(pi)
        gtr._check_fix_Q()
        gtr._eig()
        return gtr

    def _check_fix_Q(self):
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
        Wdiag = -((self.W.T*np.diagonal(self.Pi)).T).sum(axis=0)/ \
                np.diagonal(self.Pi)
        np.fill_diagonal(self.W, Wdiag)
        Q1 = self.Pi.dot(self.W)
        scale_factor = -np.sum(np.diagonal(Q1*self.Pi))
        self.W /= scale_factor
        self.mu *= scale_factor
        Q1 = self.Pi.dot(self.W)
        if (Q1.sum(axis=0) < 1e-10).sum() <  self.alphabet.shape[0]: # fix failed
            import ipdb; ipdb.set_trace()
            raise ArithmeticError("Cannot fix the diagonal of the GTR rate matrix.")


    def _eig(self):
        """
        Perform eigendecompositon of the rate matrix and stores the left- and right-
        matrices to convert the sequence profiles to the GTR matrix eigenspace
        and hence to speed-up the computations.
        """
        # eigendecomposition of the rate matrix
        eigvals, eigvecs = np.linalg.eig(self.Pi.dot(self.W))
        self.v = np.real(eigvecs)
        self.v_inv = np.linalg.inv(self.v)
        self.eigenmat = np.real(eigvals)
        return

    def prob_t(self, profile_p, profile_ch, t, rotated=False, return_log=False, ignore_gap=True):
        """
        Compute the probability of the two profiles to be separated by the time t.
        Args:
         - profile_p(np.array): parent profile of shape (L, a), where
         L - length of the sequence, a - alphabet size.

         - profile_ch(np.array): child profile of shape (L, a), where
         L - length of the sequence, a - alphabet size.

         - t (double): time (branch len), separating the profiles.

         - rotated (bool, default False): if True, assume that the supplied
         profiles are already rotated.

         - return_log(bool, default False): whether return log-probability.

        Returns:
         - prob(np.array): resulting probability.
        """
        if t < 0:
            if return_log:
                return -BIG_NUMBER
            else:
                return 0.0

        L = profile_p.shape[0]
        if L != profile_ch.shape[0]:
            raise ValueError("Sequence lengths do not match!")

        if ignore_gap:
            try:
                gap_index = np.where(self.alphabet=='-')[0][0]
            except:
                ind = np.ones(L, dtype=bool)
            else:
                ind = (profile_p.argmax(axis=1)!=gap_index)&(profile_ch.argmax(axis=1)!=gap_index)
        else:
                ind = np.ones(L, dtype=bool)

        eLambdaT = self._exp_lt(t)
        if not rotated: # we need to rotate first
            p1 = profile_p.dot(self.v) # (L x a).dot(a x a) = (L x a) - prof in eigenspace
            p2 = (self.v_inv.dot(profile_ch.T)).T # (L x a).dot(a x a) = (L x a) - prof in eigenspace
        else:
            p1 = profile_p
            p2 = profile_ch
            #prob = (profile_p*eLambdaT*profile_ch).sum(1) # sum over the alphabet

        prob = np.maximum(0,(p1*eLambdaT*p2).sum(axis=1)) # sum_i (p1_i * exp(l_i*t) * p_2_i) result = vector length L
        if return_log:
            total_prob = (np.log(prob[ind] + ttconf.TINY_NUMBER)).sum() # sum all sites
        else:
            total_prob = prob[ind].prod() # prod of all sites
        del eLambdaT, prob
        return total_prob

    def optimal_t(self, profile_p, profile_ch, rotated=False, return_log=False, ignore_gap=True):
        """
        Find the optimal distance between the two profiles
        """

        def _neg_prob(t, parent, child):
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
            return -1*self.prob_t (parent, child, t, rotated=False, return_log=True, ignore_gap=ignore_gap)

        try:
            from scipy.optimize import minimize_scalar
            opt = minimize_scalar(_neg_prob,
                    bounds=[0,ttconf.MAX_BRANCH_LENGTH],
                    method='Bounded',
                    args=(profile_p, profile_ch))
            new_len = opt["x"]
        except:
            import scipy
            print('legacy scipy', scipy.__version__)
            from scipy.optimize import fminbound
            new_len = fminbound(_neg_prob,
                    0,ttconf.MAX_BRANCH_LENGTH,
                    args=(profile_p, profile_ch))
            opt={'success':True}


        if new_len > .9 * ttconf.MAX_BRANCH_LENGTH:
            print ("WARNING: The branch length seems to be very long!")

        if opt["success"] != True:
            print ("Cannot optimize branch length, minimization failed. Return Hamming distance")
            new_len =  1 - np.all(profile_ch == profile_p, axis=1).mean()

        return  new_len

    def propagate_profile(self, profile, t, rotated=False, return_log=False):
        """
        Compute the probability of the sequence state (profile) at time (t+t0),
        given the sequence state (profile) at time t0.
        Args:
         - profile(numpy.array): sequence profile. Shape = (L, a),
         where L - sequence length, a - alphabet size.

         - t(double): time to propagate

         - rotated(bool default False): whether the supplied profile is in the
         GTR matrix eigenspace

         - return log (bool, default False): whether to return log-probability

        Returns:
         - res(np.array): profile of the sequence after time t.
         Shape = (L, a), where L - sequence length, a - alphabet size.
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
         - exp_lt(numpy.array): array of values exp(lambda(i) * t),
         where (i) - alphabet index (the eigenvalue number).
        """
        return np.exp(self.mu * t * self.eigenmat)

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
