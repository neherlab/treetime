#!/usr/local/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from .seq_utils import alphabets, profile_maps
from .gtr import GTR

def get_alphabet(a):
    if type(a)==str and a in alphabets:
        return alphabets[a]
    else:
        try:
            return np.array(a)
        except:
            raise TypeError


def JC69 (mu=1.0, alphabet="nuc", **kwargs):
    """
    Jukes-Cantor 1969 model. This model assumes equal concentrations
    of the nucleotides and equal transition rates between nucleotide states.
    For more info, see: Jukes and Cantor (1969). Evolution of Protein Molecules.
                        New York: Academic Press. pp. 21–132

    Parameters
    -----------

     mu : float
        substitution rate

     alphabet : str or character array
        specify alphabet to use.
        Available alphabets are:

         'nuc_nogap' - nucleotides only, gaps ignored

         'nuc' - nucleotide alphabet with gaps, gaps can be ignored optionally

    """
    num_chars = len(get_alphabet(alphabet))
    W, pi = np.ones((num_chars,num_chars)), np.ones(num_chars)
    gtr = GTR(alphabet=alphabet)
    gtr.assign_rates(mu=mu, pi=pi, W=W)
    return gtr

def K80(mu=1., kappa=0.1, **kwargs):
    """
    Kimura 1980 model. Assumes equal concentrations across nucleotides, but
    allows different rates between transitions and transversions. The ratio
    of the transversion/transition rates is given by kappa parameter.
    For more info, see
    Kimura (1980),  J. Mol. Evol. 16 (2): 111–120. doi:10.1007/BF01731581.

    Current implementation of the model does not account for the gaps.

    Parameters
    -----------

     mu : float
        Overall substitution rate

     kappa : float
        Ratio of transversion/transition rates
    """

    num_chars = len(alphabets['nuc_nogap'])
    pi = np.ones(len(alphabets['nuc_nogap']), dtype=float)/len(alphabets['nuc_nogap'])
    W = _create_transversion_transition_W(kappa)
    gtr = GTR(alphabet=alphabets['nuc_nogap'])
    gtr.assign_rates(mu=mu, pi=pi, W=W)
    return gtr

def F81(mu=1.0, pi=None, alphabet="nuc", **kwargs):
    """
    Felsenstein 1981 model. Assumes non-equal concentrations across nucleotides,
    but the transition rate between all states is assumed to be equal. See
    Felsenstein (1981), J. Mol. Evol. 17  (6): 368–376. doi:10.1007/BF01734359
    for details.

    Current implementation of the model does not account for the gaps (treatment of
    gaps as characters is possible if specify alphabet='nuc_gap').

    Parameters
    -----------


     mu : float
        Substitution rate

     pi : numpy.array
        Nucleotide concentrations

     alphabet : str
        Alphabet to use. POsiible values are: ['nuc', 'nuc_gap'] Default 'nuc', which discounts al gaps.
        'nuc_gap' alphabet enables treatmen of gaps as characters.
    """

    if pi is None:
        pi=0.25*np.ones(4, dtype=float)
    num_chars = len(get_alphabet(alphabet))

    pi = np.array(pi, dtype=float)

    if num_chars != len(pi) :
        pi = np.ones((num_chars, ), dtype=float)
        print ("GTR: Warning!The number of the characters in the alphabet does not match the "
            "shape of the vector of equilibrium frequencies Pi -- assuming equal frequencies for all states.")

    W = np.ones((num_chars,num_chars))
    pi /= (1.0 * np.sum(pi))
    gtr = GTR(alphabet=get_alphabet(alphabet))
    gtr.assign_rates(mu=mu, pi=pi, W=W)
    return gtr

def HKY85(mu=1.0, pi=None, kappa=0.1, **kwargs):
    """
    Hasegawa, Kishino and Yano 1985 model. Allows different concentrations of the
    nucleotides (as in F81) + distinguishes between transition/transversionsubstitutions
    (similar to K80). Link:
    Hasegawa, Kishino, Yano (1985), J. Mol. Evol. 22 (2): 160–174. doi:10.1007/BF02101694

    Current implementation of the model does not account for the gaps

    Parameters
    -----------


     mu : float
        Substitution rate

     pi : numpy.array
        Nucleotide concentrations

     kappa : float
        Ratio of transversion/transition substitution rates

    """

    if pi is None:
        pi=0.25*np.ones(4, dtype=float)
    num_chars = len(alphabets['nuc_nogap'])
    if num_chars != pi.shape[0] :
        pi = np.ones((num_chars, ), dtype=float)
        print ("GTR: Warning!The number of the characters in the alphabet does not match the "
            "shape of the vector of equilibrium frequencies Pi -- assuming equal frequencies for all states.")

    W = _create_transversion_transition_W(kappa)
    pi /= pi.sum()
    gtr = GTR(alphabet=alphabets['nuc_nogap'])
    gtr.assign_rates(mu=mu, pi=pi, W=W)
    return gtr

def T92(mu=1.0, pi_GC=0.5, kappa=0.1, **kwargs):
    """
    Tamura 1992 model. Extending Kimura  (1980) model for the case where a
    G+C-content bias exists. Link:
    Tamura K (1992),  Mol.  Biol. Evol. 9 (4): 678–687.  DOI: 10.1093/oxfordjournals.molbev.a040752

    Current implementation of the model does not account for the gaps

    Parameters
    -----------

     mu : float
        substitution rate

     pi_GC : float
        relative GC content

     kappa : float
        relative transversion/transition rate

    """

    W = _create_transversion_transition_W(kappa)
    # A C G T
    if pi_CG >=1.:
        raise ValueError("The relative CG content specified is larger than 1.0!")
    pi = np.array([(1.-pi_CG)*0.5, pi_CG*0.5, pi_CG*0.5, (1-pi_CG)*0.5])
    gtr = GTR(alphabet=alphabets['nuc_nogap'])
    gtr.assign_rates(mu=mu, pi=pi, W=W)
    return gtr

def TN93(mu=1.0, kappa1=1., kappa2=1., pi=None, **kwargs):
    """
    Tamura and Nei 1993. The model distinguishes between the two different types of
    transition: (A <-> G) is allowed to have a different rate to (C<->T).
    Transversions have the same rate. The frequencies of the nucleotides are allowed
    to be different. Link:
    Tamura, Nei (1993), MolBiol Evol. 10 (3): 512–526. DOI:10.1093/oxfordjournals.molbev.a040023

    Parameters
    -----------

     mu : float
        Substitution rate

     kappa1 : float
        relative A<-->C, A<-->T, T<-->G and G<-->C rates

     kappa2 : float
        relative C<-->T rate

    Note
    ----

     Rate of A<-->G substitution is set to one. All other rates (kappa1, kappa2)
    are specified relative to this rate

    """

    if pi is None:
        pi=0.25*np.ones(4, dtype=float)
    W = np.ones((4,4))
    W = np.array([
        [1,      kappa1, 1,      kappa1],
        [kappa1, 1,      kappa1, kappa2],
        [1,      kappa1, 1,      kappa1],
        [kappa1, kappa2, kappa1,  1]], dtype=float)

    pi /=pi.sum()
    num_chars = len(alphabets['nuc_nogap'])
    if num_chars != pi.shape[0] :
        pi = np.ones((num_chars, ), dtype=float)
        print ("GTR: Warning!The number of the characters in the alphabet does not match the "
            "shape of the vector of equilibrium frequencies Pi -- assuming equal frequencies for all states.")

    gtr = GTR(alphabet=alphabets['nuc'])
    gtr.assign_rates(mu=mu, pi=pi, W=W)
    return gtr

def _create_transversion_transition_W(kappa):
    """
    Alphabet = [A, C, G, T]
    """
    W = np.ones((4,4))
    W[0, 2]=W[1, 3]=W[2, 0]=W[3,1]=kappa
    return W

if __name__ == '__main__':
    pass

