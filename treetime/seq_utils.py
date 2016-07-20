import numpy as np

alphabets = {
            "nuc": np.array(['A', 'C', 'G', 'T', '-']),
            "aa": np.array(['-'])}

profile_maps = {'nuc':{
    'A': np.array([1, 0, 0, 0, 0], dtype='float'),
    'C': np.array([0, 1, 0, 0, 0], dtype='float'),
    'G': np.array([0, 0, 1, 0, 0], dtype='float'),
    'T': np.array([0, 0, 0, 1, 0], dtype='float'),
    '-': np.array([0, 0, 0, 0, 1], dtype='float'),
    'N': np.array([1, 1, 1, 1, 1], dtype='float'),
    'X': np.array([1, 1, 1, 1, 1], dtype='float'),
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
    }

def prepare_seq(seq):
    """
    Take the raw sequence, substitute the "overhanging" gaps with 'N' (missequenced)
    convert the sequence to the numpy array of uppercase chars.

    Args:
     - seq:  sequence as an object of SeqRecord, string or iterable

    Returns:
     - sequence(np.array): sequence as 1D numpy array of chars
    """
    try:
        sequence = ''.join(seq)
    except TypeError:
        sequence = seq
    sequence = sequence.upper()

    sequence = np.array(list(sequence))
    # substitute overhanging unsequenced tails
    sequence [:np.where(sequence != '-')[0][0]] = 'N'
    sequence [np.where(sequence != '-')[0][-1]+1:] = 'N'
    return sequence

def seq2prof(x, profile_map):
    """
    Convert the given character sequence into the profileaccording to the alphabet specified.

    Args:

     - x(numpy.array): sequence to be converted to the profile

     - profile map mapping valid characters to profiles

    Returns:

     - idx(numpy.array): profile for the character, zero array if the character not found
    """
    n_states = len(profile_map.values()[0])
    prof = np.array([profile_map[k] if k in profile_map
                    else np.ones(n_states) for k in x ])

    return prof

def prof2seq(profile, gtr, sample_from_prof=True, correct_prof=False):
    """
    Convert profile to sequence and, if requested, set the profile values (LH of
    the characters) to zeros and ones essentially converting the character
    distribution to the delta-function.

    Args:
     - profile(numpy 2D array): profile. Shape of the profile should be
     (L x a), where L - sequence length, a - alphabet size.
     - gtr (gtr.GTR) instance of teh GTR class to supply the sequence alphabet
     - correct_prof(bool, default True): whether to convert the profile to the
     delta-function

    Returns:
     -seq (numpy array of length L): sequence
     - profile(numpy 2D array of Lxa shape): the resulting profile.
    """


    if correct_prof:  # max profile value to one, others - zeros
        am = profile.argmax(axis=1)
        profile[:, :] = 0.0
        profile[np.arange(profile.shape[0]), am] = 1.0
    else: # only normalize
        profile=(profile.T/profile.sum(axis=1)).T

    if sample_from_prof:
        cumdis = profile.cumsum(axis=1).T
        randnum = np.random.random(size=cumdis.shape[1])
        seq = gtr.alphabet[np.argmax(cumdis>=randnum, axis=0)]
    else:
        seq = gtr.alphabet[profile.argmax(axis=1)]  # max LH over the alphabet

    return seq, profile
