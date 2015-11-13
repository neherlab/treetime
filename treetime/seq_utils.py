import numpy as np

alphabets = {
            "nuc": np.array(['A', 'C', 'G', 'T', '-']),
            "aa": np.array(['-'])}

full_nc_profile = {
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

def seq2prof(x, aln_type='nuc'):
    """
    Convert the given character sequence into the profileaccording to the alphabet specified.
    
    Args:
     
     - x(numpy.array): sequence to be converted to the profile
     
     - alph(str): alphabet type. Can be either 'nuc' for nucleotide alphabet, 
     or 'aa' for amino acid alphabet
    
    Returns:
     
     - idx(numpy.array): profile for the character, zero array if the character not found
    """
    
    if aln_type=='nuc':
        prof = np.array([full_nc_profile[k]
            if k in full_nc_profile else full_nc_profile['N'] for k in x])
        err = ((prof == 0.2).sum(1) != 0).sum()
        if err>0:
            print ("Seq2Profile: %d of %d characters were not identified or"
                    " not sequenced." % (err, prof.shape[0]))
    elif aln_type=='aa':
        raise NotImplementedError("Amino-acid alphabet is under development.")
    else:
        raise TypeError("Alignment type cane be either 'nuc' or 'aa'")
    return prof

def prof2seq(profile, gtr, correct_prof=True):
    """
    Convert profile to sequence and, if requested, set the profile values (LH of
    the characters) to zeros and ones essentially converting the character
    distribution to the delta-function.

    Args:
     - profile(numpy 2D array): profile. Shape of the profile should be
     (L x a), where L - sequence length, a - alphabet size.
     - gtr (gtr.GTR) instance of teh GTR class to supply the sequence alphabet
     - corerct_prof(bool, default True): whether to conver t the profile to the
     delta-function

    Returns:
     -seq (numpy array of length L): sequence
     - profile(numpy 2D array of Lxa shape): the resulting profile.
    """
    seq = gtr.alphabet[profile.argmax(axis=1)]  # max LH over the alphabet
    if correct_prof:  # max profile value to one, others - zeros
        am = profile.argmax(axis=1)
        profile[:, :] = 0.0
        profile[np.arange(profile.shape[0]), am] = 1.0
    return seq, profile
