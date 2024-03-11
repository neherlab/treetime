import numpy as np
from Bio import Seq, SeqRecord


alphabet_synonyms = {'nuc':'nuc', 'nucleotide':'nuc', 'aa':'aa', 'aminoacid':'aa',
                     'nuc_nogap':'nuc_nogap', 'nucleotide_nogap':'nuc_nogap',
                     'aa_nogap':'aa_nogap', 'aminoacid_nogap':'aa_nogap',
                     'DNA':'nuc', 'DNA_nogap':'nuc_nogap'}

alphabets = {
            "nuc":           np.array(['A', 'C', 'G', 'T', '-']),

            "nuc_nogap":np.array(['A', 'C', 'G', 'T']),

            "aa":            np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                                       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                                       'W', 'Y', '*', '-']),

            "aa_nogap": np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                                       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                                       'W', 'Y'])
            }

profile_maps = {
'nuc':{
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
    'W': np.array([1, 0, 0, 1, 0], dtype='float'),
    'K': np.array([0, 0, 1, 1, 0], dtype='float'),
    'M': np.array([1, 1, 0, 0, 0], dtype='float'),
    'D': np.array([1, 0, 1, 1, 0], dtype='float'),
    'H': np.array([1, 1, 0, 1, 0], dtype='float'),
    'B': np.array([0, 1, 1, 1, 0], dtype='float'),
    'V': np.array([1, 1, 1, 0, 0], dtype='float')
    },

'nuc_nogap':{
    'A': np.array([1, 0, 0, 0], dtype='float'),
    'C': np.array([0, 1, 0, 0], dtype='float'),
    'G': np.array([0, 0, 1, 0], dtype='float'),
    'T': np.array([0, 0, 0, 1], dtype='float'),
    '-': np.array([1, 1, 1, 1], dtype='float'), # gaps are completely ignored in distance computations
    'N': np.array([1, 1, 1, 1], dtype='float'),
    'X': np.array([1, 1, 1, 1], dtype='float'),
    'R': np.array([1, 0, 1, 0], dtype='float'),
    'Y': np.array([0, 1, 0, 1], dtype='float'),
    'S': np.array([0, 1, 1, 0], dtype='float'),
    'W': np.array([1, 0, 0, 1], dtype='float'),
    'K': np.array([0, 0, 1, 1], dtype='float'),
    'M': np.array([1, 1, 0, 0], dtype='float'),
    'D': np.array([1, 0, 1, 1], dtype='float'),
    'H': np.array([1, 1, 0, 1], dtype='float'),
    'B': np.array([0, 1, 1, 1], dtype='float'),
    'V': np.array([1, 1, 1, 0], dtype='float')
    },

'aa':{
    'A': np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Alanine         Ala
    'C': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Cysteine        Cys
    'D': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Aspartic AciD   Asp
    'E': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamic Acid   Glu
    'F': np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Phenylalanine   Phe
    'G': np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glycine         Gly
    'H': np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Histidine       His
    'I': np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Isoleucine      Ile
    'K': np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Lysine          Lys
    'L': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Leucine         Leu
    'M': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Methionine      Met
    'N': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #AsparagiNe      Asn
    'P': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Proline         Pro
    'Q': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine       Gln
    'R': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #ARginine        Arg
    'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Serine          Ser
    'T': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype='float'), #Threonine       Thr
    'V': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='float'), #Valine          Val
    'W': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype='float'), #Tryptophan      Trp
    'Y': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype='float'), #Tyrosine        Tyr
    '*': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype='float'), #stop
    '-': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float'), #gap
    'X': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype='float'), #not specified/any
    'B': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Asparagine/Aspartic Acid    Asx
    'Z': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine/Glutamic Acid     Glx
    },

'aa_nogap':{
    'A': np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Alanine         Ala
    'C': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Cysteine        Cys
    'D': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Aspartic AciD   Asp
    'E': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamic Acid   Glu
    'F': np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Phenylalanine   Phe
    'G': np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glycine         Gly
    'H': np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Histidine       His
    'I': np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Isoleucine      Ile
    'K': np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Lysine          Lys
    'L': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Leucine         Leu
    'M': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Methionine      Met
    'N': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #AsparagiNe      Asn
    'P': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Proline         Pro
    'Q': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine       Gln
    'R': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype='float'), #ARginine        Arg
    'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='float'), #Serine          Ser
    'T': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype='float'), #Threonine       Thr
    'V': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype='float'), #Valine          Val
    'W': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype='float'), #Tryptophan      Trp
    'Y': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float'), #Tyrosine        Tyr
    'X': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype='float'), #not specified/any
    'B': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Asparagine/Aspartic Acid    Asx
    'Z': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine/Glutamic Acid     Glx
    }
}


def extend_profile(gtr, aln, logger=None):
    tmp_unique_chars = []
    for seq in aln:
        tmp_unique_chars.extend(np.unique(seq))

    unique_chars = np.unique(tmp_unique_chars)
    for c in unique_chars:
        if c not in gtr.profile_map:
            gtr.profile_map[c] = np.ones(gtr.n_states)
            if logger:
                logger("WARNING: character %s is unknown. Treating it as missing information"%c,1,warn=True)


def guess_alphabet(aln):
    total=0
    nuc_count = 0
    for seq in aln:
        total += len(seq)
        for n in np.array(list('acgtACGT-N')):
            nuc_count += np.sum(seq==n)
    if nuc_count>0.9*total:
        return 'nuc'
    else:
        return 'aa'


def seq2array(seq, word_length=1, convert_upper=False, fill_overhangs=False, ambiguous='N'):
    """
    Take the raw sequence, substitute the "overhanging" gaps with 'N' (missequenced),
    and convert the sequence to the numpy array of chars.

    Parameters
    ----------
    seq : Biopython.SeqRecord, str, iterable
       Sequence as an object of SeqRecord, string or iterable
    word_length : int, optional
        1 for nucleotide or amino acids, 3 for codons etc.
    convert_upper : bool, optional
        convert the sequence to upper case
    fill_overhangs : bool
       If True, substitute the "overhanging" gaps with ambiguous character symbol
    ambiguous : char
       Specify the character for ambiguous state ('N' default for nucleotide)
    Returns
    -------
    sequence : np.array
       Sequence as 1D numpy array of chars
    """
    if isinstance(seq, str):
        seq_str = seq
    elif isinstance(seq, Seq.Seq):
        seq_str = str(seq)
    elif isinstance(seq, SeqRecord.SeqRecord):
        seq_str = str(seq.seq)
    else:
        raise TypeError("seq2array: sequence must be Bio.Seq, Bio.SeqRecord, or string. Got "+str(seq))

    if convert_upper:
        seq_str = seq_str.upper()

    if word_length==1:
        seq_array = np.array(list(seq_str))
    else:
        if len(seq_str)%word_length:
            raise ValueError("sequence length has to be multiple of word length")
        seq_array = np.array([seq_str[i*word_length:(i+1)*word_length]
                              for i in range(len(seq_str)/word_length)])

    # substitute overhanging unsequenced tails
    if fill_overhangs:
        gaps = np.where(seq_array != '-')[0]
        if len(gaps):
            seq_array[:gaps[0]] = ambiguous
            seq_array[gaps[-1]+1:] = ambiguous
        else:
            seq_array[:] = ambiguous

    return seq_array


def seq2prof(seq, profile_map):
    """
    Convert the given character sequence into the profile according to the
    alphabet specified.

    Parameters
    ----------

     seq : numpy.array
        Sequence to be converted to the profile

     profile_map : dic
        Mapping valid characters to profiles

    Returns
    -------

     idx : numpy.array
        Profile for the character. Zero array if the character not found

    """

    return np.array([profile_map[k] for k in seq])


def prof2seq(profile, gtr, sample_from_prof=False, normalize=True, rng=None):
    """
    Convert profile to sequence and normalize profile across sites.

    Parameters
    ----------

     profile : numpy 2D array
        Profile. Shape of the profile should be (L x a), where L - sequence
        length, a - alphabet size.
     gtr : gtr.GTR
        Instance of the GTR class to supply the sequence alphabet
     collapse_prof : bool
        Whether to convert the profile to the delta-function

    Returns
    -------
     seq : numpy.array
        Sequence as numpy array of length L
     prof_values :  numpy.array
        Values of the profile for the chosen sequence characters (length L)
     idx : numpy.array
        Indices chosen from profile as array of length L
    """
    if rng is None:
        rng = np.random.default_rng()
    # normalize profile such that probabilities at each site sum to one
    if normalize:
        tmp_profile, pre=normalize_profile(profile, return_offset=False)
    else:
        tmp_profile = profile

    # sample sequence according to the probabilities in the profile
    # (sampling from cumulative distribution over the different states)
    if sample_from_prof:
        cumdis = tmp_profile.cumsum(axis=1).T
        randnum = rng.random(size=cumdis.shape[1])
        idx = np.argmax(cumdis>=randnum, axis=0)
    else:
        idx = tmp_profile.argmax(axis=1)
    seq = gtr.alphabet[idx]  # max LH over the alphabet

    prof_values = tmp_profile[np.arange(tmp_profile.shape[0]), idx]

    return seq, prof_values, idx


def normalize_profile(in_profile, log=False, return_offset = True):
    """return a normalized version of a profile matrix

    Parameters
    ----------
    in_profile : np.array
        shape Lxq, will be normalized to one across each row
    log : bool, optional
        treat the input as log probabilities
    return_offset : bool, optional
        return the log of the scale factor for each row

    Returns
    -------
    tuple
        normalized profile (fresh np object) and offset (if return_offset==True)
    """
    if log:
        tmp_prefactor = in_profile.max(axis=1)
        tmp_prof = np.exp(in_profile.T - tmp_prefactor).T
    else:
        tmp_prefactor = 0.0
        tmp_prof = in_profile

    norm_vector = tmp_prof.sum(axis=1)
    return (np.einsum('ai,a->ai',tmp_prof,1.0/norm_vector),
            (np.log(norm_vector) + tmp_prefactor) if return_offset else None)

