import numpy as np
from Bio import Seq, SeqRecord


alphabet_synonyms = {'nuc':'nuc', 'nucleotide':'nuc', 'aa':'aa', 'aminoacid':'aa',
                     'nuc_nogap':'nuc_nogap', 'nucleotide_nogap':'nuc_nogap',
                     'aa_nogap':'aa_nogap', 'aminoacid_nogap':'aa_nogap',
                     'DNA':'nuc', 'DNA_nogap':'nuc_nogap'}

alphabets = {
            "nuc":           np.array(['A', 'C', 'G', 'T', '-'], dtype='S1'),

            "nuc_nogap":np.array(['A', 'C', 'G', 'T'], dtype='S1'),

            "aa":            np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                                       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                                       'W', 'Y', '*', '-'], dtype='S1'),

            "aa_nogap": np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                                       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                                       'W', 'Y'], dtype='S1')
            }

profile_maps = {
'nuc':{
    b'A': np.array([1, 0, 0, 0, 0], dtype='float'),
    b'C': np.array([0, 1, 0, 0, 0], dtype='float'),
    b'G': np.array([0, 0, 1, 0, 0], dtype='float'),
    b'T': np.array([0, 0, 0, 1, 0], dtype='float'),
    b'-': np.array([0, 0, 0, 0, 1], dtype='float'),
    b'N': np.array([1, 1, 1, 1, 1], dtype='float'),
    b'X': np.array([1, 1, 1, 1, 1], dtype='float'),
    b'R': np.array([1, 0, 1, 0, 0], dtype='float'),
    b'Y': np.array([0, 1, 0, 1, 0], dtype='float'),
    b'S': np.array([0, 1, 1, 0, 0], dtype='float'),
    b'W': np.array([1, 0, 0, 1, 0], dtype='float'),
    b'K': np.array([0, 0, 1, 1, 0], dtype='float'),
    b'M': np.array([1, 1, 0, 0, 0], dtype='float'),
    b'D': np.array([1, 0, 1, 1, 0], dtype='float'),
    b'H': np.array([1, 1, 0, 1, 0], dtype='float'),
    b'B': np.array([0, 1, 1, 1, 0], dtype='float'),
    b'V': np.array([1, 1, 1, 0, 0], dtype='float')
    },

'nuc_nogap':{
    b'A': np.array([1, 0, 0, 0], dtype='float'),
    b'C': np.array([0, 1, 0, 0], dtype='float'),
    b'G': np.array([0, 0, 1, 0], dtype='float'),
    b'T': np.array([0, 0, 0, 1], dtype='float'),
    b'-': np.array([1, 1, 1, 1], dtype='float'), # gaps are completely ignored in distance computations
    b'N': np.array([1, 1, 1, 1], dtype='float'),
    b'X': np.array([1, 1, 1, 1], dtype='float'),
    b'R': np.array([1, 0, 1, 0], dtype='float'),
    b'Y': np.array([0, 1, 0, 1], dtype='float'),
    b'S': np.array([0, 1, 1, 0], dtype='float'),
    b'W': np.array([1, 0, 0, 1], dtype='float'),
    b'K': np.array([0, 0, 1, 1], dtype='float'),
    b'M': np.array([1, 1, 0, 0], dtype='float'),
    b'D': np.array([1, 0, 1, 1], dtype='float'),
    b'H': np.array([1, 1, 0, 1], dtype='float'),
    b'B': np.array([0, 1, 1, 1], dtype='float'),
    b'V': np.array([1, 1, 1, 0], dtype='float')
    },

'aa':{
    b'A': np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Alanine         Ala
    b'C': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Cysteine        Cys
    b'D': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Aspartic AciD   Asp
    b'E': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamic Acid   Glu
    b'F': np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Phenylalanine   Phe
    b'G': np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glycine         Gly
    b'H': np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Histidine       His
    b'I': np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Isoleucine      Ile
    b'K': np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Lysine          Lys
    b'L': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Leucine         Leu
    b'M': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Methionine      Met
    b'N': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #AsparagiNe      Asn
    b'P': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Proline         Pro
    b'Q': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine       Gln
    b'R': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #ARginine        Arg
    b'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Serine          Ser
    b'T': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype='float'), #Threonine       Thr
    b'V': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='float'), #Valine          Val
    b'W': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype='float'), #Tryptophan      Trp
    b'Y': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype='float'), #Tyrosine        Tyr
    b'*': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0], dtype='float'), #stop
    b'-': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float'), #gap
    b'X': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype='float'), #not specified/any
    b'B': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Asparagine/Aspartic Acid    Asx
    b'Z': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine/Glutamic Acid     Glx
    },

'aa_nogap':{
    b'A': np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Alanine         Ala
    b'C': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Cysteine        Cys
    b'D': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Aspartic AciD   Asp
    b'E': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamic Acid   Glu
    b'F': np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Phenylalanine   Phe
    b'G': np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glycine         Gly
    b'H': np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Histidine       His
    b'I': np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Isoleucine      Ile
    b'K': np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Lysine          Lys
    b'L': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Leucine         Leu
    b'M': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Methionine      Met
    b'N': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #AsparagiNe      Asn
    b'P': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Proline         Pro
    b'Q': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine       Gln
    b'R': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype='float'), #ARginine        Arg
    b'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='float'), #Serine          Ser
    b'T': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype='float'), #Threonine       Thr
    b'V': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype='float'), #Valine          Val
    b'W': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype='float'), #Tryptophan      Trp
    b'Y': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float'), #Tyrosine        Tyr
    b'X': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype='float'), #not specified/any
    b'B': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Asparagine/Aspartic Acid    Asx
    b'Z': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine/Glutamic Acid     Glx
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
        for n in np.fromstring('acgtACGT-N', 'S1'):
            nuc_count += np.sum(seq==n)
    if nuc_count>0.9*total:
        return 'nuc'
    else:
        return 'aa'


def seq2array(seq, word_length=1, convert_upper=False, fill_overhangs=False, ambiguous=b'N'):
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
        seq_str = np.fromstring(str(seq), 'S%d'%word_length)
    elif isinstance(seq, SeqRecord.SeqRecord):
        seq_str = str(seq.seq)
    else:
        raise TypeError("seq2array: sequence must be Bio.Seq, Bio.SeqRecord, or string, got "+str(seq))

    if convert_upper:
        seq_str = seq_str.upper()

    seq_array = np.fromstring(seq_str, 'S%d'%word_length)
    # substitute overhanging unsequenced tails
    if fill_overhangs:
        gaps = np.where(seq_array != b'-')[0]
        seq_array[:gaps[0]] = ambiguous
        seq_array[gaps[-1]+1:] = ambiguous

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


def prof2seq(profile, gtr, sample_from_prof=False, normalize=True):
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

    # normalize profile such that probabilities at each site sum to one
    if normalize:
        tmp_profile, pre=normalize_profile(profile, return_offset=False)
    else:
        tmp_profile = profile

    # sample sequence according to the probabilities in the profile
    # (sampling from cumulative distribution over the different states)
    if sample_from_prof:
        cumdis = tmp_profile.cumsum(axis=1).T
        randnum = np.random.random(size=cumdis.shape[1])
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
    return (np.copy(np.einsum('ai,a->ai',tmp_prof,1.0/norm_vector)),
            (np.log(norm_vector) + tmp_prefactor) if return_offset else None)

