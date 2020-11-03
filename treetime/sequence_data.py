from __future__ import division, print_function, absolute_import
import sys
from os.path import isfile
from collections import defaultdict
import numpy as np
from Bio import SeqRecord, Seq, AlignIO, SeqIO
from treetime import config as ttconf
from treetime import MissingDataError
from .seq_utils import seq2array, guess_alphabet, alphabets

string_types = [str] if sys.version_info[0]==3 else [str, unicode]
def simple_logger(*args, **kwargs):
    print(args)

class SequenceData(object):
    """docstring for SeqData

    Attributes
    ----------
    additional_constant_sites : int
        length of the sequence without variation not included in the alignment
    aln : dict
        sequences, either sparse of full
    ambiguous : byte
        character signifying missing data
    compress : bool
        compress the alignment
    compressed_alignment : dict
        dictionary mapping sequence names to compressed sequences
    compressed_to_full_sequence_map : dict
        for each compressed position, contain a list of positions in the full alignment
    fill_overhangs : bool
        treat gaps at either end of sequence as missing data
    full_length : int
        length of the sequence
    full_to_compressed_sequence_map : np.array
        a map of each position in the full sequence to the compressed sequence
    inferred_const_sites : list
        list of positions that are constant but differ from the reference, or contain ambiguous characters
    is_sparse : bool
        whether the representation of the alignment is sparse (dict) or fill (array)
    likely_alphabet : str
        simply guess as to whether the sequence alignment is nucleotides or amino acids
    logger : callable
        function writting log messages
    multiplicity : np.array
        specifies for each column of the compressed alignment how often this pattern occurs
    nonref_positions : list
        positions where at least one sequence differs from the reference
    ref : np.array
        reference sequence (stored as np.array(dtype="S"))
    seq_multiplicity : dict
        store the multiplicity of sequence, for example read count in a deep sequencing experiment
    sequence_names : list
        list of all sequences in a fixed order
    word_length : int
        length of state (typically 1 A,C,G,T, but could be 3 for codons)
    """
    def __init__(self, aln, ref=None, logger=None, convert_upper=True,
                 sequence_length=None, compress=True, word_length=1, sequence_type=None,
                 fill_overhangs=True, seq_multiplicity=None, ambiguous=None, **kwargs):
        """construct an sequence data object

        Parameters
        ----------
        aln : Bio.Align.MultipleSeqAlignment, str
            alignment or file name
        ref : Seq, str
            sequence or file name
        logger : callable, optional
            logging function
        convert_upper : bool, optional
            convert all sequences to upper case, default true
        sequence_length : None, optional
            length of the sequence, only necessary when no alignment or ref is given
        compress : bool, optional
            compress identical alignment columns into one
        word_length : int
            length of state (typically 1 A,C,G,T, but could be 3 for codons)
        fill_overhangs : bool
            treat gaps at either end of sequence as missing data
        seq_multiplicity : dict
            store the multiplicity of sequence, for example read count in a deep sequencing experiment
        ambiguous : byte
            character signifying missing data
        **kwargs
            Description
        """
        self.logger = logger if logger else simple_logger
        self._aln = None
        self._ref = None
        self.likely_alphabet = None
        self.compressed_to_full_sequence_map = None
        self.multiplicity = None
        self.is_sparse = None
        self.convert_upper = convert_upper
        self.compress = compress
        self.seq_multiplicity = seq_multiplicity or {} # possibly a dict mapping sequences to their read cound/sample count
        self.additional_constant_sites = kwargs['additional_constant_sites'] if 'additional_constant_sites' in kwargs else 0

        # if not specified, this will be set as the alignment_length or reference length
        self._full_length = None
        self.full_length = sequence_length
        self._compressed_length = None
        self.word_length = word_length
        self.fill_overhangs = fill_overhangs
        self.ambiguous = ambiguous
        self.sequence_type = sequence_type

        self.ref = ref
        self.aln = aln


    @property
    def aln(self):
        """
        The multiple sequence alignment currently used by the TreeAnc

        :setter: Takes in alignment as MultipleSeqAlignment, str, or dict/defaultdict \
        and attaches sequences to tree nodes.
        :getter: Returns alignment as MultipleSeqAlignment or dict/defaultdict

        """
        return self._aln


    @aln.setter
    def aln(self,in_aln):
        """
        Reads in the alignment (from a dict, MultipleSeqAlignment, or file,
        as necessary), sets tree-related parameters, and attaches sequences
        to the tree nodes.

        Parameters
        ----------
        in_aln : MultipleSeqAlignment, str, dict/defaultdict
            The alignment to be read in

        """
        # load alignment from file if necessary


        from Bio.Align import MultipleSeqAlignment
        self._aln, self.is_sparse = None, None
        if in_aln is None:
            return
        elif type(in_aln) in [defaultdict, dict]:  #if input is sparse (i.e. from VCF)
            self._aln = in_aln
            self.is_sparse = True
        elif type(in_aln) in string_types and isfile(in_aln):
            if any([in_aln.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]) and (self.ref is not None):
                from .vcf_utils import read_vcf
                compress_seq = read_vcf(in_aln)
                in_aln = compress_seq['sequences']
            else:
                for fmt in ['fasta', 'phylip-relaxed', 'nexus']:
                    try:
                        in_aln=AlignIO.read(in_aln, fmt)
                    except:
                        continue

        if type(in_aln) is MultipleSeqAlignment:
            # check whether the alignment is consistent with a nucleotide alignment.
            self._aln = {s.name: seq2array(s, convert_upper=self.convert_upper,
                                           fill_overhangs=self.fill_overhangs, ambiguous=self.ambiguous)
                         for s in in_aln}
            self.check_alphabet(list(self._aln.values()))
            self.is_sparse = False
            self.logger("SequenceData: loaded alignment.",1)
        elif type(in_aln) in [dict, defaultdict]:
            self.logger("SequenceData: loaded sparse/vcf alignment.",1)
            self.check_alphabet([self.ref])
            self.is_sparse = True
            self._aln = in_aln
        else:
            raise MissingDataError("SequenceData: loading alignment failed... " + str(in_aln))

        if self.full_length:
            if self.is_sparse:
                if self.full_length!=len(self.ref):
                    self.logger("SequenceData.aln: specified sequence length doesn't match reference length, ignoring sequence length.", 1, warn=True)
                    self._full_length = len(self.ref)
            else:
                if self.full_length < in_aln.get_alignment_length():
                    raise AttributeError("SequenceData.aln: specified sequence length is smaller than alignment length!")
                elif self.full_length > in_aln.get_alignment_length():
                    self.logger("SequenceData.aln: specified sequence length doesn't match alignment length. Treating difference as constant sites.", 2, warn=True)
                    self.additional_constant_sites = max(0, self.full_length - in_aln.get_alignment_length())
        else:
            if self.is_sparse:
                self.full_length = len(self.ref)
            else:
                self.full_length = in_aln.get_alignment_length()

        self.sequence_names = list(self.aln.keys())

        self.make_compressed_alignment()


    @property
    def full_length(self):
        """length of the uncompressed sequence
        """
        return self._full_length


    @full_length.setter
    def full_length(self,L):
        """set the length of the uncompressed sequence. its inverse 'one_mutation'
        is frequently used as a general length scale. This can't be changed once
        it is set.

        Parameters
        ----------
        L : int
            length of the sequence alignment
        """
        if (not hasattr(self, '_full_length')) or self._full_length is None:
            if L:
                self._full_length = int(L)
        else:
            self.logger("Alignment: one_mutation and sequence length can only be specified once!",1)

    @property
    def compressed_length(self):
        return self._compressed_length


    @property
    def ref(self):
        """
        :setter: Sets the string reference sequence
        :getter: Returns the string reference sequence
        """
        return self._ref


    @ref.setter
    def ref(self, in_ref):
        """
        Parameters
        ----------
        in_ref : file name, str, Bio.Seq.Seq, Bio.SeqRecord.SeqRecord
            reference sequence will read and stored a byte array
        """
        read_from_file=False
        if in_ref and isfile(in_ref):
            for fmt in ['fasta', 'genbank']:
                try:
                    in_ref = SeqIO.read(in_ref, fmt)
                    self.logger("SequenceData: loaded reference sequence as %s format"%fmt,1)
                    read_from_file=True
                    break
                except:
                    continue
            if not read_from_file:
                raise TypeError('SequenceData.ref: reference sequence file %s could not be parsed, fasta and genbank formats are supported.')

        if in_ref:
            self._ref = seq2array(in_ref, fill_overhangs=False, word_length=self.word_length)
            self.full_length = self._ref.shape[0]
            self.compressed_to_full_sequence_map = None
            self.multiplicity = None


    def check_alphabet(self, seqs):
        self.likely_alphabet = guess_alphabet(seqs)

        if self.sequence_type:
            if self.likely_alphabet!=self.sequence_type:
                if self.sequence_type=='nuc':
                    self.logger("POSSIBLE ERROR: This does not look like a nucleotide alignment!", 0, warn=True)
                elif self.sequence_type=='aa':
                    self.logger("POSSIBLE ERROR: This looks like a nucleotide alignment, you indicated amino acids!", 0, warn=True)

        if self.ambiguous is None:
            self.ambiguous = 'N' if self.likely_alphabet=='nuc' else 'X'


    def make_compressed_alignment(self):
        """
        Create the compressed alignment from the full sequences. This method counts
        the multiplicity for each column of the alignment ('alignment pattern'), and
        creates the compressed alignment, where only the unique patterns are present.
        The maps from full sequence to compressed sequence and back are also stored to allow
        compressing and expanding the sequences.

        Notes
        -----
          full_to_compressed_sequence_map : (array)
             Map to reduce a sequence
          compressed_to_full_sequence_map : (dict)
             Map to restore sequence from compressed alignment
          multiplicity : (array)
            Numpy array, which stores the pattern multiplicity for each position of the compressed alignment.
          compressed_alignment : (2D numpy array)
            The compressed alignment. Shape is (N x L'), where N is number of
            sequences, L' - number of unique alignment patterns
        """

        if not self.compress: #
            self.multiplicity = np.ones(self.full_length, dtype=float)
            self.full_to_compressed_sequence_map = np.arange(self.full_length)
            self.compressed_to_full_sequence_map = {p:np.array([p]) for p in np.arange(self.full_length)}
            self._compressed_length = self._full_length
            self.compressed_alignment = self._aln
            return ttconf.SUCCESS

        self.logger("SeqData: making compressed alignment...", 1)
        # bind positions in full length sequence to that of the compressed (compressed) sequence
        self.full_to_compressed_sequence_map = np.zeros(self.full_length, dtype=int)
        # bind position in compressed sequence to the array of positions in full length sequence
        self.compressed_to_full_sequence_map = {}

        #if alignment is sparse, don't iterate over all invarible sites.
        #so pre-load alignment_patterns with the location of const sites!
        #and get the sites that we want to iterate over only!
        if self.is_sparse:
            from .vcf_utils import process_sparse_alignment
            tmp = process_sparse_alignment(self.aln, self.ref, self.ambiguous)
            compressed_aln_transpose = tmp["constant_columns"]
            alignment_patterns = tmp["constant_patterns"]
            variable_positions = tmp["variable_positions"]
            self.inferred_const_sites = tmp["constant_up_to_ambiguous"]
            self.nonref_positions = tmp["nonref_positions"]
        else: # transpose real alignment, for ease of iteration
            alignment_patterns = {}
            compressed_aln_transpose = []
            aln_transpose = np.array([self.aln[k] for k in self.sequence_names]).T
            variable_positions = np.arange(aln_transpose.shape[0])

        for pi in variable_positions:
            if self.is_sparse:
                pattern = np.array([self.aln[k][pi] if pi in self.aln[k] else self.ref[pi]
                           for k in self.sequence_names], dtype='S')
            else:
                # pylint: disable=unsubscriptable-object
                pattern = np.copy(aln_transpose[pi])

            # if the column contains only one state and ambiguous nucleotides, replace
            # those with the state in other strains right away
            unique_letters = list(np.unique(pattern))
            if len(unique_letters)==2 and self.ambiguous in unique_letters:
                other = [c for c in unique_letters if c!=self.ambiguous][0]
                #also replace in original pattern!
                pattern[pattern == self.ambiguous] = other
                unique_letters = [other]

            str_pattern = "".join(pattern.astype('U'))
            # if there is a mutation in this column, give it its private pattern
            # this is required when sampling mutations from reconstructed profiles.
            # otherwise, all mutations corresponding to the same pattern will be coupled.
            # FIXME: this could be done more efficiently
            if len(unique_letters)>1:
                str_pattern += '_%d'%pi

            # if the pattern is not yet seen,
            if str_pattern not in alignment_patterns:
                # bind the index in the compressed aln, index in sequence to the pattern string
                alignment_patterns[str_pattern] = (len(compressed_aln_transpose), [pi])
                # append this pattern to the compressed alignment
                compressed_aln_transpose.append(pattern)
            else:
                # if the pattern is already seen, append the position in the real
                # sequence to the compressed aln<->sequence_pos_indexes map
                alignment_patterns[str_pattern][1].append(pi)

        # add constant alignment column not in the alignment. We don't know where they
        # are, so just add them to the end. First, determine sequence composition.
        if self.additional_constant_sites:
            character_counts = {c:np.sum(aln_transpose==c) for c in alphabets[self.likely_alphabet+'_nogap']
                                if c not in [self.ambiguous, '-']}
            total = np.sum(list(character_counts.values()))
            additional_columns_per_character = [(c,int(np.round(self.additional_constant_sites*n/total)))
                                  for c, n in character_counts.items()]
            columns_left = self.additional_constant_sites
            pi = np.max(variable_positions)+1
            for c,n in additional_columns_per_character:
                if c==additional_columns_per_character[-1][0]:  # make sure all additions add up to the correct number to avoid rounding
                    n = columns_left
                str_pattern = c*len(self.sequence_names)
                pos_list = list(range(pi, pi+n))
                if n:
                    if str_pattern in alignment_patterns:
                        alignment_patterns[str_pattern][1].extend(pos_list)
                    else:
                        alignment_patterns[str_pattern] = (len(compressed_aln_transpose), pos_list)
                        compressed_aln_transpose.append(np.array(list(str_pattern)))
                    pi += n
                    columns_left -= n


        # count how many times each column is repeated in the real alignment
        self.multiplicity = np.zeros(len(alignment_patterns))
        for p, pos in alignment_patterns.values():
            self.multiplicity[p]=len(pos)

        # create the compressed alignment as a dictionary linking names to sequences
        tmp_compressed_alignment = np.array(compressed_aln_transpose).T
        # pylint: disable=unsubscriptable-object
        self.compressed_alignment = {k: tmp_compressed_alignment[i]
                                  for i,k in enumerate(self.sequence_names)}

        # create map to compress a sequence
        for p, pos in alignment_patterns.values():
            self.full_to_compressed_sequence_map[np.array(pos)]=p

        # create a map to reconstruct full sequence from the compressed (compressed) sequence
        for p, val in alignment_patterns.items():
            self.compressed_to_full_sequence_map[val[0]]=np.array(val[1], dtype=int)

        self.logger("SequenceData: constructed compressed alignment...", 1)
        self._compressed_length = len(self.multiplicity)
        return ttconf.SUCCESS


    def full_to_sparse_sequence(self, sequence):
        """turn a sequence into a dictionary of differences from a reference sequence

        Parameters
        ----------
        sequence : str, numpy.ndarray
            sequence to convert

        Returns
        -------
        dict
            dictionary of difference from reference
        """
        if self.ref is None:
            raise TypeError("SequenceData: sparse sequences can only be constructed when a reference sequence is defined")
        if type(sequence) is not np.ndarray:
            aseq = seq2array(sequence, fill_overhangs=False)
        else:
            aseq = sequence
        differences = np.where(self.ref!=aseq)[0]
        return {p:aseq[p] for p in differences}


    def compressed_to_sparse_sequence(self, sequence):
        """turn a compressed sequence into a list of difference from a reference

        Parameters
        ----------
        sequence : numpy.ndarray
            compressed sequence stored as array

        Returns
        -------
        dict
            dictionary of difference from reference
        """
        if self.ref is None:
            raise TypeError("SequenceData: sparse sequences can only be constructed when a reference sequence is defined")
        sparse_seq = {}
        for pos in self.nonref_positions:
            cseqLoc = self.full_to_compressed_sequence_map[pos]
            base = sequence[cseqLoc]
            if self.ref[pos] != base:
                sparse_seq[pos] = base

        return sparse_seq


    def compressed_to_full_sequence(self, sequence, include_additional_constant_sites=False, as_string=False):
        """expand a compressed sequence

        Parameters
        ----------
        sequence : np.ndarray
            compressed sequence
        include_additional_constant_sites : bool, optional
            add sites assumed constant
        as_string : bool, optional
            return a string instead of an array

        Returns
        -------
        array,str
            expanded sequence
        """
        if include_additional_constant_sites:
            L = self.full_length
        else:
            L = self.full_length - self.additional_constant_sites

        tmp_seq = sequence[self.full_to_compressed_sequence_map[:L]]
        if as_string:
            return "".join(tmp_seq.astype('U'))
        else:
            return tmp_seq

    def differences(self, seq1, seq2, seq1_compressed=True, seq2_compressed=True):
        diffs = []
        if self.is_sparse:
            if seq1_compressed: seq1 = self.compressed_to_sparse_sequence(seq1)
            if seq2_compressed: seq2 = self.compressed_to_sparse_sequence(seq2)

            for pos in set(seq1.keys()).union(seq2.keys()):
                ref_state = self.ref[pos]
                s1 = seq1.get(pos, ref_state)
                s2 = seq2.get(pos, ref_state)
                if s1!=s2:
                    diffs.append((s1,pos,s2))
        else:
            if seq1_compressed: seq1 = self.compressed_to_full_sequence(seq1)
            if seq2_compressed: seq2 = self.compressed_to_full_sequence(seq2)
            diff_pos = np.where(seq1 != seq2)[0]
            for pos in diff_pos:
                diffs.append((seq1[pos], pos, seq2[pos]))

        return sorted(diffs, key=lambda x:x[1])
