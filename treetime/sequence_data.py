from __future__ import division, print_function, absolute_import
import sys
from os.path import isfile
from collections import defaultdict
import numpy as np
from Bio import SeqRecord, Seq, AlignIO, SeqIO
from treetime import config as ttconf
from .seq_utils import seq2array, guess_alphabet

string_types = [str] if sys.version_info[0]==3 else [str, unicode]
def simple_logger(*args, **kwargs):
    print(args)

class SequenceData(object):
    """docstring for SeqData"""
    def __init__(self, aln, ref=None, logger=None, convert_upper=True,
                 sequence_length=None, reduce_alignment=True, word_length=1,
                 fill_overhangs=True, seq_multiplicity=None, ambiguous=None, **kwargs):
        self.logger = logger if logger else simple_logger
        self._aln = None
        self._ref = None
        self.likely_alphabet = None
        self.reduced_to_full_sequence_map = None
        self.pattern_multiplicity = None
        self.is_sparse = None
        self.reduce_alignment = reduce_alignment
        self.seq_multiplicity = seq_multiplicity or {} # possibly a dict mapping sequences to their read cound/sample count
        self.additional_constant_sites = kwargs['additional_constant_sites'] if 'additional_constant_sites' in kwargs else 0

        # if not specified, this will be set as the alignment_length or reference length
        self._full_length = sequence_length
        self._reduced_length = None
        self.word_length = word_length
        self.fill_overhangs = fill_overhangs
        self.ambiguous = ambiguous

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
            self.is_sparse = False
            self._aln = {s.name: seq2array(s) for s in in_aln}
            self.logger("SequenceData: loaded alignment.",1)
        elif type(in_aln) in [dict, defaultdict]:
            self.logger("SequenceData: loaded sparse/vcf alignment.",1)
            self.is_sparse = True
            self._aln = in_aln

        if self._aln is None:
            self.logger("SequenceData: loading alignment failed... " + str(in_aln),1, warn=True)
            return ttconf.ERROR

        if self.full_length:
            if self.is_sparse:
                if self.full_length!=len(self.ref):
                    self.logger("SequenceData.aln: specified sequence length doesn't match reference length, ignoring sequence length.", 1, warn=True)
                    self._full_length = len(self.ref)
            else:
                if self.full_length!= in_aln.get_alignment_length():
                    self.logger("SequenceData.aln: specified sequence length doesn't match alignment length. Treating difference as constant sites.", 2, warn=True)
                    self.additional_constant_sites = max(0, self.full_length - in_aln.get_alignment_length())
        else:
            if self.is_sparse:
                self.full_length = len(self.ref)
            else:
                self.full_length = in_aln.get_alignment_length()

        self.sequence_names = list(self.aln.keys())

        # check whether the alignment is consistent with a nucleotide alignment.
        self.likely_alphabet = guess_alphabet([self.ref] if self.is_sparse
                                               else [s for s in self.aln.values()])
        if self.ambiguous is None:
            self.ambiguous = b'N' if self.likely_alphabet=='nuc' else b'X'

        self.make_reduced_alignment()


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
            self.logger("Alignment: one_mutation and sequence length can't be reset",1)

    @property
    def reduced_length(self):
        return self._reduced_length


    @property
    def ref(self):
        """
        Get the str reference nucleotide sequence currently used by TreeAnc.
        When having read alignment in from a VCF, this is what variants map to.

        :setter: Sets the string reference sequence
        :getter: Returns the string reference sequence

        """
        # delete previous alignment if reference changes
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
            self.reduced_to_full_sequence_map = None
            self.multiplicity = None


    def make_reduced_alignment(self):
        """
        Create the reduced alignment from the full sequences attached to (some)
        tree nodes. The methods collects all sequences from the tree nodes, creates
        the alignment, counts the multiplicity for each column of the alignment
        ('alignment pattern'), and creates the reduced alignment, where only the
        unique patterns are present. The reduced alignment and the pattern multiplicity
        are sufficient for the GTR calculations and allow to save memory on profile
        instantiation.

        The maps from full sequence to reduced sequence and back are also stored to allow
        compressing and expanding the sequences.

        Notes
        -----

          full_to_reduced_sequence_map : (array)
             Map to reduce a sequence

          reduced_to_full_sequence_map : (dict)
             Map to restore sequence from reduced alignment

          multiplicity : (array)
            Numpy array, which stores the pattern multiplicity for each position of the reduced alignment.

          reduced_alignment : (2D numpy array)
            The reduced alignment. Shape is (N x L'), where N is number of
            sequences, L' - number of unique alignment patterns

          cseq : (array)
            The compressed sequence (corresponding row of the reduced alignment) attached to each node

        """

        if not self.reduce_alignment:
            self.multiplicity = np.ones(self.full_length, dtype=float)
            self.full_to_reduced_sequence_map = np.arange(self.full_length)
            self.reduced_to_full_sequence_map = {p:np.array([p]) for p in np.arange(self.full_length)}
            self._reduced_length==self._full_length
            self.reduced_alignment = self._aln
            return ttconf.SUCCESS

        self.logger("SeqData: making reduced alignment...", 1)
        # bind positions in full length sequence to that of the reduced (compressed) sequence
        self.full_to_reduced_sequence_map = np.zeros(self.full_length, dtype=int)
        # bind position in reduced sequence to the array of positions in full length sequence
        self.reduced_to_full_sequence_map = {}

        #if alignment is sparse, don't iterate over all invarible sites.
        #so pre-load alignment_patterns with the location of const sites!
        #and get the sites that we want to iterate over only!
        if self.is_sparse:
            from .vcf_utils import process_sparse_alignment
            tmp = process_sparse_alignment(self.aln, self.ref, self.ambiguous)
            reduced_aln_transpose = tmp["constant_columns"]
            alignment_patterns = tmp["constant_patterns"]
            variable_positions = tmp["variable_positions"]
            self.inferred_const_sites = tmp["constant_up_to_ambiguous"]
            self.nonref_positions = tmp["nonref_positions"]
        else: # transpose real alignment, for ease of iteration
            alignment_patterns = {}
            reduced_aln_transpose = []
            aln_transpose = np.array([self.aln[k] for k in self.sequence_names]).T
            variable_positions = np.arange(aln_transpose.shape[0])

        for pi in variable_positions:
            if self.is_sparse:
                pattern = np.array([self.aln[k][pi] if pi in self.aln[k] else self.ref[pi]
                           for k in self.sequence_names], dtype='S')
            else:
                pattern = np.copy(aln_transpose[pi])

            # if the column contains only one state and ambiguous nucleotides, replace
            # those with the state in other strains right away
            unique_letters = list(np.unique(pattern))
            if len(unique_letters)==2 and self.ambiguous in unique_letters:
                other = [c for c in unique_letters if c!=self.ambiguous][0]
                #also replace in original pattern!
                pattern[pattern == self.ambiguous] = other
                unique_letters = [other]
            # if there is a mutation in this column, give it its private pattern
            # this is required when sampling mutations from reconstructed profiles.
            # otherwise, all mutations corresponding to the same pattern will be coupled.
            # FIXME
            str_pattern = "".join(pattern.astype('U'))
            if len(unique_letters)>1:
                str_pattern += '_%d'%pi

            # if the pattern is not yet seen,
            if str_pattern not in alignment_patterns:
                # bind the index in the reduced aln, index in sequence to the pattern string
                alignment_patterns[str_pattern] = (len(reduced_aln_transpose), [pi])
                # append this pattern to the reduced alignment
                reduced_aln_transpose.append(pattern)
            else:
                # if the pattern is already seen, append the position in the real
                # sequence to the reduced aln<->sequence_pos_indexes map
                alignment_patterns[str_pattern][1].append(pi)

        # add constant alignment column not in the alignment. We don't know where they
        # are, so just add them to the end. First, determine sequence composition.
        if self.additional_constant_sites:
            character_counts = {c:np.sum(aln_transpose==c) for c in self.alphabet
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

                if str_pattern in alignment_patterns:
                    alignment_patterns[str_pattern][1].extend(pos_list)
                else:
                    alignment_patterns[str_pattern] = (len(reduced_aln_transpose), pos_list)
                    reduced_aln_transpose.append(np.array(list(str_pattern)))
                pi += n
                columns_left -= n


        # count how many times each column is repeated in the real alignment
        self.multiplicity = np.zeros(len(alignment_patterns))
        for p, pos in alignment_patterns.values():
            self.multiplicity[p]=len(pos)

        # create the reduced alignment as a dictionary linking names to sequences
        tmp_reduced_alignment = np.array(reduced_aln_transpose).T
        self.reduced_alignment = {k: tmp_reduced_alignment[i]
                                  for i,k in enumerate(self.sequence_names)}

        # create map to compress a sequence
        for p, pos in alignment_patterns.values():
            self.full_to_reduced_sequence_map[np.array(pos)]=p

        # create a map to reconstruct full sequence from the reduced (compressed) sequence
        for p, val in alignment_patterns.items():
            self.reduced_to_full_sequence_map[val[0]]=np.array(val[1], dtype=int)

        self.logger("SequenceData: constructed reduced alignment...", 1)
        self._reduced_length = len(self.multiplicity)
        return ttconf.SUCCESS


    def full_to_sparse_sequence(self, sequence):
        if self.ref is None:
            raise TypeError("SequenceData: sparse sequences can only be constructed when a reference sequence is defined")
        if type(sequence) is not np.ndarray:
            aseq = seq2array(sequence, fill_overhangs=False)
        else:
            aseq = sequence
        differences = np.where(self.ref!=aseq)[0]
        return {p:aseq[p] for p in differences}


    def reduced_to_full_sequence(self, sequence, include_additional_constant_sites=False, as_string=False):
        if include_additional_constant_sites:
            L = self.full_length
        else:
            L = self.full_length - self.additional_constant_sites

        tmp_seq = sequence[self.full_to_reduced_sequence_map[:L]]
        if as_string:
            return "".join(tmp_seq.astype('U'))
        else:
            return tmp_seq

