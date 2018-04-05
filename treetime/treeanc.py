from __future__ import print_function, division
import time
import config as ttconf
from Bio import Phylo
from Bio import AlignIO
import numpy as np
from gtr import GTR
import seq_utils
from version import tt_version as __version__
try:
    from itertools import izip
except ImportError:  #python3.x
    izip = zip


class TreeAnc(object):
    """
    Class defines simple tree object with basic interface methods: reading and
    saving from/to files, initializing leaves with sequences from the
    alignment, making ancestral state inferrence
    """

    def __init__(self, tree=None, aln=None, gtr=None, fill_overhangs=True,
                ref=None, verbose = ttconf.VERBOSE, ignore_gaps=True,
                convert_upper=True, seq_multiplicity=None, log=None, **kwargs):
        """
        TreeAnc constructor. It prepares tree, attach sequences to the leaf nodes,
        and sets some configuration parameters.

        Parameters
        ----------

         tree : str, Bio.Phylo.Tree
            Phylogenetic tree. String passed is interpreted as a filename to
            construct the Biopython tree.

         aln : str, Bio.Align.MultipleSequenceAlignment
            Sequence alignment. If a string passed, it is interpreted as the
            filename to read Biopython alignment from.

         gtr : str, GTR
            gtr model object. If string passed, it is interpreted as the type of
            the GTR model. A new GTR instance will be created for this type.
           **Note** some GTR types require additional configuration parameters.
           If the new GTR is being instantiated, these parameters are expected
           to be passed as kwargs. If nothing is passed, the default values are
           used, which might cause unexpected results.

         fill_overhangs : bool
            In some cases, the missing data on both ends of the alignment is
            filled with the gap sign('-'). As we suppose that the most
            appropriate way to deal with the missing data is to assign it to the
            "unknown" character ('N' for nucleotides, 'X' for aminoacids). If the
            parameter is set to True, the end-gaps are converted to unknown
            symbols. Otherwise, the alignment is treated as-is

         ignore_gaps: bool
            ignore gaps in branch length calculations

         verbose : int
            verbosity level as number from 0 (lowest) to 10 (highest).

         seq_multiplicity: dict
            if individual nodes in the tree correspond to multiple sampled sequences
            (i.e. read count in a deep sequencing experiment), these can be
            specified as a dictionary

        Keyword Args
        ------------

          Keyword arguments to construct GTR model

        """
        if tree is None:
            raise TypeError("TreeAnc requires a tree!")
        self.__version__ = __version__
        self.t_start = time.time()
        self.verbose = verbose
        self.log=log
        self.logger("TreeAnc: set-up",1)
        self._internal_node_count = 0
        self.use_mutation_length=False
        self.one_mutation = None
        self.fill_overhangs = fill_overhangs
        self.is_vcf = False  #this is set true when aln is set, if aln is dict
        self.seq_multiplicity = {} if seq_multiplicity is None else seq_multiplicity

        self.ignore_gaps = ignore_gaps
        self.set_gtr(gtr if gtr is not None else 'JC69', **kwargs)

        self.tree = tree
        if tree is None:
            self.logger("TreeAnc: tree loading failed! exiting",0)
            return

        # will be None if not set
        self.ref = ref

        # force all sequences to be upper case letters
        # (desired for nuc or aa, not for other discrete states)
        self.convert_upper = convert_upper

        # set alignment and attach sequences to tree on success.
        # otherwise self.aln will be None
        self.aln = aln


    def logger(self, msg, level, warn=False):
        """
        Print log message *msg* to stdout.

        Parameters
        -----------

         msg : str
            string to print on the screen

         level : int
            log-level. Only the messages with the level higher than the
            current verbose level will be shown.

         warn : bool
            warning flag. If True, the message will be displayed
            regardless of its log-level.

        """
        if level<self.verbose or (warn and level<=self.verbose):
            dt = time.time() - self.t_start
            outstr = '\n' if level<2 else ''
            outstr+=format(dt, '4.2f')+'\t'
            outstr+= level*'-'
            outstr+=msg
            try:
                log.write(outstr+'\n')
                log.flush()
            except:
                print(outstr)


####################################################################
## SET-UP
####################################################################
    @property
    def leaves_lookup(self):
        """
        Leaves lookup is the {leaf-name:leaf-node} dictionary. It enables fast
        search of a tree leaf object by its name.
        """
        return self._leaves_lookup

    @property
    def gtr(self):
        """
        Get GTR object currently used.
        """
        return self._gtr

    @gtr.setter
    def gtr(self, value):
        """
        Set a new GTR object

        Parameters
        -----------

         value :GTR
            the new GTR object
        """
        if not isinstance(value, GTR):
            raise TypeError(" GTR instance expected")
        self._gtr = value


    def set_gtr(self, in_gtr, **kwargs):
        """
        Create new GTR model if needed, and set the model as an attribute of the
        TreeAnc class

        Parameters
        -----------

         in_gtr : str, GTR
            The gtr model to be assigned. If string is passed,
            it is understood as the name of the standard GTR model, and is
            attempted to be created through GTR.standard() interface. In case
            GTR instance is passed, it is directly set as the class attribute

        Keyword Args
    ------------

         All parameters needed for the gtr creation. If none passed, defaults are assumed.
           Refer to the particular GTR models for the exact parameter values

        """
        if type(in_gtr)==str:
            self._gtr = GTR.standard(model=in_gtr, **kwargs)
            self._gtr.logger = self.logger

        elif isinstance(in_gtr, GTR):
            self._gtr = in_gtr
            self._gtr.logger=self.logger
        else:
            self.logger("TreeAnc.gtr_setter: can't interpret GTR model", 1, warn=True)
            raise TypeError("Cannot set GTR model in TreeAnc class: GTR or "
                "string expected")

        if self._gtr.ambiguous is None:
            self.fill_overhangs=False


    @property
    def tree(self):
        """
        Get reference to the phylogenetic tree currently used by the TreeAnc.
        """
        return self._tree

    @tree.setter
    def tree(self, in_tree):
        '''
        assigns a tree to the internal self._tree variable. The tree is either
        loaded from file (if in_tree is str) or assigned (if in_tree is a Phylo.tree)
        '''
        from os.path import isfile
        if isinstance(in_tree, Phylo.BaseTree.Tree):
            self._tree = in_tree
        elif type(in_tree) in [str, unicode] and isfile(in_tree):
            try:
                self._tree=Phylo.read(in_tree, 'newick')
            except:
                fmt = in_tree.split('.')[-1]
                if fmt in ['nexus', 'nex']:
                    self._tree=Phylo.read(in_tree, 'nexus')
                else:
                    self.logger('TreeAnc: could not load tree, format needs to be nexus or newick! input was '+str(in_tree),1)
                    self._tree = None
                    return
        else:
            self.logger('TreeAnc: could not load tree! input was '+str(in_tree),1)
            self._tree = None
            return

        # remove all existing sequence attributes
        for node in self._tree.find_clades():
            if hasattr(node, "sequence"):
                node.__delattr__("sequence")
            node.original_length = node.branch_length
            node.mutation_length = node.branch_length
        self.prepare_tree()


    @property
    def aln(self):
        """
        Get the multiple sequence alignment currently used by the TreeAnc
        """
        return self._aln

    @aln.setter
    def aln(self,in_aln):
        # load alignment from file if necessary
        from os.path import isfile
        from Bio.Align import MultipleSeqAlignment
        self._aln = None
        if isinstance(in_aln, MultipleSeqAlignment):
            self._aln = in_aln
        elif type(in_aln) in [str, unicode] and isfile(in_aln):
            for fmt in ['fasta', 'phylip-relaxed', 'nexus']:
                try:
                    self._aln=AlignIO.read(in_aln, 'fasta')
                    break
                except:
                    continue
        elif type(in_aln) is dict:  #if is read in from VCF file
            self._aln = in_aln
            self.is_vcf = True

        if self._aln is None:
            self.logger("TreeAnc: loading alignment failed... ",1, warn=True)
            return

        #Convert to uppercase here, rather than in _attach_sequences_to_nodes
        #(which used to do it through seq2array in seq_utils.py)
        #so that it is controlled by param convert_upper. This way for
        #mugration (ancestral reconstruction of non-sequences), you can
        #use upper- and lower case characters for discrete states!
        if (not self.is_vcf) and self.convert_upper:
            self._aln = MultipleSeqAlignment([seq.upper() for seq in self._aln])

        if hasattr(self, '_tree'):
            self._attach_sequences_to_nodes()
        else:
            self.logger("TreeAnc.aln: sequences not yet attached to tree", 3, warn=True)

    @property
    def ref(self):
        """
        Get the str reference nucleotide sequence currently used by TreeAnc
        When having read in from a VCF, this is what variants map to
        """
        return self._ref


    @ref.setter
    def ref(self, in_ref):
        self._ref = in_ref


    def _attach_sequences_to_nodes(self):
        #print ("inside attach seq to nodes")
        '''
        For each node of the tree, check whether there is a sequence available
        in the alignment and assign this sequence as a character array
        '''
        if type(self.aln) is dict:
            self.seq_len = len(self.ref)
        else:
            self.seq_len = self.aln.get_alignment_length()
        self.one_mutation = 1.0/self.seq_len

        failed_leaves= 0
        if type(self.aln) is dict:
            # if alignment is specified as difference from ref
            dic_aln = self.aln
        else:
            # if full alignment is specified
            dic_aln = {k.name: seq_utils.seq2array(k.seq, fill_overhangs=self.fill_overhangs,
                                                   ambiguous_character=self.gtr.ambiguous)
                                for k in self.aln} #

        # loop over tree,
        for l in self.tree.find_clades():
            if l.name in dic_aln:
                l.sequence= dic_aln[l.name]
                if l.name in self.seq_multiplicity:
                    l.count = self.seq_multiplicity[l.name]
                else:
                    l.count = 1.0
            elif l.is_terminal():
                self.logger("***WARNING: TreeAnc._attach_sequences_to_nodes: NO SEQUENCE FOR LEAF: %s" % l.name, 0, warn=True)
                failed_leaves += 1
                l.sequence = seq_utils.seq2array(self.gtr.ambiguous*self.seq_len, fill_overhangs=self.fill_overhangs,
                                                 ambiguous_character=self.gtr.ambiguous)
                if failed_leaves > self.tree.count_terminals() / 3:
                    self.logger("ERROR: At least 30\\% terminal nodes cannot be assigned with a sequence!\n", 0, warn=True)
                    self.logger("Are you sure the alignment belongs to the tree?", 2, warn=True)
                    break
            else: # could not assign sequence for internal node - is OK
                pass

        if failed_leaves:
            self.logger("***WARNING: TreeAnc: %d nodes don't have a matching sequence in the alignment. POSSIBLE ERROR."%failed_leaves, 0, warn=True)

        self.make_reduced_alignment()


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

        The following attributes are assigned by the method:

         - full_to_reduced_sequence_map: map to reduce a sequence
         - reduced_to_full_sequence_map: map to restore sequence from reduced alignment
         - multiplicity: numpy array, which stores the pattern multiplicity for
           each position of the reduced alignment.
         - reduced_alignment: 2D numpy array, representing the alignment. Shape is
           (N x L'), where N is number of sequences, L' - number of unique alignment patterns


        In addition, each node gets

          - cseq: compressed sequence (corresponding row of the reduced alignment)

        """

        self.logger("TreeAnc: making reduced alignment...", 1)

        from collections import defaultdict

        # bind positions in real sequence to that of the reduced (compressed) sequence
        self.full_to_reduced_sequence_map = np.zeros(self.seq_len, dtype=int)

        # bind position in reduced sequence to the array of positions in real (expanded) sequence
        self.reduced_to_full_sequence_map = {}

        #if is a dict, want to be efficient and not iterate over a bunch of const_sites
        #so pre-load alignment_patterns with the location of const sites!
        #and get the sites that we want to iterate over only!
        if type(self.aln) is dict:
            tmp_reduced_aln, alignment_patterns, positions = self.process_alignment_dict()
            seqNames = self.aln.keys() #store seqName order to put back on tree
        else:
            # transpose real alignment, for ease of iteration
            alignment_patterns = {}
            tmp_reduced_aln = []
            # NOTE the order of tree traversal must be the same as below
            # for assigning the cseq attributes to the nodes.
            seqs = [n.sequence for n in self.tree.find_clades() if hasattr(n, 'sequence')]
            if len(np.unique([len(x) for x in seqs]))>1:
                self.logger("TreeAnc: Sequences differ in in length! ABORTING",0, warn=True)
                aln_transpose = None
                return
            else:
                aln_transpose = np.array(seqs).T
                positions = range(self.seq_len)

        for pi in positions:
            if type(self.aln) is dict:
                pattern = [ self.aln[k][pi] if pi in self.aln[k].keys()
                            else self.ref[pi] for k,v in self.aln.iteritems() ]
            else:
                pattern = aln_transpose[pi]

            str_pat = "".join(pattern)
            # if the column contains only one state and ambiguous nucleotides, replace
            # those with the state in other strains right away
            unique_letters = list(np.unique(pattern))
            if hasattr(self.gtr, "ambiguous"):
                if len(unique_letters)==2 and self.gtr.ambiguous in unique_letters:
                    other = [c for c in unique_letters if c!=self.gtr.ambiguous][0]
                    str_pat = str_pat.replace(self.gtr.ambiguous, other)
                    unique_letters = [other]
            # if there is a mutation in this column, give it its private pattern
            # this is required when sampling mutations from reconstructed profiles.
            # otherwise, all mutations corresponding to the same pattern will be coupled.
            if len(unique_letters)>1:
                str_pat += '_%d'%pi

            # if the pattern is not yet seen,
            if str_pat not in alignment_patterns:
                # bind the index in the reduced aln, index in sequence to the pattern string
                alignment_patterns[str_pat] = (len(tmp_reduced_aln), [pi])
                # append this pattern to the reduced alignment
                tmp_reduced_aln.append(pattern)
            else:
                # if the pattern is already seen, append the position in the real
                # sequence to the reduced aln<->sequence_pos_indexes map
                alignment_patterns[str_pat][1].append(pi)

        # count how many times each column is repeated in the real alignment
        self.multiplicity = np.zeros(len(alignment_patterns))
        for p, pos in alignment_patterns.values():
            self.multiplicity[p]=len(pos)

        # create the reduced alignment as np array
        self.reduced_alignment = np.array(tmp_reduced_aln).T

        # create map to compress a sequence
        for p, pos in alignment_patterns.values():
            self.full_to_reduced_sequence_map[np.array(pos)]=p

        # create a map to reconstruct full sequence from the reduced (compressed) sequence
        for p, val in alignment_patterns.iteritems():
            self.reduced_to_full_sequence_map[val[0]]=np.array(val[1], dtype=int)

        # assign compressed sequences to all nodes of the tree, which have sequence assigned
        # for dict we cannot assume this is in the same order, as it does below!
        # so do it explicitly
        if type(self.aln) is dict:
            seq_reduce_align = {n:self.reduced_alignment[i]
                                for i, n in enumerate(seqNames)}
            for n in self.tree.find_clades():
                if hasattr(n, 'sequence'):
                    n.cseq = seq_reduce_align[n.name]
        else:
            # NOTE the order of tree traversal must be the same as above to catch the
            # index in the reduced alignment correctly
            seq_count = 0
            for n in self.tree.find_clades():
                if hasattr(n, 'sequence'):
                    n.cseq = self.reduced_alignment[seq_count]
                    seq_count+=1

        # sequences are overwritten during reconstruction and
        # ambiguous sites change. Keep orgininals for reference
        self.original_sequences = {n.name:n.cseq for n in self.tree.get_terminals()}

        self.logger("TreeAnc: finished reduced alignment...", 1)


    def process_alignment_dict(self):
        """
        prepare the dictionary specifying differences from a reference sequence
        to construct the reduced alignment with variable sites only. NOTE:
            - sites can be constant but different from the reference
            - sites can be constant plus a ambiguous sites

        assigns:
            - self.nonref_positions: at least one sequence is different from ref
        returns:
            - reduced_alignment_const: reduced alignment accounting for
                                       non-variable postitions
            - alignment_patterns_const:
                dict pattern -> (pos in reduced alignment, list of pos in full alignment)
            - variable_positions: list of variable positions needed to construct remaining
        """

        # number of sequences in alignment
        nseq = len(self.aln)

        from collections import defaultdict
        inv_map = defaultdict(list)
        for k,v in self.aln.iteritems():
            for pos, bs in v.iteritems():
                inv_map[pos].append(bs)

        self.nonref_positions = np.sort(inv_map.keys())
        self.inferred_const_sites = []

        ambiguous_char = self.gtr.ambiguous
        nonref_const = []
        nonref_alleles = []
        ambiguous_const = []
        variable_pos = []
        for pos, bs in inv_map.iteritems(): #loop over positions and patterns
            bases = "".join(np.unique(bs))
            if len(bs) == nseq:
                if (len(bases)<=2 and ambiguous_char in bases) or len(bases)==1:
                    # all sequences different from reference, but only one state
                    # (other than ambiguous_char) in column
                    nonref_const.append(pos)
                    nonref_alleles.append(bases.replace(ambiguous_char, ''))
                    if ambiguous_char in bases: #keep track of sites 'made constant'
                        self.inferred_const_sites.append(pos)
                else:
                    # at least two non-reference alleles
                    variable_pos.append(pos)
            else:
                # not every sequence different from reference
                if bases==ambiguous_char:
                    ambiguous_const.append(pos)
                    self.inferred_const_sites.append(pos) #keep track of sites 'made constant'
                else:
                    # at least one non ambiguous non-reference allele not in
                    # every sequence
                    variable_pos.append(pos)

        refMod = np.fromstring(self.ref, 'S1')
        # place constant non reference positions by their respective allele
        refMod[nonref_const] = nonref_alleles
        # mask variable positions
        states = self.gtr.alphabet
        # maybe states = np.unique(refMod)
        refMod[variable_pos] = '.'

        # for each base in the gtr, make constant alignment pattern and
        # assign it to all const positions in the modified reference sequence
        reduced_alignment_const = []
        alignment_patterns_const = {}
        for base in states:
            p = base*nseq
            pos = list(np.where(refMod==base)[0])
            #if the alignment doesn't have a const site of this base, don't add! (ex: no '----' site!)
            if len(pos):
                alignment_patterns_const[p] = [len(reduced_alignment_const), pos]
                reduced_alignment_const.append(list(p))


        return reduced_alignment_const, alignment_patterns_const, variable_pos


    def prepare_tree(self):
        """
        Set link to parent and net distance to root for all tree nodes.
        Should be run once the tree is read and after every tree topology or branch
        length optimizations.
        """
        if self.one_mutation is None:
            self.tree.root.branch_length = 0.001
        else:
            self.tree.root.branch_length = self.one_mutation
        self.tree.root.mutation_length = self.tree.root.branch_length
        self.tree.root.mutations = []
        self.tree.ladderize()
        self._prepare_nodes()
        self._leaves_lookup = {node.name:node for node in self.tree.get_terminals()}


    def _prepare_nodes(self):
        """
        Set auxilliary parameters to every node of the tree.
        """
        self.tree.root.up = None
        self.tree.root.bad_branch=self.tree.root.bad_branch if hasattr(self.tree.root, 'bad_branch') else False
        internal_node_count = 0
        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            internal_node_count+=1
            if clade.name is None:
                clade.name = "NODE_" + format(self._internal_node_count, '07d')
                self._internal_node_count += 1
            for c in clade.clades:
                c.bad_branch=c.bad_branch if hasattr(c, 'bad_branch') else False
                c.up = clade
        self._calc_dist2root()
        self._internal_node_count = max(internal_node_count, self._internal_node_count)

    def _calc_dist2root(self):
        """
        For each node in the tree, set its root-to-node distance as dist2root
        attribute
        """
        self.tree.root.dist2root = 0.0
        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            for c in clade.clades:
                if not hasattr(c, 'mutation_length'):
                    c.mutation_length=c.branch_length
                c.dist2root = c.up.dist2root + c.mutation_length



####################################################################
## END SET-UP
####################################################################

    def infer_gtr(self, print_raw=False, marginal=False, normalized_rate=True,
                  fixed_pi=None, pc=5.0, **kwargs):
        """
        Calculates GTR model given the multiple sequence alignment and the tree.
        It performs ancestral sequence inferrence (joint or marginal) followed by
        the branch lengths optimization. Then, the numbers of mutations are counted
        in the optimal tree and related to the time within the mutation happened.
        From this statistics, the relative state transition probabilities are inferred,
        and the transition matrix is computed.
        The result is used to construct the new GTR model of type 'custom'.
        The model is assigned to the TreeAnc and is used in the following analysis.

        Parameters
        -----------

         print_raw : bool
            Should print the inferred GTR model?

         marginal : bool
            Should use marginal sequence reconstruction?

         normalized_rate : bool
            If True, will set the mutation rate prefactor to 1.0.

         fixed_pi : np.array, None
            Provide the equilibrium character concentrations.
            If None is passed, the concentrations will be inferred from scratch.

         pc: float, 5.0
            Number of pseudo counts to use in gtr inference

        Returns
        -------

         gtr : GTR
            The inferred GTR model.
        """

        # decide which type of the Maximum-likelihood reconstruction use
        # (marginal) or (joint)
        if marginal:
            _ml_anc = self._ml_anc_marginal
        else:
            _ml_anc = self._ml_anc_joint

        self.logger("TreeAnc inferring the GTR model from the tree...", 1)
        _ml_anc(final=True, **kwargs) # call one of the reconstruction types
        alpha = list(self.gtr.alphabet)
        n=len(alpha)
        nij = np.zeros((n,n))
        Ti = np.zeros(n)

        self.logger("TreeAnc.infer_gtr: counting mutations...", 2)
        for node in self.tree.find_clades():
            if hasattr(node,'mutations'):
                for a,pos, d in node.mutations:
                    i,j = alpha.index(a), alpha.index(d)
                    nij[i,j]+=1
                    Ti[i] += 0.5*self._branch_length_to_gtr(node)
                    Ti[j] -= 0.5*self._branch_length_to_gtr(node)
                for ni,nuc in enumerate(node.cseq):
                    i = alpha.index(nuc)
                    Ti[i] += self._branch_length_to_gtr(node)*self.multiplicity[ni]
        self.logger("TreeAnc.infer_gtr: counting mutations...done", 3)
        if print_raw:
            print('alphabet:',alpha)
            print('n_ij:', nij)
            print('T_i:', Ti)
        root_state = np.array([np.sum((self.tree.root.cseq==nuc)*self.multiplicity) for nuc in alpha])

        self._gtr = GTR.infer(nij, Ti, root_state, fixed_pi=fixed_pi, pc=pc,
                              alphabet=self.gtr.alphabet, logger=self.logger,
                              prof_map = self.gtr.profile_map)
        if normalized_rate:
            self.logger("TreeAnc.infer_gtr: setting overall rate to 1.0...", 2)
            self._gtr.mu=1.0
        return self._gtr


###################################################################
### ancestral reconstruction
###################################################################
    def infer_ancestral_sequences(self,*args, **kwargs):
        """Shortcut for :meth:`reconstruct_anc`

        Reconstruct ancestral states

        Parameters
        -----------

         method : str
            Method to use. Supported values are "fitch" and "ml"

        Returns
        -------

         N_diff : int
            Number of nucleotides different from the previous
            reconstruction.  If there were no pre-set sequences, returns N*L

        """
        self.reconstruct_anc(*args,**kwargs)


    def reconstruct_anc(self, method='ml', infer_gtr=False, marginal=False, **kwargs):
        """

        Reconstruct ancestral states

        Parameters
        -----------

         method : str
            Method to use. Supported values are "fitch" and "ml"

        Returns
        -------

         N_diff : int
            Number of nucleotides different from the previous
            reconstruction.  If there were no pre-set sequences, returns N*L
        """
        self.logger("TreeAnc.infer_ancestral_sequences: method: " + method, 1)

        if method == 'ml':
            if marginal:
                _ml_anc = self._ml_anc_marginal
            else:
                _ml_anc = self._ml_anc_joint
        else:
            _ml_anc = self._fitch_anc

        if infer_gtr:
            self.infer_gtr(marginal=marginal, **kwargs)
            N_diff = _ml_anc(**kwargs)
        else:
            N_diff = _ml_anc(**kwargs)

        return N_diff


    def recover_var_ambigs(self):
        """
        Recalculates mutations using the original compressed sequence for terminal nodes
        which will recover ambiguous bases at variable sites. (See 'get_mutations')

        Once this has been run, infer_gtr and other functions which depend on self.gtr.alphabet
        will not work, as ambiguous bases are not part of that alphabet (only A, C, G, T, -).
        This is why it's left for the user to choose when to run
        """
        for node in self.tree.get_terminals():
            node.mutations = self.get_mutations(node, keep_var_ambigs=True)


    def get_mutations(self, node, keep_var_ambigs=False):
        """
        Get the mutations on a tree branch. Take compressed sequences from both sides
        of the branch (attached to the node), compute mutations between them, and
        expand these mutations to the positions in the real sequences.

        Parameters
        ----------

         node : PhyloTree.Clade
            Tree node, which is the child node attached to the branch.

         keep_var_ambigs : boolean
            If true, generates mutations based on the *original* _compressed_ sequence, which
            may include ambiguities. Note sites that only have 1 unambiguous base and ambiguous
            bases ("AAAAANN") are stripped of ambiguous bases *before* compression, so ambiguous
            bases will *not* be preserved.

        Returns
        -------

          muts : list
            List of mutations. Each mutation is represented as tuple of
            (parent_state, position, child_state).
        """

        # if ambiguous site are to be restored and node is terminal,
        # assign original sequence, else reconstructed cseq
        node_seq = node.cseq
        if keep_var_ambigs and (node.name in self.original_sequences) and node.is_terminal():
            node_seq = self.original_sequences[node.name]

        muts = []
        for p, (anc, der) in enumerate(izip(node.up.cseq, node_seq)):
            # only if the states in compressed sequences differ:
            if anc!=der:
                # expand to the positions in real sequence
                muts.extend([(anc, pos, der) for pos in self.reduced_to_full_sequence_map[p]])

        #sort by position
        return sorted(muts, key=lambda x:x[1])


    def expanded_sequence(self, node):
        """
        Get node's compressed sequence and expand it to the real sequence

        Parameters
        ----------

         node  : PhyloTree.Clade
            Tree node

        Returns
        -------

         seq : np.array
            Sequence as np.array of chars
        """
        seq = np.zeros_like(self.full_to_reduced_sequence_map, dtype='S1')
        for pos, state in enumerate(node.cseq):
            seq[self.reduced_to_full_sequence_map[pos]] = state

        return seq

    def dict_sequence(self, node, keep_var_ambigs=False):
        """
        For VCF-based TreeAnc objects, we do not want to store the entire
        sequence on every node - not space efficient! Instead, return the dict
        of mutation locations for this sequence. This is used in place of
        'expanded_sequence' for VCF-based obj throughout TreeAnc. However, users
        can still call 'expanded_sequence' if they do actually want the whole thing!

        Parameters
        ----------
         node  : PhyloTree.Clade
            Tree node

        Returns
        -------
         seq : dict
            dict where keys are position and value is the mutation

        EBH 6 Dec 2017
        """
        seq = {}

        node_seq = node.cseq
        if keep_var_ambigs and (node.name in self.original_sequences) and node.is_terminal():
            node_seq = self.original_sequences[node.name]

        for pos in self.nonref_positions:
            cseqLoc = self.full_to_reduced_sequence_map[pos]
            base = node_seq[cseqLoc]
            if self.ref[pos] != base:
                seq[pos] = base

        return seq

###################################################################
### FITCH
###################################################################
    def _fitch_anc(self, **kwargs):
        """
        Reconstruct ancestral states using Fitch's algorithm. The method requires
        sequences to be assigned to leaves. It implements the iteration from
        leaves to the root constructing the Fitch profiles for each character of
        the sequence, and then by propagating from the root to the leaves,
        reconstructs the sequences of the internal nodes.

        Keyword Args
        ------------


        Returns
        -------

         Ndiff : int
            Number of the characters that changed since the previous
            reconstruction. These changes are determined from the pre-set
            sequence attributes of the nodes. If there are no sequences available
            (i.e., no reconstruction has been made before), returns the total
            number of characters in the tree.

        """
        # set fitch profiiles to each terminal node

        for l in self.tree.get_terminals():
            l.state = [[k] for k in l.cseq]

        L = len(self.tree.get_terminals()[0].cseq)

        self.logger("TreeAnc._fitch_anc: Walking up the tree, creating the Fitch profiles",2)
        for node in self.tree.get_nonterminals(order='postorder'):
            node.state = [self._fitch_state(node, k) for k in range(L)]

        ambs = [i for i in range(L) if len(self.tree.root.state[i])>1]
        if len(ambs) > 0:
            for amb in ambs:
                self.logger("Ambiguous state of the root sequence "
                                    "in the position %d: %s, "
                                    "choosing %s" % (amb, str(self.tree.root.state[amb]),
                                                     self.tree.root.state[amb][0]), 4)
        self.tree.root.cseq = np.array([k[np.random.randint(len(k)) if len(k)>1 else 0]
                                           for k in self.tree.root.state])

        if self.is_vcf:
            self.tree.root.sequence = self.dict_sequence(self.tree.root)
        else:
            self.tree.root.sequence = self.expanded_sequence(self.tree.root)


        self.logger("TreeAnc._fitch_anc: Walking down the self.tree, generating sequences from the "
                         "Fitch profiles.", 2)
        N_diff = 0
        for node in self.tree.get_nonterminals(order='preorder'):
            if node.up != None: # not root
                sequence =  np.array([node.up.cseq[i]
                        if node.up.cseq[i] in node.state[i]
                        else node.state[i][0] for i in range(L)])
                if hasattr(node, 'sequence'):
                    N_diff += (sequence!=node.cseq).sum()
                else:
                    N_diff += L
                node.cseq = sequence
                if self.is_vcf:
                    node.sequence = self.dict_sequence(node)
                else:
                    node.sequence = self.expanded_sequence(node)
                node.mutations = self.get_mutations(node)

            node.profile = seq_utils.seq2prof(node.cseq, self.gtr.profile_map)
            del node.state # no need to store Fitch states
        self.logger("Done ancestral state reconstruction",3)
        for node in self.tree.get_terminals():
            node.profile = seq_utils.seq2prof(node.cseq, self.gtr.profile_map)
        return N_diff

    def _fitch_state(self, node, pos):
        """
        Determine the Fitch profile for a single character of the node's sequence.
        The profile is essentially the intersection between the children's
        profiles or, if the former is empty, the union of the profiles.

        Parameters
        ----------

         node : PhyloTree.Clade:
            Internal node which the profiles are to be determined

         pos : int
            Position in the node's sequence which the profiles should
            be determinedf for.

        Returns
        -------
         state : numpy.array
            Fitch profile for the character at position pos of the given node.
        """
        state = self._fitch_intersect([k.state[pos] for k in node.clades])
        if len(state) == 0:
            state = np.concatenate([k.state[pos] for k in node.clades])
        return state

    def _fitch_intersect(self, arrays):
        """
        Find the intersection of any number of 1D arrays.
        Return the sorted, unique values that are in all of the input arrays.
        Adapted from numpy.lib.arraysetops.intersect1d
        """
        def pairwise_intersect(arr1, arr2):
            s2 = set(arr2)
            b3 = [val for val in arr1 if val in s2]
            return b3

        arrays = list(arrays) # allow assignment
        N = len(arrays)
        while N > 1:
            arr1 = arrays.pop()
            arr2 = arrays.pop()
            arr = pairwise_intersect(arr1, arr2)
            arrays.append(arr)
            N = len(arrays)

        return arrays[0]



###################################################################
### Maximum Likelihood
###################################################################

    def ancestral_likelihood(self):
        """
        Calculate the likelihood of the given realization of the sequences in
        the tree

        Returns
        -------

         log_lh : float
            The tree likelihood given the sequences
        """
        log_lh = np.zeros(self.tree.root.cseq.shape[0])
        for node in self.tree.find_clades(order='postorder'):

            if node.up is None: #  root node
                # 0-1 profile
                profile = seq_utils.seq2prof(node.cseq, self.gtr.profile_map)
                # get the probabilities to observe each nucleotide
                profile *= self.gtr.Pi
                profile = profile.sum(axis=1)
                log_lh += np.log(profile) # product over all characters
                continue

            t = node.branch_length

            indices = np.array([(np.argmax(self.gtr.alphabet==a),
                        np.argmax(self.gtr.alphabet==b)) for a, b in izip(node.up.cseq, node.cseq)])

            logQt = np.log(self.gtr.expQt(t))
            lh = logQt[indices[:, 1], indices[:, 0]]
            log_lh += lh

        return log_lh

    def _branch_length_to_gtr(self, node):
        """
        Set branch lengths to either mutation lengths of given branch lengths.
        The assigend values are to be used in the following ML analysis.
        """
        if self.use_mutation_length:
            return max(ttconf.MIN_BRANCH_LENGTH*self.one_mutation, node.mutation_length)
        else:
            return max(ttconf.MIN_BRANCH_LENGTH*self.one_mutation, node.branch_length)


    def _ml_anc_marginal(self, verbose=0, store_compressed=True, final=True,
                                            sample_from_profile=False,
                                            debug=False, **kwargs):
        """
        Perform marginal ML reconstruction of the ancestral states. In contrast to
        joint reconstructions, this needs to access the probabilities rather than only
        log probabilities and is hence handled by a separate function.

        Keyword Args
        ------------

         store_lh : bool
            If True, all likelihoods will be stored for all nodes. Useful for
            testing, diagnostics and if special post-processing is required.

         verbose :int
            How verbose the output should be
        """

        tree = self.tree
        # number of nucleotides changed from prev reconstruction
        N_diff = 0

        L = self.tree.get_terminals()[0].cseq.shape[0]
        n_states = self.gtr.alphabet.shape[0]
        self.logger("TreeAnc._ml_anc_marginal: type of reconstruction: Marginal", 2)

        self.logger("Walking up the tree, computing likelihoods... ", 3)
        #  set the leaves profiles
        for leaf in tree.get_terminals():
            # in any case, set the profile
            leaf.marginal_subtree_LH = seq_utils.seq2prof(leaf.cseq, self.gtr.profile_map)
            leaf.marginal_subtree_LH_prefactor = np.zeros(L)

        # propagate leaves -->> root, set the marginal-likelihood messages
        for node in tree.get_nonterminals(order='postorder'): #leaves -> root
            # regardless of what was before, set the profile to ones
            node.marginal_subtree_LH_prefactor = np.zeros(L)
            node.marginal_subtree_LH = np.ones((L, n_states)) # we will multiply it
            for ch in node.clades:
                ch.marginal_Lx = self.gtr.propagate_profile(ch.marginal_subtree_LH,
                    self._branch_length_to_gtr(ch), return_log=False) # raw prob to transfer prob up
                node.marginal_subtree_LH *= ch.marginal_Lx
                node.marginal_subtree_LH_prefactor += ch.marginal_subtree_LH_prefactor

            pre = node.marginal_subtree_LH.sum(axis=1) #sum over nucleotide states
            node.marginal_subtree_LH = (node.marginal_subtree_LH.T/pre).T # normalize so that the sum is 1
            node.marginal_subtree_LH_prefactor += np.log(pre) # and store log-prefactor

        self.logger("Computing root node sequence and total tree likelihood...",3)
        # reconstruct the root node sequence
        tree.root.marginal_subtree_LH *= self.gtr.Pi # Msg to the root from the distant part (equ frequencies)
        pre=tree.root.marginal_subtree_LH.sum(axis=1)
        tree.root.marginal_profile = (tree.root.marginal_subtree_LH.T/pre).T
        tree.root.marginal_subtree_LH_prefactor += np.log(pre)

        # choose sequence characters from this profile.
        # treat root node differently to avoid piling up mutations on the longer branch
        if sample_from_profile=='root':
            root_sample_from_profile = True
            other_sample_from_profile = False
        elif isinstance(sample_from_profile, bool):
            root_sample_from_profile = sample_from_profile
            other_sample_from_profile = sample_from_profile

        seq, prof_vals, idxs = seq_utils.prof2seq(tree.root.marginal_profile,
                                                  self.gtr, sample_from_prof=root_sample_from_profile)

        self.tree.sequence_LH = np.log(prof_vals) + tree.root.marginal_subtree_LH_prefactor
        self.tree.sequence_marginal_LH = (self.tree.sequence_LH*self.multiplicity).sum()
        self.tree.root.cseq = seq
        if final:
            if self.is_vcf:
                self.tree.root.sequence = self.dict_sequence(self.tree.root)
            else:
                self.tree.root.sequence = self.expanded_sequence(self.tree.root)

        # need this fake msg to account for the complementary subtree when traversing tree back
        tree.root.seq_msg_from_parent = np.repeat([self.gtr.Pi], len(tree.root.cseq), axis=0)

        self.logger("Walking down the tree, computing maximum likelihood sequences...",3)
        # propagate root -->> leaves, reconstruct the internal node sequences
        # provided the upstream message + the message from the complementary subtree
        for node in tree.find_clades(order='preorder'):
            if node.up is None: # skip if node is root
                continue

            # integrate the information coming from parents with the information
            # of all children my multiplying it to the prev computed profile
            tmp_msg = np.copy(node.up.seq_msg_from_parent)
            for c in node.up.clades:
                if c != node:
                    tmp_msg*=c.marginal_Lx
            norm_vector = tmp_msg.sum(axis=1)
            tmp_msg=(tmp_msg.T/norm_vector).T
            node.seq_msg_from_parent = self.gtr.propagate_profile(tmp_msg,
                                            self._branch_length_to_gtr(node), return_log=False)
            node.marginal_profile = node.marginal_subtree_LH * node.seq_msg_from_parent

            norm_vector = node.marginal_profile.sum(axis=1)
            node.marginal_profile=(node.marginal_profile.T/norm_vector).T
            # choose sequence based maximal marginal LH.
            seq, prof_vals, idxs = seq_utils.prof2seq(node.marginal_profile, self.gtr,
                                                      sample_from_prof=other_sample_from_profile)

            if hasattr(node, 'cseq') and node.cseq is not None:
                N_diff += (seq!=node.cseq).sum()
            else:
                N_diff += L

            #assign new sequence
            node.cseq = seq
            if final:
                if self.is_vcf:
                    node.sequence = self.dict_sequence(node)
                else:
                    node.sequence = self.expanded_sequence(node)
                node.mutations = self.get_mutations(node)


        # note that the root doesn't contribute to N_diff (intended, since root sequence is often ambiguous)
        self.logger("TreeAnc._ml_anc_marginal: ...done", 3)
        if store_compressed:
            self._store_compressed_sequence_pairs()

        # do clean-up:
        if not debug:
            for node in self.tree.find_clades():
                del node.marginal_subtree_LH
                del node.marginal_subtree_LH_prefactor
                del node.seq_msg_from_parent

        return N_diff


    def _ml_anc_joint(self, verbose=0, store_compressed=True, final=True,
                                        sample_from_profile=False,
                                        debug=False, **kwargs):

        """
        Perform joint ML reconstruction of the ancestral states. In contrast to
        marginal reconstructions, this only needs to compare and multiply LH and
        can hence operate in log space.

        Keyword Args
        ------------

         store_lh : bool
            If True, all likelihoods will be stored for all nodes. Useful for
            testing, diagnostics and if special post-processing is required.

         verbose : int
            How verbose the output should be

        """
        N_diff = 0 # number of sites differ from perv reconstruction
        L = self.tree.get_terminals()[0].cseq.shape[0]
        n_states = self.gtr.alphabet.shape[0]

        self.logger("TreeAnc._ml_anc_joint: type of reconstruction: Joint", 2)

        self.logger("TreeAnc._ml_anc_joint: Walking up the tree, computing likelihoods... ", 3)
        # for the internal nodes, scan over all states j of this node, maximize the likelihood
        for node in self.tree.find_clades(order='postorder'):
            if node.up is None:
                node.joint_Cx=None # not needed for root
                continue

            # preallocate storage
            node.joint_Lx = np.zeros((L, n_states))             # likelihood array
            node.joint_Cx = np.zeros((L, n_states), dtype=int)  # max LH indices
            branch_len = self._branch_length_to_gtr(node)
            # transition matrix from parent states to the current node states.
            # denoted as Pij(i), where j - parent state, i - node state
            log_transitions = np.log(self.gtr.expQt(branch_len))

            if node.is_terminal():
                msg_from_children = np.log(np.maximum(seq_utils.seq2prof(node.cseq, self.gtr.profile_map), ttconf.TINY_NUMBER))
                msg_from_children[np.isnan(msg_from_children) | np.isinf(msg_from_children)] = -ttconf.BIG_NUMBER
            else:
                # Product (sum-Log) over all child subtree likelihoods.
                # this is prod_ch L_x(i)
                msg_from_children = np.sum(np.stack([c.joint_Lx for c in node.clades], axis=0), axis=0)

            # for every possible state of the parent node,
            # get the best state of the current node
            # and compute the likelihood of this state
            for char_i, char in enumerate(self.gtr.alphabet):
                # Pij(i) * L_ch(i) for given parent state j
                msg_to_parent = (log_transitions.T[char_i, :] + msg_from_children)
                # For this parent state, choose the best state of the current node:
                node.joint_Cx[:, char_i] = msg_to_parent.argmax(axis=1)
                # compute the likelihood of the best state of the current node
                # given the state of the parent (char_i)
                node.joint_Lx[:, char_i] = msg_to_parent.max(axis=1)

        # root node profile = likelihood of the total tree
        msg_from_children = np.sum(np.stack([c.joint_Lx for c in self.tree.root.clades], axis = 0), axis=0)
        # Pi(i) * Prod_ch Lch(i)
        self.tree.root.joint_Lx = msg_from_children + np.log(self.gtr.Pi)
        normalized_profile = (self.tree.root.joint_Lx.T - self.tree.root.joint_Lx.max(axis=1)).T

        # choose sequence characters from this profile.
        # treat root node differently to avoid piling up mutations on the longer branch
        if sample_from_profile=='root':
            root_sample_from_profile = True
        elif isinstance(sample_from_profile, bool):
            root_sample_from_profile = sample_from_profile

        seq, anc_lh_vals, idxs = seq_utils.prof2seq(np.exp(normalized_profile),
                                    self.gtr, sample_from_prof = root_sample_from_profile)

        # compute the likelihood of the most probable root sequence
        self.tree.sequence_LH = np.choose(idxs, self.tree.root.joint_Lx.T)
        self.tree.sequence_joint_LH = (self.tree.sequence_LH*self.multiplicity).sum()
        self.tree.root.cseq = seq
        self.tree.root.seq_idx = idxs
        if final:
            if self.is_vcf:
                self.tree.root.sequence = self.dict_sequence(self.tree.root)
            else:
                self.tree.root.sequence = self.expanded_sequence(self.tree.root)

        self.logger("TreeAnc._ml_anc_joint: Walking down the tree, computing maximum likelihood sequences...",3)
        # for each node, resolve the conditioning on the parent node
        for node in self.tree.find_clades(order='preorder'):

            # root node has no mutations, everything else has been alread y set
            if node.up is None:
                node.mutations = []
                continue

            # choose the value of the Cx(i), corresponding to the state of the
            # parent node i. This is the state of the current node
            node.seq_idx = np.choose(node.up.seq_idx, node.joint_Cx.T)
            # reconstruct seq, etc
            tmp_sequence = np.choose(node.seq_idx, self.gtr.alphabet)
            if hasattr(node, 'sequence') and node.cseq is not None:
                N_diff += (tmp_sequence!=node.cseq).sum()
            else:
                N_diff += L

            node.cseq = tmp_sequence
            if final:
                node.mutations = self.get_mutations(node)
                if self.is_vcf:
                    node.sequence = self.dict_sequence(node)
                else:
                    node.sequence = self.expanded_sequence(node)


        self.logger("TreeAnc._ml_anc_joint: ...done", 3)
        if store_compressed:
            self._store_compressed_sequence_pairs()

        # do clean-up
        if not debug:
            for node in self.tree.find_clades(order='preorder'):
                del node.joint_Lx
                del node.joint_Cx
                del node.seq_idx

        return N_diff


    def _store_compressed_sequence_to_node(self, node):
        """
        make a compressed representation of a pair of sequences only counting
        the number of times a particular pair of states (e.g. (A,T)) is observed
        the the aligned sequences of parent and child.

        Parameters
        -----------

         node : PhyloTree.Clade
            Tree node. **Note** because the method operates
            on the sequences on both sides of a branch, sequence reconstruction
            must be performed prior to calling this method.

        """
        seq_pairs, multiplicity = self.gtr.compress_sequence_pair(node.up.cseq,
                                              node.cseq,
                                              pattern_multiplicity = self.multiplicity,
                                              ignore_gaps = self.ignore_gaps)
        node.compressed_sequence = {'pair':seq_pairs, 'multiplicity':multiplicity}


    def _store_compressed_sequence_pairs(self):
        """
        Traverse the tree, and for each node store the compressed sequence pair.
        **Note** sequence reconstruction should be performed prior to calling
        this method.
        """
        self.logger("TreeAnc._store_compressed_sequence_pairs...",2)
        for node in self.tree.find_clades():
            if node.up is None:
                continue
            self._store_compressed_sequence_to_node(node)
        self.logger("TreeAnc._store_compressed_sequence_pairs...done",3)


###################################################################
### Branch length
###################################################################
    def optimize_branch_len(self, **kwargs):
        self.optimize_branch_length(**kwargs)

    def optimize_branch_length(self, **kwargs):
        """
        Perform optimization for the branch lengths of the whole tree or any
        subtree. **Note** this method assumes that each node stores information
        about its sequence as numpy.array object (node.sequence attribute).
        Therefore, before calling this method, sequence reconstruction with
        either of the available models must be performed.

        Keyword Args
        ------------

         verbose : int
            Output detalization

         store_old : bool
            If True, the old lenths will be saved in node._old_dist attribute.
            Useful for testing, and special post-processing.

        Returns
        -------
         None, the phylogenetic tree is modified in-place.

        """

        self.logger("TreeAnc.optimize_branch_length: running branch length optimization...",1)

        verbose = 0
        store_old_dist = False

        if 'verbose' in kwargs:
            verbose = int(kwargs['verbose'])
        if 'store_old' in kwargs:
            store_old_dist = kwargs['store_old'] == True

        for node in self.tree.find_clades(order='postorder'):
            if node.up is None: continue # this is the root
            if store_old_dist:
                node._old_length = node.branch_length

            new_len = self.optimal_branch_length(node)

            if new_len < 0:
                continue

            self.logger("Optimization results: old_len=%.4f, new_len=%.4f "
                   " Updating branch length..."%(node.branch_length, new_len), 5)

            node.branch_length = new_len
            node.mutation_length=new_len

        # as branch lengths changed, the params must be fixed
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        self._prepare_nodes()


    def optimal_branch_length(self, node):
        '''
        Calculate optimal branch length given the sequences of node and parent

        Parameters
        -----------

         node : PhyloTree.Clade
            TreeNode, attached to the branch.

        Returns
        -------

         new_len : float
            Optimal length of the given branch

        '''
        if node.up is None:
            return self.one_mutation

        parent = node.up
        if hasattr(node, 'compressed_sequence'):
            new_len = self.gtr.optimal_t_compressed(node.compressed_sequence['pair'],
                                                    node.compressed_sequence['multiplicity'])
        else:
            new_len = self.gtr.optimal_t(parent.cseq, node.cseq,
                                         pattern_multiplicity=self.multiplicity,
                                         ignore_gaps=self.ignore_gaps)
        return new_len


    def prune_short_branches(self):
        """
        If the branch length is less than the minimal value, remove the branch
        from the tree. **Requires** the ancestral sequence reconstruction
        """
        self.logger("TreeAnc.prune_short_branches: pruning short branches (max prob at zero)...", 1)
        for node in self.tree.find_clades():
            if node.up is None or node.is_terminal():
                continue

            # probability of the two seqs separated by zero time is not zero
            if self.gtr.prob_t(node.up.cseq, node.cseq, 0.0,
                               pattern_multiplicity=self.multiplicity) > 0.1:
                # re-assign the node children directly to its parent
                node.up.clades = [k for k in node.up.clades if k != node] + node.clades
                for clade in node.clades:
                    clade.up = node.up

    def optimize_sequences_and_branch_length(self,*args, **kwargs):
        """This method is a schortcut for :py:meth:`optimize_seq_and_branch_len`

        Iteratively set branch lengths and reconstruct ancestral sequences until
        the values of either former or latter do not change. The algorithm assumes
        knowing only the topology of the tree, and requires that sequences are assigned
        to all leaves of the tree. The first step is to pre-reconstruct ancestral
        states using Fitch reconstruction algorithm or ML using existing branch length
        estimates. Then, optimize branch lengths and re-do reconstruction until
        convergence using ML method.
        """
        self.optimize_seq_and_branch_len(*args,**kwargs)

    def optimize_seq_and_branch_len(self,reuse_branch_len=True,prune_short=True,
                                    max_iter=5, infer_gtr=False, **kwargs):
        """
        Iteratively set branch lengths and reconstruct ancestral sequences until
        the values of either former or latter do not change. The algorithm assumes
        knowing only the topology of the tree, and requires that sequences are assigned
        to all leaves of the tree. The first step is to pre-reconstruct ancestral
        states using Fitch reconstruction algorithm or ML using existing branch length
        estimates. Then, optimize branch lengths and re-do reconstruction until
        convergence using ML method.

        Parameters
        -----------

         reuse_branch_len : bool, default True
            If True, rely on the initial branch lenghts, and start with the
            Maximum-likelihood ancestral sequence inference using existing branch
            lengths. Otherwise, initial reconstruction of ancestral states with
            Fitch algorithm, which uses only the tree topology.

         prune_short : bool, default True
            If True, the branches with zero optimal length will be pruned from
            the tree hence creating polytomies. The polytomies could be further
            processde using resolve_polytomies from the TreeTime class.

        """
        self.logger("TreeAnc.optimize_sequences_and_branch_length: sequences...", 1)
        if reuse_branch_len:
            N_diff = self.reconstruct_anc(method='ml', infer_gtr=infer_gtr, **kwargs)
        else:
            N_diff = self.reconstruct_anc(method='fitch', infer_gtr=infer_gtr, **kwargs)

        self.optimize_branch_len(verbose=0, store_old=False)

        n = 0
        while n<max_iter:
            n += 1
            if prune_short:
                self.prune_short_branches()
            N_diff = self.reconstruct_anc(method='ml', infer_gtr=False,**kwargs)

            self.logger("TreeAnc.optimize_sequences_and_branch_length: Iteration %d."
                   " #Nuc changed since prev reconstructions: %d" %(n, N_diff), 2)

            if N_diff < 1:
                break
            self.optimize_branch_len(verbose=0, store_old=False)

        self.tree.unconstrained_sequence_LH = (self.tree.sequence_LH*self.multiplicity).sum()
        self._prepare_nodes() # fix dist2root and up-links after reconstruction
        self.logger("TreeAnc.optimize_sequences_and_branch_length: Unconstrained sequence LH:%f" % self.tree.unconstrained_sequence_LH , 2)
        return

###############################################################################
### Utility functions
###############################################################################
    def get_reconstructed_alignment(self):
        """
        Get the multiple sequence alignment including reconstructed sequences for
        the internal nodes.
        """
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        self.logger("TreeAnc.get_reconstructed_alignment ...",2)
        if not hasattr(self.tree.root, 'sequence'):
            self.logger("TreeAnc.reconstructed_alignment... reconstruction not yet done",3)
            self.reconstruct_anc('ml')

        new_aln = MultipleSeqAlignment([SeqRecord(id=n.name, seq=Seq("".join(n.sequence)), description="")
                                        for n in self.tree.find_clades()])

        return new_aln

    def get_tree_dict(self, keep_var_ambigs=False):
        """
        For VCF-based objects, returns a nested dict with all information required to
        reconstruct sequences for all nodes (terminal and internal) in the format:
        {'reference':'AGCTCGA..A',
         'sequences': { 'seq1':{4:'A', 7:'-'}, 'seq2':{100:'C'} },
         'positions': [1,4,7,10,100...],
         'inferred_const_sites': [7,100....]     <this is optional>
        }
         self.inferred_const_sites

        Reference being the reference sequence to which the variable sites are mapped;
        sequence containing a dict for each sequence with the position and base of
        mutations; and positions containing a list of all the variable positions.
        If included, inferred_const_sites is positions that were constant except
        ambiguous bases, which were converted into constant sites (ex: 'AAAN' -> 'AAAA')

        keep_var_ambigs : boolean
            If true, generates dict sequence based on the *original* _compressed_ sequence, which
            may include ambiguities. Note sites that only have 1 unambiguous base and ambiguous
            bases ("AAAAANN") are stripped of ambiguous bases *before* compression, so ambiguous
            bases will *not* be preserved.

        EBH 7 Dec 2017
        """
        if self.is_vcf:
            tree_dict = {}
            tree_dict['reference'] = self.ref
            tree_dict['positions'] = self.nonref_positions

            tree_aln = {}
            for n in self.tree.find_clades():
                if hasattr(n, 'sequence'):
                    if keep_var_ambigs: #regenerate dict to include ambig bases
                        tree_aln[n.name] = self.dict_sequence(n, keep_var_ambigs)
                    else:
                        tree_aln[n.name] = n.sequence

            tree_dict['sequences'] = tree_aln

            if len(self.inferred_const_sites) != 0:
                tree_dict['inferred_const_sites'] = self.inferred_const_sites

            return tree_dict
        else:
            raise("A dict can only be returned for trees created with VCF-input!")



if __name__=="__main__":

    from Bio import Phylo
    from StringIO import StringIO
    from Bio import Phylo,AlignIO

    tiny_tree = Phylo.read(StringIO("((A:.0060,B:.30)C:.030,D:.020)E:.004;"), 'newick')
    tiny_aln = AlignIO.read(StringIO(">A\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                                     ">B\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                                     ">C\nAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT\n"
                                     ">D\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n"
                                     ">E\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"), 'fasta')

    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']),
                       pi = np.array([0.25, 0.95, 0.005, 0.05]), W=np.ones((4,4)))

    myTree = TreeAnc(gtr=mygtr, tree = tiny_tree,
                        aln =tiny_aln, verbose = 4)

    logLH = myTree.ancestral_likelihood()
    LH = np.exp(logLH)
    print ("Net probability (for all possible realizations): " + str(np.exp(logLH).sum()))
    print (np.exp(logLH))
