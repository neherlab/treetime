from __future__ import print_function, division, absolute_import
import time, sys
import gc
from collections import defaultdict
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from Bio import AlignIO
from treetime import config as ttconf
from treetime import MissingDataError,UnknownMethodError
from .seq_utils import seq2prof, seq2array, prof2seq, normalize_profile, extend_profile
from .gtr import GTR
from .gtr_site_specific import GTR_site_specific
from .sequence_data import SequenceData

def compressed_sequence(node):
    if node.name in node.tt.data.compressed_alignment and (not node.tt.reconstructed_tip_sequences):
        return node.tt.data.compressed_alignment[node.name]
    elif hasattr(node, '_cseq'):
        return node._cseq
    elif node.is_terminal(): # node without sequence when tip-reconstruction is off.
        return None
    elif hasattr(node, '_cseq'):
        return node._cseq
    else:
        raise ValueError('Ancestral sequences are not yet inferred')

def mutations(node):
    """
    Get the mutations on a tree branch. Take compressed sequences from both sides
    of the branch (attached to the node), compute mutations between them, and
    expand these mutations to the positions in the real sequences.
    """
    if node.up is None:
        return []
    elif (not node.tt.reconstructed_tip_sequences) and node.name in node.tt.data.aln:
        return node.tt.data.differences(node.up.cseq, node.tt.data.aln[node.name], seq2_compressed=False)
    elif node.is_terminal() and (node.name not in node.tt.data.aln):
        return []
    else:
        return node.tt.data.differences(node.up.cseq, node.cseq)

string_types = [str] if sys.version_info[0]==3 else [str, unicode]
Clade.sequence = property(lambda x: x.tt.sequence(x, as_string=False))
Clade.cseq = property(compressed_sequence)
Clade.mutations = property(mutations)


class TreeAnc(object):
    """
    Class defines simple tree object with basic interface methods: reading and
    saving from/to files, initializing leaves with sequences from the
    alignment, making ancestral state inference
    """

    def __init__(self, tree=None, aln=None, gtr=None, fill_overhangs=True,
                ref=None, verbose = ttconf.VERBOSE, ignore_gaps=True,
                convert_upper=True, seq_multiplicity=None, log=None,
                 compress=True, seq_len=None,
                **kwargs):
        """
        TreeAnc constructor. It prepares the tree, attaches sequences to the leaf nodes,
        and sets some configuration parameters.

        Parameters
        ----------
        tree : str, Bio.Phylo.Tree
           Phylogenetic tree. String passed is interpreted as a filename with
           a tree in a standard format that can be parsed by the Biopython Phylo module.

        aln : str, Bio.Align.MultipleSequenceAlignment, dict
           Sequence alignment. If a string passed, it is interpreted as the
           filename to read Biopython alignment from. If a dict is given,
           this is assumed to be the output of vcf_utils.read_vcf which
           specifies for each sequence the differences from a reference

        gtr : str, GTR
           GTR model object. If string passed, it is interpreted as the type of
           the GTR model. A new GTR instance will be created for this type.

        fill_overhangs : bool
           In some cases, the missing data on both ends of the alignment is
           filled with the gap sign('-'). If set to True, the end-gaps are converted to "unknown"
           characters ('N' for nucleotides, 'X' for aminoacids). Otherwise, the alignment is treated as-is

        ref : None, optional
            Reference sequence used in VCF mode

        verbose : int
           Verbosity level as number from 0 (lowest) to 10 (highest).

        ignore_gaps : bool
           Ignore gaps in branch length calculations

        convert_upper : bool, optional
            Description

        seq_multiplicity : dict
           If individual nodes in the tree correspond to multiple sampled sequences
           (i.e. read count in a deep sequencing experiment), these can be
           specified as a dictionary. This currently only affects rooting and
           can be used to weigh individual tips by abundance or important during root search.

        compress : bool, optional
            reduce identical alignment columns to one (not useful when
            inferring site specific GTR models).

        seq_len : int, optional
            length of the sequence. this is inferred from the input alignment or the reference
            sequence in most cases but can be specified for other applications.

        **kwargs
           Keyword arguments to construct the GTR model

         .. Note::
               Some GTR types require additional configuration parameters.
               If the new GTR is being instantiated, these parameters are expected
               to be passed as kwargs. If nothing is passed, the default values are
               used, which might cause unexpected results.

        Raises
        ------
        AttributeError
            If no tree is passed in


        """
        if tree is None:
            raise TypeError("TreeAnc requires a tree!")
        self.t_start = time.time()
        self.verbose = verbose
        self.log = log
        self.ok = False
        self.data = None
        self.log_messages = set()
        self.logger("TreeAnc: set-up",1)
        self._internal_node_count = 0
        self.use_mutation_length = False
        self.ignore_gaps = ignore_gaps
        self.reconstructed_tip_sequences = False
        self.sequence_reconstruction = None

        self._tree = None
        self.tree = tree
        if tree is None:
            raise ValueError("TreeAnc: tree loading failed! exiting")

        # set up GTR model
        self._gtr = None
        self.set_gtr(gtr or 'JC69', **kwargs)

        # set alignment and attach sequences to tree on success.
        # otherwise self.data.aln will be None
        self.data = SequenceData(aln, ref=ref, logger=self.logger, compress=compress,
                                convert_upper=convert_upper, fill_overhangs=fill_overhangs, ambiguous=self.gtr.ambiguous,
                                 sequence_length=seq_len)

        if self.gtr.is_site_specific and self.data.compress:
            raise TypeError("TreeAnc: sequence compression and site specific gtr models are incompatible!" )

        if self.data.aln and self.tree:
            self._check_alignment_tree_gtr_consistency()


    def logger(self, msg, level, warn=False, only_once=False):
        """
        Print log message *msg* to stdout.

        Parameters
        -----------

         msg : str
            String to print on the screen

         level : int
            Log-level. Only the messages with a level higher than the
            current verbose level will be shown.

         warn : bool
            Warning flag. If True, the message will be displayed
            regardless of its log-level.

        """
        if only_once and msg in self.log_messages:
            return

        self.log_messages.add(msg)

        lw=80
        if level<self.verbose or (warn and level<=self.verbose):
            from textwrap import fill
            dt = time.time() - self.t_start
            outstr = '\n' if level<2 else ''
            initial_indent = format(dt, '4.2f')+'\t' + level*'-'
            subsequent_indent = " "*len(format(dt, '4.2f')) + "\t" + " "*level
            outstr += fill(msg, width=lw, initial_indent=initial_indent, subsequent_indent=subsequent_indent)
            print(outstr, file=sys.stdout)


####################################################################
## SET-UP
####################################################################
    @property
    def leaves_lookup(self):
        """
        The :code:`{leaf-name:leaf-node}` dictionary. It enables fast
        search of a tree leaf object by its name.
        """
        return self._leaves_lookup

    @property
    def gtr(self):
        """
        :setter: Sets the GTR object passed in
        :getter: Returns the current GTR object

        """
        return self._gtr

    @gtr.setter
    def gtr(self, value):
        """
        Parameters
        -----------
         value : GTR
            the new GTR object
        """
        if not isinstance(value, (GTR, GTR_site_specific)):
            raise TypeError("GTR instance expected")
        self._gtr = value


    def set_gtr(self, in_gtr, **kwargs):
        """
        Create new GTR model if needed, and set the model as an attribute of the
        TreeAnc class

        Parameters
        -----------
         in_gtr : str, GTR
            The gtr model to be assigned. If string is passed,
            it is taken as the name of a standard GTR model, and is
            attempted to be created through :code:`GTR.standard()` interface. If a
            GTR instance is passed, it is set directly .
         **kwargs
            Keyword arguments to construct the GTR model. If none are passed, defaults
            are assumed.

        """
        if isinstance(in_gtr, str):
            self._gtr = GTR.standard(model=in_gtr, **kwargs)
            self._gtr.logger = self.logger
        elif isinstance(in_gtr, (GTR, GTR_site_specific)):
            self._gtr = in_gtr
            self._gtr.logger=self.logger
        else:
            self.logger("TreeAnc.gtr_setter: can't interpret GTR model", 1, warn=True)
            raise TypeError("Cannot set GTR model in TreeAnc class: GTR or "
                "string expected")

        if self._gtr.ambiguous is None:
            self.fill_overhangs=False

    @property
    def aln(self):
        '''
        :setter: Sets the alignment
        :getter: Returns the alignment
        '''
        return self.data.aln

    @aln.setter
    def aln(self,in_aln):
        self.data.aln=in_aln
        if self.tree:
            self._check_alignment_tree_gtr_consistency()


    @property
    def tree(self):
        """
        The phylogenetic tree currently used by the TreeAnc.

        :setter: Sets the tree. Directly if passed as Phylo.Tree, or by reading from \
        file if passed as a str.
        :getter: Returns the tree as a Phylo.Tree object

        """
        return self._tree


    @tree.setter
    def tree(self, in_tree):
        '''
        assigns a tree to the internal self._tree variable. The tree is either
        loaded from file (if in_tree is str) or assigned (if in_tree is a Phylo.tree)
        '''
        from os.path import isfile
        self._tree = None
        if isinstance(in_tree, Phylo.BaseTree.Tree):
            self._tree = in_tree
        elif type(in_tree) in string_types and isfile(in_tree):
            try:
                self._tree=Phylo.read(in_tree, 'newick')
            except:
                fmt = in_tree.split('.')[-1]
                if fmt in ['nexus', 'nex']:
                    self._tree=Phylo.read(in_tree, 'nexus')
                else:
                    raise ValueError('TreeAnc: could not load tree, format needs to be nexus or newick! input was '+str(in_tree))
        else:
            raise ValueError('TreeAnc: could not load tree! input was '+str(in_tree))

        if self._tree.count_terminals()<3:
            raise ValueError('TreeAnc: tree in %s as only %d tips. Please check your tree!'%(str(in_tree), self._tree.count_terminals()))

        # remove all existing sequence attributes
        for node in self._tree.find_clades():
            node.branch_length = node.branch_length if node.branch_length else 0.0
            if hasattr(node, "_cseq"):
                node.__delattr__("_cseq")
            node.original_length = node.branch_length
            node.mutation_length = node.branch_length
        self.prepare_tree()

        if self.data:
            self._check_alignment_tree_gtr_consistency()

        return ttconf.SUCCESS


    @property
    def one_mutation(self):
        """
        Returns
        -------
        float
            inverse of the uncompressed sequene length - length scale for short branches
        """
        return 1.0/self.data.full_length if self.data.full_length else np.nan

    @one_mutation.setter
    def one_mutation(self,om):
        self.logger("TreeAnc: one_mutation can't be set",1)


    @property
    def seq_len(self):
        return self.data.full_length


    @property
    def sequence_length(self):
        return self.data.full_length


    def _check_alignment_tree_gtr_consistency(self):
        '''
        For each node of the tree, check whether there is a sequence available
        in the alignment and assign this sequence as a character array
        '''
        if len(self.tree.get_terminals()) != len(self.data.aln):
            self.logger("**WARNING: Number of tips in tree differs from number of sequences in alignment!**", 3, warn=True)
        failed_leaves= 0

        # loop over leaves and assign multiplicities of leaves (e.g. number of identical reads)
        for l in self.tree.get_terminals():
            if l.name in self.data.seq_multiplicity:
                l.count = self.data.seq_multiplicity[l.name]
            else:
                l.count = 1.0

        # loop over tree, and assign sequences
        for l in self.tree.find_clades():
            if hasattr(l, 'branch_state'): del l.branch_state
            if l.name not in self.data.compressed_alignment and l.is_terminal():
                self.logger("***WARNING: TreeAnc._attach_sequences_to_nodes: NO SEQUENCE FOR LEAF: %s" % l.name, 0, warn=True)
                failed_leaves += 1
                if failed_leaves > self.tree.count_terminals()/3:
                    raise MissingDataError("TreeAnc._check_alignment_tree_gtr_consistency: At least 30\\% terminal nodes cannot be assigned a sequence!\n"
                                           "Are you sure the alignment belongs to the tree?")
            else: # could not assign sequence for internal node - is OK
                pass

        if failed_leaves:
            self.logger("***WARNING: TreeAnc: %d nodes don't have a matching sequence in the alignment."
                        " POSSIBLE ERROR."%failed_leaves, 0, warn=True)

        # extend profile to contain additional unknown characters
        extend_profile(self.gtr, [self.data.ref] if self.data.is_sparse
                        else self.data.aln.values(), logger=self.logger)
        self.ok = True


    def prepare_tree(self):
        """
        Set link to parent and calculate distance to root for all tree nodes.
        Should be run once the tree is read and after every rerooting,
        topology change or branch length optimizations.
        """
        self.sequence_reconstruction = False
        self.tree.root.branch_length = 0.001
        self.tree.root.mutation_length = self.tree.root.branch_length
        self.tree.ladderize()
        self._prepare_nodes()
        self._leaves_lookup = {node.name:node for node in self.tree.get_terminals()}


    def _prepare_nodes(self):
        """
        Set auxilliary parameters to every node of the tree.
        """
        self.tree.root.up = None
        self.tree.root.tt = self
        self.tree.root.bad_branch=self.tree.root.bad_branch if hasattr(self.tree.root, 'bad_branch') else False

        name_set = {n.name for n in self.tree.find_clades() if n.name}
        internal_node_count = 0
        for clade in self.tree.get_nonterminals(order='preorder'): # parents first
            if clade.name is None:
                tmp = "NODE_" + format(internal_node_count, '07d')
                while tmp in name_set:
                    internal_node_count += 1
                    tmp = "NODE_" + format(internal_node_count, '07d')
                clade.name = tmp
                name_set.add(clade.name)
            internal_node_count+=1
            for c in clade.clades:
                c.up = clade
                c.tt = self

        for clade in self.tree.find_clades(order='postorder'): # children first
            if clade.is_terminal():
                clade.bad_branch = clade.bad_branch if hasattr(clade, 'bad_branch') else False
            else:
                clade.bad_branch = all([c.bad_branch for c in clade])

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
                c.dist2root = clade.dist2root + c.mutation_length



####################################################################
## END SET-UP
####################################################################


###################################################################
### ancestral reconstruction
###################################################################
    def reconstruct_anc(self,*args, **kwargs):
        """Shortcut for :py:meth:`treetime.TreeAnc.infer_ancestral_sequences`
        """
        return self.infer_ancestral_sequences(*args,**kwargs)


    def infer_ancestral_sequences(self, method='probabilistic', infer_gtr=False,
                                  marginal=False, reconstruct_tip_states=False, **kwargs):
        """Reconstruct ancestral sequences

        Parameters
        ----------
        method : str
           Method to use. Supported values are "parsimony", "fitch", "probabilistic" and "ml"
        infer_gtr : bool
           Infer a GTR model before reconstructing the sequences
        marginal : bool
           Assign sequences that are most likely after averaging over all other nodes
           instead of the jointly most likely sequences.
        reconstruct_tip_states : bool, optional
            Reconstruct sequences of terminal nodes/leaves, thereby replacing ambiguous
            characters with the inferred base/state. default: False
        **kwargs
            additional keyword arguments that are passed down to :py:meth:`TreeAnc.infer_gtr` and :py:meth:`TreeAnc._ml_anc`

        Returns
        -------
        N_diff : int
           Number of nucleotides different from the previous
           reconstruction.  If there were no pre-set sequences, returns N*L

        """
        if not self.ok:
            raise MissingDataError("TreeAnc.infer_ancestral_sequences: ERROR, sequences or tree are missing")

        self.logger("TreeAnc.infer_ancestral_sequences with method: %s, %s"%(method, 'marginal' if marginal else 'joint'), 1)

        if not reconstruct_tip_states:
            self.logger("WARNING: Previous versions of TreeTime (<0.7.0) RECONSTRUCTED sequences"
                        " of tips at positions with AMBIGUOUS bases. This resulted in"
                        " unexpected behavior is some cases and is no longer done by default."
                        " If you want to replace those ambiguous sites with their most likely state,"
                        " rerun with `reconstruct_tip_states=True` or `--reconstruct-tip-states`.", 0, warn=True, only_once=True)

        if method.lower() in ['ml', 'probabilistic']:
            if marginal:
                _ml_anc = self._ml_anc_marginal
            else:
                _ml_anc = self._ml_anc_joint
        elif method.lower() in ['fitch', 'parsimony']:
            _ml_anc = self._fitch_anc
        else:
            raise ValueError("Reconstruction method needs to be in ['ml', 'probabilistic', 'fitch', 'parsimony'], got '{}'".format(method))

        if infer_gtr:
            self.infer_gtr(marginal=marginal, **kwargs)
            N_diff = _ml_anc(reconstruct_tip_states=reconstruct_tip_states, **kwargs)
        else:
            N_diff = _ml_anc(reconstruct_tip_states=reconstruct_tip_states, **kwargs)

        return N_diff


###################################################################
### FITCH
###################################################################
    def _fitch_anc(self, **kwargs):
        """
        Reconstruct ancestral states using Fitch's algorithm. It implements
        the iteration from leaves to the root constructing the Fitch profiles for
        each character of the sequence, and then by propagating from the root
        to the leaves, reconstructs the sequences of the internal nodes.

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

        L = self.data.compressed_length

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
        self.tree.root._cseq = np.array([k[np.random.randint(len(k)) if len(k)>1 else 0]
                                           for k in self.tree.root.state])


        self.logger("TreeAnc._fitch_anc: Walking down the self.tree, generating sequences from the "
                         "Fitch profiles.", 2)
        N_diff = 0
        for node in self.tree.get_nonterminals(order='preorder'):
            if node.up != None: # not root
                sequence =  np.array([node.up._cseq[i]
                        if node.up._cseq[i] in node.state[i]
                        else node.state[i][0] for i in range(L)])

                if self.sequence_reconstruction:
                    N_diff += (sequence!=node.cseq).sum()
                else:
                    N_diff += L
                node._cseq = sequence

            del node.state # no need to store Fitch states

        self.sequence_reconstruction = 'parsimony'
        self.logger("Done ancestral state reconstruction",3)
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
    def sequence_LH(self, pos=None, full_sequence=False):
        """return the likelihood of the observed sequences given the tree

        Parameters
        ----------
        pos : int, optional
            position in the sequence, if none, the sum over all positions will be returned
        full_sequence : bool, optional
            does the position refer to the full or compressed sequence, by default compressed sequence is assumed.

        Returns
        -------
        float
            likelihood
        """
        if not hasattr(self.tree, "total_sequence_LH"):
            self.logger("TreeAnc.sequence_LH: you need to run marginal ancestral inference first!", 1)
            self.infer_ancestral_sequences(marginal=True)
        if pos is not None:
            if full_sequence:
                compressed_pos = self.data.full_to_compressed_sequence_map[pos]
            else:
                compressed_pos = pos
            return self.tree.sequence_LH[compressed_pos]
        else:
            return self.tree.total_sequence_LH


    def ancestral_likelihood(self):
        """
        Calculate the likelihood of the given realization of the sequences in
        the tree

        Returns
        -------

         log_lh : float
            The tree likelihood given the sequences
        """
        log_lh = np.zeros(self.data.multiplicity.shape[0])
        for node in self.tree.find_clades(order='postorder'):

            if node.up is None: #  root node
                # 0-1 profile
                profile = seq2prof(node.cseq, self.gtr.profile_map)
                # get the probabilities to observe each nucleotide
                profile *= self.gtr.Pi
                profile = profile.sum(axis=1)
                log_lh += np.log(profile) # product over all characters
                continue

            t = node.branch_length

            indices = np.array([(self.gtr.state_index[a], self.gtr.state_index[b])
                                 for a, b in zip(node.up.cseq, node.cseq)])

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


    def _ml_anc_marginal(self, sample_from_profile=False,
                         reconstruct_tip_states=False, debug=False, **kwargs):
        """
        Perform marginal ML reconstruction of the ancestral states. In contrast to
        joint reconstructions, this needs to access the probabilities rather than only
        log probabilities and is hence handled by a separate function.

        Parameters
        ----------
         sample_from_profile : bool or str
            assign sequences probabilistically according to the inferred probabilities
            of ancestral states instead of to their ML value. This parameter can also
            take the value 'root' in which case probabilistic sampling will happen
            at the root but at no other node.
        reconstruct_tip_states : bool, default False
            reconstruct sequence assigned to leaves, will replace ambiguous characters
            with the most likely definite character. Note that this will affect the mutations
            assigned to branches.
        """
        self.logger("TreeAnc._ml_anc_marginal: type of reconstruction: Marginal", 2)
        self.postorder_traversal_marginal()

        # choose sequence characters from this profile.
        # treat root node differently to avoid piling up mutations on the longer branch
        if sample_from_profile=='root':
            root_sample_from_profile = True
            other_sample_from_profile = False
        elif isinstance(sample_from_profile, bool):
            root_sample_from_profile = sample_from_profile
            other_sample_from_profile = sample_from_profile

        self.total_LH_and_root_sequence(sample_from_profile=root_sample_from_profile,
                                        assign_sequence=True)

        N_diff = self.preorder_traversal_marginal(reconstruct_tip_states=reconstruct_tip_states,
                                                  sample_from_profile=other_sample_from_profile,
                                                  assign_sequence=True)
        self.logger("TreeAnc._ml_anc_marginal: ...done", 3)

        self.reconstructed_tip_sequences = reconstruct_tip_states
        # do clean-up:
        if not debug:
            for node in self.tree.find_clades():
                try:
                    del node.marginal_log_Lx
                    del node.marginal_subtree_LH_prefactor
                except:
                    pass
        gc.collect()
        self.sequence_reconstruction = 'marginal'
        return N_diff


    def total_LH_and_root_sequence(self, sample_from_profile=False, assign_sequence=False):
        self.logger("Computing root node sequence and total tree likelihood...",3)
        # Msg to the root from the distant part (equ frequencies)
        if len(self.gtr.Pi.shape)==1:
            self.tree.root.marginal_outgroup_LH = np.repeat([self.gtr.Pi], self.data.compressed_length, axis=0)
        else:
            self.tree.root.marginal_outgroup_LH = np.copy(self.gtr.Pi.T)

        self.tree.root.marginal_profile, pre = normalize_profile(self.tree.root.marginal_outgroup_LH*self.tree.root.marginal_subtree_LH)
        marginal_LH_prefactor = self.tree.root.marginal_subtree_LH_prefactor + pre

        self.tree.sequence_LH = marginal_LH_prefactor
        self.tree.total_sequence_LH = (self.tree.sequence_LH*self.data.multiplicity).sum()
        self.tree.sequence_marginal_LH = self.tree.total_sequence_LH
        if assign_sequence:
            seq, prof_vals, idxs = prof2seq(self.tree.root.marginal_profile,
                                        self.gtr, sample_from_prof=sample_from_profile, normalize=False)
            self.tree.root._cseq = seq


    def postorder_traversal_marginal(self):
        L = self.data.compressed_length
        n_states = self.gtr.alphabet.shape[0]

        self.logger("Attaching sequence profiles to leafs... ", 3)
        #  set the leaves profiles. This doesn't ever need to be reassigned for leaves
        for leaf in self.tree.get_terminals():
            if not hasattr(leaf, "marginal_subtree_LH"):
                if leaf.name in self.data.compressed_alignment:
                    leaf.marginal_subtree_LH = seq2prof(self.data.compressed_alignment[leaf.name], self.gtr.profile_map)
                else:
                    leaf.marginal_subtree_LH = np.ones((L, n_states))
            if not hasattr(leaf, "marginal_subtree_LH_prefactor"):
                leaf.marginal_subtree_LH_prefactor = np.zeros(L)

        self.logger("Postorder: computing likelihoods... ", 3)
        # propagate leaves --> root, set the marginal-likelihood messages
        for node in self.tree.get_nonterminals(order='postorder'): #leaves -> root
            # regardless of what was before, set the profile to ones
            tmp_log_subtree_LH = np.zeros((L,n_states), dtype=float)
            node.marginal_subtree_LH_prefactor = np.zeros(L, dtype=float)
            for ch in node.clades:
                ch.marginal_log_Lx = self.gtr.propagate_profile(ch.marginal_subtree_LH,
                    self._branch_length_to_gtr(ch), return_log=True) # raw prob to transfer prob up
                tmp_log_subtree_LH += ch.marginal_log_Lx
                node.marginal_subtree_LH_prefactor += ch.marginal_subtree_LH_prefactor

            node.marginal_subtree_LH, offset = normalize_profile(tmp_log_subtree_LH, log=True)
            node.marginal_subtree_LH_prefactor += offset # and store log-prefactor


    def preorder_traversal_marginal(self, reconstruct_tip_states=False, sample_from_profile=False, assign_sequence=False):
        self.logger("Preorder: computing marginal profiles...",3)
        # propagate root -->> leaves, reconstruct the internal node sequences
        # provided the upstream message + the message from the complementary subtree
        N_diff = 0
        for node in self.tree.find_clades(order='preorder'):
            if node.up is None: # skip if node is root
                continue
            if hasattr(node, 'branch_state'): del node.branch_state

            # integrate the information coming from parents with the information
            # of all children my multiplying it to the prev computed profile
            node.marginal_outgroup_LH, pre = normalize_profile(np.log(np.maximum(ttconf.TINY_NUMBER, node.up.marginal_profile)) - node.marginal_log_Lx,
                                                               log=True, return_offset=False)
            if node.is_terminal() and (not reconstruct_tip_states): # skip remainder unless leaves are to be reconstructed
                continue

            tmp_msg_from_parent = self.gtr.evolve(node.marginal_outgroup_LH,
                                                 self._branch_length_to_gtr(node), return_log=False)
            node.marginal_profile, pre = normalize_profile(node.marginal_subtree_LH * tmp_msg_from_parent, return_offset=False)
            # choose sequence based maximal marginal LH.
            if assign_sequence:
                seq, prof_vals, idxs = prof2seq(node.marginal_profile, self.gtr,
                                                sample_from_prof=sample_from_profile, normalize=False)

                if self.sequence_reconstruction:
                    N_diff += (seq!=node.cseq).sum()
                else:
                    N_diff += self.data.compressed_length
                #assign new sequence
                node._cseq = seq

        return N_diff


    def _ml_anc_joint(self, sample_from_profile=False,
                      reconstruct_tip_states=False, debug=False, **kwargs):

        """
        Perform joint ML reconstruction of the ancestral states. In contrast to
        marginal reconstructions, this only needs to compare and multiply LH and
        can hence operate in log space.

        Parameters
        ----------
         sample_from_profile : str
            This parameter can take the value 'root' in which case probabilistic
            sampling will happen at the root. otherwise sequences at ALL nodes are
            set to the value that jointly optimized the likelihood.
        reconstruct_tip_states : bool, default False
            reconstruct sequence assigned to leaves, will replace ambiguous characters
            with the most likely definite character. Note that this will affect the mutations
            assigned to branches.

        """
        N_diff = 0 # number of sites differ from perv reconstruction
        L = self.data.compressed_length
        n_states = self.gtr.alphabet.shape[0]

        self.logger("TreeAnc._ml_anc_joint: type of reconstruction: Joint", 2)

        self.logger("TreeAnc._ml_anc_joint: Walking up the tree, computing likelihoods... ", 3)
        # for the internal nodes, scan over all states j of this node, maximize the likelihood
        for node in self.tree.find_clades(order='postorder'):
            if hasattr(node, 'branch_state'): del node.branch_state
            if node.up is None:
                node.joint_Cx=None # not needed for root
                continue

            branch_len = self._branch_length_to_gtr(node)
            # transition matrix from parent states to the current node states.
            # denoted as Pij(i), where j - parent state, i - node state
            log_transitions = np.log(np.maximum(ttconf.TINY_NUMBER, self.gtr.expQt(branch_len)))
            if node.is_terminal():
                if node.name in self.data.compressed_alignment:
                    tmp_prof = seq2prof(self.data.compressed_alignment[node.name], self.gtr.profile_map)
                    msg_from_children = np.log(np.maximum(tmp_prof, ttconf.TINY_NUMBER))
                else:
                    msg_from_children = np.zeros((L, n_states))

                msg_from_children[np.isnan(msg_from_children) | np.isinf(msg_from_children)] = -ttconf.BIG_NUMBER
            else:
                # Product (sum-Log) over all child subtree likelihoods.
                # this is prod_ch L_x(i)
                msg_from_children = np.sum(np.stack([c.joint_Lx for c in node.clades], axis=0), axis=0)

            # for every possible state of the parent node,
            # get the best state of the current node
            # and compute the likelihood of this state
            # preallocate storage
            node.joint_Lx = np.zeros((L, n_states))             # likelihood array
            node.joint_Cx = np.zeros((L, n_states), dtype=int)  # max LH indices
            for char_i, char in enumerate(self.gtr.alphabet):
                # Pij(i) * L_ch(i) for given parent state j
                msg_to_parent = (log_transitions[:,char_i].T + msg_from_children)
                # For this parent state, choose the best state of the current node:
                node.joint_Cx[:, char_i] = msg_to_parent.argmax(axis=1)
                # compute the likelihood of the best state of the current node
                # given the state of the parent (char_i)
                node.joint_Lx[:, char_i] = msg_to_parent.max(axis=1)

        # root node profile = likelihood of the total tree
        msg_from_children = np.sum(np.stack([c.joint_Lx for c in self.tree.root], axis = 0), axis=0)
        # Pi(i) * Prod_ch Lch(i)
        self.tree.root.joint_Lx = msg_from_children + np.log(self.gtr.Pi).T
        normalized_profile = (self.tree.root.joint_Lx.T - self.tree.root.joint_Lx.max(axis=1)).T

        # choose sequence characters from this profile.
        # treat root node differently to avoid piling up mutations on the longer branch
        if sample_from_profile=='root':
            root_sample_from_profile = True
        elif isinstance(sample_from_profile, bool):
            root_sample_from_profile = sample_from_profile

        seq, anc_lh_vals, idxs = prof2seq(np.exp(normalized_profile),
                                    self.gtr, sample_from_prof = root_sample_from_profile)

        # compute the likelihood of the most probable root sequence
        self.tree.sequence_LH = np.choose(idxs, self.tree.root.joint_Lx.T)
        self.tree.sequence_joint_LH = (self.tree.sequence_LH*self.data.multiplicity).sum()
        self.tree.root._cseq = seq
        self.tree.root.seq_idx = idxs

        self.logger("TreeAnc._ml_anc_joint: Walking down the tree, computing maximum likelihood sequences...",3)
        # for each node, resolve the conditioning on the parent node
        nodes_to_reconstruct = self.tree.get_nonterminals(order='preorder')
        if reconstruct_tip_states:
            nodes_to_reconstruct += self.tree.get_terminals()
            #TODO: Should we add tips without sequence here?

        for node in nodes_to_reconstruct:
            # root node has no mutations, everything else has been already set
            if node.up is None:
                continue

            # choose the value of the Cx(i), corresponding to the state of the
            # parent node i. This is the state of the current node
            node.seq_idx = np.choose(node.up.seq_idx, node.joint_Cx.T)
            # reconstruct seq, etc
            tmp_sequence = np.choose(node.seq_idx, self.gtr.alphabet)
            if self.sequence_reconstruction:
                N_diff += (tmp_sequence!=node.cseq).sum()
            else:
                N_diff += L

            node._cseq = tmp_sequence

        self.logger("TreeAnc._ml_anc_joint: ...done", 3)
        self.reconstructed_tip_sequences = reconstruct_tip_states

        # do clean-up
        if not debug:
            for node in self.tree.find_clades(order='preorder'):
                del node.joint_Lx
                del node.joint_Cx
                if hasattr(node, 'seq_idx'):
                    del node.seq_idx

        self.sequence_reconstruction = 'joint'
        return N_diff


###############################################################
### sequence and mutation storing
###############################################################
    def get_branch_mutation_matrix(self, node, full_sequence=False):
        """uses results from marginal ancestral inference to return a joint
        distribution of the sequence states at both ends of the branch.

        Parameters
        ----------
        node : Phylo.clade
            node of the tree
        full_sequence : bool, optional
            expand the sequence to the full sequence, if false (default)
            the there will be one mutation matrix for each column in the
            compressed alignment

        Returns
        -------
        numpy.array
            an Lxqxq stack of matrices (q=alphabet size, L (compressed)sequence length)
        """
        pp,pc = self.marginal_branch_profile(node)

        # calculate pc_i [e^Qt]_ij pp_j for each site
        expQt = self.gtr.expQt(self._branch_length_to_gtr(node)) + ttconf.SUPERTINY_NUMBER
        if len(expQt.shape)==3: # site specific model
            mut_matrix_stack = np.einsum('ai,aj,ija->aij', pc, pp, expQt)
        else:
            mut_matrix_stack = np.einsum('ai,aj,ij->aij', pc, pp, expQt)

        # normalize this distribution
        normalizer = mut_matrix_stack.sum(axis=2).sum(axis=1)
        mut_matrix_stack = np.einsum('aij,a->aij', mut_matrix_stack, 1.0/normalizer)

        # expand to full sequence if requested
        if full_sequence:
            return mut_matrix_stack[self.data.full_to_compressed_sequence_map]
        else:
            return mut_matrix_stack


    def marginal_branch_profile(self, node):
        '''
        calculate the marginal distribution of sequence states on both ends
        of the branch leading to node,

        Parameters
        ----------
        node : PhyloTree.Clade
           TreeNode, attached to the branch.

        Returns
        -------
        pp, pc : Pair of vectors (profile parent, pp) and (profile child, pc)
           that are of shape (L,n) where L is sequence length and n is alphabet size.
           note that this correspond to the compressed sequences.
        '''
        parent = node.up
        if parent is None:
            raise Exception("Branch profiles can't be calculated for the root!")
        if not hasattr(node, 'marginal_outgroup_LH'):
            raise Exception("marginal ancestral inference needs to be performed first!")

        pc = node.marginal_subtree_LH
        pp = node.marginal_outgroup_LH
        return pp, pc


    def add_branch_state(self, node):
        """add a dictionary to the node containing tuples of state pairs
        and a list of their number across the branch

        Parameters
        ----------
        node : tree.node
            attaces attribute :branch_state:
        """
        seq_pairs, multiplicity = self.gtr.state_pair(
                                       node.up.cseq, node.cseq,
                                       pattern_multiplicity = self.data.multiplicity,
                                       ignore_gaps = self.ignore_gaps)
        node.branch_state = {'pair':seq_pairs, 'multiplicity':multiplicity}


###################################################################
### Branch length optimization
###################################################################
    def optimize_branch_len(self, **kwargs):
        """Deprecated in favor of 'optimize_branch_lengths_joint'"""
        return self.optimize_branch_lengths_joint(**kwargs)

    def optimize_branch_len_joint(self, **kwargs):
        """Deprecated in favor of 'optimize_branch_lengths_joint'"""
        return self.optimize_branch_lengths_joint(**kwargs)

    def optimize_branch_lengths_joint(self, **kwargs):
        """
        Perform optimization for the branch lengths of the entire tree.
        This method only does a single path and needs to be iterated.

        **Note** this method assumes that each node stores information
        about its sequence as numpy.array object (node.sequence attribute).
        Therefore, before calling this method, sequence reconstruction with
        either of the available models must be performed.

        Parameters
        ----------
         **kwargs :
            Keyword arguments

        Keyword Args
        ------------

         store_old : bool
            If True, the old lengths will be saved in :code:`node._old_dist` attribute.
            Useful for testing, and special post-processing.


        """

        self.logger("TreeAnc.optimize_branch_length: running branch length optimization using jointML ancestral sequences",1)
        if (self.tree is None) or (self.data.aln is None):
            raise MissingDataError("TreeAnc.optimize_branch_length: ERROR, alignment or tree are missing.")

        store_old_dist = kwargs['store_old'] if 'store_old' in kwargs else False

        max_bl = 0
        for node in self.tree.find_clades(order='postorder'):
            if node.up is None: continue # this is the root
            if store_old_dist:
                node._old_length = node.branch_length

            new_len = max(0,self.optimal_branch_length(node))

            self.logger("Optimization results: old_len=%.4e, new_len=%.4e"
                        " Updating branch length..."%(node.branch_length, new_len), 5)

            node.branch_length = new_len
            node.mutation_length=new_len
            max_bl = max(max_bl, new_len)

        if max_bl>0.15:
            self.logger("TreeAnc.optimize_branch_lengths_joint: THIS TREE HAS LONG BRANCHES."
                        " \n\t ****TreeTime's JOINT IS NOT DESIGNED TO OPTIMIZE LONG BRANCHES."
                        " \n\t ****PLEASE OPTIMIZE BRANCHES USING: "
                        " \n\t ****branch_length_mode='input' or 'marginal'", 0, warn=True)

        # as branch lengths changed, the distance to root etc need to be recalculated
        self.tree.root.up = None
        self.tree.root.dist2root = 0.0
        self._prepare_nodes()
        return ttconf.SUCCESS


    def optimal_branch_length(self, node):
        '''
        Calculate optimal branch length given the sequences of node and parent

        Parameters
        ----------
        node : PhyloTree.Clade
           TreeNode, attached to the branch.

        Returns
        -------
        new_len : float
           Optimal length of the given branch

        '''
        if node.up is None:
            return self.one_mutation

        if not hasattr(node, 'branch_state'):
            self.add_branch_state(node)
        return self.gtr.optimal_t_compressed(node.branch_state['pair'],
                                    node.branch_state['multiplicity'])


    def optimal_marginal_branch_length(self, node, tol=1e-10):
        '''
        calculate the marginal distribution of sequence states on both ends
        of the branch leading to node,

        Parameters
        ----------
        node : PhyloTree.Clade
           TreeNode, attached to the branch.

        Returns
        -------
        branch_length : float
           branch length of the branch leading to the node.
           note: this can be unstable on iteration
        '''

        if node.up is None:
            return self.one_mutation
        else:
            pp, pc = self.marginal_branch_profile(node)
            return self.gtr.optimal_t_compressed((pp, pc), self.data.multiplicity, profiles=True, tol=tol)


    def optimize_tree_marginal(self, max_iter=10, infer_gtr=False, pc=1.0, damping=0.75,
                               LHtol=0.1, site_specific_gtr=False, **kwargs):
        self.infer_ancestral_sequences(marginal=True)
        oldLH = self.sequence_LH()
        self.logger("TreeAnc.optimize_tree_marginal: initial, LH=%1.2f, total branch_length %1.4f"%
                    (oldLH, self.tree.total_branch_length()), 2)
        for i in range(max_iter):
            if infer_gtr:
                self.infer_gtr(site_specific=site_specific_gtr, marginal=True, normalized_rate=True, pc=pc)
                self.infer_ancestral_sequences(marginal=True)

            old_bl = self.tree.total_branch_length()
            tol = 1e-8 + 0.01**(i+1)
            for n in self.tree.find_clades():
                if n.up is None:
                    continue
                if n.up.up is None and len(n.up.clades)==2:
                    # children of a bifurcating root!
                    n1, n2 = n.up.clades
                    total_bl = n1.branch_length+n2.branch_length
                    bl_ratio = n1.branch_length/total_bl

                    prof_c = n1.marginal_subtree_LH
                    prof_p = normalize_profile(n2.marginal_subtree_LH*self.tree.root.marginal_outgroup_LH)[0]

                    new_bl = self.gtr.optimal_t_compressed((prof_p, prof_c), self.data.multiplicity, profiles=True, tol=tol)
                    update_val = new_bl*(1-damping**(i+1)) +  total_bl*damping**(i+1)

                    n1.branch_length = update_val*bl_ratio
                    n2.branch_length = update_val*(1-bl_ratio)
                else:
                    new_val = self.optimal_marginal_branch_length(n, tol=tol)
                    update_val = new_val*(1-damping**(i+1)) + n.branch_length*damping**(i+1)
                    n.branch_length = update_val

            self.infer_ancestral_sequences(marginal=True)

            LH = self.sequence_LH()
            deltaLH = LH - oldLH
            oldLH = LH
            dbl = self.tree.total_branch_length() - old_bl
            self.logger("TreeAnc.optimize_tree_marginal: iteration %d, LH=%1.2f (%1.2f), delta branch_length=%1.4f, total branch_length %1.4f"%
                        (i, LH, deltaLH, dbl, self.tree.total_branch_length()), 2)
            if deltaLH<LHtol:
                self.logger("TreeAnc.optimize_tree_marginal: deltaLH=%f, stopping iteration."%deltaLH,1)
                break
        return ttconf.SUCCESS


    def optimize_sequences_and_branch_length(self,*args, **kwargs):
        """This method is a shortcut for :py:meth:`treetime.TreeAnc.optimize_tree`
        Deprecated in favor of 'optimize_tree'
        """
        self.logger("Deprecation warning: 'optimize_sequences_and_branch_length' will be removed and replaced by 'optimize_tree'!", 1, warn=True)
        self.optimize_tree(*args,**kwargs)

    def optimize_seq_and_branch_len(self,*args, **kwargs):
        """This method is a shortcut for :py:meth:`treetime.TreeAnc.optimize_tree`
        Deprecated in favor of 'optimize_tree'
        """
        self.logger("Deprecation warning: 'optimize_seq_and_branch_len' will be removed and replaced by 'optimize_tree'!", 1, warn=True)
        self.optimize_tree(*args,**kwargs)

    def optimize_tree(self,prune_short=True,
                      marginal_sequences=False, branch_length_mode='joint',
                      max_iter=5, infer_gtr=False, pc=1.0, **kwargs):
        """
        Iteratively set branch lengths and reconstruct ancestral sequences until
        the values of either former or latter do not change. The algorithm assumes
        knowing only the topology of the tree, and requires that sequences are assigned
        to all leaves of the tree.

        The first step is to pre-reconstruct ancestral
        states using Fitch reconstruction algorithm or ML using existing branch length
        estimates. Then, optimize branch lengths and re-do reconstruction until
        convergence using ML method.

        Parameters
        -----------
         prune_short : bool
            If True, the branches with zero optimal length will be pruned from
            the tree, creating polytomies. The polytomies could be further
            processed using :py:meth:`treetime.TreeTime.resolve_polytomies` from the TreeTime class.

         marginal_sequences : bool
            Assign sequences to their marginally most likely value, rather than
            the values that are jointly most likely across all nodes.

         branch_length_mode : str
            'joint', 'marginal', or 'input'. Branch lengths are left unchanged in case
            of 'input'. 'joint' and 'marginal' cause branch length optimization
            while setting sequences to the ML value or tracing over all possible
            internal sequence states.

         max_iter : int
            Maximal number of times sequence and branch length iteration are optimized

         infer_gtr : bool
            Infer a GTR model from the observed substitutions.

        """
        if branch_length_mode=='marginal':
            self.optimize_tree_marginal(max_iter=max_iter, infer_gtr=infer_gtr, pc=pc, **kwargs)
            if prune_short:
                self.prune_short_branches()
            return ttconf.SUCCESS
        elif branch_length_mode=='input':
            N_diff = self.reconstruct_anc(method='probabilistic', infer_gtr=infer_gtr, pc=pc,
                                          marginal=marginal_sequences, **kwargs)
            if prune_short:
                self.prune_short_branches()
            return ttconf.SUCCESS
        elif branch_length_mode!='joint':
            raise UnknownMethodError("TreeAnc.optimize_tree: `branch_length_mode` should be in ['marginal', 'joint', 'input']")

        self.logger("TreeAnc.optimize_tree: sequences...", 1)
        N_diff = self.reconstruct_anc(method='probabilistic', infer_gtr=infer_gtr, pc=pc,
                                      marginal=marginal_sequences, **kwargs)
        self.optimize_branch_len(verbose=0, store_old=False, mode=branch_length_mode)

        n = 0
        while n<max_iter:
            n += 1
            if prune_short:
                self.prune_short_branches()
            N_diff = self.reconstruct_anc(method='probabilistic', infer_gtr=False,
                                          marginal=marginal_sequences, **kwargs)

            self.logger("TreeAnc.optimize_tree: Iteration %d."
                   " #Nuc changed since prev reconstructions: %d" %(n, N_diff), 2)

            if N_diff < 1:
                break
            self.optimize_branch_lengths_joint(store_old=False)

        self.tree.unconstrained_sequence_LH = (self.tree.sequence_LH*self.data.multiplicity).sum()
        self._prepare_nodes() # fix dist2root and up-links after reconstruction
        self.logger("TreeAnc.optimize_tree: Unconstrained sequence LH:%f" % self.tree.unconstrained_sequence_LH , 2)
        return ttconf.SUCCESS


    def prune_short_branches(self):
        """
        If the branch length is less than the minimal value, remove the branch
        from the tree. **Requires** ancestral sequence reconstruction
        """
        self.logger("TreeAnc.prune_short_branches: pruning short branches (max prob at zero)...", 1)
        for node in self.tree.find_clades():
            if node.up is None or node.is_terminal():
                continue

            # probability of the two seqs separated by zero time is not zero
            if  ((node.branch_length<0.1*self.one_mutation) and
                 (self.gtr.prob_t(node.up._cseq, node._cseq, 0.0,
                                  pattern_multiplicity=self.data.multiplicity) > 0.1)):
                # re-assign the node children directly to its parent
                node.up.clades = [k for k in node.up.clades if k != node] + node.clades
                for clade in node.clades:
                    clade.up = node.up


#####################################################################
## GTR INFERENCE
#####################################################################
    def infer_gtr(self, marginal=False, site_specific=False, normalized_rate=True,
                  fixed_pi=None, pc=5.0, **kwargs):
        """
        Calculates a GTR model given the multiple sequence alignment and the tree.
        It performs ancestral sequence inferrence (joint or marginal), followed by
        the branch lengths optimization. Then, the numbers of mutations are counted
        in the optimal tree and related to the time within the mutation happened.
        From these statistics, the relative state transition probabilities are inferred,
        and the transition matrix is computed.

        The result is used to construct the new GTR model of type 'custom'.
        The model is assigned to the TreeAnc and is used in subsequent analysis.

        Parameters
        -----------

         print_raw : bool
            If True, print the inferred GTR model

         marginal : bool
            If True, use marginal sequence reconstruction

         normalized_rate : bool
            If True, sets the mutation rate prefactor to 1.0.

         fixed_pi : np.array
            Provide the equilibrium character concentrations.
            If None is passed, the concentrations will be inferred from the alignment.

         pc: float
            Number of pseudo counts to use in gtr inference

        Returns
        -------

         gtr : GTR
            The inferred GTR model
        """
        if site_specific and self.data.compress:
            raise TypeError("TreeAnc.infer_gtr(): sequence compression and site specific GTR models are incompatible!" )

        if not self.ok:
            raise MissingDataError("TreeAnc.infer_gtr: ERROR, sequences or tree are missing", 0)

        # if ancestral sequences are not in place, reconstruct them
        if marginal and self.sequence_reconstruction!='marginal':
            self._ml_anc_marginal(**kwargs)
        elif not self.sequence_reconstruction:
            self._ml_anc_joint(**kwargs)

        n = self.gtr.n_states
        L = len(self.tree.root._cseq)
        # matrix of mutations n_{ij}: i = derived state, j=ancestral state
        n_ija = np.zeros((n,n,L))
        T_ia = np.zeros((n,L))

        self.logger("TreeAnc.infer_gtr: counting mutations...", 2)
        for node in self.tree.get_nonterminals():
            for c in node:
                if marginal:
                    mut_stack = np.transpose(self.get_branch_mutation_matrix(c, full_sequence=False), (1,2,0))
                    T_ia += 0.5*self._branch_length_to_gtr(c) * mut_stack.sum(axis=0) * self.data.multiplicity
                    T_ia += 0.5*self._branch_length_to_gtr(c) * mut_stack.sum(axis=1) * self.data.multiplicity
                    n_ija += mut_stack * self.data.multiplicity
                else:
                    for a,pos, d in c.mutations:
                        try:
                            i,j = self.gtr.state_index[d], self.gtr.state_index[a]
                        except:
                            # ambiguous positions
                            continue
                        cpos = self.data.full_to_compressed_sequence_map[pos]
                        n_ija[i,j,cpos]+=1
                        T_ia[j,cpos] += 0.5*self._branch_length_to_gtr(c)
                        T_ia[i,cpos] -= 0.5*self._branch_length_to_gtr(c)

                    for i, nuc in enumerate(self.gtr.alphabet):
                        cseq = c.cseq
                        if cseq is not None:
                            ind = cseq==nuc
                            T_ia[i,ind] += self._branch_length_to_gtr(c)*self.data.multiplicity[ind]

        self.logger("TreeAnc.infer_gtr: counting mutations...done", 3)

        if site_specific:
            if marginal:
                root_state = self.tree.root.marginal_profile.T
            else:
                root_state = seq2prof(self.tree.root.cseq, self.gtr.profile_map).T
            self._gtr = GTR_site_specific.infer(n_ija, T_ia, pc=pc,
                                root_state=root_state, logger=self.logger,
                                alphabet=self.gtr.alphabet, prof_map=self.gtr.profile_map)
        else:
            root_state = np.array([np.sum((self.tree.root.cseq==nuc)*self.data.multiplicity)
                                   for nuc in self.gtr.alphabet])
            n_ij = n_ija.sum(axis=-1)
            self._gtr = GTR.infer(n_ij, T_ia.sum(axis=-1), root_state, fixed_pi=fixed_pi, pc=pc,
                                  alphabet=self.gtr.alphabet, logger=self.logger,
                                  prof_map = self.gtr.profile_map)

        if normalized_rate:
            self.logger("TreeAnc.infer_gtr: setting overall rate to 1.0...", 2)
            if site_specific:
                self._gtr.mu /= self._gtr.average_rate().mean()
            else:
                self._gtr.mu=1.0
        return self._gtr


    def infer_gtr_iterative(self, max_iter=10, site_specific=False, LHtol=0.1,
                            pc=1.0, normalized_rate=False):
        """infer GTR model by iteratively estimating ancestral sequences and the GTR model

        Parameters
        ----------
        max_iter : int, optional
            maximal number of iterations
        site_specific : bool, optional
            use a site specific model
        LHtol : float, optional
            stop iteration when LH improvement falls below this cutoff
        pc : float, optional
            pseudocount to use
        normalized_rate : bool, optional
            set the overall rate to 1 (makes sense when optimizing branch lengths as well)

        Returns
        -------
        str
            success/failure code
        """
        self.infer_ancestral_sequences(marginal=True)
        old_p = np.copy(self.gtr.Pi)
        old_LH = self.sequence_LH()

        for i in range(max_iter):
            self.infer_gtr(site_specific=site_specific, marginal=True,
                           normalized_rate=normalized_rate, pc=pc)
            self.infer_ancestral_sequences(marginal=True)

            dp = np.abs(self.gtr.Pi - old_p).mean() if self.gtr.Pi.shape==old_p.shape else np.nan

            deltaLH = self.sequence_LH() - old_LH

            old_p = np.copy(self.gtr.Pi)
            old_LH = self.sequence_LH()

            self.logger("TreeAnc.infer_gtr_iterative: iteration %d, LH=%1.2f (%1.2f), deltaP=%1.4f"%
                        (i, old_LH, deltaLH, dp), 2)
            if deltaLH<LHtol:
                self.logger("TreeAnc.infer_gtr_iterative: deltaLH=%f, stopping iteration."%deltaLH,1)
                break
        return ttconf.SUCCESS

    def optimize_gtr_rate(self):
        """Estimate the overal rate of the GTR model by optimizing the full
        likelihood of the sequence data.
        """
        from scipy.optimize import minimize_scalar
        def cost_func(sqrt_mu):
            self.gtr.mu = sqrt_mu**2
            self.postorder_traversal_marginal()
            self.total_LH_and_root_sequence(sample_from_profile=False,
                                            assign_sequence=False)
            return -self.sequence_LH()

        old_mu = self.gtr.mu
        try:
            sol = minimize_scalar(cost_func,bracket=[0.01*np.sqrt(old_mu), np.sqrt(old_mu),100*np.sqrt(old_mu)])
        except:
            self.gtr.mu=old_mu
            self.logger('treeanc:optimize_gtr_rate: optimization failed, continuing with previous mu',1,warn=True)
            return

        if sol['success']:
            self.gtr.mu = sol['x']**2
            self.logger('treeanc:optimize_gtr_rate: optimization successful. Overall rate estimated to be %f'%self.gtr.mu,1)
        else:
            self.gtr.mu=old_mu
            self.logger('treeanc:optimize_gtr_rate: optimization failed, continuing with previous mu',1,warn=True)


###############################################################################
### Utility functions
###############################################################################
    def get_reconstructed_alignment(self, reconstruct_tip_states=False):
        """
        Get the multiple sequence alignment, including reconstructed sequences for
        the internal nodes.

        Parameters
        ----------
        reconstructed_tip_sequences : bool, optional
            return reconstructed sequences of terminal nodes. If these have not
            been reconstructed yet, this will trigger a rerun of `infer_ancestral_sequences`

        Returns
        -------
        new_aln : MultipleSeqAlignment
           Alignment including sequences of all internal nodes

        """
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        self.logger("TreeAnc.get_reconstructed_alignment ...",2)
        if (not self.sequence_reconstruction) or (reconstruct_tip_states != self.reconstructed_tip_sequences):
            self.logger("TreeAnc.reconstructed_alignment... reconstruction not yet done",3)
            self.infer_ancestral_sequences(reconstruct_tip_states=reconstruct_tip_states)

        if self.data.is_sparse:
            new_aln = {'sequences': {n.name: self.data.compressed_to_sparse_sequence(n.cseq)
                       for n in self.tree.find_clades()}}
            new_aln['reference'] = self.data.ref
            new_aln['positions'] = self.data.nonref_positions
            new_aln['inferred_const_sites'] = self.data.inferred_const_sites
        else:
            new_aln = MultipleSeqAlignment([SeqRecord(id=n.name,
                                              seq=Seq(self.sequence(n, reconstructed=reconstruct_tip_states,
                                                      as_string=True, compressed=False)), description="")
                                        for n in self.tree.find_clades()])

        return new_aln


    def sequence(self, node, reconstructed=False, as_string=True, compressed=False):
        """return the sequence of a node.

        Parameters
        ----------
        node : Phylo.node, str
            node in tree
        reconstructed : bool, optional
            return the reconstructed sequence also for terminal nodes. this will replace
            ambiguous sites with the most likely sequence state.
        as_string : bool, optional
            return the sequence as character array rather than contiguous string
        compressed : bool, optional
            return the a sequence where unique alignment patterns are reduced to
            one alignment column each

        Returns
        -------
        str or np.array
            sequence of node
        """
        if type(node)==str:
            if node in self.leaves_lookup:
                nodes = self.leaves_lookup
            else:
                raise ValueError("TreeAnc.sequence accepts strings are argument only when the node is terminal and present in the leave lookup table")

        if reconstructed and not self.reconstructed_tip_sequences:
            raise ValueError("TreeAnc.sequence can only return reconstructed terminal nodes if TreeAnc.infer_ancestral_sequences was run with this the flag `reconstruct_tip_states`.")

        if compressed:
            if (not reconstructed) and (node.name in self.data.compressed_alignment):
                tmp_seq = self.data.compressed_alignment[node.name]
            else:
                tmp_seq = node.cseq
        else:
            if (not reconstructed) and (node.name in self.data.aln):
                tmp_seq = self.data.aln[node.name]
            elif node.cseq is not None:
                tmp_seq = self.data.compressed_to_full_sequence(node.cseq, as_string=False)
            else:
                tmp_seq = np.array([self.gtr.ambiguous or 'N']*self.data.compressed_length)

        return "".join(tmp_seq) if as_string else np.copy(tmp_seq)


    def get_tree_dict(self, keep_var_ambigs=False):
        return self.get_reconstructed_alignment()


    def recover_var_ambigs(self):
        self.logger("TreeAnc: recover_var_ambigs: calls to recover_var_ambigs are no longer necessary since tip states are not inferred unless explicitly specified using `reconstruct_tip_states=True`.", 0, warn=True)
        if self.reconstructed_tip_sequences:
            self.logger("Your code reconstructed tip states, please change the call of ancestral inference in your code",0, warn=True)
        else:
            self.logger("Your analysis did not reconstructed tip states, you can remove the call of `recover_var_ambigs`",0, warn=True)
