from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node
from lib import SeqPartition, SeqInfoParsimony, find_char_ranges, VarPos, RangeCollection_intersection, RangeCollection, Mut
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map

NON_CHAR = '.'
VARIABLE = '~'
FILL_CHAR = ' '

def fitch_backwards(graph: Graph):
  n_seq_partitions = len(graph.partitions)
  alphabets = [''.join(p.gtr.alphabet) for p in graph.partitions]
  alphabets_gapN = [a+'-N' for a in alphabets]
  def fitch_backwards_node(node: GraphNodeBackward):
    node.payload.fitch = []
    if not node.is_leaf: # leaf nodes have a full sequence attached
      node.payload.full_seq = []

    for si in range(n_seq_partitions):
      if node.is_leaf: # process sequences on leaves
        seq = node.payload.full_seq[si]
        # code the missing regions and gap regions along with the ambiguous nucleotides
        seq_info = SeqInfoParsimony(unknown=find_char_ranges(seq, 'N'), gaps=find_char_ranges(seq, '-'),
                           non_char=RangeCollection(find_char_ranges(seq, 'N').ranges + find_char_ranges(seq, '-').ranges),
                           variable={pos:VarPos(graph.partitions[si].profile(nuc), None)
                                    for pos, nuc in enumerate(seq) if nuc not in alphabets_gapN[si]})
        node.payload.fitch.append(seq_info)
      else: # process internal nodes
        # init the local representation with gaps, unknowns
        # need to account for parts of the sequence transmitted along edges
        seq_info = SeqInfoParsimony(gaps=RangeCollection_intersection([c.fitch[si].gaps for c,e in node.children]),
                           unknown=RangeCollection_intersection([c.fitch[si].unknown for c,e in node.children]),
                           non_char=RangeCollection_intersection([c.fitch[si].non_char for c,e in node.children]),
                           variable={})
        # define dummy sequence
        full_seq = np.array([FILL_CHAR]*graph.partitions[0].length)
        # fill indeterminate positions
        for r in seq_info.non_char:
          for pos in range(r.start, r.end):
            full_seq[pos]=NON_CHAR

        # process all positions where the children are variable (there are probably better ways to do this)
        # need to account for parts of the sequence transmitted along edges
        variable_pos = np.unique(np.concatenate([np.array(list(c.fitch[si].variable.keys()), dtype=int) for c,e in node.children]))
        for pos in variable_pos:
          # 2D float array
          # stack all ingroup profiles of children, use the determinate profile if position is not variable.
          child_profiles = []
          for c,e in node.children:
            if e.transmission and (not e.transmission.contains(pos)):
              # transmission field is not currently used
              continue
            if c.fitch[si].non_char.contains(pos):
              # this position does not have character state information
              continue
            if pos in c.fitch[si].variable:
              child_profiles.append(c.fitch[si].variable[pos].profile)
            else:
              child_profiles.append(graph.partitions[si].profile(c.full_seq[si][pos]))

          isect = np.prod(child_profiles, axis=0)
          # if we save the states of the children for each position that is variable in the node, we would not need the full_seq on the forward pass
          if isect.sum()==1:
            full_seq[pos]=alphabets[si][np.argmax(isect)]
          elif isect.sum()>1:
            seq_info.variable[pos]=VarPos(isect, None)
            full_seq[pos]=VARIABLE
          else:
            child_union = np.sum(child_profiles, axis=0)
            seq_info.variable[pos]=VarPos(np.array([1.0 if x>0 else 0.0 for x in child_union]), None)
            full_seq[pos]=VARIABLE

        # process all positions where the children are fixed or completely unknown in some children
        # this is ridiculously slow, even though we are only comparing nucleotides in the children
        # it should be possible to make this much faster.
        for pos, (nuc, child_states) in enumerate(zip(full_seq, zip(*[c.full_seq[si] for c,e in node.children]))):
          if nuc!=FILL_CHAR: # already touched this position
            continue
          determined_states = set([x for x in child_states if x in alphabets[si]])
          if len(determined_states)==1:
            full_seq[pos]=determined_states.pop()
          elif len(determined_states)>1:
            # if we save the states of the children for each position that is variable in the node, we would not need the full_seq on the forward pass
            seq_info.variable[pos] = VarPos(np.sum([graph.partitions[si].profile(x) for x in determined_states], axis=0), None)
            full_seq[pos]=VARIABLE
          else:
            import ipdb; ipdb.set_trace()

        node.payload.fitch.append(seq_info)
        node.payload.full_seq.append(full_seq)

  graph.par_iter_backward(fitch_backwards_node)

def add_mutation(pos, pnuc, cnuc, composition):
  composition[pnuc] -= 1
  composition[cnuc] += 1
  return Mut(pnuc, pos, cnuc)

def fitch_forward(graph: Graph):
  n_seq_partitions = len(graph.partitions)
  alphabets = [''.join(p.gtr.alphabet) for p in graph.partitions]

  def fitch_forward_node(node: GraphNodeForward):
    if node.is_root:
      for si, (seq_info, full_seq) in enumerate(zip(node.payload.fitch, node.payload.full_seq)):
        for pos, p in seq_info.variable.items():
          p.state = alphabets[si][np.argmax(p.profile)]
          full_seq[pos] = p.state
        seq_info.fixed_composition = {s:0 for s in alphabets[si]}
        for s in full_seq:
          if s in alphabets[si]: #there could be positions that are gap or N everywhere, should be over complement of `non_char`
            seq_info.fixed_composition[s] += 1
    else:
      # only deal with one parent for now
      parent, edge = node.parents[0]
      edge.muts = []
      # loop over sequence partitions
      for si, (seq_info, full_seq) in enumerate(zip(node.payload.fitch, node.payload.full_seq)):
        pseq=parent.full_seq[si]
        seq_info.fixed_composition = {s:k for s,k in parent.fitch[si].fixed_composition.items()}

        # fill in the indeterminate positions.
        for r in seq_info.non_char:
          for pos in range(r.start, r.end):
            full_seq[pos]=pseq[pos]

        muts = []
        # for each variable position, pick a state or a mutation
        for pos, p in seq_info.variable.items():
          pnuc=pseq[pos]
          # check whether parent is in child profile (sum>0 --> parent state is in profile)
          if np.sum(graph.partitions[si].profile(pnuc)*p.profile):
            full_seq[pos] = pnuc
          else:
            cnuc = alphabets[si][np.argmax(p.profile)]
            full_seq[pos] = cnuc
            muts.append(add_mutation(pos, pnuc, cnuc, seq_info.fixed_composition))
        for pos, pvar in parent.fitch[si].variable.items():
          if pos in seq_info.variable: continue
          # NOTE: access to full_seq would not be necessary if we had saved the
          # child state of variable positions in the backward pass
          if pvar.state!=full_seq[pos]:
            muts.append(add_mutation(pos, pvar.state, full_seq[pos], seq_info.fixed_composition))
        for pos in seq_info.variable:
          # saving the reference state at variable positions for probabilistic inference
          seq_info.variable[pos].state = full_seq[pos]
        edge.muts.append(muts)
      # print(edge.muts)


  def clean_up_node(node):
    for full_seq, seq_info in zip(node.payload.full_seq, node.payload.fitch):
      if not node.is_leaf: #delete the variable position everywhere instead of leaves
        seq_info.variable = {}
      for r in seq_info.non_char: # remove the undetermined counts from the counts of fixed positions
        for pos in range(r.start, r.end):
          seq_info.fixed_composition[full_seq[pos]] -= 1
      for pos, p in seq_info.variable.items():
        seq_info.fixed_composition[p.state] -= 1

    if not node.is_root:
      node.payload.full_seq = [[] for n in range(n_seq_partitions)]

  graph.par_iter_forward(fitch_forward_node)
  graph.par_iter_forward(clean_up_node)


def reconstruct_sequence(graph: Graph, node: Node):
  '''
  reconstruct a specific sequence using the full sequence at the root and implanting mutations
  '''
  path_to_root = [node]
  p = node
  while not p.is_root():
    assert len(p._inbound)==1
    p = graph.get_node(graph.get_edge(p._inbound[0]).source())
    path_to_root.append(p)

  seqs = [np.copy(s) for s in path_to_root[-1].payload().full_seq]
  # walk back from root to node
  for n in path_to_root[::-1][1:]:
    muts=graph.get_edge(n._inbound[0]).payload().muts
    for s, mset  in zip(seqs, muts):
      for m in mset:
        s[m.pos] = m.qry
  for partition, s, seq_rep in zip(graph.partitions, seqs, node.payload().fitch):
    for r in seq_rep.gaps:
      for pos in range(r.start, r.end):
        s[pos]='-'
    for r in seq_rep.unknown:
      for pos in range(r.start, r.end):
        s[pos]=partition.gtr.ambiguous
    for pos, p in seq_rep.variable.items():
      s[pos]=partition.code(tuple(p.profile))
  return seqs


if __name__=="__main__":
  fname_nwk = 'data/ebola/ebola.nwk'
  fname_seq = 'data/ebola/ebola_dna.fasta'
  fname_nwk = 'test_scripts/data/tree.nwk'
  fname_seq = 'test_scripts/data/sequences.fasta'
  with open(fname_nwk) as fh:
    nwkstr = fh.read()
  G = graph_from_nwk_str(nwk_string=nwkstr, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  name_to_key={n.payload().name: n.key() for n in G.get_leaves()}
  aln = AlignIO.read(fname_seq, 'fasta')
  aln_dict = {s.id:np.fromiter(str(s.seq).upper(), 'U1') for s in aln}

  gtr = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
  G.partitions = [SeqPartition(gtr=gtr, length=aln.get_alignment_length(), profile_map=profile_map)]


  for sname, seq in aln_dict.items():
    G.get_node(name_to_key[sname]).payload().full_seq = [np.copy(seq)]

  fitch_backwards(G)
  fitch_forward(G)

  # check output
  for leaf in G.leaves:
    node = G.get_node(leaf)
    nname = node.payload().name
    rec_seq = reconstruct_sequence(G, node)[0]
    print(nname, np.sum(rec_seq!=aln_dict[nname]))


'''
# NOTES
## Gaps and insertions
I have preliminarily included an InDel struct on the edges. This should allow us to treat indels
properly in ancestral reconstructions and other computations (currently, they are treated like missing info).
It should then also allow to to remove the gap and unknown fields from the nodes since these could be
accumulated from indels along branches like mutations (though it might be useful to keep them). Properly
accounting for indels would be a significant advance since few tree-builders do (if any). The arithmetic would
require a lot of 'interval/range collection' operations.



'''
# how to avoid carrying along the sequence

# - pre-compute nucleotide composition
# - positions become variable only if there is a mutation somewhere
# - save the parsimony state along with the variable profiles
# - when collecting variable positions from children, infer the parent parsimony state from mutations

