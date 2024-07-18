from Bio import AlignIO
from typing import Dict
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node, SparseSeqNode, SparseSeqInfo, SparseSeqDis, SparseSeqEdge
from lib import SeqPartition, SeqInfoParsimony, find_char_ranges, VarPos, RangeCollection_intersection, RangeCollection, Mut
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map
from itertools import chain
import time

NON_CHAR = '.'
VARIABLE = '~'
FILL_CHAR = ' '

def seq_info_from_array(seq: np.array, profile: Dict[str,np.array], alphabet_gapN: str) -> SparseSeqInfo:
    seq_dis = SparseSeqDis(variable={pos:VarPos(profile(nuc), None) for pos, nuc in enumerate(seq) if nuc not in alphabet_gapN})
    # code the missing regions and gap regions along with the ambiguous nucleotides
    unknown = find_char_ranges(seq, 'N')
    gaps = find_char_ranges(seq, '-')
    return(SparseSeqInfo(unknown=unknown, gaps=gaps,
                        non_char=RangeCollection(unknown.ranges + gaps.ranges),
                        distribution=seq_dis, sequence=seq))

def attach_seqs_to_graph(graph, aln_list, gtr_list):
  name_to_key = {}
  for n in graph.get_leaves():
    n.payload().sparse_sequences = []
    name_to_key[n.payload().name] = n.key()

  graph.partitions = []
  for aln, gtr in zip(aln_list, gtr_list):
    L = len(aln[list(aln.keys())[0]]) # stupid way to get alignment length
    partition = SeqPartition(gtr=gtr, length=L, profile_map=profile_map)
    graph.partitions.append(partition)
    alphabet_gapN = ''.join(gtr.alphabet)+'-N'
    profile = partition.profile

    for leaf in graph.get_leaves():
      leaf_name = leaf.payload().name
      leaf.payload().sparse_sequences.append(seq_info_from_array(np.fromiter(str(aln[leaf_name]).upper(), 'U1'), profile, alphabet_gapN))

def init_sparse_sequences(graph, aln_list, gtr_list):
  t0 = time.time()
  print(f"elapsed: {time.time()-t0:.3f}")
  attach_seqs_to_graph(graph, aln_list, gtr_list)
  print(f"elapsed: {time.time()-t0:.3f}")
  fitch_backwards(graph)
  print(f"elapsed: {time.time()-t0:.3f}")
  fitch_forward(graph)
  print(f"elapsed: {time.time()-t0:.3f}")


def fitch_backwards(graph: Graph):
  n_partitions = len(graph.partitions)
  def fitch_backwards_node(node: GraphNodeBackward):
    if node.is_leaf: # leaves have been dealt with in init
      return

    # init internal node and edges
    node.payload.sparse_sequences = []
    for c,e in node.children:
      # not clear where the transmission info on the edge is going to come from
      e.sparse_sequences = [SparseSeqEdge() for _ in range(n_partitions)]

    for si in range(n_partitions):
      # short hands
      child_seqs = [c.sparse_sequences[si] for c,e in node.children]
      child_edges = [e.sparse_sequences[si] for c,e in node.children]
      gtr = graph.partitions[si].gtr
      L = graph.partitions[si].length
      profile = graph.partitions[si].profile

      # init the local representation with gaps, unknowns
      # need to account for parts of the sequence transmitted along edges
      seq_dis = SparseSeqDis(variable={})
      seq_info = SparseSeqInfo(gaps=RangeCollection_intersection([cseq.gaps for cseq in child_seqs]),
                          unknown=RangeCollection_intersection([cseq.unknown for cseq in child_seqs]),
                          non_char=RangeCollection_intersection([cseq.non_char for cseq in child_seqs]))
      # define dummy sequence
      full_seq = np.array([FILL_CHAR]*L)
      # fill indeterminate positions
      for r in seq_info.non_char:
        for pos in range(r.start, r.end):
          full_seq[pos]=NON_CHAR

      # process all positions where the children are variable (there are probably better ways to do this)
      # need to account for parts of the sequence transmitted along edges
      variable_pos = sorted(set(chain.from_iterable([cseq.distribution.variable.keys() for cseq in child_seqs])))
      for pos in variable_pos:
        # 2D float array
        # stack all ingroup profiles of children, use the determinate profile if position is not variable.
        child_profiles = []
        for cseq,e in zip(child_seqs, child_edges):
          if e.transmission and (not e.transmission.contains(pos)):
            # transmission field is not currently used
            continue
          if cseq.non_char.contains(pos):
            # this position does not have character state information
            continue
          if pos in cseq.distribution.variable:
            child_profiles.append(cseq.distribution.variable[pos].profile)
          else:
            child_profiles.append(profile(cseq.sequence[pos]))

        isect = np.prod(child_profiles, axis=0)
        # if we save the states of the children for each position that is variable in the node, we would not need the full_seq on the forward pass
        if isect.sum()==1:
          full_seq[pos]=gtr.alphabet[np.argmax(isect)]
        elif isect.sum()>1:
          seq_dis.variable[pos]=VarPos(isect, None)
          full_seq[pos]=VARIABLE
        else:
          child_union = np.sum(child_profiles, axis=0)
          seq_dis.variable[pos]=VarPos(np.array([1.0 if x>0 else 0.0 for x in child_union]), None)
          full_seq[pos]=VARIABLE

      # process all positions where the children are fixed or completely unknown in some children
      # this is ridiculously slow, even though we are only comparing nucleotides in the children
      # it should be possible to make this much faster.
      for pos, (nuc, child_states) in enumerate(zip(full_seq, zip(*[cseq.sequence for cseq in child_seqs]))):
        if nuc!=FILL_CHAR: # already touched this position
          continue
        determined_states = set([x for x in child_states if x in gtr.alphabet])
        if len(determined_states)==1:
          full_seq[pos]=determined_states.pop()
        elif len(determined_states)>1:
          # if we save the states of the children for each position that is variable in the node, we would not need the full_seq on the forward pass
          seq_dis.variable[pos] = VarPos(np.sum([profile(x) for x in determined_states], axis=0), None)
          full_seq[pos]=VARIABLE
        else:
          import ipdb; ipdb.set_trace()

      seq_info.distribution=seq_dis
      seq_info.sequence=full_seq
      node.payload.sparse_sequences.append(seq_info)

  graph.par_iter_backward(fitch_backwards_node)

def add_mutation(pos, pnuc, cnuc, composition):
  composition[pnuc] -= 1
  composition[cnuc] += 1
  return Mut(pnuc, pos, cnuc)

def fitch_forward(graph: Graph):
  n_partitions = len(graph.partitions)

  def fitch_forward_node(node: GraphNodeForward):
    for si in range(n_partitions):
      gtr = graph.partitions[si].gtr
      profile = graph.partitions[si].profile
      seq_info = node.payload.sparse_sequences[si]

      if node.is_root:
        for pos, p in seq_info.distribution.variable.items():
            p.state = gtr.alphabet[np.argmax(p.profile)]
            seq_info.sequence[pos] = p.state

        seq_info.composition = {s:0 for s in gtr.alphabet}
        for s in seq_info.sequence:
          if s in gtr.alphabet: #there could be positions that are gap or N everywhere, should be over complement of `non_char`
            seq_info.composition[s] += 1
      else:
        # short hands
        pseq = node.parents[0][0].sparse_sequences[si]
        pedge = node.parents[0][1].sparse_sequences[si]
        seq_info = node.payload.sparse_sequences[si]
        pedge.muts = []

        # copy parent seqinfo
        seq_info.composition = {s:k for s,k in pseq.composition.items()}

        # fill in the indeterminate positions.
        for r in seq_info.non_char:
          for pos in range(r.start, r.end):
            seq_info.sequence[pos]=pseq.sequence[pos]

        # for each variable position, pick a state or a mutation
        for pos, p in seq_info.distribution.variable.items():
          pnuc=pseq.sequence[pos]
          # check whether parent is in child profile (sum>0 --> parent state is in profile)
          if np.sum(profile(pnuc)*p.profile):
            seq_info.sequence[pos] = pnuc
          else:
            cnuc = gtr.alphabet[np.argmax(p.profile)]
            seq_info.sequence[pos] = cnuc
            pedge.muts.append(add_mutation(pos, pnuc, cnuc, seq_info.composition))

          p.state = seq_info.sequence[pos]

        for pos, pvar in pseq.distribution.variable.items():
          if pos in seq_info.distribution.variable: continue
          # NOTE: access to full_seq would not be necessary if we had saved the
          # child state of variable positions in the backward pass
          node_nuc = seq_info.sequence[pos]
          if pvar.state!=node_nuc:
            pedge.muts.append(add_mutation(pos, pvar.state, node_nuc, seq_info.composition))

  def clean_up_node(node):
    for seq_info in node.payload.sparse_sequences:
      if not node.is_leaf: #delete the variable position everywhere instead of leaves
        seq_info.distribution.variable = {}
      for r in seq_info.non_char: # remove the undetermined counts from the counts of fixed positions
        for pos in range(r.start, r.end):
          seq_info.composition[seq_info.sequence[pos]] -= 1
      seq_info.distribution.fixed_counts = {k:v for k,v in seq_info.composition.items()}

      for pos, p in seq_info.distribution.variable.items():
        seq_info.distribution.fixed_counts[p.state] -= 1

      if not node.is_root:
        seq_info.sequence = []

  graph.par_iter_forward(fitch_forward_node)
  graph.par_iter_forward(clean_up_node)


def reconstruct_sequence(graph: Graph, node: Node):
  '''
  reconstruct a specific sequence using the full sequence at the root and implanting mutations.
  Note that this works only for trees with one root.
  '''
  path_to_root = [node]
  p = node
  while not p.is_root():
    assert len(p._inbound)==1
    p = graph.get_node(graph.get_edge(p._inbound[0]).source())
    path_to_root.append(p)

  seqs = []
  for si in range(len(graph.partitions)):
    seq = np.copy(path_to_root[-1].payload().sparse_sequences[si].sequence)
    # walk back from root to node
    for n in path_to_root[::-1][1:]:
      for m in graph.get_edge(n._inbound[0]).payload().sparse_sequences[si].muts:
        seq[m.pos] = m.qry

    seq_info = node.payload().sparse_sequences[si]
    for r in seq_info.gaps:
      for pos in range(r.start, r.end):
        seq[pos]='-'
    for r in seq_info.unknown:
      for pos in range(r.start, r.end):
        seq[pos]=graph.partitions[si].gtr.ambiguous
    for pos, p in seq_info.distribution.variable.items():
      seq[pos]=graph.partitions[si].code(tuple(p.profile))
    seqs.append(seq)

  return seqs

def tests():
  aln = {"root":"ACAGCCATGTATTG",
         "AB":"ACATCCCTGTA-TG",
         "A":"ACATCGCCNNA--G",
         "B":"GCATCCCTGTA-NG",
         "CD":"CCGGCCATGTATTG",
         "C":"CCGGCGATGTRTTG",
         "D":"TCGGCCGTGTRTTG"}

  tree = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;"
  profile = lambda x: profile_map[x]

  # initializion of leaves
  seq_info = seq_info_from_array(np.fromiter(aln['B'], 'U1'), profile, 'ACGT-N')
  assert str(seq_info.gaps.ranges) == "[Range(start=11, end=12)]"
  assert str(seq_info.unknown.ranges) == "[Range(start=12, end=13)]"

  seq_info = seq_info_from_array(np.fromiter(aln['D'], 'U1'), profile, 'ACGT-N')
  assert str(seq_info.gaps.ranges) == "[]"
  assert str(seq_info.distribution.variable) == "{10: VarPos(profile=array([1., 0., 1., 0.]), state=None)}"

  G = graph_from_nwk_str(nwk_string=tree, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  gtr = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
  attach_seqs_to_graph(G, [aln], [gtr])
  fitch_backwards(G)
  seq_info = G.get_one_root().payload().sparse_sequences[0]
  assert tuple(sorted(seq_info.distribution.variable.keys())) == (0, 2, 3, 5, 6)
  assert "".join(seq_info.sequence) == '~C~~C~~TGTATTG'

  fitch_forward(G)
  seq_info = G.get_one_root().payload().sparse_sequences[0]
  # note that this is contingent on picking the first of ACGT when there are multiple options
  assert "".join(seq_info.sequence) == aln['root']
  for n in G.get_nodes():
    rec_seq = ''.join(reconstruct_sequence(G, n)[0])
    assert rec_seq==aln[n.payload().name]

  muts_expected = [
        "['C6G', 'T8C']",
        "['A1G']",
        "['G4T', 'A7C']",
        "['C6G']",
        "['C1T', 'A7G']",
        "['A1C', 'A3G']",
        ]
  for mexp,e in zip(muts_expected,G.get_edges()):
    muts = str(sorted(e.payload().sparse_sequences[0].muts))
    assert muts==mexp

if __name__=="__main__":
  fname_nwk = 'data/ebola/ebola.nwk'
  fname_seq = 'data/ebola/ebola_dna.fasta'
  fname_nwk = 'test_scripts/data/tree.nwk'
  fname_seq = 'test_scripts/data/sequences.fasta'
  with open(fname_nwk) as fh:
    nwkstr = fh.read()
  G = graph_from_nwk_str(nwk_string=nwkstr, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)

  aln = {seq.id: str(seq.seq).upper() for seq in AlignIO.read(fname_seq, 'fasta')}
  gtr = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
  init_sparse_sequences(G, [aln], [gtr])

  failed = []
  t0=time.time()
  # check output
  for leaf in G.leaves:
    node = G.get_node(leaf)
    nname = node.payload().name
    rec_seq = ''.join(reconstruct_sequence(G, node)[0])
    if rec_seq!=aln[nname]:
      failed.append((nname,rec_seq, aln[nname]))

  if failed:
    print(failed)
  else:
    print("all reconstructed successfully")
  print(f"reconstructed in {time.time()-t0:.3f}s")

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

