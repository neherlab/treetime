from dataclasses import dataclass
from Bio import AlignIO
import numpy as np
from typing import Optional, List, Dict
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward
from lib import SeqPartition, SeqInfo, find_char_ranges, VarPos, RangeCollection_intersection, RangeCollection
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map, default_profile

NON_CHAR = '.'
VARIABLE = '~'
FILL_CHAR = ' '

def fitch_backwards(graph: Graph):
  n_seq_partitions = len(graph.partitions)

  def fitch_backwards_node(node: GraphNodeBackward):
    print(node.payload.name)
    if node.is_leaf: # process sequences on leaves
      node.payload.seq_info_ingroup = []
      for si,seq in enumerate(node.payload.full_seq):
        # code the missing regions and gap regions along with the ambiguous nucleotides
        seq_info = SeqInfo(unknown=find_char_ranges(seq, 'N'), gaps=find_char_ranges(seq, '-'),
                           non_char=RangeCollection(find_char_ranges(seq, 'N').ranges + find_char_ranges(seq, '-').ranges),
                           variable={pos:VarPos(graph.partitions[si].profile(nuc), None)
                                    for pos, nuc in enumerate(seq) if nuc not in 'ACTG-N'})
        node.payload.seq_info_ingroup.append(seq_info)
    else:
      # process internal nodes, again for each sequence partition
      node.payload.seq_info_ingroup = []
      node.payload.full_seq = []
      for si in range(n_seq_partitions):
        # init the local representation with gaps, unknowns
        # need to account for parts of the sequence transmitted along edges
        seq_info = SeqInfo(gaps=RangeCollection_intersection([c.seq_info_ingroup[si].gaps for c,e in node.children]),
                           unknown=RangeCollection_intersection([c.seq_info_ingroup[si].unknown for c,e in node.children]),
                           non_char=RangeCollection_intersection([c.seq_info_ingroup[si].non_char for c,e in node.children]),
                           variable={})
        # define dummy sequence
        full_seq = np.array([FILL_CHAR]*graph.partitions[0].length)
        # fill indeterminate positions
        for r in seq_info.non_char:
          for pos in range(r.start, r.end):
            full_seq[pos]=NON_CHAR

        # process all positions where the children are variable (there are probably better ways to )
        # need to account for parts of the sequence transmitted along edges
        variable_pos = np.unique(np.concatenate([np.array(list(c.seq_info_ingroup[si].variable.keys()), dtype=int) for c,e in node.children]))
        for pos in variable_pos:
          # 2D float array
          # stack all ingroup profiles of children, use the determinate profile if position is not variable.
          child_profiles = []
          for c,e in node.children:
            if e.transmission and (not e.transmission.contains(pos)):
              # transmission field is not currently used
              continue
            if pos in c.seq_info_ingroup[si].variable:
              child_profiles.append(c.seq_info_ingroup[si].variable[pos].profile)
            else:
              child_profiles.append(graph.partitions[si].profile(c.full_seq[si][pos]))

          isect = np.prod(child_profiles, axis=0)
          if isect.sum()==1:
            full_seq[pos]=graph.partitions[si].gtr.alphabet[np.argmax(isect)]
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

          determined_states = set([x for x in child_states if x in graph.partitions[si].gtr.alphabet])
          if len(determined_states)==1:
            full_seq[pos]=determined_states.pop()
          elif len(determined_states)>1:
            seq_info.variable[pos] = VarPos(np.sum([graph.partitions[si].profile(x) for x in determined_states], axis=0), None)
            full_seq[pos]=VARIABLE
          else:
            import ipdb; ipdb.set_trace()

        node.payload.seq_info_ingroup.append(seq_info)
        node.payload.full_seq.append(full_seq)

  graph.par_iter_backward(fitch_backwards_node)


if __name__=="__main__":
#   fname_nwk = 'data/ebola/ebola.nwk'
#   fname_seq = 'data/ebola/ebola_dna.fasta'
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

  # generate_sparse_sequence_representation(G)

  # # for n in G.nodes:
  # #   print(n.payload().name, n.payload().seq[0].state_composition)
  # ingroup_profiles(G, [gtr])
  # print(calculate_root_state(G, [gtr]))

  # # check output
  # for leaf in G.leaves:
  #   node = G.get_node(leaf)
  #   nname = node.payload().name
  #   rec_seq = reconstruct_sequence(G, node)[0]
  #   print(nname, np.sum(rec_seq!=aln[nname]))


# how to avoid carrying along the sequence

# - pre-compute nucleotide composition
# - positions become variable only if there is a mutation somewhere
# - save the parsimony state along with the variable profiles
# - when collecting variable positions from children, infer the parent parsimony state from mutations

