from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node
from lib import SeqPartition, SeqInfoLh, VarPos, RangeCollection_intersection, RangeCollection, Mut
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map
from fitch import fitch_backwards, fitch_forward

def get_variable_states(node, seq_index):
  variable_pos = {}
  child_states = [{} for c, e in node.children]
  for ci, (c,e) in enumerate(node.children):
    # go over all mutations and get reference state
    for m in e.muts[seq_index]:
      variable_pos[m.pos] = m.ref  # this might be set multiple times, but the reference state should always be the same
      child_states[ci][m.pos] = m.qry

    # go over child variable position and get reference state
    for pos, p in c.seq_info_ingroup[seq_index].variable.items():
      if pos not in variable_pos:
        variable_pos[pos] = p.state

  return variable_pos, child_states

def ingroup_profiles(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.partitions]
  gtrs = [p.gtr for p in graph.partitions]

  eps=1e-6
  def calculate_ingroup_node(node: GraphNodeBackward) -> None:
    node.payload.seq_info_ingroup = []
    for si,fitch_seq_info in enumerate(node.payload.fitch):
      if node.is_leaf:
        # this is mostly a copy (or ref here) of the fitch state.
        node.payload.seq_info_ingroup.append(SeqInfoLh(unknown=fitch_seq_info.unknown, gaps=fitch_seq_info.gaps, non_char=fitch_seq_info.non_char,
                                             fixed_composition=fitch_seq_info.fixed_composition, variable=fitch_seq_info.variable,
                                             fixed={state:graph.partitions[si].profile(state) for state in alphabets[si]}))
      else: # internal nodes
        # get all variable positions, the reference state, and the child states at these positions
        variable_pos, child_states = get_variable_states(node, si)
        seq_info = SeqInfoLh(unknown=fitch_seq_info.unknown, gaps=fitch_seq_info.gaps, non_char=fitch_seq_info.non_char,
                             fixed_composition={k:v for k,v in fitch_seq_info.fixed_composition.items()}, variable={},
                             fixed={state:graph.partitions[si].profile(state) for state in alphabets[si]})
        expQT = [gtrs[si].expQt(e.branch_length or 0.0).T for c,e in node.children]

        for c,e in node.children:
          seq_info.logLH += c.seq_info_ingroup[si].logLH

        # go over all putatively variable positions
        for pos, state in variable_pos.items():
          # collect the profiles of children to multiply
          child_profiles = []
          for ci,(c,e) in enumerate(node.children):
            if e.transmission and (not e.transmission.contains(pos)):
              # transmission field is not currently used
              continue
            # need to check whether position is determined
            if c.seq_info_ingroup[si].non_char.contains(pos):
              # this position does not have character state information
              continue
            if pos in c.seq_info_ingroup[si].variable:
              v = c.seq_info_ingroup[si].variable[pos].profile
            elif pos in child_states[ci]:
              v = c.seq_info_ingroup[si].fixed[child_states[ci][pos]]
            else:
              v = c.seq_info_ingroup[si].fixed[state]

            child_profiles.append(expQT[ci].dot(v))

            # calculate new profile and likelihood contribution
            vec = np.prod(child_profiles, axis=0)
            vec_norm = vec.sum()
            seq_info.logLH  += np.log(vec_norm)

          # add position to variable states if the subleading states have a probability exceeding eps
          if vec.max()<(1-eps)*vec_norm:
            seq_info.variable[pos] = VarPos(vec/vec_norm, state)

            # this position is now accounted for, hence we can subtract it from the count of fixed nucs
            seq_info.fixed_composition[state] -= 1

        # collect contribution from the fixed sites
        for state in alphabets[si]:
          # indeterminate parts in some children are not handled correctly here.
          # they should not contribute to the product. This will require some additional
          # handling or could be handled by treating these positions as variable
          child_profiles = []
          for ci,(c,e) in enumerate(node.children):
            child_profiles.append(expQT[ci].dot(c.seq_info_ingroup[si].fixed[state]))
          vec = np.prod(child_profiles, axis=0)
          vec_norm = vec.sum()

          seq_info.logLH += seq_info.fixed_composition[state]*np.log(vec_norm)
          seq_info.fixed[state] = vec/vec_norm

        node.payload.seq_info_ingroup.append(seq_info)

  G.par_iter_backward(calculate_ingroup_node)

def calculate_root_state(graph: Graph):
  root_node = graph.get_roots()[0].payload()
  logLH = 0
  for si, seq_info in enumerate(root_node.seq_info_ingroup):
    gtr = graph.partitions[si].gtr
    seq_profile = SeqInfoLh(unknown=seq_info.unknown, gaps=seq_info.gaps, non_char=seq_info.non_char,
                            variable={pos:v for pos, v in seq_info.variable.items()},
                            fixed={pos:v for pos, v in seq_info.fixed.items()},
                            fixed_composition={pos:v for pos, v in seq_info.fixed_composition.items()})
    seq_profile.logLH += seq_info.logLH
    for pos, p in seq_info.variable.items():
      vec = p.profile*gtr.Pi
      vec_norm = vec.sum()
      seq_profile.logLH += np.log(vec_norm)
      seq_profile.variable[pos] = VarPos(vec/vec_norm, p.state)
    for state in gtr.alphabet:
      vec = seq_info.fixed[state]*gtr.Pi
      vec_norm = vec.sum()
      seq_profile.logLH += np.log(vec_norm)*seq_info.fixed_composition[state]
      seq_info.fixed[state] = vec/vec_norm
    logLH += seq_profile.logLH
  return logLH



if __name__=="__main__":
  # fname_nwk = 'data/ebola/ebola.nwk'
  # fname_seq = 'data/ebola/ebola_dna.fasta'
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

  ingroup_profiles(G)
  print(calculate_root_state(G))

