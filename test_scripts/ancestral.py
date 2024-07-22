from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node
from lib import SeqPartition, SparseSeqDis, VarPos, RangeCollection_intersection, RangeCollection, Mut
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map
from fitch import init_sparse_sequences

def get_variable_states_children(child_seqs, child_edges):
  variable_pos = {}
  child_states = [{} for c in child_edges]
  for ci, edge in enumerate(child_edges):
    # go over all mutations and get reference state
    for m in edge.muts:
      variable_pos[m.pos] = m.ref  # this might be set multiple times, but the reference state should always be the same
      child_states[ci][m.pos] = m.qry

  for ci, seq_info in enumerate(child_seqs):
    # go over child variable position and get reference state
    for pos, p in seq_info.msg_to_parents.variable.items():
      if pos not in variable_pos:
        variable_pos[pos] = p.state

  return variable_pos, child_states

def get_variable_states_parents(node_seq_info, parents, parent_edges):
  variable_pos = {pos: p.state for pos, p in node_seq_info.variable.items()}
  parent_states = [{} for p in parents]
  for pi, pseq in enumerate(parents):
    # go over all mutations and get reference state
    for m in parent_edges[pi].muts:
      variable_pos[m.pos] = m.qry  # this might be set multiple times, but the reference state should always be the same
      parent_states[pi][m.pos] = m.ref

    # go over child variable position and get reference state
    for pos, p in pseq.variable.items():
      if pos not in variable_pos:
        variable_pos[pos] = p.state

  return variable_pos, parent_states


def propagate(origin, expQt, seq_dis, variable_pos, child_states, non_char, transmission):
  message = SparseSeqDis(origin=origin, fixed_counts={k:v for k,v in seq_dis.fixed_counts.items()}, variable={}, fixed={})
  for pos, state in variable_pos.items():
    if transmission and (not transmission.contains(pos)):
      # transmission field is not currently used
      continue
    if non_char.contains(pos):
      continue

    if pos in seq_dis.variable:
      v = seq_dis.variable[pos].profile
    elif pos in child_states:
      v = seq_dis.fixed[child_states[pos]]
      message.fixed_counts[child_states[pos]] -= 1
    else:
      v = seq_dis.fixed[state]
      message.fixed_counts[state] -= 1

    message.variable[pos] = VarPos(expQt.dot(v), state)

  for s, p in seq_dis.fixed.items():
    message.fixed[s] = expQt.dot(p)

  return message

def combine_messages(seq_dis, messages, variable_pos, eps, alphabet):
  # go over all putatively variable positions
  for pos, state in variable_pos.items():
    # collect the profiles of children to multiply
    msg_profiles = []
    for msg in messages:
      if pos in msg.variable:
        msg_profiles.append(msg.variable[pos].profile)

    # calculate new profile and likelihood contribution
    try:
      vec = np.prod(msg_profiles, axis=0)
    except:
      print(msg_profiles)
      import ipdb; ipdb.set_trace()
    vec_norm = vec.sum()

    # add position to variable states if the subleading states have a probability exceeding eps
    if vec.max()<(1-eps)*vec_norm:
      if len(vec.shape)>1:
        import ipdb; ipdb.set_trace()
      seq_dis.variable[pos] = VarPos(vec/vec_norm, state)
      seq_dis.logLH  += np.log(vec_norm)

      # this position is now accounted for, hence we can subtract it from the count of fixed nucs
      seq_dis.fixed_counts[state] -= 1

  # collect contribution from the fixed sites
  for state in alphabet:
    # indeterminate parts in some children are not handled correctly here.
    # they should not contribute to the product. This will require some additional
    # handling or could be handled by treating these positions as variable
    msg_profiles = []
    for msg in messages:
      msg_profiles.append(msg.fixed[state])
    vec = np.prod(msg_profiles, axis=0)
    vec_norm = vec.sum()

    seq_dis.logLH += seq_dis.fixed_counts[state]*np.log(vec_norm)
    seq_dis.fixed[state] = vec/vec_norm


def sparse_ingroup_profiles(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.partitions]
  gtrs = [p.gtr for p in graph.partitions]

  eps=1e-6
  def calculate_ingroup_node(node: GraphNodeBackward) -> None:

    for si,seq_info in enumerate(node.payload.sparse_sequences):
      if node.is_leaf:
        # this is mostly a copy (or ref here) of the fitch state.
        seq_info.msg_to_parent = SparseSeqDis(fixed_counts=seq_info.state.composition, variable=seq_info.state.distribution.variable, origin='input',
                                                fixed={state:graph.partitions[si].profile(state) for state in alphabets[si]})
      else: # internal nodes
        # get all variable positions, the reference state, and the child states at these positions
        child_expQt = [gtrs[si].expQt(e.branch_length or 0.0).T for c,e in node.children]
        child_seqs = [c.sparse_sequences[si] for c,e in node.children]
        child_edges = [e.sparse_sequences[si] for c,e in node.children]
        variable_pos, child_states = get_variable_states_children(child_seqs, child_edges)

        seq_dis = SparseSeqDis(fixed_counts={k:v for k,v in seq_info.state.composition.items()})

        seq_info.msgs_from_children = []
        for ci, (c,e) in enumerate(node.children):
          seq_dis.logLH += child_seqs[ci].msg_to_parents.logLH
          seq_info.msgs_from_children.append(propagate(origin=c.name, expQt=child_expQt[ci],
                                                       seq_dis=child_seqs[ci].msg_to_parent, variable_pos=variable_pos,
                                                       child_states = child_states[ci], non_char=child_seqs[ci].state.non_char, transmission=None))

        combine_messages(seq_dis, seq_info.msgs_from_children, variable_pos, eps, alphabets[si])

        seq_info.msg_to_parent=seq_dis

  graph.par_iter_backward(calculate_ingroup_node)


def outgroup_profiles(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.partitions]
  gtrs = [p.gtr for p in graph.partitions]

  eps=1e-6
  def calculate_outgroup_node(node: GraphNodeForward) -> None:
    for si, seq_info in enumerate(node.payload.sparse_sequences):
      variable_pos, parent_states = get_variable_states_parents(seq_info.msg_to_parents,
                                                                [p.sparse_sequences[si].profile for p,e in node.parents],
                                                                [e.sparse_sequences[si] for p,e in node.parents])
      for p, e in node.parents:
        pseq_info = p.sparse_sequences[si]
        messages = pseq_info.msgs_from_parents + [m for m in pseq_info.msgs_from_children if m.origin!=node.payload.name]
        # gaps, unknown, etc should be the from the combined messages
        seq_dis = SparseSeqDis(fixed_counts={k:v for k,v in pseq_info.profile.fixed_counts.items()})

        combine_messages(seq_dis, messages, variable_pos=variable_pos, eps=eps, alphabet=alphabets[si])

        seq_info.msgs_from_parents.append(propagate(origin=p.name, expQt=gtrs[si].expQt(e.branch_length or 0.0),
                                                     seq_dis=seq_dis, variable_pos=variable_pos, child_states=parent_states,
                                                     non_char=pseq_info.state.non_char, transmission=None))

      # gaps, unknown, etc should be the from the combined messages
      seq_dis = SparseSeqDis(fixed_counts={k:v for k,v in seq_info.state.composition.items()})
      combine_messages(seq_dis=seq_dis, messages=seq_info.msgs_from_parents + [seq_info.msg_to_parent],
                       variable_pos=variable_pos, eps=eps, alphabet=alphabets[si])
      seq_info.profile = seq_dis

  graph.par_iter_forward(calculate_outgroup_node)

def calculate_root_state(graph: Graph):
  root_node = graph.get_roots()[0].payload()
  logLH = 0
  for si, seq_info in enumerate(root_node.sparse_sequences):
    gtr = graph.partitions[si].gtr
    seq_profile = SparseSeqDis(fixed_counts={pos:v for pos, v in seq_info.msg_to_parent.fixed_counts.items()},
                               logLH=seq_info.msg_to_parents.logLH)

    for pos, p in seq_info.msg_to_parent.variable.items():
      vec = p.profile*gtr.Pi
      vec_norm = vec.sum()
      seq_profile.logLH += np.log(vec_norm)
      seq_profile.variable[pos] = VarPos(vec/vec_norm, p.state)
    for state, p in seq_info.msg_to_parent.fixed.items():
      vec = p*gtr.Pi
      vec_norm = vec.sum()
      seq_profile.logLH += np.log(vec_norm)*seq_info.msg_to_parent.fixed_counts[state]
      seq_profile.fixed[state] = vec/vec_norm
    logLH += seq_profile.logLH

    seq_info.profile=seq_profile
  return logLH



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

  sparse_ingroup_profiles(G)
  print(calculate_root_state(G))
  outgroup_profiles(G)

