from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward
from lib import SparseSeqDis, VarPos
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map
from fitch import init_sequences_sparse

eps=1e-6

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


def propagate(expQt, seq_dis, variable_pos, child_states, non_char, transmission):
  message = SparseSeqDis(fixed_counts={k:v for k,v in seq_dis.fixed_counts.items()}, variable={}, fixed={})
  for pos, state in variable_pos.items():
    if transmission and (not transmission.contains(pos)):
      # transmission field is not currently used
      continue
    if non_char.contains(pos):
      continue

    if pos in seq_dis.variable:
      v = seq_dis.variable[pos].dis
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

def combine_messages(seq_dis, messages, variable_pos, eps, alphabet, gtr_weight=None):
  # go over all putatively variable positions
  for pos, state in variable_pos.items():
    # collect the profiles of children to multiply
    msg_dis = [] if gtr_weight is None else [gtr_weight]
    for msg in messages:
      if pos in msg.variable:
        msg_dis.append(msg.variable[pos].dis)

    # calculate new profile and likelihood contribution
    vec = np.prod(msg_dis, axis=0)
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
    msg_dis = [] if gtr_weight is None else [gtr_weight]
    for msg in messages:
      msg_dis.append(msg.fixed[state])
    vec = np.prod(msg_dis, axis=0)
    vec_norm = vec.sum()

    seq_dis.logLH += seq_dis.fixed_counts[state]*np.log(vec_norm)
    seq_dis.fixed[state] = vec/vec_norm


def ingroup_profiles_sparse(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.sparse_partitions]
  gtrs = [p.gtr for p in graph.sparse_partitions]

  def calculate_ingroup_node(node: GraphNodeBackward) -> None:

    for si,seq_info in enumerate(node.payload.sparse_sequences):
      if node.is_leaf:
        # this is mostly a copy (or ref here) of the fitch state.
        seq_info.msg_to_parents = SparseSeqDis(fixed_counts=seq_info.seq.composition, variable=seq_info.seq.fitch.variable,
                                                fixed={state:graph.sparse_partitions[si].profile(state) for state in alphabets[si]})
      else: # internal nodes
        child_expQt = [gtrs[si].expQt(e.branch_length or 0.0).T for c,e in node.children]
        child_seqs =  [c.sparse_sequences[si] for c,e in node.children]
        child_edges = [e.sparse_sequences[si] for c,e in node.children]

        # get all variable positions, the reference state, and the child states at these positions
        variable_pos, child_states = get_variable_states_children(child_seqs, child_edges)

        seq_dis = SparseSeqDis(fixed_counts={k:v for k,v in seq_info.seq.composition.items()})

        seq_info.msgs_from_children = {}
        for ci, (c,e) in enumerate(node.children):
          seq_dis.logLH += child_seqs[ci].msg_to_parents.logLH
          seq_info.msgs_from_children[c.name] = propagate(expQt=child_expQt[ci],
                                                       seq_dis=child_seqs[ci].msg_to_parents, variable_pos=variable_pos,
                                                       child_states = child_states[ci], non_char=child_seqs[ci].seq.non_char,
                                                       transmission=None)

        combine_messages(seq_dis, list(seq_info.msgs_from_children.values()), variable_pos, eps, alphabets[si])
        seq_info.msg_to_parents=seq_dis

  graph.par_iter_backward(calculate_ingroup_node)


def outgroup_profiles_sparse(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.sparse_partitions]
  gtrs = [p.gtr for p in graph.sparse_partitions]

  def calculate_outgroup_node(node: GraphNodeForward) -> None:
    if node.is_root:
      return
    for si, seq_info in enumerate(node.payload.sparse_sequences):
      variable_pos, parent_states = get_variable_states_parents(seq_info.msg_to_parents,
                                                                [p.sparse_sequences[si].profile for p,e in node.parents],
                                                                [e.sparse_sequences[si] for p,e in node.parents])
      msgs_from_parents = []
      for p, e in node.parents:
        pseq_info = p.sparse_sequences[si]
        # gaps, unknown, etc should be the from the combined messages
        seq_dis = SparseSeqDis(fixed_counts={k:v for k,v in pseq_info.profile.fixed_counts.items()})
        msgs_from_parents.append(propagate(expQt=gtrs[si].expQt(e.branch_length or 0.0),
                                           seq_dis=pseq_info.msgs_to_children[node.payload.name],
                                           variable_pos=variable_pos, child_states=parent_states,
                                           non_char=pseq_info.seq.non_char, transmission=None))

      # gaps, unknown, etc should be the from the combined messages
      seq_dis = SparseSeqDis(fixed_counts={k:v for k,v in seq_info.seq.composition.items()})
      combine_messages(seq_dis=seq_dis, messages=msgs_from_parents + [seq_info.msg_to_parents],
                       variable_pos=variable_pos, eps=eps, alphabet=alphabets[si])
      seq_info.profile = seq_dis

      # precalculate messages to children that summarize into from their siblings (from the backward pass) and the parent
      # note that this could be replaced by "subtracting/dividing" the profile by the msg from the respective child
      # in some circumstances, this might be more efficient
      seq_info.msgs_to_children = {}
      for c in seq_info.msgs_from_children:
        msgs = [m for k,m in seq_info.msgs_from_children.items() if k!=c]
        seq_dis = SparseSeqDis(fixed_counts={pos:v for pos, v in seq_info.msg_to_parents.fixed_counts.items()})
        combine_messages(seq_dis, msgs, variable_pos=variable_pos, eps=eps, alphabet=gtr.alphabet)
        seq_info.msgs_to_children[c] = seq_dis

  graph.par_iter_forward(calculate_outgroup_node)

def calculate_root_state_sparse(graph: Graph):

  for root_node in  graph.get_roots():
    logLH = 0
    for si, seq_info in enumerate(root_node.payload().sparse_sequences):
      gtr = graph.sparse_partitions[si].gtr
      seq_profile = SparseSeqDis(fixed_counts={pos:v for pos, v in seq_info.msg_to_parents.fixed_counts.items()},
                                logLH=seq_info.msg_to_parents.logLH)

      # multiply the info from the tree with the GTR equilibrium probabilities (variable and fixed)
      for pos, p in seq_info.msg_to_parents.variable.items():
        vec = p.dis*gtr.Pi
        vec_norm = vec.sum()
        seq_profile.logLH += np.log(vec_norm)
        seq_profile.variable[pos] = VarPos(vec/vec_norm, p.state)
      for state, p in seq_info.msg_to_parents.fixed.items():
        vec = p*gtr.Pi
        vec_norm = vec.sum()
        seq_profile.logLH += np.log(vec_norm)*seq_info.msg_to_parents.fixed_counts[state]
        seq_profile.fixed[state] = vec/vec_norm

      seq_info.profile=seq_profile
      logLH += seq_profile.logLH

      # calculate messages to children
      variable_pos = {pos:p.state for pos, p in seq_info.msg_to_parents.variable.items()}
      seq_info.msgs_to_children = {}
      for c,e in graph.children_of(root_node.key()):
        cname = c.payload().name
        msgs = [m for k,m in seq_info.msgs_from_children.items() if k!=cname]
        seq_dis = SparseSeqDis(fixed_counts={pos:v for pos, v in seq_info.msg_to_parents.fixed_counts.items()})
        combine_messages(seq_dis, msgs, variable_pos=variable_pos, eps=eps, alphabet=gtr.alphabet,
                        gtr_weight=gtr.Pi)
        seq_info.msgs_to_children[cname] = seq_dis

    return logLH

def ancestral_sparse(graph: Graph) -> float:
  ingroup_profiles_sparse(graph)
  logLH = calculate_root_state_sparse(graph)
  outgroup_profiles_sparse(graph)
  return logLH

def tests():
  aln = {"root":"ACAGCCATGTATTG--",
         "AB":"ACATCCCTGTA-TG--",
         "A":"ACATCGCCNNA--GAC",
         "B":"GCATCCCTGTA-NG--",
         "CD":"CCGGCCATGTATTG--",
         "C":"CCGGCGATGTRTTG--",
         "D":"TCGGCCGTGTRTTG--"}

  tree = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;"

  G = graph_from_nwk_str(nwk_string=tree, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  gtr = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
  init_sequences_sparse(G, [aln], [gtr])

  ingroup_profiles_sparse(G)
  seq_info_root = G.get_one_root().payload().sparse_sequences[0]

  assert tuple(sorted(seq_info_root.msg_to_parents.variable.keys()))==(0,2,3,5,6,7,10)
  assert tuple([round(x,8) for x in seq_info_root.msg_to_parents.variable[0].dis])==(0.34485164, 0.17637237, 0.22492433, 0.25385166)

  calculate_root_state_sparse(G)
  assert tuple([round(x,8) for x in seq_info_root.profile.variable[0].dis])==(0.28212327, 0.21643546, 0.13800802, 0.36343326)
  assert np.abs(seq_info_root.profile.fixed['G']-np.array([1.76723056e-04, 2.65084585e-04, 9.99248927e-01, 3.09265349e-04])).sum()<1e-6

  outgroup_profiles_sparse(G)
  node_AB = G.get_node(G.nodes[1].key()).payload().sparse_sequences[0]
  assert np.abs(node_AB.profile.variable[0].dis-np.array([0.51275208, 0.09128506, 0.24647255, 0.14949031])).sum()<1e-6

  # redo the tests with two sequence sets using different GTRs
  G2 = graph_from_nwk_str(nwk_string=tree, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  gtr1 = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
  gtr2 = GTR.custom(pi=[0.4, 0.15, 0.12, 0.33], alphabet='nuc_nogap')
  init_sequences_sparse(G2, [aln, aln], [gtr1, gtr2])

  ingroup_profiles_sparse(G2)
  seq_info_root_1 = G2.get_one_root().payload().sparse_sequences[0]
  seq_info_root_2 = G2.get_one_root().payload().sparse_sequences[1]

  assert tuple(sorted(seq_info_root_1.msg_to_parents.variable.keys()))==(0,2,3,5,6,7,10)
  assert tuple(sorted(seq_info_root_2.msg_to_parents.variable.keys()))==(0,2,3,5,6,7,10)
  assert tuple([round(x,8) for x in seq_info_root_1.msg_to_parents.variable[0].dis])==(0.34485164, 0.17637237, 0.22492433, 0.25385166)
  assert tuple([round(x,8) for x in seq_info_root_2.msg_to_parents.variable[0].dis])==(0.17085068, 0.31568135, 0.26076566, 0.25270231)

  calculate_root_state_sparse(G2)
  assert tuple([round(x,8) for x in seq_info_root_1.profile.variable[0].dis])==(0.28212327, 0.21643546, 0.13800802, 0.36343326)
  assert tuple([round(x,8) for x in seq_info_root_2.profile.variable[0].dis])==(0.29664652, 0.20554302, 0.13582953, 0.36198093)
  assert np.abs(seq_info_root_2.profile.fixed['G']-np.array([2.78247843e-04, 1.04342941e-04, 9.99387855e-01, 2.29554470e-04])).sum()<1e-6

  outgroup_profiles_sparse(G2)
  node_AB = G2.get_node(G2.nodes[1].key()).payload().sparse_sequences[1]
  assert np.abs(node_AB.profile.variable[0].dis-np.array([0.52331521, 0.08336271, 0.24488808, 0.148434])).sum()<1e-6



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
  init_sequences_sparse(G, [aln], [gtr])

  print("LogLH", ancestral_sparse(G))

  seq_info = G.get_roots()[0].payload().sparse_sequences[0]