
from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node
from lib import RangeCollection_intersection, RangeCollection, Mut, SeqPartition, DenseSeqNode, DenseSeqEdge, DenseSeqInfo, DenseSeqDis, find_char_ranges
from treetime import GTR
from treetime.seq_utils import seq2prof
from profile_map import profile_map
from typing import List, Dict
from payload import NodePayload, EdgePayload

profile_map['-'] = profile_map['N']

def attach_seqs_to_graph_dense(graph: Graph, aln_list: List[Dict[str,np.array]], gtr_list: List[GTR]) -> None:
  name_to_key = {}
  for n in graph.get_leaves():
    n.payload().dense_sequences = []
    name_to_key[n.payload().name] = n.key()
  for edge in graph.get_edges():
    edge.payload().dense_sequences = []

  graph.dense_partitions = []
  for aln, gtr in zip(aln_list, gtr_list):
    L = len(aln[list(aln.keys())[0]]) # stupid way to get alignment length
    partition = SeqPartition(gtr=gtr, length=L, profile_map=profile_map)
    graph.dense_partitions.append(partition)

    for leaf in graph.get_leaves():
      leaf_name = leaf.payload().name
      seq_array = np.fromiter(aln[leaf_name], 'U1')
      seq_node = DenseSeqNode(seq=DenseSeqInfo(sequence=seq_array, gaps=find_char_ranges(seq_array, '-')),
                              msg_to_parents=DenseSeqDis(dis=seq2prof(seq_array, profile_map=profile_map)))
      leaf.payload().dense_sequences.append(seq_node)

    for edge in graph.get_edges():
      edge.payload().dense_sequences.append(DenseSeqEdge())

def init_sequences_dense(graph: Graph, aln_list: List[Dict[str,np.array]], gtr_list: List[GTR]) -> None:
  attach_seqs_to_graph_dense(graph, aln_list, gtr_list)

def combine_dense_messages(msgs):
  seq_dis = DenseSeqDis(logLH=np.sum([m.logLH for m in msgs.values()]))
  prod_dis = np.prod([m.dis for m in msgs.values()], axis=0)
  norm = prod_dis.sum(axis=1)
  seq_dis.logLH += np.sum(np.log(norm))
  seq_dis.dis = (prod_dis.T/norm).T
  return seq_dis

def ingroup_profiles_dense(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.dense_partitions]
  gtrs = [p.gtr for p in graph.dense_partitions]
  seq_len = [p.length for p in graph.dense_partitions]

  def calculate_ingroup_node(node: GraphNodeBackward) -> None:
    if node.is_leaf:
      return # the msg to parents for leaves was put in init phase

    node.payload.dense_sequences=[]
    for si,gtr in enumerate(gtrs):
      child_expQt = [gtr.expQt(e.branch_length or 0.0) for c,e in node.children]
      child_seqs =  [c.dense_sequences[si] for c,e in node.children]
      child_edges = [e.dense_sequences[si] for c,e in node.children]
      gap_intersection = RangeCollection_intersection([cseq.seq.gaps for cseq in child_seqs])
      seq_info = DenseSeqInfo(gaps=gap_intersection)

      msgs_from_children = {}
      for ci, (c,e) in enumerate(node.children):
        seq_dis = DenseSeqDis(dis = np.ones((seq_len[si], len(alphabets[si]))))
        seq_dis.logLH += child_seqs[ci].msg_to_parents.logLH
        if child_edges[ci].transmission:
          for r in child_edges[ci].transmission:
            seq_dis.dis[r.start:r.end] *= child_seqs[ci].msg_to_parents.dis[r.start:r.end].dot(child_expQt[ci])
        else:
          seq_dis.dis *= child_seqs[ci].msg_to_parents.dis.dot(child_expQt[ci])
        msgs_from_children[c.name] = seq_dis

      msg_to_parents = combine_dense_messages(msgs_from_children)
      node.payload.dense_sequences.append(DenseSeqNode(seq=seq_info, msgs_from_children=msgs_from_children,
                                                       msg_to_parents=msg_to_parents))

  graph.par_iter_backward(calculate_ingroup_node)

def outgroup_profiles_dense(graph: Graph):
  alphabets = [''.join(p.gtr.alphabet) for p in graph.dense_partitions]
  gtrs = [p.gtr for p in graph.dense_partitions]
  seq_len = [p.length for p in graph.dense_partitions]

  def calculate_outgroup_node(node: GraphNodeBackward) -> None:
    if node.is_root:
      return
    for si,seq_info in enumerate(node.payload.dense_sequences):
      msgs_from_parents = {}
      for p,e in node.parents:
        eQt = gtrs[si].expQt(e.branch_length or 0.0).T
        pseq = p.dense_sequences[si]
        seq_dis = DenseSeqDis(dis = np.ones((seq_len[si], len(alphabets[si]))),
                              logLH=pseq.msgs_to_children[node.payload.name].logLH)
        if e.dense_sequences[si].transmission:
          for r in e.dense_sequences[si].transmission:
            seq_dis.dis[r.start:r.end] *= pseq.msgs_to_children[node.payload.name].dis[r.start:r.end].dot(eQt)
        else:
          # could make this a copy
          seq_dis.dis *= pseq.msgs_to_children[node.payload.name].dis.dot(eQt)
        msgs_from_parents[p.name] = seq_dis
      msgs_from_parents['children'] = seq_info.msg_to_parents
      seq_info.profile = combine_dense_messages(msgs_from_parents)

      seq_info.msgs_to_children = {}
      for cname in seq_info.msgs_from_children:
        seq_info.msgs_to_children[cname] = DenseSeqDis(dis = seq_info.profile.dis/seq_info.msgs_from_children[cname].dis,
                              logLH=seq_info.profile.logLH-seq_info.msgs_from_children[cname].logLH)


  graph.par_iter_forward(calculate_outgroup_node)

def calculate_root_state_dense(graph: Graph):
  logLH=0
  for root_node in graph.get_roots():
    for si,seq_info in enumerate(root_node.payload().dense_sequences):
      gtr = graph.dense_partitions[si].gtr
      seq_dis = DenseSeqDis(logLH=seq_info.msg_to_parents.logLH)
      seq_dis.dis = seq_info.msg_to_parents.dis*gtr.Pi
      norm = seq_dis.dis.sum(axis=1)
      seq_dis.logLH += np.sum(np.log(norm))
      seq_dis.dis = (seq_dis.dis.T/norm).T

      seq_info.profile = seq_dis
      seq_info.msgs_to_children = {}
      for cname in seq_info.msgs_from_children:
        seq_info.msgs_to_children[cname] = DenseSeqDis(dis = seq_info.profile.dis/seq_info.msgs_from_children[cname].dis,
                              logLH=seq_info.profile.logLH-seq_info.msgs_from_children[cname].logLH)
      logLH+=seq_info.profile.logLH
  return logLH

def ancestral_dense(graph: Graph) -> float:
  ingroup_profiles_dense(graph)
  logLH = calculate_root_state_dense(graph)
  outgroup_profiles_dense(graph)
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
  init_sequences_dense(G, [aln], [gtr])

  print("LogLH", ancestral_dense(G))
