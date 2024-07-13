from dataclasses import dataclass
from Bio import AlignIO
import numpy as np
from typing import Optional, List, Dict
from lib import Graph, AutoRepr, Mut, graph_from_nwk_str, Node, GraphNodeBackward, GraphNodeForward
from treetime import GTR

# definition and function to handle mixed sites
#                        "A,C,G,T"
default_profile = np.array([1,1,1,1], dtype=float)
profiles = {'A':np.array([1,0,0,0], dtype=float),
            'C':np.array([0,1,0,0], dtype=float),
            'G':np.array([0,0,1,0], dtype=float),
            'T':np.array([0,0,0,1], dtype=float),
            'R':np.array([1,0,1,0], dtype=float),
            'Y':np.array([0,1,0,1], dtype=float),
            'S':np.array([0,1,1,0], dtype=float),
            'W':np.array([1,0,0,1], dtype=float),
            'K':np.array([0,0,1,1], dtype=float),
            'M':np.array([1,1,0,0], dtype=float),
            'B':np.array([0,1,1,1], dtype=float),
            'D':np.array([1,0,1,1], dtype=float),
            'H':np.array([1,1,0,1], dtype=float),
            'V':np.array([1,1,1,0], dtype=float),
            'N':default_profile,
            '.':default_profile}

reverse_profile = {tuple(p): nuc for nuc, p in profiles.items()}

n_seq_partitions = 1
alphabet='ACGT'

@dataclass
class Range(AutoRepr):
  start: int
  end: int

@dataclass
class VarPos(AutoRepr):
  profile: np.array
  state: str

@dataclass
class Seq(AutoRepr):
  gaps: List[Range]
  unknown: List[Range]
  indeterminate: List[Range]
  ingroupLH: List[float]
  variable_ingroup:  Optional[Dict[int, VarPos]]
  variable_outgroup: Optional[Dict[int, VarPos]]
  variable_profile:  Optional[Dict[int, VarPos]]
  fixed_ingroup:     Optional[Dict[str, np.array]]
  fixed_outgroup:    Optional[Dict[str, np.array]]
  fixed_profile:     Optional[Dict[str, np.array]]
  state_composition: Optional[Dict[str, int]]


@dataclass
class NodePayload(AutoRepr):
  name: Optional[str]
  date: Optional[float]
  full_seq: Optional[List[str]]
  seq: Optional[List[Seq]]
  seq_len: Optional[List[int]]


@dataclass
class EdgePayload():
  branchlength: Optional[float]
  muts: List[List[Mut]]

def find_char_ranges(seq, char):
  ranges = []
  start = None
  for pos, nuc in enumerate(seq):
    if nuc != char and (start is not None):
      ranges.append(Range(start, pos))
      start = None
    elif nuc == char and start is None:
      start = pos
  if start:
    ranges.append(Range(start, len(seq)))
  return ranges

def range_intersection(range_sets: List[List[Range]]):
  '''
  Note, this assumes sorted ranges
  '''
  if any([len(r)==0 for r in range_sets]):
    return []

  current_ranges = [Range(r.start, r.end) for r in range_sets[0]]
  for next_ranges in range_sets[1:]:
    new_ranges = []
    ri1, ri2 = 0, 0
    r1 = current_ranges[ri1]
    r2 = next_ranges[ri2]
    while ri1<len(current_ranges) and ri2<len(next_ranges):
      if r2.start>r1.end:
        ri1 += 1
        if ri1<len(current_ranges):
          r1 = current_ranges[ri1]
      elif r1.start>r2.end:
        ri2 += 1
        if ri2<len(next_ranges):
          r2 = next_ranges[ri2]
      else:
        new_ranges.append(Range(max(r1.start, r2.start), min(r1.end, r2.end)))
        if r1.end<r2.end:
          ri1 += 1
          if ri1<len(current_ranges):
            r1 = current_ranges[ri1]
        else:
          ri2 += 1
          if ri2<len(next_ranges):
            r2 = next_ranges[ri2]
    current_ranges = new_ranges
  return current_ranges

def contains(ranges: List[Range], pos: int) -> bool:
  for r in ranges:
    if r.start<=pos and pos<r.end: return True
  return False

def generate_sparse_sequence_representation(graph):
  '''
  Starting from graph where full sequences are attached at the leaves, determine
  the root sequence and the mutations necessary to generate all sequences from the root
  while saving the parts of the sequence that are indeterminate (via gaps or Ns) as ranges.
  Letters that are considered indeterminate are probably best abstracted.
  '''
  def fitch_backwards(node: GraphNodeBackward):
    if node.is_leaf: # process sequences on leaves
      for si in range(n_seq_partitions):
        seq = node.payload.full_seq[si]
        # code the missing regions and gap regions along with the ambiguous nucleotides
        seq_rep = Seq(gaps=find_char_ranges(seq, '-'),
                      unknown=find_char_ranges(seq, 'N'),
                      indeterminate=sorted(find_char_ranges(seq, '-')+find_char_ranges(seq, 'N'), key=lambda x:x.start), ingroupLH=0.0,
                      variable_ingroup={pos:VarPos(profiles.get(nuc, default_profile), None) for pos, nuc in enumerate(seq) if nuc not in 'ACTG-N'},
                      variable_outgroup={}, variable_profile={}, fixed_ingroup={}, fixed_outgroup={}, fixed_profile={}, state_composition={})
        node.payload.seq.append(seq_rep)
    else:
      # process internal nodes, again for each sequence partition
      for si in range(n_seq_partitions):
        # init the local representation with gaps, unknowns
        seq_rep = Seq(gaps=range_intersection([c.seq[si].gaps for c,e in node.children]),
                      unknown=range_intersection([c.seq[si].unknown for c,e in node.children]),
                      indeterminate=range_intersection([c.seq[si].indeterminate for c,e in node.children]),
                      ingroupLH=0.0, variable_ingroup={}, variable_outgroup={}, variable_profile={},
                      fixed_ingroup={}, fixed_outgroup={}, fixed_profile={}, state_composition={})
        # define dummy sequence
        full_seq = np.array(['?']*node.children[0][0].seq_len[si])
        variable_ingroup = {}
        # fill indeterminate positions
        for r in seq_rep.indeterminate:
          for pos in range(r.start, r.end):
            full_seq[pos]='.'

        # process all positions where the children are variable
        variable_pos = np.unique(np.concatenate([np.array(list(c.seq[si].variable_ingroup.keys()), dtype=int) for c,e in node.children]))
        for pos in variable_pos:
          # 2D float array
          # stack all ingroup profiles of children, use the determinate profile if position is not variable.
          child_profiles = np.array([c.seq[si].variable_ingroup[pos].profile
                                      if pos in c.seq[si].variable_ingroup
                                      else profiles[c.full_seq[si][pos]]
                                      for c,e in node.children])

          isect = child_profiles.prod(axis=0)
          if isect.sum()==1:
            full_seq[pos]=alphabet[np.argmax(isect)]
          elif isect.sum()>1:
            variable_ingroup[pos]=VarPos(isect, None)
            full_seq[pos]='~'
          else:
            child_union = np.sum(child_profiles, axis=0)
            variable_ingroup[pos]=VarPos(np.array([1.0 if x>0 else 0.0 for x in child_union]), None)
            full_seq[pos]='~'

        # process all positions where the children are fixed or completely unknown in some children
        # this is ridiculously slow, even though we are only comparing nucleotides in the children
        for pos, (nuc, child_states) in enumerate(zip(full_seq, zip(*[c.full_seq[si] for c,e in node.children]))):
          if nuc!='?':
            continue

          determined_states = set([x for x in child_states if x in alphabet])
          if len(determined_states)==1:
            full_seq[pos]=determined_states.pop()
          elif len(determined_states)>1:
            variable_ingroup[pos] = VarPos(np.sum([profiles[x] for x in determined_states], axis=0), None)
            full_seq[pos]='~'
          else:
            import ipdb; ipdb.set_trace()


        seq_rep.variable_ingroup.update(variable_ingroup)
        node.payload.seq.append(seq_rep)
        node.payload.full_seq.append(full_seq)
        node.payload.seq_len.append(len(full_seq))

  def fitch_forward(node: GraphNodeForward):
    if node.is_root:
      for seq_rep, full_seq in zip(node.payload.seq, node.payload.full_seq):
        for pos, p in seq_rep.variable_ingroup.items():
          full_seq[pos] = alphabet[np.argmax(p.profile)]
        seq_rep.state_composition = {s:0 for s in alphabet}
        for s in full_seq:
          seq_rep.state_composition[s] += 1
    else:
      # only deal with one parent for now
      parent, edge = node.parents[0]
      # loop over sequence partitions
      for si, (seq_rep, full_seq) in enumerate(zip(node.payload.seq, node.payload.full_seq)):
        pseq=parent.full_seq[si]
        seq_rep.state_composition = {s:k for s,k in parent.seq[si].state_composition.items()}

        # fill in the indeterminate positions.
        for r in seq_rep.indeterminate:
          for pos in range(r.start, r.end):
            full_seq[pos]=pseq[pos]
            seq_rep.state_composition[pseq[pos]] -= 1

        muts = []
        # for each variable position, pick a state or a mutation
        for pos, p in seq_rep.variable_ingroup.items():
          pnuc=pseq[pos]
          if np.sum(profiles[pnuc]*p.profile):
            full_seq[pos] = pnuc
          else:
            cnuc = alphabet[np.argmax(p.profile)]
            full_seq[pos] = cnuc
            # could be factored out
            muts.append(Mut(pnuc, pos, cnuc))
            seq_rep.state_composition[pnuc] -= 1
            seq_rep.state_composition[cnuc] += 1
        for pos in parent.seq[si].variable_ingroup:
          if pseq[pos]!=full_seq[pos]:
            # could be factored out
            muts.append(Mut(pseq[pos], pos, full_seq[pos]))
            seq_rep.state_composition[pseq[pos]] -= 1
            seq_rep.state_composition[full_seq[pos]] += 1
        for pos in seq_rep.variable_ingroup:
          seq_rep.variable_ingroup[pos].state = full_seq[pos]
        edge.muts.append(muts)
      # print(edge.muts)

  def clean_up(node):
    if not node.is_root:
      node.payload.full_seq = [[] for n in range(n_seq_partitions)]
    if not node.is_leaf:
      for seq_rep in node.payload.seq:
        seq_rep.variable_ingroup = {}

  graph.par_iter_backward(visitor=fitch_backwards)

  graph.par_iter_forward(visitor=fitch_forward)

  graph.par_iter_forward(visitor=clean_up)


def reconstruct_sequence(G: Graph, node: Node):
  '''
  reconstruct a specific sequence using the full sequence at the root and implanting mutations
  '''
  path_to_root = [node]
  p = node
  while not p.is_root():
    assert len(p._inbound)==1
    p = G.get_node(G.get_edge(p._inbound[0]).source())
    path_to_root.append(p)

  seqs = [np.copy(s) for s in path_to_root[-1].payload().full_seq]
  for n in path_to_root[::-1][1:]:
    muts=G.get_edge(n._inbound[0]).payload().muts
    for s, mset  in zip(seqs, muts):
      for m in mset:
        s[m.pos] = m.qry
  for s, seq_rep in zip(seqs, node.payload().seq):
    for r in seq_rep.gaps:
      for pos in range(r.start, r.end):
        s[pos]='-'
    for r in seq_rep.unknown:
      for pos in range(r.start, r.end):
        s[pos]='N'
    for pos, p in seq_rep.variable_ingroup.items():
      s[pos]=reverse_profile[tuple(p.profile)]
  return seqs


def ingroup_profiles(G: Graph, gtrs: List[GTR]):
  eps=1e-6
  def calculate_ingroup(node: GraphNodeBackward) -> None:
    print(node.payload.name)
    if node.is_leaf:
      for seq_rep in node.payload.seq:
        seq_rep.fixed_ingroup = {state: profiles[state] for state in alphabet}
      return

    # get all variable positions and the reference state
    variable_pos = [{} for _ in range(n_seq_partitions)]
    child_states = [[{} for _ in range(n_seq_partitions)] for c, e in node.children]
    for ci, (c,e) in enumerate(node.children):
      # go over all mutations and get reference state
      for si, mutset in enumerate(e.muts):
        for m in mutset:
          variable_pos[si][m.pos] = m.ref
          child_states[ci][si][m.pos] = m.qry

      # go over child variable position and get reference state
      for si, seq_rep in enumerate(c.seq):
        for pos, p in seq_rep.variable_ingroup.items():
          if pos not in variable_pos[si]:
            variable_pos[si][pos] = p.state

    expQT = [[gtr.expQt(e.branchlength or 0.0) for c,e in node.children] for gtr in gtrs]
    for si, seq_rep in enumerate(node.payload.seq):
      seq_rep.ingroupLH = 0.0
      for c,e in node.children:
        seq_rep.ingroupLH += c.seq[si].ingroupLH

      variable_ingroup = {}
      for pos, state in variable_pos[si].items():
        child_profiles = []
        for ci,(c,e) in enumerate(node.children):
          # need to check whether position is determined
          if not contains(c.seq[si].indeterminate, pos):
            if pos in c.seq[si].variable_ingroup:
              v = c.seq[si].variable_ingroup[pos].profile
            elif pos in child_states[ci][si]:
              v = c.seq[si].fixed_ingroup[child_states[ci][si][pos]]
            else:
              v=c.seq[si].fixed_ingroup[state]
            child_profiles.append(expQT[si][ci].dot(v))

        vec = np.prod(child_profiles, axis=0)
        vec_norm = vec.sum()
        seq_rep.ingroupLH  += np.log(vec_norm)

        # add position to variable states if the subleading states have a probability exceeding eps
        if vec.max()<(1-eps)*vec_norm:
          variable_ingroup[pos] = VarPos(vec/vec_norm, state)

          # this position is accounted for, hence we can subtract it from the count of fixed nucs
          if state in 'ACGT':
              seq_rep.state_composition[state] -= 1

      seq_rep.variable_ingroup = variable_ingroup
      # collect contribution from the inert sites
      for state in alphabet:
        child_profiles = []
        for ci,(c,e) in enumerate(node.children):
          child_profiles.append(expQT[si][ci].dot(c.seq[si].fixed_ingroup[state]))
        vec = np.prod(child_profiles, axis=0)
        vec_norm = vec.sum()

        seq_rep.ingroupLH  += seq_rep.state_composition[state]*np.log(vec_norm)
        seq_rep.fixed_ingroup[state] = vec/vec_norm

  G.par_iter_backward(calculate_ingroup)

def calculate_root_state(G: Graph, gtrs: List[GTR]):
  root_node = G.get_roots()[0].payload()
  logLH = 0
  for si, seq_rep in enumerate(root_node.seq):
    logLH += seq_rep.ingroupLH
    for pos, p in seq_rep.variable_ingroup.items():
      vec = p.profile*gtrs[si].Pi
      vec_norm = vec.sum()
      logLH += np.log(vec_norm)
      seq_rep.variable_profile[pos] = VarPos(vec/vec_norm, p.state)
      print(pos, vec/vec_norm)
    for state in alphabet:
      vec = seq_rep.fixed_ingroup[state]*gtrs[si].Pi
      vec_norm = vec.sum()
      logLH += np.log(vec_norm)*seq_rep.state_composition[state]
      seq_rep.fixed_profile[state] = vec/vec_norm
  return logLH

if __name__=="__main__":
  fname_nwk = 'data/ebola/ebola.nwk'
  fname_seq = 'data/ebola/ebola_dna.fasta'
  fname_nwk = 'test_scripts/data/tree.nwk'
  fname_seq = 'test_scripts/data/sequences.fasta'
  with open(fname_nwk) as fh:
    nwkstr = fh.read()
  G = graph_from_nwk_str(nwk_string=nwkstr, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  name_to_key={n.payload().name: n.key() for n in G.get_leaves()}
  aln = {s.id:np.fromiter(str(s.seq).upper(), 'U1') for s in AlignIO.read(fname_seq, 'fasta')}

  for sname, seq in aln.items():
    G.get_node(name_to_key[sname]).payload().full_seq.append(np.copy(seq))
    G.get_node(name_to_key[sname]).payload().seq_len.append(len(seq))

  generate_sparse_sequence_representation(G)
  gtr = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')

  ingroup_profiles(G, [gtr])
  print(calculate_root_state(G, [gtr]))

  # check output
  for leaf in G.leaves:
    node = G.get_node(leaf)
    nname = node.payload().name
    rec_seq = reconstruct_sequence(G, node)[0]
    print(nname, np.sum(rec_seq!=aln[nname]))


# how to avoid carrying along the sequence

# - pre-compute nucleotide composition
# - positions become variable only if there is a mutation somewhere
# - save the parsimony state along with the variable profiles
# - when collecting variable positions from children, infer the parent parsimony state from mutations

