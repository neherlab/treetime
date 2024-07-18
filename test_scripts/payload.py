from typing import Optional, List, Dict, Tuple
from lib import Graph, AutoRepr, Mut, InDel, Node, SeqInfoLh, SeqInfoParsimony, RangeCollection, SparseSeqNode, SparseSeqEdge
from dataclasses import dataclass, field
import numpy as np

@dataclass
class NodePayload(AutoRepr):
  name: str
  date: Optional[float] = None
  sparse_sequences: List[SparseSeqNode] = field(default_factory=list)
  # dense_sequences:  List[DenseSeqNode] = []
  # discrete_traits:  List[DiscreteTraitNode] = []

  # full_seq: Optional[List[str]] = None

  # fitch: Optional[List[SeqInfoParsimony]] = None

  # # there is one for each sequence partition. Maybe make plural obvious?
  # seq_info: Optional[List[SeqInfoLh]] = None
  # seq_info_ingroup:  Optional[List[SeqInfoLh]] = None  #these are calculated for the node from its children
  # seq_info_outgroup: Optional[List[SeqInfoLh]] = None  #these are calculated for each child, i.e. the outbound edges
  # child_messages: Optional[List[List[Tuple[str, SeqInfoLh]]]] = None # messages received from children. Could move to edge
  # parent_messages: Optional[List[List[Tuple[str, SeqInfoLh]]]] = None # messages received from parents. Could move to edge

@dataclass
class EdgePayload():
  '''
  - length of the branch
  - mutations (character changes in the parsimony based description of the sequence)
  - indels (stretches of sequence that are lost of gained)
  '''
  branch_length: float = 0.0

  sparse_sequences: List[SparseSeqEdge] = field(default_factory=list)
  # dense_sequences:  List[DenseSeqEdge] = []
  # discrete_traits: List[DiscreteTraitEdge] = []


