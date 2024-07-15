from typing import Optional, List, Dict
from lib import Graph, AutoRepr, Mut, InDel, Node, SequenceInfo
from dataclasses import dataclass


@dataclass
class NodePayload(AutoRepr):
  name: Optional[str]
  date: Optional[float]
  full_seq: Optional[List[str]]
  SeqInfo: List[SequenceInfo]
  SeqInfo_ingroup:  Optional[List[List[SequenceInfo]]]  #these are calculated for each parent, i.e. the inbound edges
  SeqInfo_outgroup: Optional[List[List[SequenceInfo]]] #these are calculated for each child, i.e. the outbound edges


@dataclass
class EdgePayload():
  '''
  - length of the branch
  - mutations (character changes in the parsimony based description of the sequence)
  - indels (stretches of sequence that are lost of gained)
  '''
  branch_length: Optional[float]
  muts: List[Mut]
  indels: List[InDel]


