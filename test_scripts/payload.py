from typing import Optional, List, Dict, Tuple
from lib import Graph, AutoRepr, Mut, InDel, Node, SeqInfoLh, SeqInfoParsimony, RangeCollection
from dataclasses import dataclass


@dataclass
class NodePayload(AutoRepr):
  name: str
  date: Optional[float] = None
  full_seq: Optional[List[str]] = None
  fitch: Optional[List[SeqInfoParsimony]] = None

  # there is one for each sequence partition. Maybe make plural obvious?
  seq_info: Optional[List[SeqInfoLh]] = None
  seq_info_ingroup:  Optional[List[SeqInfoLh]] = None  #these are calculated for the node from its children
  seq_info_outgroup: Optional[List[SeqInfoLh]] = None  #these are calculated for each child, i.e. the outbound edges
  child_messages: Optional[List[List[Tuple[str, SeqInfoLh]]]] = None # messages received from children.
  parent_messages: Optional[List[List[Tuple[str, SeqInfoLh]]]] = None # messages received from parents.

@dataclass
class EdgePayload():
  '''
  - length of the branch
  - mutations (character changes in the parsimony based description of the sequence)
  - indels (stretches of sequence that are lost of gained)
  '''
  branch_length: float = 0.0
  muts:   Optional[List[Mut]] = None
  indels: Optional[List[InDel]] = None
  transmission: RangeCollection = None


