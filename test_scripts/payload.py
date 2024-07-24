from typing import Optional, List
from lib import AutoRepr, SparseSeqNode, SparseSeqEdge, DenseSeqEdge, DenseSeqNode
from dataclasses import dataclass, field

@dataclass
class NodePayload(AutoRepr):
  name: str
  date: Optional[float] = None
  sparse_sequences: List[SparseSeqNode] = field(default_factory=list)
  dense_sequences:  List[DenseSeqNode] = field(default_factory=list)
  # discrete_traits:  List[DiscreteTraitNode] = field(default_factory=list)

@dataclass
class EdgePayload():
  '''
  - length of the branch
  - mutations (character changes in the parsimony based description of the sequence)
  - indels (stretches of sequence that are lost of gained)
  '''
  branch_length: float = 0.0

  sparse_sequences: List[SparseSeqEdge] = field(default_factory=list)
  dense_sequences:  List[DenseSeqEdge] = field(default_factory=list)
  # discrete_traits: List[DiscreteTraitEdge] = field(default_factory=list)


