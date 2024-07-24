from typing import Optional, List, Dict
from dataclasses import dataclass, field
from lib import AutoRepr, RangeCollection, InDel
import numpy as np


@dataclass
class DenseSeqInfo(AutoRepr):
  '''
  This class is meant to contain information related to the parts of a sequence necessary for likelihood calculation.
  This information includes
  - gaps (ranges)
  '''
  gaps: RangeCollection = field(default_factory=RangeCollection)  # save gaps in the sequence as ranges
  sequence: Optional[np.array] = None

@dataclass
class DenseSeqDis(AutoRepr):
  logLH: float = 0.0
  dis: np.array = field(default_factory=lambda: np.array([]))

@dataclass
class DenseSeqNode(AutoRepr):
  seq: DenseSeqInfo
  profile: DenseSeqDis = field(default_factory=DenseSeqDis)
  msg_to_parents:    DenseSeqDis = field(default_factory=DenseSeqDis)      # there might be multiple parents, but all parents only see info from children
  msgs_to_children:   Dict[str, DenseSeqDis] = field(default_factory=dict)
  msgs_from_children: Dict[str, DenseSeqDis] = field(default_factory=dict)


@dataclass
class DenseSeqEdge():
  indels: List[InDel] = field(default_factory=list)
  transmission: Optional[RangeCollection] = None



