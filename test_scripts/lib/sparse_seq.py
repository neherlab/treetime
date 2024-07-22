from typing import Optional, List, Dict
from dataclasses import dataclass, field
from .mut import Mut, InDel
from .range import RangeCollection, Range
from .str import AutoRepr
import numpy as np

@dataclass
class VarPos(AutoRepr):
  profile: np.array  # array of floats of size 'alphabet'
  state: str

@dataclass
class Deletion(AutoRepr):
  deleted: int = 0  # number of times deletion is observed
  ins: int = 0   # or not
  alt: str = ''

@dataclass
class SparseSeqDis(AutoRepr):
  '''
  - variable (probability vector for each variable position collecting information from children)
  - fixed (probability vector for the state of fixed positions based on information from children)
  - logLH (total_LH)
  '''
  origin: str = None # needed to filter sibling messages
  variable: Dict[int, VarPos] = field(default_factory=dict)
  variable_indel: Dict[Range, Deletion] = field(default_factory=dict)
  fixed: Dict[str, np.array] = field(default_factory=dict)
  fixed_counts: Dict[str, int] = field(default_factory=dict)
  logLH: float = 0.0

@dataclass
class SparseSeqInfo(AutoRepr):
  '''
  This class is meant to contain information related to the parts of a sequence necessary for likelihood calculation.
  This information includes
  - unknown (ranges)
  - gaps (ranges)
  - non_char (ranges): any position that does not evolve according to the substition model, i.e. gap or N
  '''
  unknown: RangeCollection
  gaps: RangeCollection
  non_char: RangeCollection
  composition: Dict[str, int] = field(default_factory=dict) # count of all characters in the region that is not `non_char`
  sequence: Optional[np.array] = None
  distribution: SparseSeqDis = field(default_factory=SparseSeqDis)

@dataclass
class SparseSeqNode(AutoRepr):
  state: SparseSeqInfo
  msgs_to_parents:    SparseSeqDis       # there might be multiple parents, but all parents only see info from children
  msgs_from_parents:  List[SparseSeqDis] = field(default_factory=list) # lists correspond to lists of children or parents
  msgs_to_children:   List[SparseSeqDis] = field(default_factory=list)
  msgs_from_children: List[SparseSeqDis] = field(default_factory=list)


@dataclass
class SparseSeqEdge():
  muts:   List[Mut] = field(default_factory=list)
  indels: List[InDel] = field(default_factory=list)
  transmission: Optional[RangeCollection] = None


