from typing import Optional, List, Dict
from dataclasses import dataclass, field
from .mut import Mut, InDel
from .ranges import RangeCollection, Range
from .str import AutoRepr
import numpy as np

@dataclass
class VarPos(AutoRepr):
  dis: np.array  # array of floats of size 'alphabet'
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
  variable: Dict[int, VarPos] = field(default_factory=dict)
  variable_indel: Dict[Range, Deletion] = field(default_factory=dict)
  fixed: Dict[str, np.array] = field(default_factory=dict)
  fixed_counts: Dict[str, int] = field(default_factory=dict)
  logLH: float = 0.0

@dataclass
class FitchVar(AutoRepr):
  '''
  - variable (probability vector for each variable position collecting information from children)
  '''
  variable: Dict[int, VarPos] = field(default_factory=dict)
  variable_indel: Dict[Range, Deletion] = field(default_factory=dict)


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
  fitch: FitchVar = field(default_factory=FitchVar)

@dataclass
class SparseSeqNode(AutoRepr):
  seq: SparseSeqInfo
  profile: SparseSeqDis = field(default_factory=SparseSeqDis)
  msg_to_parents:    SparseSeqDis = field(default_factory=SparseSeqDis)      # there might be multiple parents, but all parents only see info from children
  msgs_to_children:   Dict[str, SparseSeqDis] = field(default_factory=dict)
  msgs_from_children: Dict[str, SparseSeqDis] = field(default_factory=dict)


@dataclass
class SparseSeqEdge():
  muts:   List[Mut] = field(default_factory=list)
  indels: List[InDel] = field(default_factory=list)
  transmission: Optional[RangeCollection] = None


