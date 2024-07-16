from dataclasses import dataclass
from .str import AutoRepr
from typing import List, Dict, Optional
from .range import RangeCollection
import numpy as np
from treetime import GTR

class SeqPartition(AutoRepr):
  def __init__(self, gtr: GTR, length: int, profile_map: Dict[str, np.array]):
    self.gtr=gtr
    self.length=length
    self._profile_map=profile_map
    self._reverse_profile_map = {tuple(p): nuc for nuc, p in self._profile_map.items()}

  def profile(self, s: str) -> np.array:
    # need to handle wrong input
    return self._profile_map[s]

  def code(self, p) -> str:
    # will and should only work for input profiles
    return self._reverse_profile_map[tuple(p)]


@dataclass
class VarPos(AutoRepr):
  profile: np.array  # array of floats of size 'alphabet'
  state: str

@dataclass
class SeqInfoLh(AutoRepr):
  '''
  This class is meant to contain information related to the parts of a sequence necessary for likelihood calculation.
  This information includes
  - unknown (ranges)
  - gaps (ranges)
  - non_nuc (ranges): any position that does not evolve according to the substition model, i.e. gap or N
  - variable (probability vector for each variable position collecting information from children)
  - fixed (probability vector for the state of fixed positions based on information from children)
  - fixed_composition (this is important for likelihood calculations)
  '''
  unknown: RangeCollection
  gaps: RangeCollection
  non_char: RangeCollection
  variable: Dict[int, VarPos]
  fixed: Optional[Dict[str, np.array]] = None
  fixed_composition: Optional[Dict[str, int]] = None

@dataclass
class SeqInfoParsimony(AutoRepr):
  '''
  This class is meant to contain information related to the parts of a sequence necessary for likelihood calculation.
  This information includes
  - unknown (ranges)
  - gaps (ranges)
  - non_nuc (ranges): any position that does not evolve according to the substition model, i.e. gap or N
  - variable (0/1 vector for each variable position collecting information from children)
  - fixed_composition (this is important for likelihood calculations)
  '''
  unknown: RangeCollection
  gaps: RangeCollection
  non_char: RangeCollection
  variable: Dict[int, VarPos]
  fixed_composition: Optional[Dict[str, int]] = None


