from dataclasses import dataclass
from .str import AutoRepr
from typing import List, Dict
from .mut import Mut, InDel
from .range import RangeCollection
import numpy as np

@dataclass
class VarPos(AutoRepr):
  profile: np.array  # array of floats of size 'alphabet'
  state: str

@dataclass
class SequenceInfo(AutoRepr):
    '''
    This class is meant to contain information related to the parts of a sequence necessary for likelihood calculation.
    This information includes
    - unknown (ranges)
    - variable (probability vector for each variable position collecting information from children)
    - fixed (probability vector for the state of fixed positions based on information from children)
    - fixed_composition (this is important for likelihood calculations)
    - gaps (ranges)
    '''
    unknown: RangeCollection
    gaps: RangeCollection
    variable: Dict[int, VarPos]
    fixed: Dict[str, np.array]
    fixed_composition: Dict[str, int]



