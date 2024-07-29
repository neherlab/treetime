from dataclasses import dataclass, field
from .str import AutoRepr
from typing import Dict, Optional
import numpy as np

@dataclass
class ClockSet(AutoRepr):
    t_sum: float = 0.0
    tsq_sum: float = 0.0
    d_sum: float = 0.0
    dsq_sum: float = 0.0
    dt_sum: float = 0.0
    norm: float = 0.0

    def add(self, other):
        self.t_sum += other.t_sum
        self.tsq_sum += other.tsq_sum
        self.d_sum += other.d_sum
        self.dsq_sum += other.dsq_sum
        self.dt_sum += other.dt_sum
        self.norm += other.norm

    def subtract(self, other):
        self.t_sum -= other.t_sum
        self.tsq_sum -= other.tsq_sum
        self.d_sum -= other.d_sum
        self.dsq_sum -= other.dsq_sum
        self.dt_sum -= other.dt_sum
        self.norm -= other.norm


@dataclass
class ClockData(AutoRepr):
    date: Optional[float] = None
    total: ClockSet = field(default_factory=ClockSet)
    to_parent: ClockSet = field(default_factory=ClockSet)
    to_children: Dict[str, ClockSet] =  field(default_factory=dict)
    from_children: Dict[str, ClockSet] =  field(default_factory=dict)

@dataclass
class ClockModel(AutoRepr):
    rate: float
    intercept: float
    hessian: np.array
    chisq: float = 0.0
