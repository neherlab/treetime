from dataclasses import dataclass
from typing import List, Optional

@dataclass
class Mut:
  ref: str
  pos: int
  qry: str

  def __repr__(self) -> str:
    return f"'{self.ref}{self.pos + 1}{self.qry}'"

  @classmethod
  def from_str(cls, mut_str):
    if len(mut_str) < 3:
      raise ValueError(f"Invalid mutation: '{mut_str}'")

    ref = mut_str[0]

    try:
      pos = int(mut_str[1:-1]) - 1
    except ValueError:
      raise ValueError(f"Invalid mutation: '{mut_str}': position must be an integer")

    qry = mut_str[-1]

    self = cls(ref, pos, qry)

    if not self.is_mut():
      raise ValueError(f"Invalid mutation: '{mut_str}': query and ref must be different")

    return self

  def is_mut(self):
    return self.qry != self.ref

  def is_del(self):
    return self.qry == '-'

  def __gt__(self, other):
    if self.pos>other.pos: return True
    elif self.pos==other.pos:
      if self.ref>other.ref:
        return True
      if self.ref==other.ref:
        return self.qry>self.ref
    return False

@dataclass
class InDel:
  pos: int
  length: int
  seq: Optional[list[str]]
  deletion: bool = True # deletion if True, insertion if False
