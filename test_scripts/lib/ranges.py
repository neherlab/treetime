from dataclasses import dataclass
from .str import AutoRepr
from typing import List

@dataclass
class Range(AutoRepr):
  start: int
  end: int

  def __hash__(self) -> int:
    return tuple.__hash__((self.start, self.end))

# there are probably existing interval collections available.
class RangeCollection(AutoRepr):
  def __init__(self, ranges: List[Range]) -> None:
    self.ranges = sorted(ranges, key=lambda x:x.start)

  def contains(self, pos: int) -> bool:
    for r in self.ranges:
      if r.start<=pos and pos<r.end:
        return True
      elif r.start>pos:
        return False
    return False

  def add(self, r: Range):
    # should check overlap
    self.ranges.append(r)
    self.ranges.sort(key=lambda x:x.start)

  def __len__(self) -> int:
    return len(self.ranges)

  def __getitem__(self, i) -> Range:
    return self.ranges[i]

  # could add an interator method
  # also a position iterator (for r in ranges: for p in range(r.start, r.end):)


def RangeCollection_complement(ranges: RangeCollection, global_start: int, global_end: int) -> RangeCollection:
  if len(ranges)==0:
    return RangeCollection([Range(global_start, global_end)])

  new_ranges = []
  # add part before the first range
  if ranges[0].start>global_start:
    new_ranges.append(Range(global_start, ranges[0].start))

  for ri in range(len(ranges)-1):
    r1, r2 = ranges[ri], ranges[ri+1]

    if r2.start<global_start or r1.end>global_end:
      # outside the range
      continue
    elif r1.end<global_start:
      # starts before the global_start
      if r2.start<global_end:
        #ends in the range
        new_ranges.append(Range(global_start, r2.start))
      else:
        new_ranges.append(Range(global_start, global_end))
    else:
      # starts inside the range
      if r2.start<global_end and r1.end<r2.start:
        #ends in the range
        new_ranges.append(Range(r1.end, r2.start))
      else:
        new_ranges.append(Range(r1.end, global_end))
        # can stop here since ranges are sorted.
        return RangeCollection(new_ranges)

  # add final complement extending to the end of the range
  if ranges[-1].end<global_end:
    new_ranges.append(Range(ranges[-1].end, global_end))

  return RangeCollection(new_ranges)

def RangeCollection_intersection(range_collections: List[RangeCollection]) -> RangeCollection:
  '''
  Note, this assumes sorted ranges
  '''
  if any([len(r)==0 for r in range_collections]):
    return RangeCollection([])

  current_ranges = RangeCollection(range_collections[0])
  for next_ranges in range_collections[1:]:
    new_ranges = []
    # index walker variables
    ri1, ri2 = 0, 0
    r1 = current_ranges[ri1]
    r2 = next_ranges[ri2]
    while ri1<len(current_ranges) and ri2<len(next_ranges):
      if r2.start>r1.end:
        ri1 += 1
        if ri1<len(current_ranges):
          r1 = current_ranges[ri1]
      elif r1.start>r2.end:
        ri2 += 1
        if ri2<len(next_ranges):
          r2 = next_ranges[ri2]
      else:
        new_range = Range(max(r1.start, r2.start), min(r1.end, r2.end))
        if new_range.start<new_range.end:
          new_ranges.append(new_range)
        if r1.end<r2.end:
          ri1 += 1
          if ri1<len(current_ranges):
            r1 = current_ranges[ri1]
        else:
          ri2 += 1
          if ri2<len(next_ranges):
            r2 = next_ranges[ri2]
    current_ranges = RangeCollection(new_ranges)

  return current_ranges

def RangeCollection_difference(rc1, rc2) -> RangeCollection:
  if len(rc1.ranges)==0:
    return rc1
  # restrict the complement to the range of rc1 -- could be entire sequence as well
  rc2_comp = RangeCollection_complement(rc2, global_start=rc1.ranges[0].start,
                                             global_end=rc1.ranges[-1].end)
  return RangeCollection_intersection([rc1, rc2_comp])


def find_char_ranges(seq, char):
  ranges = []
  start = None
  for pos, nuc in enumerate(seq):
    if nuc != char and (start is not None):
      ranges.append(Range(start, pos))
      start = None
    elif nuc == char and start is None:
      start = pos
  if start:
    ranges.append(Range(start, len(seq)))
  return RangeCollection(ranges)

