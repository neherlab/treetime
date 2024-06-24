import unittest
from typing import Optional

from ..str import AutoRepr
from ..nwk import graph_to_nwk_str, graph_from_nwk_str


class NodePayload(AutoRepr):
  def __init__(self, name: Optional[str]):
    self.name: Optional[str] = name


class EdgePayload(AutoRepr):
  def __init__(self, weight: Optional[float]):
    self.weight: Optional[float] = weight


class TestNewickConversion(unittest.TestCase):

  def test_conversion_roundtrip(self):
    nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;"
    graph_from_nwk = graph_from_nwk_str(nwk, NodePayload, EdgePayload)
    restored = graph_to_nwk_str(graph_from_nwk)
    self.assertEqual(nwk, restored)


if __name__ == "__main__":
  unittest.main()
