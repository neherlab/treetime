import unittest

from .. import Graph
from ..json import graph_to_json, graph_from_json


class TestGraphSerialization(unittest.TestCase):
  def setUp(self):
    self.graph = Graph()
    self.a = self.graph.add_node("A")
    self.b = self.graph.add_node("B")
    self.node_c = self.graph.add_node("C")
    self.node_d = self.graph.add_node("D")
    self.node_e = self.graph.add_node("E")
    self.graph.add_edge(self.a, self.b, "A->B")
    self.graph.add_edge(self.a, self.node_c, "A->C")
    self.graph.add_edge(self.b, self.node_d, "B->D")
    self.graph.add_edge(self.node_c, self.node_e, "C->E")
    self.graph.build()

  def test_graph_roundtrip(self):
    json_str = graph_to_json(self.graph)

    new_graph = graph_from_json(json_str)

    self.assertEqual(len(new_graph.nodes), len(self.graph.nodes))
    self.assertEqual(len(new_graph.edges), len(self.graph.edges))
    self.assertEqual(len(new_graph.roots), len(self.graph.roots))
    self.assertEqual(len(new_graph.leaves), len(self.graph.leaves))

    self.assertEqual(new_graph.get_node_payloads(), self.graph.get_node_payloads())
