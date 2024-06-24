import unittest
from .. import Graph, GraphNodeForward, GraphNodeBackward


class TestGraph(unittest.TestCase):
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

  def test_par_iter_forward(self):
    result = []

    def explorer(node: GraphNodeForward):
      result.append(node.payload)

    self.graph.par_iter_forward(explorer)
    expected = ["A", "B", "D", "C", "E"]
    self.assertEqual(result, expected)

  def test_par_iter_backward(self):
    result = []

    def explorer(node: GraphNodeBackward):
      result.append(node.payload)

    self.graph.par_iter_backward(explorer)
    expected = ["D", "B", "E", "C", "A"]
    self.assertEqual(result, expected)


if __name__ == '__main__':
  unittest.main()
