from dataclasses import dataclass
from typing import Optional, List
from lib import Graph, GraphNodeForward, AutoRepr, GraphNodeBackward, Mut


@dataclass
class NodePayload(AutoRepr):
  name: Optional[str]
  muts: List[Mut]


@dataclass
class EdgePayload(AutoRepr):
  name: Optional[str]
  weight: Optional[float]


def main():
  graph = Graph()
  a = graph.add_node(NodePayload("A", [Mut.from_str('A123C')]))
  b = graph.add_node(NodePayload("B", [Mut.from_str('A123C')]))
  c = graph.add_node(NodePayload("C", [Mut.from_str('A123C')]))
  d = graph.add_node(NodePayload("D", [Mut.from_str('A123C')]))
  e = graph.add_node(NodePayload("E", [Mut.from_str('A123C')]))
  graph.add_edge(a, b, EdgePayload("A->B", 0.1))
  graph.add_edge(a, c, EdgePayload("A->C", 0.2))
  graph.add_edge(b, d, EdgePayload("B->D", 0.3))
  graph.add_edge(c, e, EdgePayload("C->E", 0.4))
  graph.build()

  result = []

  def explorer(node: GraphNodeForward):
    result.append(node.payload)

  graph.par_iter_forward(explorer)

  print(result)

  def explorer(node: GraphNodeBackward):
    result.append(node.payload)

  graph.par_iter_backward(explorer)

  print(result)


if __name__ == '__main__':
  main()
