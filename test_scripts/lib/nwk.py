from typing import Optional, Callable
from io import StringIO
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from . import Graph, N, E, GraphNodeKey, Node, Edge

CreateNodePayloadFunction = Callable[[Optional[str]], N]
CreateEdgePayloadFunction = Callable[[Optional[float]], E]


def graph_from_nwk_str(
  nwk_string: str,
  node_payload_factory: CreateNodePayloadFunction,
  edge_payload_factory: CreateEdgePayloadFunction,
) -> Graph[N, E]:
  tree = Phylo.read(StringIO(nwk_string), "newick")

  graph = Graph()

  def recurse(tree_node: Clade) -> GraphNodeKey:
    node_key = graph.add_node(node_payload_factory(name=tree_node.name))
    for child in tree_node.clades:
      child_key = recurse(child)
      graph.add_edge(node_key, child_key, edge_payload_factory(child.branch_length))
    return node_key

  root = tree.root
  recurse(root)
  graph.build()

  return graph


def graph_to_nwk_str(graph: Graph[N, E]) -> str:
  def recurse(node: Node[N], edge: Optional[Edge[E]]) -> Clade:
    name = node.payload().name
    branch_length = edge.payload().weight if edge is not None else None
    children = graph.children_of(node.key())
    clades = []
    if children:
      clades = [recurse(child, child_edge) for child, child_edge in children]
    return Clade(name=name, branch_length=branch_length, clades=clades)

  root = graph.get_one_root()
  tree = recurse(root, None)

  with StringIO() as f:
    Phylo.write(tree, f, "newick")
    return f.getvalue().strip()
