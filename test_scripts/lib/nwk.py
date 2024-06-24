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

  def add_node(tree_node: Clade) -> GraphNodeKey:
    node_key = graph.add_node(node_payload_factory(tree_node.name))
    for child in tree_node.clades:
      child_key = add_node(child)
      graph.add_edge(node_key, child_key, edge_payload_factory(tree_node.branch_length))
    return node_key

  root = tree.root
  add_node(root)
  graph.build()

  return graph


def graph_to_nwk_string(graph: Graph[N, E]) -> str:
  def build_phylo_tree(node: Node[N], edge: Optional[Edge[E]]) -> Clade:
    if node.is_leaf():
      name = node.payload().name
      branch_length = edge.payload().weight
      return Clade(name=name, branch_length=branch_length)
    else:
      children = graph.children_of(node.key())
      clades = [build_phylo_tree(child, edge) for child, edge in children]
      return Clade(clades=clades)

  root = graph.get_one_root()
  tree = build_phylo_tree(root, None)

  with StringIO() as f:
    Phylo.write(tree, f, "newick")
    return f.getvalue().strip()
