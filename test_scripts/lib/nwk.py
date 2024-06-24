from typing import Optional, Callable
from io import StringIO
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from . import Graph, N, E, GraphNodeKey


def graph_from_nwk_str(
  nwk_string: str,
  node_payload_factory: Callable[[Optional[str], Optional[float]], N],
  edge_payload_factory: Callable[[Optional[float]], E],
) -> Graph[N, E]:
  tree = Phylo.read(StringIO(nwk_string), "newick")

  graph = Graph()

  def add_node(tree_node: Clade) -> GraphNodeKey:
    node_key = graph.add_node(node_payload_factory(tree_node.name, tree_node.branch_length))
    for child in tree_node.clades:
      child_key = add_node(child)
      graph.add_edge(node_key, child_key, edge_payload_factory(tree_node.branch_length))
    return node_key

  root = tree.root
  add_node(root)
  graph.build()

  return graph


def graph_to_nwk_string(graph: Graph[N, E]) -> str:
  def build_phylo_tree(node_key: GraphNodeKey) -> Clade:
    node = graph.get_node(node_key)
    if node.is_leaf():
      name = node.payload().name
      branch_length = node.payload().weight
      return Clade(name=name, branch_length=branch_length)
    else:
      children = graph.child_keys_of(node_key)
      clades = [build_phylo_tree(child_key) for child_key in children]
      return Clade(clades=clades)

  root_key = graph.get_one_root().key()
  tree = build_phylo_tree(root_key)

  with StringIO() as f:
    Phylo.write(tree, f, "newick")
    return f.getvalue().strip()
