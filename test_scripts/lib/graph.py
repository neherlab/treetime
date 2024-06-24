from typing import List, Tuple, Callable, Any, Optional, TypeVar, Generic

from .key import GraphEdgeKey, GraphNodeKey
from .node import Node, N
from .edge import Edge, E


class GraphNodeForward:
  """
  Subset of node data which is safe to access during forward parallel traversal
  """

  def __init__(self, graph: 'Graph', node: Node[N]):
    self.is_root: bool = node.is_root()
    self.is_leaf: bool = node.is_leaf()
    self.key: GraphNodeKey = node.key()
    self.payload: N = node.payload()
    self.parents: List[Tuple[N, E]] = [
      (parent.payload(), edge.payload())
      for parent, edge in graph.parents_of(self.key)
    ]


class GraphNodeBackward:
  """
  Subset of node data which is safe to access during backward parallel traversal
  """

  def __init__(self, graph: 'Graph', node: Node[N]):
    self.is_root: bool = node.is_root()
    self.is_leaf: bool = node.is_leaf()
    self.key: GraphNodeKey = node.key()
    self.payload: N = node.payload()
    self.children: List[Tuple[N, E]] = [
      (child.payload(), edge.payload())
      for child, edge in graph.children_of(self.key)
    ]


class GraphNodeSafe:
  """
  Subset of node data which is safe to access during unordered parallel iteration.
  """

  def __init__(self, graph: 'Graph', node: Node[N]):
    self.is_root: bool = node.is_root()
    self.is_leaf: bool = node.is_leaf()
    self.key: GraphNodeKey = node.key()
    self.payload: N = node.payload()
    self.children: List[Tuple[Node[N], Edge[E]]] = graph.children_of(self.key)
    self.parents: List[Tuple[Node[N], Edge[E]]] = graph.parents_of(self.key)


T = TypeVar('T')


class Graph(Generic[N, E]):
  def __init__(self):
    self.nodes: List[Node[N]] = []
    self.edges: List[Edge[E]] = []
    self.roots: List[GraphNodeKey] = []
    self.leaves: List[GraphNodeKey] = []

  def parents_of(self, node_key: GraphNodeKey) -> List[Tuple[Node[N], Edge[E]]]:
    node: Node[N] = self.get_node(node_key)
    parents_with_edges = []
    for edge_key in node.inbound():
      edge: Edge[E] = self.get_edge(edge_key)
      if edge:
        parent_node = self.get_node(edge.source())
        if parent_node:
          parents_with_edges.append((parent_node, edge))

    return parents_with_edges

  def parent_keys_of(self, node_key: GraphNodeKey) -> List[GraphNodeKey]:
    node: Node[N] = self.get_node(node_key)
    parent_keys = []
    for edge_key in node.inbound():
      edge: Edge[E] = self.get_edge(edge_key)
      if edge:
        parent_keys.append(edge.source())

    return parent_keys

  def one_parent_of(self, node_key: GraphNodeKey) -> Node[N]:
    parents = self.parents_of(node_key)
    if len(parents) > 1:
      raise NotImplemented(f"Only one parent per node is currently supported, but found: {len(parents)}")
    return parents[0][0]

  def children_of(self, node_key: GraphNodeKey) -> List[Tuple[Node[N], Edge[E]]]:
    node: Node[N] = self.get_node(node_key)
    children_with_edges = []
    for edge_key in node.outbound():
      edge = self.get_edge(edge_key)
      if edge:
        child_node = self.get_node(edge.target())
        if child_node:
          children_with_edges.append((child_node, edge))

    return children_with_edges

  def child_keys_of(self, node_key: GraphNodeKey) -> List[GraphNodeKey]:
    node: Node[N] = self.get_node(node_key)
    child_keys = []
    for edge_key in node.outbound():
      edge: Edge[E] = self.get_edge(edge_key)
      if edge:
        child_keys.append(edge.target())

    return child_keys

  def get_node(self, node_key: GraphNodeKey) -> Node[N]:
    index = node_key.index
    if index < 0 or index >= len(self.nodes):
      raise KeyError(f"When looking for node '{index}': index out of range")
    return self.nodes[index]

  def get_edge(self, edge_key: GraphEdgeKey) -> Edge[E]:
    index = edge_key.index
    if index < 0 or index >= len(self.edges):
      raise KeyError(f"When looking for edge '{index}': index out of range")
    return self.edges[index]

  def for_each(self, func: Callable[['GraphNodeSafe'], None]) -> None:
    for node in self.nodes:
      func(GraphNodeSafe(self, node))

  def map(self, func: Callable[['GraphNodeSafe'], N]) -> List[T]:
    return [func(GraphNodeSafe(self, node)) for node in self.nodes]

  def filter_map(self, func: Callable[['GraphNodeSafe'], Optional[T]]) -> List[T]:
    return [func(GraphNodeSafe(self, node)) for node in self.nodes if func(GraphNodeSafe(self, node)) is not None]

  def num_nodes(self) -> int:
    return len(self.nodes)

  def num_roots(self) -> int:
    return len(self.roots)

  def num_leaves(self) -> int:
    return len(self.leaves)

  def get_nodes(self) -> List[Node[N]]:
    return self.nodes

  def get_node_payloads(self) -> List[N]:
    return [node.payload() for node in self.nodes]

  def get_roots(self) -> List[Node[N]]:
    return [self.get_node(idx) for idx in self.roots]

  def get_one_root(self) -> Node[N]:
    roots = self.get_roots()
    if len(roots) > 1:
      raise NotImplemented(f"Only one root is currently supported, but found: {len(roots)}")
    return roots[0]

  def get_leaves(self) -> List[Node[N]]:
    return [self.get_node(idx) for idx in self.leaves]

  def get_edges(self) -> List[Edge[E]]:
    return self.edges

  def add_node(self, node_payload: N) -> GraphNodeKey:
    node_key = GraphNodeKey(len(self.nodes))
    self.nodes.append(Node[N](node_key, node_payload))
    return node_key

  def add_edge(self, source_key: GraphNodeKey, target_key: GraphNodeKey, edge_payload: Any) -> None:
    if source_key == target_key:
      raise ValueError("Cannot connect node to itself")
    source_node = self.get_node(source_key)
    target_node = self.get_node(target_key)
    if not source_node or not target_node:
      raise ValueError("Source or target node not found")
    edge_key = GraphEdgeKey(len(self.edges))
    new_edge = Edge[E](edge_key, source_key, target_key, edge_payload)
    if any(key == target_key for key in source_node.outbound()):
      raise ValueError("Nodes are already connected")
    self.edges.append(new_edge)
    source_node.add_outbound(edge_key)
    target_node.add_inbound(edge_key)

  def build(self) -> None:
    self.roots = [node.key() for node in self.nodes if node.is_root()]
    self.leaves = [node.key() for node in self.nodes if node.is_leaf()]

  def par_iter_forward(self, visitor: Callable[['GraphNodeForward'], None]) -> None:
    """
    Emulates interface of parallel forward traversal (from roots to leaves).

    Should not expect any particular order - there is no order in parallel traversal. The only guarantee is that
    all nodes are visited after all of their parents have been already visited.
    """

    # The fact that this is a preorder traversal is a coincidence, not to be relied on
    return self._iter_depth_first_preorder(lambda node: visitor(GraphNodeForward(self, node)))

  def par_iter_backward(self, visitor: Callable[['GraphNodeBackward'], None]) -> None:
    """
    Emulates interface of parallel backward traversal (from leaves to roots).

    Should not expect any particular order - there is no order in parallel traversal. The only guarantee is that
    all nodes are visited after all of their children have been already visited.
    """

    # The fact that this is a postorder traversal is a coincidence, not to be relied on
    return self._iter_depth_first_postorder(lambda node: visitor(GraphNodeBackward(self, node)))

  def _iter_depth_first_preorder(self, visitor: Callable[[Node], None]) -> None:
    """
    Should not use these. But if **really** needed, then here they are.
    """
    stack: List[Node[N]] = self.get_roots()
    visited = set()
    while stack:
      node = stack.pop()
      if node.key() not in visited:
        visitor(node)
        visited.add(node.key())
        children = reversed(list(child for child, _ in self.children_of(node.key())))
        stack.extend(children)

  def _iter_depth_first_postorder(self, visitor: Callable[[Node], None]) -> None:
    """
    Should not use these. But if **really** needed, then here they are.
    """
    stack: List[Tuple[Node[N], bool]] = [(root, False) for root in self.get_roots()]
    visited = set()
    while stack:
      node, visited_flag = stack.pop()
      if node.key() in visited:
        continue

      if visited_flag:
        visitor(node)
        visited.add(node.key())
      else:
        stack.append((node, True))
        children = reversed(list(child for child, _ in self.children_of(node.key())))
        stack.extend((child, False) for child in children if child.key() not in visited)
