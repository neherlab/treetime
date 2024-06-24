from typing import TypeVar, Generic, List

from .key import GraphNodeKey, GraphEdgeKey

N = TypeVar('N')


class Node(Generic[N]):
  def __init__(self, key: GraphNodeKey, payload: N):
    self._key: GraphNodeKey = key
    self._payload: N = payload
    self._inbound: List[GraphEdgeKey] = []
    self._outbound: List[GraphEdgeKey] = []

  def key(self) -> GraphNodeKey:
    return self._key

  def payload(self) -> N:
    return self._payload

  def is_root(self) -> bool:
    return len(self._inbound) == 0

  def is_leaf(self) -> bool:
    return len(self._outbound) == 0

  def inbound(self) -> List[GraphEdgeKey]:
    return self._inbound

  def outbound(self) -> List[GraphEdgeKey]:
    return self._outbound

  def add_inbound(self, edge_key: GraphEdgeKey) -> None:
    self._inbound.append(edge_key)

  def add_outbound(self, edge_key: GraphEdgeKey) -> None:
    self._outbound.append(edge_key)
