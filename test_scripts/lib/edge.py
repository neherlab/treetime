from typing import TypeVar, Generic

from .key import GraphEdgeKey, GraphNodeKey
from .str import AutoRepr

E = TypeVar('E')


class Edge(Generic[E], AutoRepr):
  def __init__(self, key: GraphEdgeKey, source: GraphNodeKey, target: GraphNodeKey, payload: E):
    self._key: GraphEdgeKey = key
    self._source: GraphNodeKey = source
    self._target: GraphNodeKey = target
    self._payload: E = payload

  def key(self) -> GraphEdgeKey:
    return self._key

  def source(self) -> GraphNodeKey:
    return self._source

  def target(self) -> GraphNodeKey:
    return self._target

  def payload(self) -> E:
    return self._payload
