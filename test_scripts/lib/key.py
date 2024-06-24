class GraphNodeKey:
  def __init__(self, index: int):
    self.index: int = index

  def to_json(self):
    return self.index

  def __repr__(self):
    return self.index.__repr__()


class GraphEdgeKey:
  def __init__(self, index: int):
    self.index: int = index

  def to_json(self):
    return self.index

  def __repr__(self):
    return self.index.__repr__()
