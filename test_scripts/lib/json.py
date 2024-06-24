import jsonpickle

from test_scripts.lib import Graph


def graph_to_json(self: Graph):
  return jsonpickle.encode(self, warn=True, indent=2)


def graph_from_json(json_str: str) -> Graph:
  return jsonpickle.decode(json_str)
