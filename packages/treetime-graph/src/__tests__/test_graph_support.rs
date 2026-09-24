#[cfg(test)]
pub(crate) mod tests {
  use crate::edge::GraphEdgeKey;
  use crate::graph::Graph;
  use crate::node::GraphNodeKey;
  use eyre::Report;
  use itertools::Itertools;
  use std::collections::BTreeMap;
  use treetime_utils::make_report;

  pub(crate) struct NamedGraph {
    pub(crate) graph: Graph,
    keys: BTreeMap<&'static str, GraphNodeKey>,
  }

  impl NamedGraph {
    pub(crate) fn new(nodes: &[&'static str], edges: &[(&'static str, &'static str)]) -> Result<Self, Report> {
      let mut graph = Graph::new();
      let keys: BTreeMap<_, _> = nodes.iter().map(|name| (*name, graph.add_node())).collect();
      for (source, target) in edges {
        graph.add_edge(keys[source], keys[target])?;
      }
      graph.build()?;
      Ok(Self { graph, keys })
    }

    pub(crate) fn key(&self, name: &str) -> GraphNodeKey {
      self.keys[name]
    }

    pub(crate) fn name(&self, key: GraphNodeKey) -> &'static str {
      self
        .keys
        .iter()
        .find_map(|(name, node_key)| (*node_key == key).then_some(*name))
        .unwrap_or_else(|| panic!("Node {key} has no name"))
    }

    pub(crate) fn edge(&self, source: &str, target: &str) -> Result<GraphEdgeKey, Report> {
      let (source_key, target_key) = (self.key(source), self.key(target));
      self
        .graph
        .get_edges()
        .find(|edge| edge.source() == source_key && edge.target() == target_key)
        .map(|edge| edge.key())
        .ok_or_else(|| make_report!("Edge {source}->{target} not found"))
    }

    pub(crate) fn children(&self, name: &str) -> Vec<&'static str> {
      let node = self.graph.get_node(self.key(name)).expect("node exists");
      self
        .graph
        .children_keys_of(node)
        .map(|(child_key, _)| self.name(child_key))
        .collect_vec()
    }
  }
}
