#[cfg(test)]
pub(crate) mod tests {
  pub(crate) mod generators {
    use super::helpers::DagSpec;
    use itertools::Itertools;
    use proptest::prelude::Just;
    use proptest::sample::subsequence;
    use proptest::strategy::Strategy;

    pub(crate) fn gen_dag(max_nodes: usize) -> impl Strategy<Value = DagSpec> {
      (1..max_nodes)
        .prop_flat_map(|count| {
          let parents = (0..count)
            .map(|node| subsequence((0..node).collect_vec(), 0..=node.min(3)))
            .collect_vec();
          (parents, Just((0..count).collect_vec()).prop_shuffle())
        })
        .prop_flat_map(|(parents, keys)| {
          let edges = parents
            .iter()
            .enumerate()
            .flat_map(|(child, parents)| parents.iter().map(move |&parent| (parent, child)))
            .collect_vec();
          (Just(keys), Just(edges).prop_shuffle())
        })
        .prop_map(|(keys, edges)| DagSpec { keys, edges })
    }
  }

  pub(crate) mod helpers {
    use crate::edge::GraphEdgeKey;
    use crate::graph::Graph;
    use crate::key_map::NodeMap;
    use crate::node::GraphNodeKey;
    use crate::node_names::NodeNames;
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use std::iter::repeat_with;
    use treetime_utils::make_report;

    pub(crate) fn names_on_star(input: NodeMap<Option<String>>) -> Result<NodeNames, Report> {
      let size = input.keys().map(|key| key.as_usize() + 1).max().unwrap_or(1);
      let graph = graph_from_parents(&vec![0; size])?;
      Ok(NodeNames::from_input(input, &graph)?.names)
    }

    pub(crate) fn graph_from_parents(parents: &[usize]) -> Result<Graph, Report> {
      let mut graph = Graph::new();
      let keys = repeat_with(|| graph.add_node()).take(parents.len()).collect_vec();
      for (child, &parent) in parents.iter().enumerate().skip(1) {
        graph.add_edge(keys[parent], keys[child])?;
      }
      graph.build()?;
      Ok(graph)
    }

    pub(crate) fn graph_from_dag(dag: &DagSpec) -> Result<Graph, Report> {
      let mut graph = Graph::new();
      let created = repeat_with(|| graph.add_node()).take(dag.keys.len()).collect_vec();
      let key = |node: usize| created[dag.keys[node]];
      for &(parent, child) in &dag.edges {
        graph.add_edge(key(parent), key(child))?;
      }
      graph.build()?;
      Ok(graph)
    }

    #[derive(Clone, Debug)]
    pub(crate) struct DagSpec {
      pub(crate) keys: Vec<usize>,
      pub(crate) edges: Vec<(usize, usize)>,
    }

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
}
