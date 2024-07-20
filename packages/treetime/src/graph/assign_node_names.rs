use crate::graph::edge::GraphEdge;
use crate::graph::graph::{Graph, GraphNodeForward};
use crate::graph::node::{GraphNode, Named};
use std::collections::HashSet;

pub fn assign_node_names<N: GraphNode + Named, E: GraphEdge, D: Sync + Send>(graph: &Graph<N, E, D>) {
  let mut names = graph
    .get_node_payloads()
    .map(|node| node.read().name().to_owned())
    .collect::<HashSet<String>>();

  let mut internal_node_counter = 0;

  graph.iter_depth_first_preorder_forward(
    |GraphNodeForward {
       key,
       mut payload,
       parents,
       is_leaf,
       is_root,
       ..
     }| {
      if payload.name().is_empty() {
        let mut name = format!("NODE_{internal_node_counter:07}");
        while names.contains(&name) {
          internal_node_counter += 1;
          name = format!("NODE_{internal_node_counter:07}");
        }
        names.insert(name.clone());
        payload.set_name(&name);
        internal_node_counter += 1;
      }
    },
  );
}
