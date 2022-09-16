use crate::graph::edge::GraphEdge;
use crate::graph::graph::{Graph, GraphNodeForward};
use crate::graph::node::GraphNode;
use itertools::Itertools;
use std::collections::HashSet;

pub fn assign_node_names<N: GraphNode, E: GraphEdge>(graph: &mut Graph<N, E>) {
  let mut names = graph
    .get_node_payloads()
    .map(|node| node.read().name().to_owned())
    .collect::<HashSet<String>>();

  let mut node_counter = 0_i64;
  graph.iter_depth_first_preorder_forward(
    |GraphNodeForward {
       key,
       payload,
       parents,
       is_leaf,
       ..
     }| {
      if payload.name().is_empty() {
        let mut name = format!("NODE_{node_counter:07}");
        while names.contains(&name) {
          node_counter += 1;
          name = format!("NODE_{node_counter:07}");
        }
        names.insert(name.clone());
        payload.set_name(&name);
      }
    },
  );
}
