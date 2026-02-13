use crate::edge::GraphEdge;
use crate::graph::Graph;
use crate::graph_traverse::GraphNodeForward;
use crate::node::{GraphNode, Named};
use std::collections::BTreeSet;

pub fn assign_node_names<N: GraphNode + Named, E: GraphEdge, D: Sync + Send>(graph: &Graph<N, E, D>) {
  let mut names = graph
    .get_node_payloads()
    .map(|node| {
      node
        .read_arc()
        .name()
        .map_or_else(|| "Unknown".to_owned(), |name| name.as_ref().to_owned())
    })
    .collect::<BTreeSet<String>>();

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
      if payload.name().is_none_or(|name| name.as_ref().is_empty()) {
        let mut name = format!("NODE_{internal_node_counter:07}");
        while names.contains(&name) {
          internal_node_counter += 1;
          name = format!("NODE_{internal_node_counter:07}");
        }
        payload.set_name(Some(&name));
        names.insert(name);
        internal_node_counter += 1;
      }
    },
  );
}
