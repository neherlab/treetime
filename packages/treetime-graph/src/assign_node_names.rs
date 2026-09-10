use crate::edge::GraphEdge;
use crate::graph::Graph;
use crate::graph_traverse::GraphNodeForward;
use crate::node::{GraphNode, GraphNodeKey, Named};
use crate::value_maps::node_names;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};

/// Assign synthetic `NODE_{counter:07}` names to unnamed internal nodes in DFS-preorder, and return
/// the resulting node-keyed name map.
///
/// The returned map is `node_names(graph)` taken after the write, so it holds exactly the
/// `Option<String>` each consumer would read off the payload at this program point: the parsed name
/// for named nodes and the freshly assigned synthetic name for internals. Threading this map lets
/// name consumers read the value instead of the payload; the payload write stays transitional so
/// not-yet-migrated readers and the parse path keep the same names.
pub fn assign_node_names<N: GraphNode + Named, E: GraphEdge, D: Sync + Send>(
  graph: &Graph<N, E, D>,
) -> Result<BTreeMap<GraphNodeKey, Option<String>>, Report> {
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
      Ok(())
    },
  )?;

  Ok(node_names(graph))
}
