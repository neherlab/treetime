use crate::make_error;
use crate::representation::partition::discrete::PartitionDiscrete;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::payload::discrete::{DiscreteEdgeData, DiscreteNodeData};
use eyre::Report;
use itertools::Itertools;
use parking_lot::Mutex;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::sync::mutex::extract_parallel_error;

pub fn run_discrete_marginal<N, E>(graph: &Graph<N, E, ()>, partition: &mut PartitionDiscrete) -> Result<f64, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let n_states = partition.n_states();

  // Initialize leaf nodes from traits
  for leaf in graph.get_leaves() {
    let leaf_key = leaf.read_arc().key();
    partition
      .nodes
      .entry(leaf_key)
      .or_insert_with(|| DiscreteNodeData::missing(n_states));
  }

  // Backward pass: postorder traversal (leaves to root)
  discrete_marginal_backward(graph, partition)?;

  // Forward pass: preorder traversal (root to leaves)
  discrete_marginal_forward(graph, partition)?;

  // Return total log likelihood from root
  let root = graph.get_exactly_one_root()?;
  let root_key = root.read_arc().key();
  let log_lh = partition.get_log_lh(root_key);

  Ok(log_lh)
}

fn discrete_marginal_backward<N, E>(graph: &Graph<N, E, ()>, partition: &mut PartitionDiscrete) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));
  let partition = Arc::new(Mutex::new(partition));

  graph.par_iter_breadth_first_backward(|node| {
    let mut partition = partition.lock();
    if let Err(e) = partition.process_node_backward(&node) {
      let mut guard = error.lock();
      if guard.is_none() {
        *guard = Some(e);
      }
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });

  extract_parallel_error(error)
}

fn discrete_marginal_forward<N, E>(graph: &Graph<N, E, ()>, partition: &mut PartitionDiscrete) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));
  let partition = Arc::new(Mutex::new(partition));

  graph.par_iter_breadth_first_forward(|node| {
    let mut partition = partition.lock();
    if let Err(e) = partition.process_node_forward(graph, &node) {
      let mut guard = error.lock();
      if guard.is_none() {
        *guard = Some(e);
      }
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });

  extract_parallel_error(error)
}

pub fn attach_traits<N, E>(
  partition: &mut PartitionDiscrete,
  graph: &Graph<N, E, ()>,
  traits: &BTreeMap<String, String>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let n_states = partition.n_states();
  validate_trait_names(graph, traits)?;

  for leaf in graph.get_leaves() {
    let leaf_guard = leaf.read_arc();
    let leaf_key = leaf_guard.key();
    let leaf_payload = leaf_guard.payload().read_arc();

    let leaf_name = leaf_payload.name().map(|n| n.as_ref().to_owned()).unwrap_or_default();

    let node_data = if let Some(trait_value) = traits.get(&leaf_name) {
      if let Some(index) = partition.states.get_index(trait_value) {
        DiscreteNodeData::from_observed(index, n_states)
      } else {
        DiscreteNodeData::missing(n_states)
      }
    } else {
      DiscreteNodeData::missing(n_states)
    };

    partition.nodes.insert(leaf_key, node_data);
  }

  // Initialize edges
  for edge in graph.get_edges() {
    let edge_key = edge.read_arc().key();
    partition.edges.insert(edge_key, DiscreteEdgeData::default());
  }

  Ok(())
}

fn validate_trait_names<N, E>(graph: &Graph<N, E, ()>, traits: &BTreeMap<String, String>) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let leaf_names = graph
    .get_leaves()
    .iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      let payload = leaf.payload().read_arc();
      payload.name().map(|name| name.as_ref().to_owned()).unwrap_or_default()
    })
    .collect_vec();

  let missing_in_metadata = leaf_names
    .iter()
    .filter(|leaf_name| !traits.contains_key(*leaf_name))
    .cloned()
    .collect_vec();
  if !missing_in_metadata.is_empty() {
    return make_error!(
      "Mugration: tree leaves missing from metadata: {}",
      missing_in_metadata.iter().join(", ")
    );
  }

  let missing_in_tree = traits
    .keys()
    .filter(|trait_name| !leaf_names.iter().any(|leaf_name| leaf_name == *trait_name))
    .cloned()
    .collect_vec();
  if !missing_in_tree.is_empty() {
    return make_error!(
      "Mugration: metadata names missing from tree leaves: {}",
      missing_in_tree.iter().join(", ")
    );
  }

  Ok(())
}
