use crate::substitution_model::{NUCLEOTIDES, SubstitutionModel};
use eyre::Report;
use itertools::Itertools;
use rand::Rng;
use rand_distr::{Distribution, Exp1};
use std::iter::repeat_with;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_utils::{make_error, make_internal_error};

pub fn evolve_site(
  state: usize,
  duration: f64,
  model: &SubstitutionModel,
  rng: &mut impl Rng,
) -> Result<SiteEvolution, Report> {
  if state >= NUCLEOTIDES.len() {
    return make_error!(
      "State {state} is not one of the {} nucleotide states",
      NUCLEOTIDES.len()
    );
  }
  validate_duration(duration, "Duration")?;
  Ok(jump_process(state, duration, model, rng))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SiteEvolution {
  pub state: usize,
  pub jumps: usize,
}

pub fn evolve_sequences(
  graph: &Graph,
  durations: &EdgeMap<f64>,
  model: &SubstitutionModel,
  length: usize,
  rng: &mut impl Rng,
) -> Result<SimulatedSequences, Report> {
  if length == 0 {
    return make_error!("Sequence length must be positive, got 0");
  }
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  for edge in graph.get_edges() {
    let key = edge.key();
    let Some(&duration) = durations.get(key) else {
      return make_error!("Edge {key} has no duration");
    };
    validate_edge_duration(key, duration)?;
  }

  let mut states: NodeMap<Vec<usize>> = NodeMap::new();
  let mut substitutions = EdgeMap::new();
  for &node in schedule.preorder() {
    let node_states = match schedule.parents(node) {
      [] => repeat_with(|| model.sample_state(rng)).take(length).collect_vec(),
      [(parent, edge)] => {
        let duration = durations[edge];
        let evolutions = states[parent]
          .iter()
          .map(|&state| jump_process(state, duration, model, rng))
          .collect_vec();
        substitutions.insert(*edge, evolutions.iter().map(|evolution| evolution.jumps).sum());
        evolutions.iter().map(|evolution| evolution.state).collect_vec()
      },
      parents => {
        return make_internal_error!("Node {node} has {} parents in a tree", parents.len());
      },
    };
    states.insert(node, node_states);
  }

  let sequences = states
    .into_iter()
    .map(|(node, node_states)| (node, node_states.into_iter().map(|state| NUCLEOTIDES[state]).collect()))
    .collect();
  Ok(SimulatedSequences {
    sequences,
    substitutions,
  })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SimulatedSequences {
  pub sequences: NodeMap<Seq>,
  pub substitutions: EdgeMap<usize>,
}

impl SimulatedSequences {
  pub fn leaf_records(&self, graph: &Graph, names: &NodeNames) -> Vec<AlignmentRecord> {
    leaf_records(graph, names, &self.sequences)
  }
}

pub fn leaf_records(graph: &Graph, names: &NodeNames, sequences: &NodeMap<Seq>) -> Vec<AlignmentRecord> {
  graph
    .get_leaves()
    .map(|leaf| AlignmentRecord {
      name: names[leaf.key()].to_owned(),
      seq: sequences[&leaf.key()].clone(),
    })
    .collect()
}

pub(crate) fn validate_edge_duration(edge: GraphEdgeKey, duration: f64) -> Result<(), Report> {
  validate_duration(duration, &format!("Duration of edge {edge}"))
}

fn validate_duration(duration: f64, what: &str) -> Result<(), Report> {
  if !(duration.is_finite() && duration >= 0.0) {
    return make_error!("{what} must be finite and non-negative, got {duration}");
  }
  Ok(())
}

fn jump_process(state: usize, duration: f64, model: &SubstitutionModel, rng: &mut impl Rng) -> SiteEvolution {
  let mut state = state;
  let mut jumps = 0;
  let mut elapsed = 0.0;
  while let Some(jump) = model.jump(state) {
    let exp1: f64 = Exp1.sample(rng);
    elapsed += exp1 / jump.exit_rate;
    if elapsed > duration {
      break;
    }
    state = jump.targets.sample(rng);
    jumps += 1;
  }
  SiteEvolution { state, jumps }
}
