use crate::clock::clock_model::ClockModel;
use crate::clock::clock_state::ClockState;
use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::coalescent::CoalescentModel;
use crate::optimize::gather::{
  gather_timetree_edge_contributions, gather_timetree_edge_indel_counts, timetree_total_sequence_length,
};
use crate::optimize::indel::estimate_indel_rate;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::timetree::timetree_state::TimetreeState;
use crate::timetree::utils::initialize_node_divergences;
use eyre::Report;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

/// Target resolution of every *stored* timetree time-distribution grid (design D3, proposal Part D).
///
/// This is a minimum point count, not a fixed one: the mass re-window
/// ([`rewindow_to_mass`](treetime_distribution::rewindow_to_mass)) resamples the mass domain at
/// spacing `min(mass_width / (GRID_POINTS - 1), source_dx)`, so a stored grid holds at least this many
/// points and is never coarser than the distribution it was re-windowed from. Changing this one
/// constant changes the baseline stored-grid resolution everywhere (construction, backward fold,
/// forward refinement).
pub const GRID_POINTS: usize = 300;

/// Per-side tail mass fraction trimmed when sizing a soft grid edge by probability mass (design D4,
/// proposal Axis 2).
///
/// A soft side's grid edge is placed so that `EPS` of the total mass lies beyond it; the already-fitted
/// tail law carries that mass rather than discarding it, so accuracy is bounded by tail-law fidelity,
/// not by `EPS`. A grid therefore holds at least `1 - 2*EPS` of the mass across its two soft sides.
pub const EPS: f64 = 5e-4;

/// Infer node times.
///
/// `coalescent` is the prior imposed on those times, or `None` for a run that carries no
/// coalescent prior. It is supplied rather than derived here because its lineage counts must be
/// held fixed across passes; see [`CoalescentModel`].
///
/// `branch_lengths` and `names` are the per-edge branch-length and per-node name maps the caller
/// snapshots from the current tree just before this pass. The branch-distribution construction reads
/// each edge's length from `branch_lengths`, and the forward pass reads each node's label from
/// `names`. The caller re-snapshots them after any length or
/// topology change so each pass sees the current tree.
#[allow(clippy::too_many_arguments)]
pub fn run_timetree(
  graph: &mut Graph,
  constraints: &DateConstraints,
  partitions: &[PartitionTimetree],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  clock_model: &ClockModel,
  coalescent: Option<&CoalescentModel>,
  no_indels: bool,
  mut state: TimetreeState,
  clock_state: &mut ClockState,
) -> Result<TimetreeState, Report> {
  info!("# Running timetree inference");

  info!("## Calculating divergence distances");
  initialize_node_divergences(graph, clock_state, branch_lengths, names)?;

  // Rebuild the state's maps for the current topology, carrying every value-resident date field
  // forward. Times, distributions, bad-branch flags, and time lengths all live on the state
  // (their producers write them straight into it); the fixed date constraints stay in the durable
  // `constraints` the passes borrow.
  state.reseed_from_values(graph);

  info!("## Using clock model");
  let clock_rate = clock_model.clock_rate();
  info!("**Clock rate:** {clock_rate:.6e}");

  if !partitions.is_empty() {
    info!("## Computing branch distributions from partitions");
    compute_branch_distributions_marginal_mode(graph, partitions, branch_lengths, clock_rate, no_indels, &mut state)?;
  } else {
    info!("## Creating branch distributions from input lengths");
    create_branch_distributions_input_mode(graph, branch_lengths, clock_rate, &mut state)?;
  }

  info!("## Propagating distributions backward");
  propagate_distributions_backward(graph, constraints, coalescent, &mut state)?;

  info!("## Propagating distributions forward");
  propagate_distributions_forward(graph, constraints, names, &mut state)?;

  info!("# Timetree inference completed");
  Ok(state)
}

/// Weight given to the freshly inferred value when a commit damps against the previous one.
///
/// Each pass re-infers every node time at once, so an undamped write-back oscillates: measured
/// on `data/ebola/20` with the coalescent lineage counts held fixed, the per-round maximum time
/// change repeated bit-identically from round 2 onward. Blending halves is the same remedy
/// [`apply_damping`](crate::optimize::iteration::apply_damping) applies in the branch-length
/// optimization loop. It changes the path, not the answer: `b = (1 - f) b + f b` for any `f`.
pub const CLOCK_BRANCH_LENGTH_DAMPING: f64 = 0.5;

/// Commit each edge's clock-constrained branch length, `clock_rate * gamma * (t_child - t_parent)`.
///
/// This is the clock-constrained length the marginal reconstruction propagates sequence profiles
/// along afterwards (see [`timetree_branch_lengths`]). It is the constrained M-step of
/// the refinement loop: `branch_length` stays the free ML or input estimate, while profiles move
/// along lengths the inferred times imply. v0 gets this for free by running the whole timetree in
/// divergence units, where `branch_length = clock_length` is dimensionally a no-op; v1 works in
/// calendar time, so the conversion is explicit.
///
/// `damping` is the weight given to the newly computed value, the rest carried over from the
/// previous commit. Pass 1.0 wherever the previous value describes a different tree -- the first
/// commit of a run, and the one after polytomy resolution has re-parented edges -- and
/// [`CLOCK_BRANCH_LENGTH_DAMPING`] inside the refinement loop.
///
/// An edge whose duration is negative gets zero. The forward pass clamps internal children to
/// their parent but leaves observed leaf dates alone, so a leaf dated before its parent reaches
/// here; that is a real inconsistency in the input or the fit, and it is reported rather than
/// silently floored. Edges whose endpoints are not both dated are left untouched.
pub fn commit_clock_branch_lengths(
  graph: &Graph,
  clock_rate: f64,
  damping: f64,
  clock_branch_lengths: &mut BTreeMap<GraphEdgeKey, f64>,
  state: &TimetreeState,
) {
  let node_time = |key| state.nodes.get(&key).and_then(|node| node.time);

  // The committed value blends against the previous one held in the routed map, so the fold builds
  // the new values in parallel first (each reading its own previous value from the shared map) and
  // inserts them serially, without a lock.
  let previous_lengths: &BTreeMap<GraphEdgeKey, f64> = clock_branch_lengths;
  let committed: Vec<(GraphEdgeKey, f64, bool)> = graph
    .get_edges()
    .par_iter()
    .filter_map(|edge_ref| {
      let edge_ref = edge_ref.read_arc();
      let key = edge_ref.key();
      let (Some(parent_time), Some(child_time)) = (node_time(edge_ref.source()), node_time(edge_ref.target())) else {
        return None;
      };

      let duration = child_time - parent_time;
      let fresh = clock_rate * state.edge(key).gamma * duration.max(0.0);
      let value = match previous_lengths.get(&key) {
        Some(previous) => (1.0 - damping) * previous + damping * fresh,
        None => fresh,
      };

      Some((key, value, duration < 0.0))
    })
    .collect();

  let mut inverted = 0_usize;
  for (key, value, is_inverted) in committed {
    clock_branch_lengths.insert(key, value);
    inverted += usize::from(is_inverted);
  }

  if inverted > 0 {
    warn!(
      "Timetree: {inverted} branch(es) run backwards in time, i.e. the child is dated before its \
       parent. Their clock branch lengths were committed as zero. This is expected only where an \
       observed leaf date conflicts with the fitted clock, since the forward pass clamps internal \
       nodes to their parent but leaves leaf dates as given."
    );
  }
}

fn compute_branch_distributions_marginal_mode(
  graph: &Graph,
  partitions: &[PartitionTimetree],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_rate: f64,
  no_indels: bool,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  let total_sites = timetree_total_sequence_length(partitions);
  let one_mutation = 1.0 / total_sites as f64;

  let indel_counts = gather_timetree_edge_indel_counts(graph, partitions);
  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, &indel_counts, branch_lengths)
  };
  let contributions = gather_timetree_edge_contributions(graph, partitions)?;

  info!(
    "Computing branch distributions from {} partition(s) with {} total sites",
    partitions.len(),
    total_sites
  );
  debug!("One mutation = {one_mutation:.6e} substitutions/site");
  debug!("Indel rate = {indel_rate:.6e} indels/(site*time)");

  // Compute each edge's branch-length distribution in parallel, reading its relaxed-clock rate from
  // the value state, and carry the distribution out to insert into the value serially: the value's
  // per-edge map cannot be written from parallel workers without a lock. The immutable reborrow ends
  // at the `collect`, so the serial inserts can take a mutable borrow.
  let edge_states: &TimetreeState = state;
  let distributions: Vec<(GraphEdgeKey, Option<f64>, Arc<Distribution<NegLog>>)> = graph
    .get_edges()
    .par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Option<f64>, Arc<Distribution<NegLog>>), Report> {
        let edge_key = edge_ref.read_arc().key();
        let branch_length = branch_lengths[&edge_key].unwrap_or(one_mutation);
        let gamma = edge_states.edge(edge_key).gamma;

        debug!("Edge {edge_key:?}: input branch_length = {branch_length:.6e}, gamma = {gamma:.4}");

        let contributions = &contributions[&edge_key];
        let indel_count: usize = if no_indels { 0 } else { indel_counts[&edge_key] };
        let distribution = compute_branch_length_distribution(
          contributions,
          indel_count,
          indel_rate,
          branch_length,
          one_mutation,
          GRID_POINTS,
          clock_rate,
          gamma,
        )?;

        let time_length = distribution.likely_time();
        if let Some(likely_time) = time_length {
          debug!("Edge {edge_key:?}: distribution peak at time = {likely_time:.6e}");
        }

        Ok((edge_key, time_length, distribution))
      },
    )
    .collect::<Result<Vec<_>, Report>>()?;

  for (key, time_length, distribution) in distributions {
    let entry = state.edge_mut(key);
    entry.time_length = time_length;
    entry.branch_length_distribution = Some(distribution);
  }
  Ok(())
}

pub(super) fn create_branch_distributions_input_mode(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_rate: f64,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  // Build each edge's point branch-length distribution in parallel, reading its relaxed-clock rate
  // from the value state, and carry the distribution out to insert into the value serially. An edge
  // with neither a branch length nor a time length is left untouched: its previous value-side
  // distribution and time length carry over unchanged. The immutable reborrow ends at the `collect`,
  // so the serial inserts can take a mutable borrow.
  let edge_states: &TimetreeState = state;
  let distributions: Vec<(GraphEdgeKey, f64, Arc<Distribution<NegLog>>)> = graph
    .get_edges()
    .par_iter()
    .filter_map(|edge_ref| {
      let key = edge_ref.read_arc().key();
      // TODO: this is wrong. The branch length distribution should be a gamma distribution with branch_length/one_mutation
      // as the shape parameter. n_mut = branch_length/one_mutation --> P(dt) = (mu*dt)^n_mut * exp(-mu*dt) / n_mut!
      let time_duration = if let Some(branch_length) = branch_lengths[&key] {
        // Convert branch length (substitutions/site) to time duration (years)
        // gamma > 1 means faster evolution, so same substitutions correspond to shorter time
        let effective_clock_rate = clock_rate * edge_states.edge(key).gamma;
        Some(branch_length / effective_clock_rate)
      } else {
        edge_states.edge(key).time_length
      };

      time_duration.map(|time_duration| {
        // Negative-log ordinate `0` is the `NegLog` multiplicative identity (probability 1).
        let distribution = Arc::new(Distribution::point(time_duration, 0.0));
        (key, time_duration, distribution)
      })
    })
    .collect();

  for (key, time_duration, distribution) in distributions {
    let entry = state.edge_mut(key);
    entry.time_length = Some(time_duration);
    entry.branch_length_distribution = Some(distribution);
  }

  Ok(())
}

/// The per-edge branch length each post-commit marginal reconstruction propagates sequence profiles
/// along, keyed by edge, sourcing the clock-constrained length from the routed commit map.
///
/// The timetree counterpart of [`profile_branch_lengths`](crate::ancestral::marginal::profile_branch_lengths):
/// the value is the committed clock length held in `clock_branch_lengths` when the edge has one, and
/// the edge's own branch length otherwise. Used only after a commit, where the clock length has been
/// established; before the first commit the two collectors agree, because no clock length exists yet.
pub fn timetree_branch_lengths(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> BTreeMap<GraphEdgeKey, f64> {
  graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      let branch_length = clock_branch_lengths
        .get(&key)
        .copied()
        .or_else(|| branch_lengths[&key])
        .unwrap_or(0.0);
      (key, branch_length)
    })
    .collect()
}
