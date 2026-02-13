use crate::commands::timetree::coalescent::contributions::compute_node_contributions;
use crate::commands::timetree::coalescent::events::collect_tree_events;
use crate::commands::timetree::coalescent::integration::compute_integral_merger_rate;
use crate::commands::timetree::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::commands::timetree::timetree_traits::TimetreeNode;
use crate::distribution::distribution::{Distribution, DistributionNegLog};
use eyre::Report;
use indexmap::IndexMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};

/// Computes Kingman coalescent prior contributions for all nodes in the phylogenetic tree.
///
/// Returns distributions that encode coalescent likelihood contributions for each node,
/// to be multiplied with node time distributions during backward pass optimization.
///
/// # Meaning
///
/// The Kingman coalescent model provides a probabilistic framework for assessing whether
/// the timing and structure of a phylogenetic tree are consistent with a given population
/// history. Not all tree topologies and divergence times are equally likely - trees with
/// many lineages coalescing simultaneously are less probable than gradual coalescence,
/// especially in large populations.
///
/// This function computes how likely each node's divergence time is under the coalescent
/// model, given the population size history Tc(t) and the number of concurrent lineages k(t).
/// These likelihood contributions act as priors that guide the time tree inference toward
/// more biologically plausible configurations by penalizing unlikely coalescence patterns.
///
/// # Notation and Terms
///
/// - `t` - time (negative values for past, zero at present)
/// - `k(t)` - number of concurrent lineages at time t
/// - `Tc(t)` - coalescence time scale (effective population size) at time t, controls merger rate
/// - `κ(t)` - branch merger rate: rate for one branch to merge with any other = (k(t)-1)/(2*Tc(t))
/// - `λ(t)` - total merger rate: rate for any merger to occur = k(t)*(k(t)-1)/(2*Tc(t))
/// - `I(t)` - cumulative merger rate: I(t) = ∫₀ᵗ κ(t') dt'
///
/// # Kingman Coalescent Probability Density
///
/// The coalescent contributions differ between leaf and internal nodes:
///
/// - Leaf nodes: exp(I(t))
///   Represents survival probability that the lineage existed without coalescing
///   from present to sampling time. Important for:
///   - Uncertain/range dates: helps infer actual sampling time
///   - Precise dates: encodes survival cost that propagates to ancestors
///
/// - Internal nodes with m children (m≥2): λ(t)^(m-1) · exp(-I(t))
///   - λ(t)^(m-1): probability density of m-way merger at time t
///   - exp(-I(t)): probability of no merger before time t
///
/// Note: The sign convention follows Python v0 which stores -I(t) in neg-log space.
///
/// # Returns
///
/// Map from node keys to distributions representing coalescent prior contributions (NegLog space).
/// Each distribution should be multiplied with the node's time distribution.
pub fn compute_coalescent_contributions<N, E, D>(
  graph: &Graph<N, E, D>,
  tc: &Distribution,
) -> Result<IndexMap<GraphNodeKey, Arc<DistributionNegLog>>, Report>
where
  N: GraphNode + TimetreeNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let (present_time, events_calendar) = collect_tree_events(graph)?;
  let events_tbp: Vec<_> = events_calendar
    .iter()
    .map(|(t, delta)| (present_time - *t, *delta))
    .collect();

  let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
  let integral_merger_rate = compute_integral_merger_rate(tc, &lineage_counts)?;
  compute_node_contributions(graph, &integral_merger_rate, tc, &lineage_counts, present_time)
}
