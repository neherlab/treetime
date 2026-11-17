use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::gtr::partition_model::PartitionModel;
use crate::partition::marginal::mixture::site_log_lh::{INLINE_CATEGORIES, SparseMixture, SparseSiteLogLh};
use crate::partition::marginal::posterior_states::most_likely_state;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalPasses, MarginalUpdate};
use crate::partition::marginal::sparse::forward::compute_ml_subs_for_nodes;
use crate::partition::marginal::sparse::message::EPS;
use crate::partition::marginal::sparse::partition::{PartitionMarginalSparse, SparseMarginalEdges};
use crate::partition::storage::sparse::{
  FixedProfiles, INLINE_STATES, SiteProfile, SiteScales, SparseEdgeBackward, SparseEdgeForward, SparseNodeState,
  SparseSeqDistribution,
};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::{EitherOrBoth, Itertools, izip};
use ndarray::aview1;
use smallvec::{SmallVec, smallvec};
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, AsciiTable, LogLh};
use treetime_utils::{make_internal_error, make_internal_report};

const GROUP_SPREAD_LIMIT: f64 = EPS / 2.0;

pub(crate) fn sparse_mixture_update(
  partition: &PartitionMarginalSparse,
  model: &PartitionModel,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  input: &NodeMap<SparseNodeState>,
) -> Result<(SparseMixtureUpdate, SparseMixtureState), Report> {
  require_site_scales(partition)?;
  let root_key = graph.root_key()?;
  let mut categories = vec![];
  let mut site_log_lh = vec![];
  for gtr in model.category_gtrs()? {
    let MarginalUpdate {
      mut node_states,
      mut edges,
      ..
    } = partition.marginal_update(&gtr, graph, branch_lengths, input)?;
    let root = node_states
      .get_mut(root_key)
      .ok_or_else(|| make_internal_report!("Rate category has no state at the root {root_key}"))?;
    site_log_lh.push(SparseSiteLogLh::from_pass(root, &edges.backward, &partition.alphabet));
    root.profile.scales = SiteScales::default();
    for edge in edges.backward.values_mut() {
      edge.msg_to_parent.scales = SiteScales::default();
      edge.msg_from_child.scales = SiteScales::default();
    }
    categories.push(SparseCategory {
      gtr,
      node_states,
      edges,
    });
  }
  let mixture = SparseMixture::new(
    &site_log_lh,
    model.gamma().weights(),
    &partition.obs_nodes[&root_key].composition,
    &partition.alphabet,
  )?;
  let node_states = mix_sparse_node_states(partition, &categories, &mixture)?;
  let Some(first) = categories.first() else {
    return make_internal_error!("A rate mixture needs at least one category");
  };
  let estimates = mixed_estimates(partition, graph, &node_states, &first.edges.backward)?;
  let update = MarginalUpdate {
    node_states,
    edges: MarginalEdges {
      backward: EdgeMap::new(),
      forward: EdgeMap::new(),
      estimates,
    },
    log_lh: mixture.log_lh,
  };
  Ok((
    update,
    SparseMixtureState {
      categories,
      mixture,
      branch_lengths: branch_lengths.clone(),
      endpoints: graph
        .get_edges()
        .map(|edge| Ok((edge.key(), graph.edge_endpoints(edge.key())?)))
        .collect::<Result<_, Report>>()?,
    },
  ))
}

pub(crate) fn sparse_mixture_log_lh(
  partition: &PartitionMarginalSparse,
  model: &PartitionModel,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  input: &NodeMap<SparseNodeState>,
) -> Result<LogLh, Report> {
  require_site_scales(partition)?;
  let root_key = graph.root_key()?;
  let input = PartitionMarginalSparse::fresh_backward_input(input);
  let site_log_lh = model
    .category_gtrs()?
    .iter()
    .map(|gtr| {
      let MarginalBackward { node_states, backward } =
        partition.marginal_backward(gtr, graph, branch_lengths, &input)?;
      Ok(SparseSiteLogLh::from_pass(
        &node_states[&root_key],
        &backward,
        &partition.alphabet,
      ))
    })
    .try_collect::<_, Vec<_>, Report>()?;
  Ok(
    SparseMixture::new(
      &site_log_lh,
      model.gamma().weights(),
      &partition.obs_nodes[&root_key].composition,
      &partition.alphabet,
    )?
    .log_lh,
  )
}

pub(crate) type SparseMixtureUpdate = MarginalUpdate<SparseNodeState, SparseEdgeBackward, SparseEdgeForward, Vec<Sub>>;

#[derive(Clone, Debug)]
pub(crate) struct SparseMixtureState {
  pub(crate) categories: Vec<SparseCategory>,
  pub(crate) mixture: SparseMixture,
  pub(crate) branch_lengths: EdgeMap<f64>,
  pub(crate) endpoints: EdgeMap<(GraphNodeKey, GraphNodeKey)>,
}

#[derive(Clone, Debug)]
pub struct SparseCategory {
  pub(crate) gtr: GTR,
  pub(crate) node_states: NodeMap<SparseNodeState>,
  pub(crate) edges: SparseMarginalEdges,
}

fn require_site_scales(partition: &PartitionMarginalSparse) -> Result<(), Report> {
  if partition.track_site_scales {
    Ok(())
  } else {
    make_internal_error!("A sparse rate mixture needs per-site scales, but the partition does not track them")
  }
}

fn mix_sparse_node_states(
  partition: &PartitionMarginalSparse,
  categories: &[SparseCategory],
  mixture: &SparseMixture,
) -> Result<NodeMap<SparseNodeState>, Report> {
  let Some(first) = categories.first() else {
    return make_internal_error!("A rate mixture needs at least one category");
  };
  first
    .node_states
    .keys()
    .map(|key| {
      let states = categories
        .iter()
        .map(|category| {
          category
            .node_states
            .get(key)
            .ok_or_else(|| make_internal_report!("A rate category has no state at node {key}"))
        })
        .collect::<Result<SmallVec<[&SparseNodeState; 8]>, Report>>()?;
      let composition = &partition.obs_nodes[&key].composition;
      Ok((key, mix_sparse_node(&partition.alphabet, composition, &states, mixture)))
    })
    .collect()
}

fn mix_sparse_node(
  alphabet: &Alphabet,
  composition: &Composition,
  states: &[&SparseNodeState],
  mixture: &SparseMixture,
) -> SparseNodeState {
  let first = states[0];
  let n_states = alphabet.n_canonical();
  let groups: SmallVec<[MixedGroup; INLINE_STATES]> = first
    .profile
    .fixed
    .iter()
    .map(|entry| {
      let weights = mixture.group(entry.state).map(|group| group.posterior.as_slice());
      let profiles = states
        .iter()
        .filter_map(|state| state.profile.fixed.get(entry.state).map(|fixed| fixed.dis.as_slice()))
        .collect::<SmallVec<[&[f64]; INLINE_CATEGORIES]>>();
      let profiles_by_category = states
        .iter()
        .map(|state| state.profile.fixed.get(entry.state).map(|fixed| fixed.dis.as_slice()));
      MixedGroup {
        state: entry.state,
        dis: weighted_profile(n_states, weights, profiles_by_category).unwrap_or_else(|| entry.dis.clone()),
        spread: profile_spread(n_states, &profiles),
      }
    })
    .collect();
  let group = |state: AsciiChar| groups.iter().find(|group| group.state == state);

  let positions = states
    .iter()
    .map(|state| state.profile.variable.keys().copied())
    .kmerge()
    .dedup()
    .merge_join_by(mixture.listed(), |pos, (listed, ..)| pos.cmp(listed))
    .map(|position| match position {
      EitherOrBoth::Left(pos) => (pos, true, None),
      EitherOrBoth::Right((pos, _, site)) => (pos, false, Some(site)),
      EitherOrBoth::Both(pos, (_, _, site)) => (pos, true, Some(site)),
    });
  let mut variable = VarPosMap::default();
  let mut labeled = AsciiTable::<usize>::default();
  for (pos, in_variable, listed) in positions {
    let state = first.sequence[pos];
    if !in_variable && group(state).is_some_and(|group| group.spread <= GROUP_SPREAD_LIMIT) {
      continue;
    }
    let weights = listed.map_or_else(|| mixture.unlisted_posterior(state), |site| &site.posterior);
    let profiles = states.iter().map(|node| {
      node.profile.variable.get(pos).map(|var| var.dis).or_else(|| {
        alphabet
          .is_canonical(state)
          .then(|| node.profile.fixed.get(state).map(|fixed| fixed.dis.as_slice()))
          .flatten()
      })
    });
    let Some(mixed) = weighted_profile(n_states, Some(weights), profiles) else {
      continue;
    };
    let differs = group(state).is_none_or(|group| izip!(&mixed, &group.dis).any(|(a, b)| (a - b).abs() > EPS));
    if !in_variable && !differs {
      continue;
    }
    let label = if alphabet.is_canonical(state) {
      state
    } else {
      states
        .iter()
        .find_map(|node| node.profile.variable.get(pos).map(|var| var.state))
        .unwrap_or_else(|| alphabet.char(most_likely_state(&aview1(&mixed)).unwrap_or(0)))
    };
    labeled[label] += 1;
    variable.insert(pos, label, &mixed);
  }

  let mut fixed = FixedProfiles::default();
  for MixedGroup { state, dis, .. } in groups {
    let count = composition.get(state).unwrap_or(0).saturating_sub(labeled[state]);
    fixed.insert(state, count, dis);
  }
  SparseNodeState {
    sequence: Arc::clone(&first.sequence),
    profile: SparseSeqDistribution {
      variable,
      fixed,
      log_lh: mixture.log_lh,
      scales: SiteScales::default(),
    },
    globally_missing: Arc::clone(&first.globally_missing),
  }
}

struct MixedGroup {
  state: AsciiChar,
  dis: SiteProfile,
  spread: f64,
}

fn profile_spread(n_states: usize, profiles: &[&[f64]]) -> f64 {
  (0..n_states)
    .map(|i| {
      let (min, max) = profiles
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), profile| {
          (min.min(profile[i]), max.max(profile[i]))
        });
      max - min
    })
    .fold(0.0, f64::max)
}

fn weighted_profile<'a>(
  n_states: usize,
  weights: Option<&[f64]>,
  profiles: impl Iterator<Item = Option<&'a [f64]>>,
) -> Option<SiteProfile> {
  let mut mixed: SiteProfile = smallvec![0.0; n_states];
  let mut total = 0.0;
  for (k, profile) in profiles.enumerate() {
    let Some(profile) = profile else {
      continue;
    };
    let weight = weights.map_or(1.0, |weights| weights[k]);
    for (value, p) in izip!(mixed.iter_mut(), profile) {
      *value += weight * p;
    }
    total += weight;
  }
  if total <= 0.0 {
    return None;
  }
  for value in &mut mixed {
    *value /= total;
  }
  Some(mixed)
}

fn mixed_estimates(
  partition: &PartitionMarginalSparse,
  graph: &Graph,
  node_states: &NodeMap<SparseNodeState>,
  backward: &EdgeMap<SparseEdgeBackward>,
) -> Result<EdgeMap<Vec<Sub>>, Report> {
  graph
    .get_edges()
    .map(|edge| {
      let edge_key = edge.key();
      let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
      let subs = compute_ml_subs_for_nodes(
        &partition.alphabet,
        &node_states[&parent_key],
        &partition.obs_nodes[&parent_key],
        &node_states[&child_key],
        &partition.obs_nodes[&child_key],
        &partition.obs_edges[&edge_key],
        &backward[&edge_key],
      )?;
      Ok((edge_key, subs))
    })
    .collect()
}
