use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::mixture::dense::DenseMixtureState;
use crate::partition::marginal::mixture::site_log_lh::{CategoryValues, MixedSite};
use crate::partition::marginal::mixture::sparse::SparseMixtureState;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::optimize;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::optimize::sparse::GroupedSiteContributions;
use crate::partition::storage::sparse::{SparseNodeState, SparseSeqDistribution};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, Array2, ArrayView1, Axis, aview1, concatenate};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::key_map::NodeMap;
use treetime_primitives::AsciiChar;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::interval::range::range_contains;
use treetime_utils::make_internal_report;

pub(crate) fn dense_mixture_contribution(
  partition: &PartitionMarginalDense,
  state: &DenseMixtureState,
  weights: &Array1<f64>,
  edge_key: GraphEdgeKey,
) -> Result<OptimizationContribution, Report> {
  let branch_length = state.branch_lengths[&edge_key].max(partition.inputs.min_branch_length);
  let blocks = state
    .categories
    .iter()
    .map(|category| {
      let backward = &category.backward[&edge_key];
      let forward = &category.forward[&edge_key];
      backward.msg_to_parent.dis.dot(category.gtr.v()) * forward.msg_to_child.dis.dot(&category.gtr.v_inv().t())
    })
    .collect_vec();
  let rate_eigvals = state
    .categories
    .iter()
    .map(|category| category.gtr.rate_eigvals())
    .collect_vec();
  let site_log_lh = &state.mixture.site_log_lh;
  let n_sites = site_log_lh.nrows();
  let mut scales = Array2::<f64>::zeros((n_sites, blocks.len()));
  for (k, (block, eigvals, &weight)) in izip!(&blocks, &rate_eigvals, weights).enumerate() {
    let decay = eigvals.mapv(|l| (l * branch_length).exp());
    for (site, row) in block.outer_iter().enumerate() {
      scales[[site, k]] = scale_exponent(weight, site_log_lh[[site, k]], row, &decay);
    }
  }
  for mut row in scales.outer_iter_mut() {
    let max = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    row.mapv_inplace(|exponent| if max.is_finite() { (exponent - max).exp() } else { 0.0 });
  }
  let scaled = izip!(&blocks, scales.columns())
    .map(|(block, scale)| block * &scale.insert_axis(Axis(1)))
    .collect_vec();
  let coefficients = concatenate(Axis(1), &scaled.iter().map(|block| block.view()).collect_vec())?;
  let rate_eigvals = concatenate(Axis(0), &rate_eigvals.iter().map(|values| values.view()).collect_vec())?;
  Ok(OptimizationContribution::Dense(
    optimize::dense::PartitionContribution::mixture(coefficients, rate_eigvals),
  ))
}

pub(crate) fn sparse_mixture_contribution(
  partition: &PartitionMarginalSparse,
  node_states: &NodeMap<SparseNodeState>,
  state: &SparseMixtureState,
  weights: &Array1<f64>,
  edge_key: GraphEdgeKey,
) -> Result<OptimizationContribution, Report> {
  let (parent_key, child_key) = state.endpoints[&edge_key];
  let parent_sequence = &node_states[&parent_key].sequence;
  let child_sequence = &node_states[&child_key].sequence;
  let alphabet = &partition.alphabet;
  let branch_length = fix_branch_length(partition.total_length, state.branch_lengths[&edge_key]);
  let categories = state
    .categories
    .iter()
    .map(|category| EdgeCategory {
      gtr: &category.gtr,
      to_child: &category.edges.forward[&edge_key].msg_to_child,
      to_parent: &category.edges.backward[&edge_key].msg_to_parent,
      decay: category.gtr.rate_eigvals().mapv(|l| (l * branch_length).exp()),
    })
    .collect_vec();
  let Some(first) = state.categories.first() else {
    return Err(make_internal_report!("A rate mixture needs at least one category"));
  };
  let child_no_data = &first.edges.backward[&edge_key].child_no_data;

  let mut fitch_states = BTreeMap::new();
  for sub in partition.obs_edges[&edge_key].fitch_subs() {
    fitch_states.entry(sub.pos()).or_insert_with(|| (sub.reff(), sub.qry()));
  }
  let variable = categories
    .iter()
    .flat_map(|category| {
      category
        .to_child
        .variable
        .keys()
        .chain(category.to_parent.variable.keys())
        .copied()
    })
    .chain(partition.obs_edges[&edge_key].fitch_subs().iter().map(Sub::pos))
    .collect::<BTreeSet<_>>();
  let listed = state
    .mixture
    .listed()
    .map(|(pos, _)| pos)
    .filter(|&pos| {
      !variable.contains(&pos)
        && !range_contains(child_no_data, pos)
        && alphabet.is_canonical(parent_sequence[pos])
        && alphabet.is_canonical(child_sequence[pos])
    })
    .collect_vec();

  let mut sites = GroupedSiteContributions::with_capacity(variable.len() + listed.len() + alphabet.n_canonical());
  let mut taken_from_groups = BTreeMap::<AsciiChar, usize>::new();
  for &pos in variable.iter().chain(&listed) {
    let (parent_state, child_state) = fitch_states.get(&pos).copied().unwrap_or_else(|| {
      (
        variable_label(&categories, pos, true).unwrap_or(parent_sequence[pos]),
        variable_label(&categories, pos, false).unwrap_or(child_sequence[pos]),
      )
    });
    if !categories[0].to_child.variable.contains_key(pos)
      && !range_contains(child_no_data, pos)
      && alphabet.is_canonical(parent_state)
    {
      *taken_from_groups.entry(parent_state).or_default() += 1;
    }
    let site = state
      .mixture
      .site(pos, parent_state)
      .ok_or_else(|| make_internal_report!("Position {pos} has no mixture weights"))?;
    let rows = categories
      .iter()
      .map(|category| {
        let parent = category.to_child.variable.get(pos).map(|var| var.dis).or_else(|| {
          category
            .to_child
            .fixed
            .get(parent_state)
            .map(|fixed| fixed.dis.as_slice())
        });
        let child = category.to_parent.variable.get(pos).map(|var| var.dis).or_else(|| {
          category
            .to_parent
            .fixed
            .get(child_state)
            .map(|fixed| fixed.dis.as_slice())
        });
        Some(category.row(parent?, child?))
      })
      .collect::<Option<Vec<_>>>();
    if let Some(rows) = rows {
      sites.add(1.0, mixed_row(&categories, &rows, site, weights)?);
    }
  }
  for fixed in categories[0].to_child.fixed.iter() {
    let taken = taken_from_groups.get(&fixed.state).copied().unwrap_or(0);
    let count = fixed.count.saturating_sub(taken);
    if count == 0 {
      continue;
    }
    let Some(site) = state.mixture.group(fixed.state) else {
      continue;
    };
    let rows = categories
      .iter()
      .map(|category| {
        let parent = category.to_child.fixed.get(fixed.state)?;
        let child = category.to_parent.fixed.get(fixed.state)?;
        Some(category.row(&parent.dis, &child.dis))
      })
      .collect::<Option<Vec<_>>>();
    if let Some(rows) = rows {
      sites.add(count.to_f64_exact(), mixed_row(&categories, &rows, site, weights)?);
    }
  }
  let rate_eigvals = concatenate(
    Axis(0),
    &categories
      .iter()
      .map(|category| category.gtr.rate_eigvals())
      .collect_vec()
      .iter()
      .map(|values| values.view())
      .collect_vec(),
  )?;
  Ok(OptimizationContribution::Sparse(
    optimize::sparse::PartitionContribution::mixture(sites.into_vec(), rate_eigvals),
  ))
}

struct EdgeCategory<'a> {
  gtr: &'a GTR,
  to_child: &'a SparseSeqDistribution,
  to_parent: &'a SparseSeqDistribution,
  decay: Array1<f64>,
}

impl EdgeCategory<'_> {
  fn row(&self, parent: &[f64], child: &[f64]) -> Array1<f64> {
    aview1(child).dot(self.gtr.v()) * aview1(parent).dot(&self.gtr.v_inv().t())
  }
}

fn variable_label(categories: &[EdgeCategory<'_>], pos: usize, parent_side: bool) -> Option<AsciiChar> {
  categories.iter().find_map(|category| {
    let (own, other) = if parent_side {
      (category.to_child, category.to_parent)
    } else {
      (category.to_parent, category.to_child)
    };
    own
      .variable
      .get(pos)
      .or_else(|| other.variable.get(pos))
      .map(|var| var.state)
  })
}

fn mixed_row(
  categories: &[EdgeCategory<'_>],
  rows: &[Array1<f64>],
  site: &MixedSite,
  weights: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  let exponents: CategoryValues = izip!(categories, rows, &site.site_log_lh, weights)
    .map(|(category, row, &site_log_lh, &weight)| scale_exponent(weight, site_log_lh, row.view(), &category.decay))
    .collect();
  let max = exponents.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let scaled = izip!(rows, &exponents)
    .map(|(row, &exponent)| {
      let scale = if max.is_finite() { (exponent - max).exp() } else { 0.0 };
      row * scale
    })
    .collect_vec();
  Ok(concatenate(
    Axis(0),
    &scaled.iter().map(|values| values.view()).collect_vec(),
  )?)
}

#[expect(
  clippy::needless_pass_by_value,
  reason = "array views are Copy; clippy misses the Copy impl of ndarray 0.17 views"
)]
fn scale_exponent(weight: f64, site_log_lh: f64, row: ArrayView1<'_, f64>, decay: &Array1<f64>) -> f64 {
  if site_log_lh == f64::NEG_INFINITY {
    return f64::NEG_INFINITY;
  }
  let site_lh = row.dot(decay).max(f64::MIN_POSITIVE);
  weight.ln() + site_log_lh - site_lh.ln()
}
