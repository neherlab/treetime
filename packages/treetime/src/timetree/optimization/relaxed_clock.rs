use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::reroot::EdgeMergeInfo;
use treetime_graph::schedule::GraphSchedule;

const DEFAULT_RELAX: f64 = 1.0;

const UNIT_GAMMA: f64 = 1.0;

pub(crate) fn apply_relaxed_clock(
  graph: &Graph,
  branch_lengths: &EdgeMap<Option<f64>>,
  params: &[f64],
  one_mutation: f64,
  clock_rate: f64,
  date_branch_lengths: &EdgeMap<Option<f64>>,
) -> Result<EdgeMap<f64>, Report> {
  let RelaxedClockPrior { slack, coupling } = RelaxedClockPrior::of(params);

  let c = 1.0 / one_mutation;

  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;

  let mut coeffs: NodeMap<RelaxedClockCoeffs> = NodeMap::new();

  for &key in schedule.postorder() {
    let mut node_coeffs = RelaxedClockCoeffs::default();

    let (opt_len, act_len) = match schedule.parents(key).first() {
      Some(&(_, edge_key)) => {
        let opt_len = branch_lengths[&edge_key].unwrap_or(0.0);
        let act_len = date_branch_lengths[&edge_key].map_or(opt_len, |duration| duration * clock_rate);
        (opt_len, act_len)
      },
      None => (one_mutation, one_mutation),
    };

    let denom = opt_len + one_mutation;
    node_coeffs.k2 = slack + c * act_len * act_len / denom;
    node_coeffs.k1 = -2.0 * (c * act_len * opt_len / denom + slack);

    for (child_key, _) in schedule.children(key) {
      if let Some(child_coeffs) = coeffs.get(child_key) {
        let denom = coupling + child_coeffs.k2;
        if denom.abs() > 1e-10 {
          let ratio = coupling / denom;
          node_coeffs.k2 += coupling * (1.0 - ratio).powi(2) + child_coeffs.k2 * ratio.powi(2);
          node_coeffs.k1 += coupling * (1.0 - ratio) * child_coeffs.k1 / denom
            - coupling * child_coeffs.k1 * child_coeffs.k2 / denom.powi(2)
            + coupling * child_coeffs.k1 / denom;
        }
      }
    }

    coeffs.insert(key, node_coeffs);
  }

  let mut gammas: NodeMap<f64> = NodeMap::new();

  for &key in schedule.preorder() {
    let node_coeffs = coeffs.get(key).cloned().unwrap_or_default();

    let gamma = match schedule.parents(key).first() {
      None => {
        if node_coeffs.k2.abs() > 1e-10 {
          (-0.5 * node_coeffs.k1 / node_coeffs.k2).max(0.1)
        } else {
          1.0
        }
      },
      Some((parent_key, _)) => {
        let parent_gamma = gammas.get(parent_key).copied().unwrap_or(1.0);

        let denom = coupling + node_coeffs.k2;
        if denom.abs() > 1e-10 {
          ((coupling * parent_gamma - 0.5 * node_coeffs.k1) / denom).max(0.1)
        } else {
          1.0
        }
      },
    };

    gammas.insert(key, gamma);
  }

  let edge_gammas = graph
    .get_edges()
    .map(|edge| (edge.key(), gammas[&edge.target()]))
    .collect();
  Ok(edge_gammas)
}

pub(crate) struct RelaxedClockPrior {
  pub slack: f64,
  pub coupling: f64,
}

impl RelaxedClockPrior {
  pub(crate) fn of(params: &[f64]) -> Self {
    Self {
      slack: params.first().copied().unwrap_or(DEFAULT_RELAX),
      coupling: params.get(1).copied().unwrap_or(DEFAULT_RELAX),
    }
  }
}

#[derive(Clone, Default)]
struct RelaxedClockCoeffs {
  k1: f64,
  k2: f64,
}

#[must_use]
pub(crate) fn unit_gammas(graph: &Graph) -> EdgeMap<f64> {
  graph.get_edges().map(|edge| (edge.key(), UNIT_GAMMA)).collect()
}

#[must_use]
pub(crate) fn reconcile_gammas(graph: &Graph, mut gammas: EdgeMap<f64>, merges: &[EdgeMergeInfo]) -> EdgeMap<f64> {
  for merge in merges {
    gammas.remove(merge.parent_edge_key);
    if let Some(child_gamma) = gammas.remove(merge.child_edge_key) {
      gammas.insert(merge.merged_edge_key, child_gamma);
    }
  }
  graph
    .get_edges()
    .map(|edge| (edge.key(), gammas.get(edge.key()).copied().unwrap_or(UNIT_GAMMA)))
    .collect()
}
