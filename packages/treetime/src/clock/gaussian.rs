//! Gaussian message passing on a dated tree at a given clock rate.
//!
//! A branch with `d` substitutions (branch length times sequence length) has an elapsed time approximated by
//! `N(nu / eta, nu / eta^2)` with `nu = d + 1`, the moment match of its Gamma(`d + 1`, `eta`) likelihood under a flat
//! prior on time. The `+1` is the hazard of the event that ends a branch (a split, or sampling at a tip), so `eta` is
//! the rate of the combined process: `eta = mu (D + M + N) / D` for the substitution rate `mu`, with `D` the
//! substitutions, `M` the internal and `N` the terminal branches. A message about a node's time `t` is
//! `(S, R, W)` with mean `S - R / eta` and variance `W / eta^2`; `S`, `R`, `W` do not depend on `eta`. A message
//! crosses a branch toward the parent as `(S, R + nu, W + nu)` and toward the child as `(S, R - nu, W + nu)`; messages
//! about the same time combine by summing their natural parameters `(1/W, S/W, R/W)`.
//!
//! The combination of messages `i` at a node leaves `-(a eta^2 - 2 b eta + c) / 2` in the log-likelihood, with
//! `a = sum S_i^2 / W_i - S^2 / W`, `b = sum S_i R_i / W_i - S R / W`, `c = sum R_i^2 / W_i - R^2 / W`. With the root
//! time optimized, `log L(eta) = kappa ln eta - A eta^2 / 2 + B eta - C / 2`, `kappa` the number of dated tips and `A`,
//! `B`, `C` the sums over nodes. Following `gaussian_timetree.py` and its derivation (neherlab/gaussian-timetree).

use crate::branch_lengths::branch_length_or_zero;
use crate::timetree::inference::result::NodeTimes;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Node};
use treetime_graph::pass::{GraphMapOutputs, GraphPass, GraphPassNodeOutput};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_report};

/// Precision below which a combination of messages carries no information
const MIN_PRECISION: f64 = 1e-300;

/// Time of a node as a Gaussian
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct NodeTimeGaussian {
  pub mean: f64,
  pub std: f64,
}

/// Rate and root uncertainty and the node-time marginals of the Gaussian approximation at a given rate
#[derive(Clone, Debug)]
pub(crate) struct GaussianClock {
  /// Standard deviation of the substitution rate, per site and unit time
  /// `None` when the information matrix of rate and root time is not positive definite
  pub rate_std: Option<f64>,
  pub root_time_std: Option<f64>,
  /// Marginal time of every internal node that has a dated tip below it
  pub node_times: BTreeMap<GraphNodeKey, NodeTimeGaussian>,
}

#[derive(Clone, Copy, Debug)]
struct Message {
  s: f64,
  r: f64,
  w: f64,
}

impl Message {
  fn toward_parent(self, nu: f64) -> Self {
    Self {
      r: self.r + nu,
      w: self.w + nu,
      ..self
    }
  }

  fn toward_child(self, nu: f64) -> Self {
    Self {
      r: self.r - nu,
      w: self.w + nu,
      ..self
    }
  }

  fn natural(self) -> Natural {
    Natural {
      p: 1.0 / self.w,
      ps: self.s / self.w,
      pr: self.r / self.w,
    }
  }
}

/// Natural parameters `(1/W, S/W, R/W)` of a message, which add when messages combine
#[derive(Clone, Copy, Debug, Default)]
struct Natural {
  p: f64,
  ps: f64,
  pr: f64,
}

impl Natural {
  fn message(self) -> Option<Message> {
    (self.p > MIN_PRECISION).then(|| Message {
      s: self.ps / self.p,
      r: self.pr / self.p,
      w: 1.0 / self.p,
    })
  }

  fn plus(self, other: Self) -> Self {
    Self {
      p: self.p + other.p,
      ps: self.ps + other.ps,
      pr: self.pr + other.pr,
    }
  }

  fn minus(self, other: Self) -> Self {
    Self {
      p: self.p - other.p,
      ps: self.ps - other.ps,
      pr: self.pr - other.pr,
    }
  }
}

/// Combination of the messages at one node, and the coefficient `a` it adds to the log-likelihood; `b` and `c` shift
/// the rate optimum but not its curvature, so they are not needed at a given rate
fn combine(messages: &[Message]) -> Option<(Message, f64)> {
  let total = messages.iter().fold(Natural::default(), |acc, m| acc.plus(m.natural()));
  let combined = total.message()?;
  let a = messages.iter().map(|m| m.s * m.s / m.w).sum::<f64>() - combined.s * combined.s / combined.w;
  Some((combined, a))
}

#[derive(Clone, Copy, Debug, Default)]
struct Below {
  /// Combined message from the children, about the node's own time
  down: Option<Message>,
  a: f64,
}

/// Gaussian timetree at substitution rate `clock_rate` (per site and unit time). A tip takes part when it has a date in
/// `tip_dates` and is not flagged in `excluded`; the others send no message.
pub(crate) fn gaussian_clock(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  sequence_length: f64,
  tip_dates: &NodeTimes,
  excluded: &BTreeMap<GraphNodeKey, bool>,
  clock_rate: f64,
) -> Result<GaussianClock, Report> {
  if !(sequence_length > 0.0 && clock_rate > 0.0) {
    return make_error!(
      "Gaussian clock needs a positive sequence length and clock rate, got {sequence_length} and {clock_rate}"
    );
  }
  let pass = GraphPass::new(graph)?;
  // a tip's own message: its date, known exactly
  let tips: BTreeMap<GraphNodeKey, Option<Message>> = graph
    .get_nodes()
    .map(|n| {
      let included = n.is_leaf() && !excluded.get(&n.key()).copied().unwrap_or(false);
      let date = tip_dates.get(&n.key()).copied().flatten().filter(|_| included);
      (n.key(), date.map(|s| Message { s, r: 0.0, w: 0.0 }))
    })
    .collect();
  let nu: BTreeMap<GraphEdgeKey, f64> = graph
    .get_edges()
    .map(|e| {
      (
        e.key(),
        branch_length_or_zero(branch_lengths, e.key()) * sequence_length + 1.0,
      )
    })
    .collect();

  let GraphMapOutputs {
    nodes: below,
    edges: up,
  } = backward(&pass, &tips, &nu)?;
  let total = forward(&pass, &tips, &nu, &below, &up)?;

  let mut substitutions = 0.0;
  let mut internal_branches = 0_usize;
  let mut tip_branches = 0_usize;
  for edge in graph.get_edges() {
    if up.get(&edge.key()).copied().flatten().is_none() {
      continue;
    }
    substitutions += nu[&edge.key()] - 1.0;
    if graph.get_node(edge.target()).is_some_and(Node::is_leaf) {
      tip_branches += 1;
    } else {
      internal_branches += 1;
    }
  }
  if substitutions <= 0.0 {
    return make_error!("Gaussian clock needs substitutions on the dated part of the tree");
  }
  let events = substitutions + (internal_branches + tip_branches).to_f64_exact();
  let eta = clock_rate * sequence_length * events / substitutions;

  let root_key = graph.root_key()?;
  let Some(root) = below[&root_key].down else {
    return make_error!("Gaussian clock needs dated tips below the root");
  };
  let a = below.values().map(|node| node.a).sum();
  let kappa = tip_branches.to_f64_exact();
  let stds = rate_and_root_std(a, kappa, root, eta);

  let node_times = graph
    .get_nodes()
    .filter(|n| !n.is_leaf())
    .filter_map(|n| {
      let m = total[&n.key()].message()?;
      Some((
        n.key(),
        NodeTimeGaussian {
          mean: m.s - m.r / eta,
          std: m.w.sqrt() / eta,
        },
      ))
    })
    .collect();

  Ok(GaussianClock {
    rate_std: stds.map(|(eta_std, _)| eta_std * substitutions / (events * sequence_length)),
    root_time_std: stds.map(|(_, root_std)| root_std),
    node_times,
  })
}

/// Messages toward the root: per node the combination from its children, per edge the message it carries upward
fn backward(
  pass: &GraphPass,
  tips: &BTreeMap<GraphNodeKey, Option<Message>>,
  nu: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<GraphMapOutputs<Below, Option<Message>>, Report> {
  pass.map_backward(
    tips,
    nu,
    |key| Err(make_internal_report!("Gaussian clock is missing node {key}")),
    |context| {
      let (node, own) = if context.is_leaf {
        (Below::default(), *context.input)
      } else {
        let messages: Vec<Message> = context
          .children
          .iter()
          .filter_map(|child| child.edge.copied().flatten())
          .collect();
        let combined = combine(&messages);
        let node = Below {
          down: combined.map(|(m, _)| m),
          a: combined.map_or(0.0, |(_, a)| a),
        };
        (node, node.down)
      };
      let parent_message = context.parent_edge.map(|(_, &nu)| own.map(|m| m.toward_parent(nu)));
      Ok(GraphPassNodeOutput { node, parent_message })
    },
  )
}

/// Messages toward the tips, combined with `below` into the natural parameters of each node's marginal: the outside
/// message of a node is its parent's total without the node's own upward message, carried across the edge
fn forward(
  pass: &GraphPass,
  tips: &BTreeMap<GraphNodeKey, Option<Message>>,
  nu: &BTreeMap<GraphEdgeKey, f64>,
  below: &BTreeMap<GraphNodeKey, Below>,
  up: &BTreeMap<GraphEdgeKey, Option<Message>>,
) -> Result<BTreeMap<GraphNodeKey, Natural>, Report> {
  let GraphMapOutputs { nodes, .. } = pass.map_forward(
    tips,
    nu,
    |key| Err(make_internal_report!("Gaussian clock is missing node {key}")),
    |context| {
      let outside = match (context.parent_edge, context.parent) {
        (Some((edge, &nu)), Some(parent_total)) => {
          let own_up = up
            .get(&edge)
            .copied()
            .flatten()
            .map_or_else(Natural::default, Message::natural);
          let parent_total: Natural = *parent_total;
          parent_total.minus(own_up).message().map(|m| m.toward_child(nu))
        },
        _ => None,
      };
      let node = [below[&context.key].down, outside]
        .into_iter()
        .flatten()
        .fold(Natural::default(), |acc, m| acc.plus(m.natural()));
      Ok(GraphPassNodeOutput::<Natural, ()> {
        node,
        parent_message: None,
      })
    },
  )?;
  Ok(nodes)
}

/// Standard deviations of `eta` and the root time from the observed information of the fixed-root log-likelihood
/// `kappa ln eta - A' eta^2 / 2 + B' eta - C' / 2` with `A' = A + (t_r - S_r)^2 / W_r`, at `t_r = S_r - R_r / eta`:
///
/// `I = [[kappa / eta^2 + A + R_r^2 / (eta^2 W_r), -R_r / W_r], [-R_r / W_r, eta^2 / W_r]]`, `cov = I^-1`.
fn rate_and_root_std(a: f64, kappa: f64, root: Message, eta: f64) -> Option<(f64, f64)> {
  let i_rate = kappa / eta.powi(2) + a + root.r.powi(2) / (eta.powi(2) * root.w);
  let i_cross = -root.r / root.w;
  let i_root = eta.powi(2) / root.w;
  let det = i_rate * i_root - i_cross.powi(2);
  (det > 0.0 && i_rate > 0.0).then(|| ((i_root / det).sqrt(), (i_rate / det).sqrt()))
}
