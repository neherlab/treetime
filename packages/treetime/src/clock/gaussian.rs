use crate::branch_lengths::branch_length_or_zero;
use crate::timetree::inference::result::NodeTimes;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Node};
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::NodeEdgeMaps;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;

const MIN_PRECISION: f64 = 1e-300;

pub(crate) struct GaussianTree<'g> {
  schedule: GraphSchedule<'g>,
  internal_nodes: Vec<GraphNodeKey>,
  sequence_length: f64,
  nu: BTreeMap<GraphEdgeKey, f64>,
  below: BTreeMap<GraphNodeKey, Below>,
  up: BTreeMap<GraphEdgeKey, Option<Message>>,
  root: Message,
  substitutions: f64,
  events: f64,
  kappa: f64,
  a: f64,
  b: f64,
}

impl<'g> GaussianTree<'g> {
  pub(crate) fn new(
    graph: &'g Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    sequence_length: usize,
    tip_dates: &NodeTimes,
    excluded: &BTreeMap<GraphNodeKey, bool>,
  ) -> Result<Self, Report> {
    if sequence_length == 0 {
      return make_error!("Gaussian clock needs a positive sequence length");
    }
    let sequence_length = sequence_length.to_f64_exact();
    let schedule = GraphSchedule::new(graph)?;
    schedule.require_tree()?;
    let root = schedule.roots()[0];
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
    let NodeEdgeMaps {
      nodes: below,
      edges: up,
    } = backward(&schedule, tips, &nu)?;

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
    let Some(root) = below[&root].down else {
      return make_error!("Gaussian clock needs dated tips below the root");
    };
    let [a, b] = below
      .values()
      .fold([0.0; 2], |acc, node| [acc[0] + node.ab[0], acc[1] + node.ab[1]]);
    Ok(Self {
      schedule,
      internal_nodes: graph.get_nodes().filter(|n| !n.is_leaf()).map(Node::key).collect(),
      sequence_length,
      nu,
      below,
      up,
      root,
      substitutions,
      events: substitutions + (internal_branches + tip_branches).to_f64_exact(),
      kappa: tip_branches.to_f64_exact(),
      a,
      b,
    })
  }

  pub(crate) fn clock(&self, clock_rate: f64) -> Result<GaussianClock, Report> {
    if clock_rate.is_nan() || clock_rate <= 0.0 {
      return make_error!("Gaussian clock needs a positive clock rate, got {clock_rate}");
    }
    let eta = self.eta(clock_rate);
    let total = forward(&self.schedule, &self.nu, &self.below, &self.up)?;
    let stds = rate_and_root_std(self.a, self.kappa, self.root, eta);

    let node_times = self
      .internal_nodes
      .iter()
      .filter_map(|key| {
        let m = total[key].message()?;
        Some((
          *key,
          NodeTimeGaussian {
            mean: m.s - m.r / eta,
            std: m.w.sqrt() / eta,
          },
        ))
      })
      .collect();

    Ok(GaussianClock {
      rate_std: stds.map(|(eta_std, _)| self.mu(eta_std)),
      root_time_std: stds.map(|(_, root_std)| root_std),
      node_times,
    })
  }

  pub(crate) fn rate(&self) -> Result<f64, Report> {
    let (a, b, kappa) = (self.a, self.b, self.kappa);
    if a.is_nan() || a <= 0.0 {
      return make_error!("Gaussian clock rate is undefined: the dates of the tips do not vary along the tree");
    }
    let eta = (b + (b * b + 4.0 * a * kappa).sqrt()) / (2.0 * a);
    Ok(self.mu(eta))
  }

  fn eta(&self, mu: f64) -> f64 {
    mu * self.sequence_length * self.events / self.substitutions
  }

  fn mu(&self, eta: f64) -> f64 {
    eta * self.substitutions / (self.events * self.sequence_length)
  }
}

#[derive(Clone, Debug)]
pub(crate) struct GaussianClock {
  pub rate_std: Option<f64>,
  pub root_time_std: Option<f64>,
  pub node_times: BTreeMap<GraphNodeKey, NodeTimeGaussian>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct NodeTimeGaussian {
  pub mean: f64,
  pub std: f64,
}

fn backward(
  schedule: &GraphSchedule<'_>,
  tips: BTreeMap<GraphNodeKey, Option<Message>>,
  nu: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<NodeEdgeMaps<Below, Option<Message>>, Report> {
  let outputs = schedule.scan_up(tips, |key, tip, from_children| {
    let (node, own) = if schedule.children(key).is_empty() {
      (Below::default(), tip)
    } else {
      let messages = from_children.into_iter().flatten().collect_vec();
      let combined = combine(&messages);
      let node = Below {
        down: combined.map(|(m, _)| m),
        ab: combined.map_or([0.0; 2], |(_, ab)| ab),
      };
      (node, node.down)
    };
    let to_parents = schedule
      .parents(key)
      .iter()
      .map(|&(_, edge_key)| own.map(|m| m.toward_parent(nu[&edge_key])))
      .collect_vec();
    let to_parent = to_parents.first().copied().flatten();
    Ok(((node, to_parents), to_parent))
  })?;
  outputs.split_parent_edges()
}

fn combine(messages: &[Message]) -> Option<(Message, [f64; 2])> {
  let total = messages.iter().fold(Natural::default(), |acc, m| acc.plus(m.natural()));
  let combined = total.message()?;
  let a = messages.iter().map(|m| m.s * m.s / m.w).sum::<f64>() - combined.s * combined.s / combined.w;
  let b = messages.iter().map(|m| m.s * m.r / m.w).sum::<f64>() - combined.s * combined.r / combined.w;
  Some((combined, [a, b]))
}

fn forward(
  schedule: &GraphSchedule<'_>,
  nu: &BTreeMap<GraphEdgeKey, f64>,
  below: &BTreeMap<GraphNodeKey, Below>,
  up: &BTreeMap<GraphEdgeKey, Option<Message>>,
) -> Result<BTreeMap<GraphNodeKey, Natural>, Report> {
  let down = below.iter().map(|(&key, node)| (key, node.down));
  let totals = schedule.scan_down(down, |key, down, parent_totals: Vec<Natural>| {
    let outside = match (schedule.parents(key).first(), parent_totals.first()) {
      (Some(&(_, edge_key)), Some(&parent_total)) => {
        let own_up = up
          .get(&edge_key)
          .copied()
          .flatten()
          .map_or_else(Natural::default, Message::natural);
        parent_total
          .minus(own_up)
          .message()
          .map(|m| m.toward_child(nu[&edge_key]))
      },
      _ => None,
    };
    let node = [down, outside]
      .into_iter()
      .flatten()
      .fold(Natural::default(), |acc, m| acc.plus(m.natural()));
    Ok((node, node))
  })?;
  Ok(totals.into_map())
}

fn rate_and_root_std(a: f64, kappa: f64, root: Message, eta: f64) -> Option<(f64, f64)> {
  let i_rate = kappa / eta.powi(2) + a + root.r.powi(2) / (eta.powi(2) * root.w);
  let i_cross = -root.r / root.w;
  let i_root = eta.powi(2) / root.w;
  let det = i_rate * i_root - i_cross.powi(2);
  (det > 0.0 && i_rate > 0.0).then(|| ((i_root / det).sqrt(), (i_rate / det).sqrt()))
}

#[derive(Clone, Copy, Debug, Default)]
struct Below {
  down: Option<Message>,
  ab: [f64; 2],
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
