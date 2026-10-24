use crate::alphabet::alphabet::Alphabet;
use crate::build::mat::Mat;
use crate::build::pipeline::{BuildParams, BuildStats, InitialTree, SampleOrder, SearchMode, place_samples};
use crate::build::query::Query;
use crate::build::spr::spr_round;
use eyre::Report;
use itertools::Itertools;
use rand::seq::SliceRandom;
use std::collections::{BTreeMap, VecDeque};
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::cast::usize_to_u64;
use treetime_utils::sync::random::get_random_number_generator;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RebuildParams {
  pub radius: usize,
  pub max_leaves: usize,
  pub orders: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RebuildStats {
  pub tried: usize,
  pub improved: usize,
  pub gain: usize,
  pub too_large: usize,
}

struct Ball {
  top: GraphNodeKey,
  members: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
  interior: Vec<GraphNodeKey>,
}

impl Ball {
  fn around(mat: &Mat, center: GraphNodeKey, params: &RebuildParams) -> Option<Self> {
    let mut distance: BTreeMap<GraphNodeKey, usize> = BTreeMap::from([(center, 0)]);
    let mut queue = VecDeque::from([center]);
    while let Some(node) = queue.pop_front() {
      let next = distance[&node] + 1;
      if next > params.radius {
        continue;
      }
      for (neighbor, _) in mat.children(node).chain(mat.parent(node)) {
        if let std::collections::btree_map::Entry::Vacant(entry) = distance.entry(neighbor) {
          entry.insert(next);
          queue.push_back(neighbor);
        }
      }
      if distance.len() > 4 * params.max_leaves {
        return None;
      }
    }
    let top = *distance.keys().find(|&&node| {
      mat
        .parent(node)
        .is_none_or(|(parent, _)| !distance.contains_key(&parent))
    })?;
    let members = distance.keys().copied().filter(|&node| node != top).collect_vec();
    let (leaves, interior): (Vec<_>, Vec<_>) = members
      .iter()
      .partition(|&&node| mat.is_leaf(node) || distance[&node] == params.radius);
    (leaves.len() <= params.max_leaves).then_some(Self {
      top,
      members,
      leaves,
      interior,
    })
  }

  fn cost(&self, mat: &Mat) -> usize {
    self
      .members
      .iter()
      .filter_map(|&node| mat.parent(node))
      .map(|(_, edge)| mat.edge_subs(edge).len())
      .sum()
  }
}

pub(crate) fn rebuild_round(
  mat: &mut Mat,
  seqs: &[Seq],
  alphabet: &Alphabet,
  params: &RebuildParams,
) -> Result<RebuildStats, Report> {
  let mut stats = RebuildStats::default();
  for center in mat.preorder() {
    if !mat.has_node(center) || mat.is_leaf(center) {
      continue;
    }
    let Some(ball) = Ball::around(mat, center, params) else {
      stats.too_large += 1;
      continue;
    };
    if ball.interior.is_empty() || ball.leaves.len() < 3 {
      continue;
    }
    stats.tried += 1;
    let before = ball.cost(mat);
    let local = best_local_tree(mat, &ball, seqs, alphabet, params.orders, center)?;
    if local.parsimony() < before {
      stats.improved += 1;
      stats.gain += before - local.parsimony();
      let old_states = ball
        .leaves
        .iter()
        .filter(|&&leaf| !mat.is_leaf(leaf))
        .map(|&leaf| (leaf, mat.node_state(leaf)))
        .collect_vec();
      mat.splice(ball.top, &ball.members, &ball.interior, &local, &ball.leaves)?;
      let allows = |sample: usize, pos: usize, state| {
        let observed = seqs[sample][pos];
        alphabet.is_unknown(observed)
          || alphabet.is_gap(observed)
          || alphabet
            .canonical_states(observed)
            .is_ok_and(|states| states.is_empty() || states.contains(state))
      };
      for (leaf, old_state) in &old_states {
        mat.restate(*leaf, old_state, allows)?;
      }
    }
  }
  Ok(stats)
}

fn best_local_tree(
  mat: &Mat,
  ball: &Ball,
  seqs: &[Seq],
  alphabet: &Alphabet,
  orders: usize,
  center: GraphNodeKey,
) -> Result<Mat, Report> {
  let top_state = mat.node_state(ball.top);
  let mut frame = mat.reference().clone();
  for &(pos, state) in &top_state {
    frame[pos] = state;
  }
  let queries = ball
    .leaves
    .iter()
    .map(|&leaf| {
      mat.sample_of(leaf).map_or_else(
        || {
          let state = mat.node_state(leaf);
          Query::from_state_and_sets(
            &mat.relative_state(&state, &top_state),
            &mat.plurality_sets(leaf, &state),
            &frame,
          )
        },
        |sample| Query::from_seq(&seqs[sample], &frame, alphabet),
      )
    })
    .collect_vec();
  let params = BuildParams {
    search: SearchMode::Exhaustive,
    batch_size: 1,
    order: SampleOrder::Input,
    spr_rounds: 0,
    spr_slack: None,
    fitch_refresh: false,
    compare_exhaustive: false,
    rebuild: None,
    reprune: None,
    scaffold: None,
    initial: InitialTree::Placement,
  };
  let mut rng = get_random_number_generator(usize_to_u64(center.as_usize()));
  let mut order = (0..queries.len()).collect_vec();
  let mut best: Option<Mat> = None;
  for attempt in 0..orders.max(1) {
    if attempt > 0 {
      order.shuffle(&mut rng);
    }
    let mut local = Mat::new(frame.clone());
    place_samples(&mut local, &queries, &order, &params, &mut BuildStats::default())?;
    loop {
      let before = local.parsimony();
      spr_round(&mut local, &queries, None)?;
      if local.parsimony() >= before {
        break;
      }
    }
    if best.as_ref().is_none_or(|best| local.parsimony() < best.parsimony()) {
      best = Some(local);
    }
  }
  best.ok_or_else(|| crate::make_internal_report!("Local rebuild produced no tree"))
}
