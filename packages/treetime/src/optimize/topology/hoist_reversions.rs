use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{SparseNodeObs, SparseNodeState};
use crate::seq::indel::{InDel, compose_indels, sort_indels};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::izip;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AsciiChar;

pub fn count_child_reversions(
  sparse: &[PartitionMarginalSparse],
  parent_edge_key: GraphEdgeKey,
  sibling_edge_key: Option<GraphEdgeKey>,
  child_edge_key: GraphEdgeKey,
) -> usize {
  sparse
    .iter()
    .map(|partition| {
      let parent_subs: &[Sub] = partition
        .obs_edges
        .get(&parent_edge_key)
        .map_or(&[], |e| e.fitch_subs());
      let child_subs: &[Sub] = partition.obs_edges.get(&child_edge_key).map_or(&[], |e| e.fitch_subs());
      match sibling_edge_key {
        Some(sibling_edge_key) => {
          let sibling_subs: &[Sub] = partition
            .obs_edges
            .get(&sibling_edge_key)
            .map_or(&[], |e| e.fitch_subs());
          let (augmented, _) = augment_parent_with_sibling(parent_subs, sibling_subs);
          count_reversions(&augmented, child_subs)
        },
        None => count_reversions(parent_subs, child_subs),
      }
    })
    .sum()
}

fn augment_parent_with_sibling(parent_subs: &[Sub], sibling_subs: &[Sub]) -> (Vec<Sub>, BTreeSet<usize>) {
  let parent_positions: BTreeSet<usize> = parent_subs.iter().map(Sub::pos).collect();
  let mut augmented = parent_subs.to_vec();
  let mut sibling_sourced = BTreeSet::new();
  for sibling_sub in sibling_subs {
    if parent_positions.contains(&sibling_sub.pos()) {
      continue;
    }
    let mut inverted = sibling_sub.clone();
    inverted.invert();
    augmented.push(inverted);
    sibling_sourced.insert(sibling_sub.pos());
  }
  augmented.sort_by_key(Sub::pos);
  (augmented, sibling_sourced)
}

pub fn slide_bifurcating_root_for_child(
  sparse: &mut [PartitionMarginalSparse],
  node_states: &mut [BTreeMap<GraphNodeKey, SparseNodeState>],
  root_key: GraphNodeKey,
  parent_edge_key: GraphEdgeKey,
  sibling_edge_key: GraphEdgeKey,
  child_edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  for (partition, node_states) in izip!(sparse.iter_mut(), node_states.iter_mut()) {
    let parent_subs = partition
      .obs_edges
      .get(&parent_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let sibling_subs = partition
      .obs_edges
      .get(&sibling_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let child_by_pos: BTreeMap<usize, (AsciiChar, AsciiChar)> = partition
      .obs_edges
      .get(&child_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec())
      .iter()
      .map(|c| (c.pos(), (c.reff(), c.qry())))
      .collect();
    let parent_positions: BTreeSet<usize> = parent_subs.iter().map(Sub::pos).collect();

    let mut hoisted_parent = parent_subs;
    let mut remaining_sibling = Vec::new();
    let mut slid_any = false;
    for sibling_sub in sibling_subs {
      let reverted_by_child = child_by_pos
        .get(&sibling_sub.pos())
        .is_some_and(|&(reff, qry)| reff == sibling_sub.reff() && qry == sibling_sub.qry());
      if !parent_positions.contains(&sibling_sub.pos()) && reverted_by_child {
        let pos = sibling_sub.pos();
        partition.root_sequence[pos] = sibling_sub.qry();
        if let Some(root_node) = node_states.get_mut(&root_key) {
          if pos < root_node.sequence.len() {
            root_node.sequence[pos] = sibling_sub.qry();
          }
        }
        let mut inverted = sibling_sub.clone();
        inverted.invert();
        hoisted_parent.push(inverted);
        slid_any = true;
      } else {
        remaining_sibling.push(sibling_sub);
      }
    }

    if !slid_any {
      continue;
    }
    hoisted_parent.sort_by_key(Sub::pos);

    if let Some(sibling_edge) = partition.obs_edges.get_mut(&sibling_edge_key) {
      sibling_edge.set_fitch_subs(remaining_sibling);
    }
    partition
      .obs_edges
      .entry(parent_edge_key)
      .or_default()
      .set_fitch_subs(hoisted_parent);
  }
  Ok(())
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn hoist_reverting_child(
  graph: &mut Graph,
  sparse: &mut [PartitionMarginalSparse],
  parent_edge_key: GraphEdgeKey,
  child_edge_key: GraphEdgeKey,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<GraphNodeKey, Report> {
  let u_key = graph.get_source_node_key(parent_edge_key)?;

  let mut splits = Vec::with_capacity(sparse.len());
  let mut total_parent_subs = 0_usize;
  let mut total_hoisted_subs = 0_usize;
  for family in sparse.iter() {
    let obs_edges = &family.obs_edges;
    let parent_subs = obs_edges
      .get(&parent_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let child_subs = obs_edges
      .get(&child_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let parent_indels = obs_edges.get(&parent_edge_key).map_or(Vec::new(), |e| e.indels.clone());
    let child_indels = obs_edges.get(&child_edge_key).map_or(Vec::new(), |e| e.indels.clone());

    let sub_split = split_subs(&parent_subs, &child_subs)?;
    let indel_split = split_indels(&parent_indels, &child_indels);

    total_parent_subs += parent_subs.len();
    total_hoisted_subs += sub_split.hoisted.len();

    splits.push(EdgeSplit {
      hoisted: sub_split.hoisted,
      kept: sub_split.kept,
      composed: sub_split.composed,
      indels: indel_split,
    });
  }

  let bl_uv = branch_lengths[&parent_edge_key].unwrap_or(0.0);
  let bl_vc = branch_lengths[&child_edge_key].unwrap_or(0.0);
  let bl_un = if total_parent_subs > 0 {
    bl_uv * (total_hoisted_subs as f64) / (total_parent_subs as f64)
  } else {
    0.0
  };
  let bl_nv = bl_uv - bl_un;
  let bl_nc = bl_nv + bl_vc;

  let n_key = graph.add_node();
  let un_edge_key = graph.add_edge(u_key, n_key)?;
  branch_lengths.insert(un_edge_key, Some(bl_un));
  graph.reparent_edge(parent_edge_key, n_key)?;
  graph.reparent_edge(child_edge_key, n_key)?;
  branch_lengths.insert(parent_edge_key, Some(bl_nv));
  branch_lengths.insert(child_edge_key, Some(bl_nc));

  for (partition, split) in sparse.iter_mut().zip(splits) {
    let mut node_n = SparseNodeObs::empty(&partition.alphabet);
    node_n.composition = partition.obs_nodes[&u_key].composition.clone();
    partition.obs_nodes.entry(n_key).or_insert(node_n);

    let un_edge = partition.obs_edges.entry(un_edge_key).or_default();
    un_edge.set_fitch_subs(split.hoisted);
    un_edge.indels = split.indels.hoisted;

    let nv_edge = partition.obs_edges.entry(parent_edge_key).or_default();
    nv_edge.set_fitch_subs(split.kept);
    nv_edge.indels = split.indels.kept;

    let nc_edge = partition.obs_edges.entry(child_edge_key).or_default();
    nc_edge.set_fitch_subs(split.composed);
    nc_edge.indels = split.indels.composed;
  }

  Ok(n_key)
}

struct SubSplit {
  hoisted: Vec<Sub>,
  kept: Vec<Sub>,
  composed: Vec<Sub>,
}

fn split_subs(parent_subs: &[Sub], child_subs: &[Sub]) -> Result<SubSplit, Report> {
  debug_assert!(
    parent_subs.is_sorted_by(|a, b| a.pos() < b.pos()),
    "parent_subs not sorted by unique position"
  );
  debug_assert!(
    child_subs.is_sorted_by(|a, b| a.pos() < b.pos()),
    "child_subs not sorted by unique position"
  );

  let mut hoisted = Vec::new();
  let mut kept = Vec::new();
  let mut composed = Vec::new();
  let mut pi = 0;
  let mut ci = 0;

  while pi < parent_subs.len() && ci < child_subs.len() {
    let ps = &parent_subs[pi];
    let cs = &child_subs[ci];
    match ps.pos().cmp(&cs.pos()) {
      Ordering::Less => {
        hoisted.push(ps.clone());
        pi += 1;
      },
      Ordering::Greater => {
        composed.push(cs.clone());
        ci += 1;
      },
      Ordering::Equal => {
        debug_assert_eq!(
          ps.qry(),
          cs.reff(),
          "Substitution chain broken at position {}: parent produces {} but child expects {}",
          ps.pos(),
          ps.qry(),
          cs.reff()
        );
        kept.push(ps.clone());
        if ps.reff() != cs.qry() {
          composed.push(Sub::new(ps.reff(), ps.pos(), cs.qry())?);
        }
        pi += 1;
        ci += 1;
      },
    }
  }

  hoisted.extend_from_slice(&parent_subs[pi..]);
  composed.extend_from_slice(&child_subs[ci..]);

  Ok(SubSplit {
    hoisted,
    kept,
    composed,
  })
}

fn count_reversions(parent_subs: &[Sub], child_subs: &[Sub]) -> usize {
  debug_assert!(
    parent_subs.is_sorted_by(|a, b| a.pos() < b.pos()),
    "parent_subs not sorted by unique position"
  );
  debug_assert!(
    child_subs.is_sorted_by(|a, b| a.pos() < b.pos()),
    "child_subs not sorted by unique position"
  );

  let mut count = 0;
  let mut pi = 0;
  let mut ci = 0;
  while pi < parent_subs.len() && ci < child_subs.len() {
    let ps = &parent_subs[pi];
    let cs = &child_subs[ci];
    match ps.pos().cmp(&cs.pos()) {
      Ordering::Less => pi += 1,
      Ordering::Greater => ci += 1,
      Ordering::Equal => {
        if ps.reff() == cs.qry() {
          count += 1;
        }
        pi += 1;
        ci += 1;
      },
    }
  }
  count
}

struct IndelSplit {
  hoisted: Vec<InDel>,
  kept: Vec<InDel>,
  composed: Vec<InDel>,
}

fn split_indels(parent_indels: &[InDel], child_indels: &[InDel]) -> IndelSplit {
  if indels_interact(parent_indels, child_indels) {
    let mut parent = parent_indels.to_vec();
    let mut child = child_indels.to_vec();
    sort_indels(&mut parent);
    sort_indels(&mut child);
    let composed = compose_indels(&parent, &child);
    IndelSplit {
      hoisted: Vec::new(),
      kept: parent,
      composed,
    }
  } else {
    IndelSplit {
      hoisted: parent_indels.to_vec(),
      kept: Vec::new(),
      composed: child_indels.to_vec(),
    }
  }
}

fn indels_interact(parent_indels: &[InDel], child_indels: &[InDel]) -> bool {
  parent_indels.iter().any(|p| {
    child_indels
      .iter()
      .any(|c| ranges_overlap_or_adjacent(p.range, c.range))
  })
}

fn ranges_overlap_or_adjacent((a_lo, a_hi): (usize, usize), (b_lo, b_hi): (usize, usize)) -> bool {
  a_lo <= b_hi && b_lo <= a_hi
}

struct EdgeSplit {
  hoisted: Vec<Sub>,
  kept: Vec<Sub>,
  composed: Vec<Sub>,
  indels: IndelSplit,
}
