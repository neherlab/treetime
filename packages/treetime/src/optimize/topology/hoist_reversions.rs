use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::dense::{DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use crate::partition::storage::sparse::SparseNodePartition;
use crate::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::seq::indel::{InDel, compose_indels, sort_indels};
use crate::seq::mutation::Sub;
use eyre::Report;
use parking_lot::RwLock;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AsciiChar;

/// Count the substitution reversions a child edge applies to its parent edge, summed
/// across sparse partitions.
///
/// A position is a reversion when the parent edge carries $a \to b$ and the child edge
/// carries the exact inverse $b \to a$, so composing the two edges cancels the change.
/// This is the gain $\lvert R\rvert$ of hoisting that child (see [`hoist_reverting_child`]);
/// each reversion is one mutation the move removes.
///
/// When `sibling_edge_key` is `Some`, the node sits directly under a bifurcating (degree-2)
/// root and the parent set is augmented with the inverted sibling substitutions at
/// parent-empty positions (see [`augment_parent_with_sibling`]), so a reversion split across
/// the root is counted. The augmentation is scoring only; [`slide_bifurcating_root_for_child`]
/// materializes it before the hoist.
pub(crate) fn count_child_reversions(
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  parent_edge_key: GraphEdgeKey,
  sibling_edge_key: Option<GraphEdgeKey>,
  child_edge_key: GraphEdgeKey,
) -> usize {
  sparse
    .iter()
    .map(|partition| {
      let partition = partition.read_arc();
      let parent_subs: &[Sub] = partition.edges.get(&parent_edge_key).map_or(&[], |e| e.fitch_subs());
      let child_subs: &[Sub] = partition.edges.get(&child_edge_key).map_or(&[], |e| e.fitch_subs());
      match sibling_edge_key {
        Some(sibling_edge_key) => {
          let sibling_subs: &[Sub] = partition.edges.get(&sibling_edge_key).map_or(&[], |e| e.fitch_subs());
          let (augmented, _) = augment_parent_with_sibling(parent_subs, sibling_subs);
          count_reversions(&augmented, child_subs)
        },
        None => count_reversions(parent_subs, child_subs),
      }
    })
    .sum()
}

/// Build the parent-edge substitution set as seen across a bifurcating root.
///
/// A degree-2 root is a reversible pass-through: the edge above a child `v` is one half of the
/// single unrooted edge `v--s` to the sibling `s`. At a site where `v` equals the root but `s`
/// differs, the distinguishing substitution sits on the sibling edge (`root_char -> s_char`),
/// not on `v`'s parent edge. The effective parent substitution there is its inverse
/// (`s_char -> v_char`, with `v_char` the root char), which a child of `v` can revert.
///
/// Returns the augmented, position-sorted parent substitutions and the set of positions
/// contributed by the sibling. A site carried by both edges is skipped rather than composed:
/// under Fitch a degree-2 root state is always one of its two children's states, so no site is
/// ever non-empty on both root edges, and skipping keeps the augmentation exactly invertible by
/// the slide, which only rewrites parent-empty sites.
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

/// Re-root a bifurcating root toward the sibling at the positions the chosen child reverts.
///
/// For each sparse partition, a site qualifies when the sibling edge carries a substitution
/// there, the parent edge is empty there, and `child_edge_key` carries the same substitution
/// (the child changes the root char to the sibling char). Such a site is re-rooted onto the
/// sibling's state: the root char becomes `s_char`, the sibling substitution is dropped, and its
/// inverse (`s_char -> v_char`) is added to the parent edge. The reversion is then local to the
/// parent edge and [`hoist_reverting_child`] removes it unchanged.
///
/// The re-rooting preserves every leaf and internal MAP sequence and the total mutation count;
/// it only moves one substitution from the sibling half of the root edge to the parent half. It
/// is therefore substitution-count neutral on its own, so the caller MUST follow it immediately
/// with the hoist that consumes the exposed reversion. Restricting to a degree-2 root guarantees
/// both root edges are rewritten consistently, because the root has exactly these two children.
///
/// The root's clamped MAP sequence (`root_sequence` and the root node's `seq.sequence`) is moved
/// in lock-step, so the marginal pass that reads the clamp on the next optimizer iteration stays
/// consistent with the rewritten edges.
pub(crate) fn slide_bifurcating_root_for_child(
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  root_key: GraphNodeKey,
  parent_edge_key: GraphEdgeKey,
  sibling_edge_key: GraphEdgeKey,
  child_edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  for partition in sparse {
    let mut partition = partition.write_arc();

    let parent_subs = partition
      .edges
      .get(&parent_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let sibling_subs = partition
      .edges
      .get(&sibling_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let child_by_pos: BTreeMap<usize, (AsciiChar, AsciiChar)> = partition
      .edges
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
        // Keep the root node's clamped MAP sequence in step with `root_sequence`. The forward
        // pass only clamps the root MAP when it is empty, so a bare `root_sequence` edit would
        // leave the populated root MAP stale for the next iteration's marginal reconstruction.
        if let Some(root_node) = partition.nodes.get_mut(&root_key) {
          if pos < root_node.seq.sequence.len() {
            root_node.seq.sequence[pos] = sibling_sub.qry();
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

    if let Some(sibling_edge) = partition.edges.get_mut(&sibling_edge_key) {
      sibling_edge.set_fitch_subs(remaining_sibling);
    }
    partition
      .edges
      .entry(parent_edge_key)
      .or_default()
      .set_fitch_subs(hoisted_parent);
  }
  Ok(())
}

/// Insert a new node $N$ between $u$ and $v$ that groups $v$ and one reverting child $c$,
/// hoisting the non-reverted substitutions above $N$.
///
/// Let $e_p = u \to v$ carry the parent substitutions $M_v$ and $e_c = v \to c$ carry $M_c$.
/// Each position of $M_v$ falls into one of three disjoint sets relative to $M_c$:
///
/// - $T$: positions untouched by the child. Hoisted onto $u \to N$.
/// - $H$: chained positions ($a \to b$ then $b \to d$, $d \neq a$). Kept on $N \to v$; the
///   composed $a \to d$ moves to $N \to c$.
/// - $R$: reverted positions ($a \to b$ then $b \to a$). Kept on $N \to v$; nothing on
///   $N \to c$, because the composition cancels. This is where the move removes one mutation
///   per reverted position.
///
/// The resulting edges carry:
///
/// | edge | substitutions |
/// | --- | --- |
/// | $u \to N$ | $T$ |
/// | $N \to v$ | $H \cup R$ (the original $M_v$ entries at child-shared positions) |
/// | $N \to c$ | $\mathrm{compose}(M_v, M_c)$ at $M_c$ positions ($H' \cup D$) |
///
/// where $D$ are the child's own positions absent from $M_v$. The net substitution change is
/// $\Delta = -\lvert R\rvert$: no mutation is ever added.
///
/// Branch lengths are split in proportion to substitution count so that root-to-$v$ and
/// root-to-$c$ distances are preserved exactly; the next optimizer iteration re-fits them.
///
/// Indels use an all-or-nothing rule (see [`split_indels`]) that is always distance
/// preserving. The substitution gain is unaffected by the indel handling.
///
/// The relocated edges $e_p$ and $e_c$ keep their edge keys via
/// [`Graph::reparent_edge`](treetime_graph::graph::Graph::reparent_edge), so their partition
/// entries stay valid and only their substitution and indel content is rewritten. Only the
/// fresh $u \to N$ edge and the node $N$ are new keys to register.
///
/// Returns the key of the new node $N$.
pub(crate) fn hoist_reverting_child(
  graph: &mut GraphAncestral,
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense: &[Arc<RwLock<PartitionMarginalDense>>],
  parent_edge_key: GraphEdgeKey,
  child_edge_key: GraphEdgeKey,
) -> Result<GraphNodeKey, Report> {
  let u_key = graph.get_source_node_key(parent_edge_key)?;

  // Compute the per-partition split before touching the graph, so the reads see the
  // pre-move edge state.
  let mut splits = Vec::with_capacity(sparse.len());
  let mut total_parent_subs = 0_usize;
  let mut total_hoisted_subs = 0_usize;
  for partition in sparse {
    let partition = partition.read_arc();
    let parent_subs = partition
      .edges
      .get(&parent_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let child_subs = partition
      .edges
      .get(&child_edge_key)
      .map_or(Vec::new(), |e| e.fitch_subs().to_vec());
    let parent_indels = partition
      .edges
      .get(&parent_edge_key)
      .map_or(Vec::new(), |e| e.indels.clone());
    let child_indels = partition
      .edges
      .get(&child_edge_key)
      .map_or(Vec::new(), |e| e.indels.clone());

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

  // Distance-preserving branch-length split, proportional to substitution count. The move
  // only fires when R is non-empty, so `total_parent_subs >= 1` and the ratio is well defined.
  let bl_uv = edge_branch_length(graph, parent_edge_key);
  let bl_vc = edge_branch_length(graph, child_edge_key);
  let bl_un = if total_parent_subs > 0 {
    bl_uv * (total_hoisted_subs as f64) / (total_parent_subs as f64)
  } else {
    0.0
  };
  let bl_nv = bl_uv - bl_un;
  let bl_nc = bl_nv + bl_vc;

  // Graph surgery: add N, connect u -> N, then relocate the two existing edges under N.
  let n_key = graph.add_node(NodeAncestral::default());
  let un_edge_key = graph.add_edge(
    u_key,
    n_key,
    EdgeAncestral {
      branch_length: Some(bl_un),
    },
  )?;
  graph.reparent_edge(parent_edge_key, n_key)?; // e_p becomes N -> v
  graph.reparent_edge(child_edge_key, n_key)?; // e_c becomes N -> c
  set_edge_branch_length(graph, parent_edge_key, bl_nv);
  set_edge_branch_length(graph, child_edge_key, bl_nc);

  // Sparse partition bookkeeping. The relocated edges keep their keys, so only their content
  // is rewritten; the u -> N edge and node N are inserted fresh.
  for (partition, split) in sparse.iter().zip(splits) {
    let mut partition = partition.write_arc();

    let mut node_n = SparseNodePartition::empty(&partition.alphabet);
    node_n.seq.composition = partition.nodes[&u_key].seq.composition.clone();
    partition.nodes.entry(n_key).or_insert(node_n);

    let un_edge = partition.edges.entry(un_edge_key).or_default();
    un_edge.set_fitch_subs(split.hoisted);
    un_edge.indels = split.indels.hoisted;

    let nv_edge = partition.edges.entry(parent_edge_key).or_default();
    nv_edge.set_fitch_subs(split.kept);
    nv_edge.indels = split.indels.kept;

    let nc_edge = partition.edges.entry(child_edge_key).or_default();
    nc_edge.set_fitch_subs(split.composed);
    nc_edge.indels = split.indels.composed;
  }

  // Dense partitions carry no mutation lists, but they key node and edge state by graph key
  // and must learn about the new node and edge. `apply_reroot` registers reroot-created keys
  // the same way; mirroring it keeps dense state consistent when both families coexist.
  for partition in dense {
    let mut partition = partition.write_arc();
    partition.data.nodes.entry(n_key).or_insert_with(|| DenseNodePartition {
      seq: DenseSeqInfo::default(),
      profile: DenseSeqDistribution::default(),
    });
    partition.data.edges.entry(un_edge_key).or_default();
  }

  Ok(n_key)
}

/// One partition's substitution split of $M_v$ against $M_c$ (see [`hoist_reverting_child`]).
struct SubSplit {
  /// $T$: parent positions untouched by the child. Goes to $u \to N$.
  hoisted: Vec<Sub>,
  /// $H \cup R$: original parent entries at positions the child also touches. Goes to $N \to v$.
  kept: Vec<Sub>,
  /// $H' \cup D$: $\mathrm{compose}(M_v, M_c)$ restricted to child positions. Goes to $N \to c$.
  composed: Vec<Sub>,
}

/// Split parent-edge substitutions against one child edge into the hoist's three edges.
///
/// Both inputs are position-sorted with at most one entry per position (the fitch-subs
/// invariant). The single merge-walk keeps every output position-sorted.
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
        hoisted.push(ps.clone()); // T: parent-only position
        pi += 1;
      },
      Ordering::Greater => {
        composed.push(cs.clone()); // D: child-only position
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
        kept.push(ps.clone()); // H or R: keep original parent entry on N -> v
        if ps.reff() != cs.qry() {
          composed.push(Sub::new(ps.reff(), ps.pos(), cs.qry())?); // H': net a -> d on N -> c
        }
        // ps.reff() == cs.qry(): reversion, cancels, nothing on N -> c
        pi += 1;
        ci += 1;
      },
    }
  }

  hoisted.extend_from_slice(&parent_subs[pi..]); // remaining parent-only -> T
  composed.extend_from_slice(&child_subs[ci..]); // remaining child-only -> D

  Ok(SubSplit {
    hoisted,
    kept,
    composed,
  })
}

/// Count reversions between one parent edge and one child edge (single partition).
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

/// One partition's indel split for the hoist.
struct IndelSplit {
  hoisted: Vec<InDel>,
  kept: Vec<InDel>,
  composed: Vec<InDel>,
}

/// Split parent-edge indels against a child edge, all-or-nothing.
///
/// `compose_indels` merges overlapping and adjacent ranges rather than being position-keyed,
/// so an exact three-way split is not always definable. When no child indel overlaps or is
/// adjacent to any parent indel, the parent indels move cleanly above $N$ and $N \to c$ carries
/// only the child's own indels. Otherwise the parent indels stay on $N \to v$ and $N \to c$
/// carries the full composition. Both branches preserve the root-to-$v$ and root-to-$c$ indel
/// content exactly.
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

/// Whether any parent indel overlaps or is adjacent to any child indel.
///
/// Adjacency matters because `compose_indels` merges touching deletions, so a clean hoist is
/// only possible when the ranges are strictly separated.
fn indels_interact(parent_indels: &[InDel], child_indels: &[InDel]) -> bool {
  parent_indels.iter().any(|p| {
    child_indels
      .iter()
      .any(|c| ranges_overlap_or_adjacent(p.range, c.range))
  })
}

/// Whether two half-open ranges overlap or touch at an endpoint.
fn ranges_overlap_or_adjacent((a_lo, a_hi): (usize, usize), (b_lo, b_hi): (usize, usize)) -> bool {
  a_lo <= b_hi && b_lo <= a_hi
}

/// Combined substitution and indel split for one partition.
struct EdgeSplit {
  hoisted: Vec<Sub>,
  kept: Vec<Sub>,
  composed: Vec<Sub>,
  indels: IndelSplit,
}

fn edge_branch_length(graph: &GraphAncestral, edge_key: GraphEdgeKey) -> f64 {
  graph
    .get_edge(edge_key)
    .and_then(|edge| edge.read_arc().payload().read_arc().branch_length())
    .unwrap_or(0.0)
}

fn set_edge_branch_length(graph: &GraphAncestral, edge_key: GraphEdgeKey, branch_length: f64) {
  if let Some(edge) = graph.get_edge(edge_key) {
    edge
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(branch_length));
  }
}
