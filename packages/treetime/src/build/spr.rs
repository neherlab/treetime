use crate::build::insert::{Insertee, insert};
use crate::build::mat::Mat;
use crate::build::placement::{Target, place_exhaustive, place_greedy};
use crate::build::query::Query;
use eyre::Report;
use std::borrow::Cow;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SprStats {
  pub tried: usize,
  pub moves: usize,
  pub gain: usize,
  pub visited: usize,
  pub parsimony_before: usize,
  pub parsimony_after: usize,
}

pub(crate) fn spr_round(mat: &mut Mat, queries: &[Query], slack: Option<usize>) -> Result<SprStats, Report> {
  let mut stats = SprStats::default();
  for node in mat.preorder().into_iter().rev() {
    if !mat.has_node(node) {
      continue;
    }
    let Some((parent, _)) = mat.parent(node) else {
      continue;
    };
    if parent == mat.root() && mat.n_children(parent) == 1 {
      continue;
    }
    let query = match mat.sample_of(node) {
      Some(sample) => Cow::Borrowed(&queries[sample]),
      None => Cow::Owned(Query::from_diffs(mat.node_state(node))),
    };
    let before = mat.parsimony();
    let detached = mat.detach(node)?;
    let prefer = match detached.reattach_hint {
      [Some(_), Some(child)] => mat.parent(child).map(|(_, merged)| Target::Edge(merged)),
      [Some(parent), None] => Some(Target::Node(parent)),
      _ => None,
    };
    let placement = match slack {
      Some(slack) => {
        let starts = detached.reattach_hint.iter().flatten().copied().collect::<Vec<_>>();
        place_greedy(mat, &query, &starts, slack, prefer)?
      },
      None => place_exhaustive(mat, &query, prefer, None, true)?,
    };
    let inserted = insert(mat, &query, placement.target, Insertee::Subtree(node))?;
    if inserted.cost != placement.cost {
      return Err(crate::make_internal_report!(
        "SPR of node {node}: the placement scored {} but its insertion below {} costs {}",
        placement.cost,
        inserted.attached_to,
        inserted.cost
      ));
    }
    stats.tried += 1;
    stats.visited += placement.visited;
    let after = mat.parsimony();
    if after < before {
      stats.moves += 1;
      stats.gain += before - after;
    } else if after > before {
      return Err(crate::make_internal_report!(
        "SPR of node {node} raised the parsimony score from {before} to {after}"
      ));
    }
  }
  Ok(stats)
}
