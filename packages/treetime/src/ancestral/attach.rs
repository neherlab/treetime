use crate::alphabet::alphabet::Alphabet;
use crate::progress::LogSink;
use crate::progress_warn;
use crate::seq::alignment::NodeSeqInput;
use crate::{make_error, make_internal_report};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlphabetLike, Seq, seq};
use treetime_utils::cast::ToF64Exact;

pub fn complete_alignment_for_leaves(
  graph: &Graph,
  nodes: &mut BTreeMap<GraphNodeKey, NodeSeqInput>,
  alignment_length: usize,
  alphabet: &Alphabet,
  limit: MissingLeafLimit,
  source: &str,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let n_leaves = graph.num_leaves();
  let missing = graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .filter(|key| nodes[key].seq.is_none())
    .collect_vec();

  let n_missing = missing.len();
  if n_missing > 0 {
    for key in &missing {
      let name = node_name_or_key(*key, nodes[key].name.as_deref());
      progress_warn!(
        log,
        "No sequence found for leaf '{name}' in {source}; treating it as fully ambiguous (missing data)."
      );
    }
    progress_warn!(
      log,
      "{n_missing} of {n_leaves} tips have no matching sequence in {source} and are treated as fully ambiguous."
    );
  }

  if limit.applies() && n_missing.to_f64_exact() > n_leaves.to_f64_exact() / 3.0 {
    let message = format!(
      "More than one third of terminal nodes ({n_missing} of {n_leaves}) cannot be assigned a sequence. \
       Are you sure the alignment belongs to the tree?"
    );
    return match limit {
      MissingLeafLimit::Enforced => make_error!("{message}"),
      MissingLeafLimit::Overridable { .. } => {
        make_error!("{message} Pass --ignore-missing-alns to proceed, treating missing tips as fully ambiguous.")
      },
    };
  }

  let unknown = alphabet.unknown();
  for key in missing {
    let node = nodes
      .get_mut(&key)
      .ok_or_else(|| make_internal_report!("Leaf {key} has no sequence input"))?;
    node.seq = Some(seq![unknown; alignment_length]);
  }

  Ok(())
}

pub fn sanitize_to_alphabet(seq: &Seq, alphabet: &Alphabet) -> (Seq, usize) {
  let unknown = alphabet.unknown();
  let mut changed = 0_usize;
  let sanitized = seq
    .iter()
    .map(|&c| {
      if alphabet.contains(c) {
        c
      } else {
        changed += 1;
        unknown
      }
    })
    .collect();
  (sanitized, changed)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingLeafLimit {
  Enforced,
  Overridable { ignored: bool },
}

impl MissingLeafLimit {
  fn applies(self) -> bool {
    !matches!(self, Self::Overridable { ignored: true })
  }
}
