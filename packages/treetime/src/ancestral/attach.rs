use crate::alphabet::alphabet::Alphabet;
use crate::seq::alignment::get_common_length;
use crate::{make_error, make_report};
use eyre::Report;
use log::warn;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlignmentRecord, AlphabetLike, Seq, seq};

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn complete_alignment_for_leaves(
  graph: &Graph,
  mut sequences: Vec<AlignmentRecord>,
  alphabet: &Alphabet,
  ignore_missing_alns: bool,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<Vec<AlignmentRecord>, Report> {
  let alignment_length = get_common_length(&sequences)?;

  let present: BTreeSet<String> = sequences.iter().map(|record| record.name.clone()).collect();

  let mut missing = Vec::new();
  let mut n_leaves = 0_usize;
  for leaf in graph.get_leaves() {
    n_leaves += 1;
    let name = names[&leaf.key()]
      .clone()
      .ok_or_else(|| {
        make_report!("Expected all leaf nodes to have names, so they can be matched to their sequences. Found a leaf node with no name.")
      })?;
    if !present.contains(&name) {
      missing.push(name);
    }
  }
  drop(present);

  let n_missing = missing.len();
  if n_missing > 0 {
    for name in &missing {
      warn!("No sequence found for leaf '{name}'; treating it as fully ambiguous (missing data).");
    }
    warn!(
      "{n_missing} of {n_leaves} tips have no matching sequence in the alignment and are treated as fully ambiguous."
    );
  }

  if !ignore_missing_alns && (n_missing as f64) > (n_leaves as f64) / 3.0 {
    return make_error!(
      "At least one third of terminal nodes ({n_missing} of {n_leaves}) cannot be assigned a sequence. \
       Are you sure the alignment belongs to the tree? \
       Pass --ignore-missing-alns to proceed, treating missing tips as fully ambiguous."
    );
  }

  for name in missing {
    sequences.push(AlignmentRecord {
      name,
      seq: seq![alphabet.unknown(); alignment_length],
    });
  }

  Ok(sequences)
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
