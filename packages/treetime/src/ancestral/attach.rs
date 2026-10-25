use crate::alphabet::alphabet::Alphabet;
use crate::progress::LogSink;
use crate::progress_warn;
use crate::seq::alignment::NodeSeqInput;
use crate::{make_error, make_internal_report};
use eyre::Report;
use itertools::Itertools;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_primitives::{AlphabetLike, Seq, seq};
use treetime_utils::cast::ToF64Exact;

pub fn complete_alignment_for_leaves(
  graph: &Graph,
  nodes: &mut NodeMap<NodeSeqInput>,
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
  warn_leaves_without_data(
    &missing.iter().map(|key| nodes[key].name.clone()).collect_vec(),
    n_leaves,
    LeafData::Sequence,
    source,
    log,
  );

  if !limit.ignored && n_missing.to_f64_exact() > n_leaves.to_f64_exact() / 3.0 {
    return make_error!(
      "More than one third of terminal nodes ({n_missing} of {n_leaves}) cannot be assigned a sequence. \
       Are you sure the alignment belongs to the tree? \
       Pass --allow-missing-sequences to proceed, treating missing tips as fully ambiguous."
    );
  }

  let unknown = alphabet.unknown();
  for key in missing {
    let node = nodes
      .get_mut(key)
      .ok_or_else(|| make_internal_report!("Leaf {key} has no sequence input"))?;
    node.seq = Some(seq![unknown; alignment_length]);
  }

  Ok(())
}

pub fn warn_leaves_without_data(
  missing_names: &[String],
  n_leaves: usize,
  data: LeafData,
  source: &str,
  log: &dyn LogSink,
) {
  let noun = data.noun();
  for name in missing_names {
    progress_warn!(
      log,
      "No {noun} found for leaf '{name}' in {source}; treating it as fully ambiguous (missing data)."
    );
  }
  let n_missing = missing_names.len();
  if n_missing > 0 {
    progress_warn!(
      log,
      "{n_missing} of {n_leaves} tips have no matching {noun} in {source} and are treated as fully ambiguous."
    );
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LeafData {
  Sequence,
  TraitValue,
}

impl LeafData {
  const fn noun(self) -> &'static str {
    match self {
      Self::Sequence => "sequence",
      Self::TraitValue => "trait value",
    }
  }
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
pub struct MissingLeafLimit {
  pub ignored: bool,
}
