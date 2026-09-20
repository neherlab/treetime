use crate::alphabet::alphabet::Alphabet;
use crate::seq::alignment::get_common_length;
use crate::{make_error, make_report};
use eyre::Report;
use log::warn;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlignmentRecord, AlphabetLike, Seq, seq};

#[allow(clippy::as_conversions, reason = "count/index numeric cast is exact for the domain range")]
/// Complete an alignment so every tree leaf has a sequence, matching v0 missing-tip semantics.
///
/// Tips absent from the alignment are treated as fully ambiguous: each gets an all-`unknown`
/// sequence of the alignment length, which is equivalent to v0's uniform leaf profile. The same
/// treatment applies to nucleotide and amino-acid partitions, so the missing-data policy lives in
/// one shared place rather than in each reconstruction backend.
///
/// If more than one third of the tips are missing, the function errors, unless `ignore_missing_alns`
/// is set. Extra alignment records that match no leaf are left in place and ignored by attachment
/// (which matches by name).
///
/// Reference: TreeTime v0 `TreeAnc._check_alignment_tree_gtr_consistency`
/// (`packages/legacy/treetime/treetime/treeanc.py:419-440`): per-leaf warning, abort when
/// `failed_leaves > tree.count_terminals()/3` unless `ignore_missing_alns`, missing leaves later
/// assigned a uniform (all-states-equal) profile at reconstruction time.
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

  // v0 uses strict `>` with true (float) division, so exactly one third missing does not abort.
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

/// Map characters not present in `alphabet` to its unknown state, returning the sanitized sequence
/// and the number of characters changed.
///
/// Amino-acid translation files are read with the stop-inclusive `Aa` alphabet so that a stop codon
/// `*` is accepted. When reconstruction runs over an empirical 20-amino-acid model (no stop), the
/// stop codon and any other out-of-alphabet character must be folded into the unknown state `X`
/// before attachment, matching augur's amino-acid sequence correction (`_make_seq_corrector('aa')`).
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
