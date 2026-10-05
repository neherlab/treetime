use crate::annotated_graph::AnnotatedTreeView;
use crate::mutation_filter::UnknownBridge;
use eyre::{Report, WrapErr, eyre};
use std::collections::BTreeSet;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_write_str};
use treetime_io::usher_mat::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};
use treetime_primitives::AsciiChar;
use treetime_utils::{make_error, make_internal_error};

#[derive(Debug)]
#[expect(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct MatOutput {
  pub(crate) tree: UsherTree,
  pub(crate) gaps: MatGapCounts,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[expect(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct MatGapCounts {
  pub(crate) deletions: usize,
  pub(crate) insertions: usize,
  pub(crate) substitutions: usize,
}

impl MatGapCounts {
  pub(crate) fn warning(self) -> Option<String> {
    let mut clauses = vec![];
    if self.deletions > 0 {
      clauses.push(format!("wrote {} deletion(s) as missing data (N)", self.deletions));
    }
    if self.insertions > 0 || self.substitutions > 0 {
      clauses.push(format!(
        "left out {} insertion(s) and {} substitution(s) in alignment columns where the root sequence, which is the MAT reference, has a gap",
        self.insertions, self.substitutions
      ));
    }
    (!clauses.is_empty()).then(|| format!("UShER MAT has no gap state: {}", clauses.join("; ")))
  }

  fn count(&mut self, mutations: &[Mutation], reference_gaps: &BTreeSet<usize>) {
    for mutation in mutations
      .iter()
      .filter(|mutation| mutation.track == MutationTrack::Nucleotide)
    {
      match &mutation.event {
        MutationEvent::Deletion(_) => self.deletions += 1,
        MutationEvent::Insertion(segment) => {
          if reference_gaps.range(segment.range.0..segment.range.1).next().is_some() {
            self.insertions += 1;
          }
        },
        MutationEvent::Substitution(substitution) => {
          if reference_gaps.contains(&substitution.pos()) {
            self.substitutions += 1;
          }
        },
      }
    }
  }
}

pub(crate) fn mat_tree(tree: &AnnotatedTreeView<'_>) -> Result<MatOutput, Report> {
  let graph = tree.graph();
  let sequences = graph.sequences.as_ref();
  let reference = sequences.map(|sequences| sequences.root_sequence.as_str());
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let reference_gaps: BTreeSet<usize> = reference.map_or_else(BTreeSet::new, |reference| {
    reference
      .bytes()
      .enumerate()
      .filter(|&(_, state)| state == u8::from(alphabet.gap()))
      .map(|(pos, _)| pos)
      .collect()
  });
  let mut missing_data = UnknownBridge::new(alphabet.unknown());
  let mut gaps = MatGapCounts::default();
  let order = tree.tree().preorder();
  let mut node_mutations = Vec::with_capacity(order.len());
  let mut condensed_nodes = Vec::with_capacity(order.len());
  let mut metadata = Vec::with_capacity(order.len());
  for &key in order {
    let name = graph.names[&key].clone().unwrap_or_default();
    let mutations = match (sequences, tree.tree().parent(key)) {
      (Some(sequences), Some((parent_key, edge_key))) => {
        let mutations = &sequences.edge_mutations[&edge_key];
        gaps.count(mutations, &reference_gaps);
        let mutations = gaps_as_missing_data(mutations, alphabet.unknown())?;
        missing_data.bridge_edge(parent_key, key, tree.tree().children(key).len(), mutations)?
      },
      _ => vec![],
    };
    let mut mutations = mutations
      .iter()
      .filter(|mutation| !is_in_reference_gap(mutation, &reference_gaps))
      .map(|mutation| mat_mutation(mutation, reference, &alphabet, &name))
      .collect::<Result<Vec<_>, _>>()?;
    mutations.sort_by_key(|mutation| mutation.position);
    node_mutations.push(UsherMutationList { mutation: mutations });
    condensed_nodes.push(UsherTreeNode {
      node_name: name,
      condensed_leaves: vec![],
    });
    metadata.push(UsherMetadata {
      clade_annotations: vec![],
    });
  }
  Ok(MatOutput {
    tree: UsherTree {
      newick: nwk_write_str(
        tree.tree(),
        graph.names,
        graph.tree_branch_lengths(),
        &NwkWriteOptions::default(),
        &NwkNodeComments::new(),
      )?,
      node_mutations,
      condensed_nodes,
      metadata,
    },
    gaps,
  })
}

pub(crate) fn mat_mutation(
  mutation: &Mutation,
  reference: Option<&str>,
  alphabet: &Alphabet,
  node_name: &str,
) -> Result<UsherMutation, Report> {
  if mutation.track != MutationTrack::Nucleotide {
    return make_internal_error!(
      "Node '{node_name}' has an amino-acid mutation, but UShER MAT stores nucleotide mutations only"
    );
  }
  let MutationEvent::Substitution(substitution) = &mutation.event else {
    return make_internal_error!(
      "Node '{node_name}' has an insertion or deletion that was not converted to missing data"
    );
  };
  let reference = reference.ok_or_else(|| {
    eyre!("Node '{node_name}' has nucleotide mutations, but UShER MAT requires a root nucleotide reference")
  })?;
  let position = substitution
    .pos()
    .checked_add(1)
    .ok_or_else(|| eyre!("Node '{node_name}' mutation coordinate overflow"))?;
  let position = i32::try_from(position).wrap_err_with(|| {
    format!("Node '{node_name}' mutation position {position} exceeds the UShER MAT i32 coordinate range")
  })?;
  let reference_state = reference.as_bytes().get(substitution.pos()).copied().ok_or_else(|| {
    eyre!(
      "Node '{node_name}' mutation position {} is outside the root nucleotide reference of length {}",
      substitution.pos() + 1,
      reference.len()
    )
  })?;
  let reference_state = AsciiChar::try_new(reference_state)?;
  Ok(UsherMutation {
    position,
    ref_nuc: mat_nucleotide(reference_state, node_name, "root reference")?,
    par_nuc: mat_nucleotide(substitution.reff(), node_name, "parent")?,
    mut_nuc: mat_nucleotide_states(substitution.qry(), alphabet, node_name)?,
    chromosome: String::new(),
  })
}

fn gaps_as_missing_data(mutations: &[Mutation], unknown: AsciiChar) -> Result<Vec<Mutation>, Report> {
  let mut converted = Vec::with_capacity(mutations.len());
  for mutation in mutations {
    let (segment, is_deletion) = match (&mutation.track, &mutation.event) {
      (MutationTrack::Nucleotide, MutationEvent::Deletion(segment)) => (segment, true),
      (MutationTrack::Nucleotide, MutationEvent::Insertion(segment)) => (segment, false),
      _ => {
        converted.push(mutation.clone());
        continue;
      },
    };
    for (pos, &state) in (segment.range.0..segment.range.1).zip(segment.sequence.iter()) {
      if state == unknown {
        continue;
      }
      let substitution = if is_deletion {
        Sub::new(state, pos, unknown)?
      } else {
        Sub::new(unknown, pos, state)?
      };
      converted.push(Mutation::substitution(MutationTrack::Nucleotide, substitution));
    }
  }
  Ok(converted)
}

fn is_in_reference_gap(mutation: &Mutation, reference_gaps: &BTreeSet<usize>) -> bool {
  matches!(
    (&mutation.track, &mutation.event),
    (MutationTrack::Nucleotide, MutationEvent::Substitution(substitution)) if reference_gaps.contains(&substitution.pos())
  )
}

fn mat_nucleotide_states(nucleotide: AsciiChar, alphabet: &Alphabet, node_name: &str) -> Result<Vec<i32>, Report> {
  let rejected = || {
    format!(
      "Node '{node_name}' has child nucleotide '{nucleotide}', but UShER MAT accepts only A, C, G, T, IUPAC ambiguity codes, or N"
    )
  };
  let states = alphabet.canonical_states(nucleotide).wrap_err_with(rejected)?;
  if states.is_empty() {
    return make_internal_error!(
      "Node '{node_name}' has child nucleotide '{nucleotide}', which stands for no nucleotide"
    );
  }
  states
    .iter()
    .map(|state| mat_nucleotide(state, node_name, "child"))
    .collect()
}

fn mat_nucleotide(nucleotide: AsciiChar, node_name: &str, role: &str) -> Result<i32, Report> {
  match char::from(nucleotide).to_ascii_uppercase() {
    'A' => Ok(0),
    'C' => Ok(1),
    'G' => Ok(2),
    'T' => Ok(3),
    state => {
      make_error!("Node '{node_name}' has {role} nucleotide '{state}', but UShER MAT accepts only A, C, G, or T")
    },
  }
}
