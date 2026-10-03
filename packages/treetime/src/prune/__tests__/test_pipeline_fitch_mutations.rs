#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::cancel::NoopCancel;
  use crate::progress::NoopProgress;
  use crate::prune::pipeline::{PruneInput, PruneOutput, PruneParams, run};
  use crate::seq::mutation::{Mutation, MutationTrack, Sub};
  use crate::test_utils::find_edge_key;
  use eyre::Report;
  use helpers::*;
  use itertools::Itertools;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AlignmentRecord, AsciiChar, Seq};

  #[test]
  fn test_prune_pipeline_edge_mutations_are_hand_derived_fitch_substitutions() -> Result<(), Report> {
    let output = run_prune_empty()?;
    let actual = edge_mutations_by_name(&output)?;
    let expected = btreemap! {
      ("root", "AB") => vec![sub(b'G', 3, b'T')?],
      ("root", "C")  => vec![sub(b'A', 0, b'T')?],
      ("root", "D")  => vec![],
      ("AB", "A")    => vec![],
      ("AB", "B")    => vec![sub(b'C', 1, b'G')?],
    };
    assert_eq!((5, expected), (output.graph.get_edges().count(), actual));
    Ok(())
  }

  #[test]
  fn test_prune_pipeline_root_sequence_is_hand_derived_fitch_root() -> Result<(), Report> {
    let output = run_prune_empty()?;
    assert_eq!(Seq::try_from_str("ACGG")?, output.partitions[0].fitch_root_sequence());
    Ok(())
  }

  #[test]
  fn test_prune_pipeline_mutation_count_equals_fitch_parsimony_score() -> Result<(), Report> {
    let output = run_prune_empty()?;
    let total: usize = edge_mutations_by_name(&output)?.values().map(Vec::len).sum();
    assert_eq!(3, total);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run_prune_empty() -> Result<PruneOutput, Report> {
      let parsed = nwk_read_str("(((A:0.1,B:0.1)AB:0.1,C:0.1)ABC:0.1,D:0.1)root;")?;
      let names = parsed.names();
      let sequences = vec![
        record("A", "ACGT")?,
        record("B", "AGGT")?,
        record("C", "TCGG")?,
        record("D", "ACGG")?,
      ];
      let params = PruneParams {
        prune_short: None,
        prune_empty: true,
        merge_shared_mutations: false,
        node_names: btreeset! {},
      };
      let input = PruneInput {
        graph: parsed.graph,
        alphabet: Alphabet::new(AlphabetName::Nuc)?,
        sequences: Some(sequences),
        branch_lengths: parsed.branch_lengths,
      };
      run(&params, input, &names, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())
    }

    pub(super) fn edge_mutations_by_name(
      output: &PruneOutput,
    ) -> Result<BTreeMap<(&'static str, &'static str), Vec<Mutation>>, Report> {
      let edges = [("root", "AB"), ("root", "C"), ("root", "D"), ("AB", "A"), ("AB", "B")];
      edges
        .into_iter()
        .map(|(source, target)| {
          let edge_key = find_edge_key(&output.graph, &output.names, source, target)
            .ok_or_else(|| eyre::eyre!("edge {source}->{target} missing after pruning"))?;
          let mutations = output.partitions[0].edge_fitch_mutations(edge_key, &MutationTrack::Nucleotide)?;
          Ok(((source, target), mutations))
        })
        .try_collect()
    }

    pub(super) fn sub(reff: u8, pos: usize, qry: u8) -> Result<Mutation, Report> {
      Ok(Mutation::substitution(
        MutationTrack::Nucleotide,
        Sub::new(AsciiChar::try_new(reff)?, pos, AsciiChar::try_new(qry)?)?,
      ))
    }

    fn record(name: &str, sequence: &str) -> Result<AlignmentRecord, Report> {
      Ok(AlignmentRecord {
        name: name.to_owned(),
        seq: Seq::try_from_str(sequence)?,
      })
    }
  }
}
