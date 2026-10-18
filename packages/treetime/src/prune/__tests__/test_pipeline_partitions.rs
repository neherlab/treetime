#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{FILE_1, FILE_2, internal_node_names, run_prune_empty};
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;

  #[rustfmt::skip]
  #[rstest]
  #[case::first_only( &[FILE_1],         btreeset! {"root", "AB"})]
  #[case::second_only(&[FILE_2],         btreeset! {"root", "ABC"})]
  #[case::both(       &[FILE_1, FILE_2], btreeset! {"root", "ABC", "AB"})]
  #[trace]
  fn test_prune_empty_collapses_a_branch_only_when_every_alignment_has_no_substitution_on_it(
    #[case] files: &[&[(&str, &str)]],
    #[case] expected: BTreeSet<&str>,
  ) -> Result<(), Report> {
    let output = run_prune_empty(files)?;
    let numbers = output.results.iter().map(|result| result.number).collect::<Vec<_>>();
    assert_eq!(
      (expected, (1..=files.len()).collect::<Vec<_>>()),
      (internal_node_names(&output), numbers)
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::progress::NoopProgress;
    use crate::prune::pipeline::{PruneInput, PruneOutput, PruneParams, run};
    use crate::test_utils::{leaf_seq_inputs, sequence_source};
    use eyre::Report;
    use maplit::btreeset;
    use std::collections::BTreeSet;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};

    const TREE: &str = "((((A:0.1,B:0.1)AB:0.1,C:0.1)ABC:0.1,D:0.1)ABCD:0.1,E:0.1)root;";

    pub(super) const FILE_1: &[(&str, &str)] = &[("A", "T"), ("B", "T"), ("C", "A"), ("D", "A"), ("E", "A")];

    pub(super) const FILE_2: &[(&str, &str)] = &[("A", "G"), ("B", "G"), ("C", "G"), ("D", "C"), ("E", "C")];

    pub(super) fn run_prune_empty(files: &[&[(&str, &str)]]) -> Result<PruneOutput, Report> {
      let parsed = nwk_read(TREE.as_bytes())?;
      let names = parsed.names();
      let sources = files
        .iter()
        .map(|file| {
          let records = file
            .iter()
            .map(|(name, seq)| {
              Ok(AlignmentRecord {
                name: (*name).to_owned(),
                seq: Seq::try_from_str(seq)?,
              })
            })
            .collect::<Result<Vec<_>, Report>>()?;
          sequence_source(Alphabet::default(), leaf_seq_inputs(&parsed.graph, &names, records))
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let params = PruneParams {
        prune_short: None,
        prune_empty: true,
        merge_shared_mutations: false,
        node_keys: btreeset! {},
      };
      let input = PruneInput {
        graph: parsed.graph,
        sources,
        branch_lengths: parsed.branch_lengths,
      };
      run(&params, input, &names, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())
    }

    pub(super) fn internal_node_names(output: &PruneOutput) -> BTreeSet<&str> {
      output
        .graph
        .get_internal_nodes()
        .filter_map(|node| output.names[&node.key()].as_deref())
        .collect()
    }
  }
}
