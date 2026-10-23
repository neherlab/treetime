#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::test_utils::{TestDataset, split_at_gap_free_column};
  use eyre::Report;
  use helpers::{
    FILE_1, FILE_2, TIED_FILE_1, TIED_FILE_2, clades, internal_node_names, records, run_prune, run_prune_empty,
  };
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

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
    let numbers = output.alignments.iter().map(|result| result.number).collect::<Vec<_>>();
    assert_eq!(
      (expected, (1..=files.len()).collect::<Vec<_>>()),
      (internal_node_names(&output), numbers)
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zika( TestDataset::Zika20)]
  #[case::flu(  TestDataset::Flu20)]
  #[case::ebola(TestDataset::Ebola20)]
  #[trace]
  fn test_prune_two_column_blocks_give_the_tree_of_one_alignment(#[case] dataset: TestDataset) -> Result<(), Report> {
    let whole = dataset.alignment()?;
    let [first, second] = split_at_gap_free_column(&whole)?;
    let whole = run_prune(dataset.tree(), &[whole])?;
    let split = run_prune(dataset.tree(), &[first, second])?;
    assert_eq!(clades(&whole)?, clades(&split)?);
    pretty_assert_map_abs_diff_eq!(
      branch_lengths_or_zero(&whole.branch_lengths),
      branch_lengths_or_zero(&split.branch_lengths),
      epsilon = 1e-8
    );
    Ok(())
  }

  #[test]
  fn test_prune_merges_the_same_of_two_tied_sibling_groups_in_either_alignment_order() -> Result<(), Report> {
    let tree = "(A:0.1,B:0.1,C:0.1,D:0.1)root;";
    let forward = run_prune(tree, &[records(TIED_FILE_1)?, records(TIED_FILE_2)?])?;
    let backward = run_prune(tree, &[records(TIED_FILE_2)?, records(TIED_FILE_1)?])?;
    assert_eq!(
      (
        btreeset! {btreeset! {"A", "B"}, btreeset! {"A", "B", "C", "D"}},
        clades(&forward)?
      ),
      (clades(&forward)?, clades(&backward)?)
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::partition::resolve::{PartitionLayout, resolve_partition_layouts};
    use crate::progress::NoopProgress;
    use crate::prune::pipeline::{PruneInput, PruneOutput, PruneParams, run};
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{leaf_seq_inputs, sequence_source};
    use eyre::Report;
    use maplit::btreeset;
    use std::collections::BTreeSet;
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};

    const TREE: &str = "((((A:0.1,B:0.1)AB:0.1,C:0.1)ABC:0.1,D:0.1)ABCD:0.1,E:0.1)root;";

    pub(super) const FILE_1: &[(&str, &str)] = &[("A", "T"), ("B", "T"), ("C", "A"), ("D", "A"), ("E", "A")];

    pub(super) const FILE_2: &[(&str, &str)] = &[("A", "G"), ("B", "G"), ("C", "G"), ("D", "C"), ("E", "C")];

    pub(super) const TIED_FILE_1: &[(&str, &str)] = &[("A", "T"), ("B", "T"), ("C", "A"), ("D", "A")];

    pub(super) const TIED_FILE_2: &[(&str, &str)] = &[("A", "G"), ("B", "C"), ("C", "G"), ("D", "C")];

    pub(super) fn records(file: &[(&str, &str)]) -> Result<Vec<AlignmentRecord>, Report> {
      file
        .iter()
        .map(|(name, seq)| {
          Ok(AlignmentRecord {
            name: (*name).to_owned(),
            seq: Seq::try_from_str(seq)?,
          })
        })
        .collect()
    }

    pub(super) fn run_prune(tree: &str, files: &[Vec<AlignmentRecord>]) -> Result<PruneOutput, Report> {
      let parsed = nwk_read(tree.as_bytes())?;
      let names = parsed.names.clone();
      let sources = files
        .iter()
        .map(|records| {
          sequence_source(
            Alphabet::default(),
            leaf_seq_inputs(&parsed.graph, &names, records.clone()),
          )
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let params = PruneParams {
        prune_short: None,
        prune_empty: true,
        merge_shared_mutations: true,
        node_keys: btreeset! {},
      };
      let input = PruneInput {
        graph: parsed.graph,
        names,
        layouts: whole_layouts(&sources)?,
        alignments: sources,
        branch_lengths: parsed.branch_lengths,
      };
      run(&params, input, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())
    }

    pub(super) fn clades(output: &PruneOutput) -> Result<BTreeSet<BTreeSet<&str>>, Report> {
      let mut ids = TipIds::default();
      let keys = graph_clade_keys(&output.graph, &output.names, &mut ids)?;
      Ok(
        output
          .graph
          .get_internal_nodes()
          .map(|node| keys[&node.key()].iter().map(|&id| ids.name(id)).collect())
          .collect(),
      )
    }

    pub(super) fn run_prune_empty(files: &[&[(&str, &str)]]) -> Result<PruneOutput, Report> {
      let parsed = nwk_read(TREE.as_bytes())?;
      let names = parsed.names.clone();
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
        names,
        layouts: whole_layouts(&sources)?,
        alignments: sources,
        branch_lengths: parsed.branch_lengths,
      };
      run(&params, input, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())
    }

    fn whole_layouts(sources: &[SequenceSource]) -> Result<Vec<PartitionLayout>, Report> {
      let lengths = sources.iter().map(|source| source.length).collect::<Vec<_>>();
      Ok(resolve_partition_layouts(&[], &lengths)?.layouts)
    }

    pub(super) fn internal_node_names(output: &PruneOutput) -> BTreeSet<&str> {
      output
        .graph
        .get_internal_nodes()
        .map(|node| &output.names[node.key()])
        .collect()
    }
  }
}
