#[cfg(test)]
mod tests {
  use crate::timetree::params::BranchLengthMode;
  use eyre::Report;
  use helpers::{Sequences, run_topology};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const ZERO_LENGTH_CHAIN: &str = "((((A:0.1,B:0.1)X:0.0,C:0.1)Y:0.0,D:0.1)Z:0.05,(E:0.1,F:0.1)W:0.1)root;";

  const ZERO_LENGTH_CHAIN_FASTA: &str = indoc! {"
    >A
    AGGTACGTACGTACGTACGTTCGTACGTACGTACGTACGT
    >B
    ACGTAGGTACGTACGTACGTTCGTACGTACGTACGTACGT
    >C
    ACGTACGTAGGTACGTACGTTCGTACGTACGTACGTACGT
    >D
    ACGTACGTACGTAGGTACGTTCGTACGTACGTACGTACGT
    >E
    ACGTACGTACGTACGTACGTACGTCCGTAAGTACGTACGT
    >F
    ACGTACGTACGTACGTACGTACGTCCGTACGTAAGTACGT
  "};

  const SHORT_EDGE: &str = "(((A:0.1,B:0.1,G:0.1)X:0.001,C:0.1,D:0.1)Y:0.05,(E:0.1,F:0.1)W:0.1)root;";

  #[rustfmt::skip]
  #[rstest]
  #[case::marginal_sparse(BranchLengthMode::Marginal, Some(false), Sequences::Alignment, "((D,C,A,B)Z,(E,F)W)root;"        )]
  #[case::marginal_dense( BranchLengthMode::Marginal, Some(true),  Sequences::Alignment, "((D,C,A,B)Z,(E,F)W)root;"        )]
  #[case::input_sparse(   BranchLengthMode::Input,    Some(false), Sequences::Alignment, "((D,C,A,B)Z,(E,F)W)root;"        )]
  #[case::input_dense(    BranchLengthMode::Input,    Some(true),  Sequences::Alignment, "((D,C,A,B)Z,(E,F)W)root;"        )]
  #[case::input_length(   BranchLengthMode::Input,    None,        Sequences::Length,    "((((A,B)X,C)Y,D)Z,(E,F)W)root;" )]
  #[trace]
  fn test_pipeline_short_branches_collapse_a_zero_length_chain_with_an_alignment(
    #[case] mode: BranchLengthMode,
    #[case] dense: Option<bool>,
    #[case] sequences: Sequences,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let actual = run_topology(ZERO_LENGTH_CHAIN, ZERO_LENGTH_CHAIN_FASTA, mode, dense, sequences, &[], true)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_change(   helpers::SHORT_EDGE_NO_CHANGE,    "((C,D,A,B,G)Y,(E,F)W)root;"   )]
  #[case::substitution(helpers::SHORT_EDGE_SUBSTITUTION, "(((A,B,G)X,C,D)Y,(E,F)W)root;")]
  #[case::gap_change(  helpers::SHORT_EDGE_GAP_CHANGE,   "(((A,B,G)X,C,D)Y,(E,F)W)root;")]
  #[trace]
  fn test_pipeline_short_branches_keep_an_edge_below_a_tenth_of_a_mutation_that_changes_the_sequence(
    #[case] fasta: &str,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let actual = run_topology(SHORT_EDGE, fasta, BranchLengthMode::Input, Some(false), Sequences::Alignment, &[], true)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::marginal(BranchLengthMode::Marginal)]
  #[case::input(   BranchLengthMode::Input)]
  #[trace]
  fn test_pipeline_short_branches_keep_a_dated_internal_node(#[case] mode: BranchLengthMode) -> Result<(), Report> {
    let actual = run_topology(
      ZERO_LENGTH_CHAIN,
      ZERO_LENGTH_CHAIN_FASTA,
      mode,
      Some(false),
      Sequences::Alignment,
      &[("Y", 2005.0)],
      true,
    )?;

    assert_eq!("(((C,A,B)Y,D)Z,(E,F)W)root;", actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::marginal(BranchLengthMode::Marginal)]
  #[case::input(   BranchLengthMode::Input)]
  #[trace]
  fn test_pipeline_short_branches_stay_when_the_collapse_is_off(#[case] mode: BranchLengthMode) -> Result<(), Report> {
    let actual = run_topology(
      ZERO_LENGTH_CHAIN,
      ZERO_LENGTH_CHAIN_FASTA,
      mode,
      Some(false),
      Sequences::Alignment,
      &[],
      false,
    )?;

    assert_eq!("((((A,B)X,C)Y,D)Z,(E,F)W)root;", actual);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::partition::split::SequenceInput;
    use crate::progress::NoopProgress;
    use crate::test_utils::{dates_by_node, leaf_seq_inputs, marginal_timetree_params, sequence_source};
    use crate::timetree::params::{BranchLengthMode, TimetreeParams};
    use crate::timetree::pipeline::{self, SequenceLength, TimetreeInput, TimetreeSequences};
    use eyre::Report;
    use indoc::indoc;
    use itertools::Itertools;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::dates_csv::DateConstraint;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
    use treetime_primitives::AlignmentRecord;
    use treetime_utils::o;

    pub(super) const SHORT_EDGE_NO_CHANGE: &str = indoc! {"
      >A
      AGGTACGTACGTACGTACGTTCGTACGTACGTACGTACGT
      >B
      ACGTAGGTACGTACGTACGTTCGTACGTACGTACGTACGT
      >G
      ACGTACGTACGTACGTAGGTTCGTACGTACGTACGTACGT
      >C
      ACGTACGTAGGTACGTACGTTCGTACGTACGTACGTACGT
      >D
      ACGTACGTACGTAGGTACGTTCGTACGTACGTACGTACGT
      >E
      ACGTACGTACGTACGTACGTACGTCCGTAAGTACGTACGT
      >F
      ACGTACGTACGTACGTACGTACGTCCGTACGTAAGTACGT
    "};

    pub(super) const SHORT_EDGE_SUBSTITUTION: &str = indoc! {"
      >A
      AGGTACGTACGTACGTACGTTCGTACGTACGTACGTGCGT
      >B
      ACGTAGGTACGTACGTACGTTCGTACGTACGTACGTGCGT
      >G
      ACGTACGTACGTACGTAGGTTCGTACGTACGTACGTGCGT
      >C
      ACGTACGTAGGTACGTACGTTCGTACGTACGTACGTACGT
      >D
      ACGTACGTACGTAGGTACGTTCGTACGTACGTACGTACGT
      >E
      ACGTACGTACGTACGTACGTACGTCCGTAAGTACGTACGT
      >F
      ACGTACGTACGTACGTACGTACGTCCGTACGTAAGTACGT
    "};

    pub(super) const SHORT_EDGE_GAP_CHANGE: &str = indoc! {"
      >A
      AGGTACGTACGTACGTACGTTCGTACGTACGTACGTA--T
      >B
      ACGTAGGTACGTACGTACGTTCGTACGTACGTACGTA--T
      >G
      ACGTACGTACGTACGTAGGTTCGTACGTACGTACGTA--T
      >C
      ACGTACGTAGGTACGTACGTTCGTACGTACGTACGTACGT
      >D
      ACGTACGTACGTAGGTACGTTCGTACGTACGTACGTACGT
      >E
      ACGTACGTACGTACGTACGTACGTCCGTAAGTACGTACGT
      >F
      ACGTACGTACGTACGTACGTACGTCCGTACGTAAGTACGT
    "};

    const SEQUENCE_LENGTH: usize = 40;

    const CLOCK_RATE: f64 = 1e-3;

    const TIP_DATES: [(&str, f64); 7] = [
      ("A", 2010.0),
      ("B", 2011.0),
      ("C", 2012.0),
      ("D", 2013.0),
      ("E", 2014.0),
      ("F", 2015.0),
      ("G", 2016.0),
    ];

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Sequences {
      Alignment,
      Length,
    }

    pub(super) fn run_topology(
      newick: &str,
      fasta: &str,
      mode: BranchLengthMode,
      dense: Option<bool>,
      sequences: Sequences,
      internal_dates: &[(&str, f64)],
      collapse_short_branches: bool,
    ) -> Result<String, Report> {
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let sequences = match sequences {
        Sequences::Alignment => {
          let alphabet = Alphabet::default();
          let records = fasta_read(fasta.as_bytes(), &alphabet)?
            .into_iter()
            .map(AlignmentRecord::from)
            .collect_vec();
          let source = sequence_source(alphabet, leaf_seq_inputs(&nwk_parsed.graph, &names, records))?;
          TimetreeSequences::Alignments(SequenceInput::whole_alignments(
            vec![source],
            &GtrModel::JC69(JC69Params::default()),
          )?)
        },
        Sequences::Length => TimetreeSequences::Length(SequenceLength {
          length: Some(SEQUENCE_LENGTH),
          model: GtrModel::JC69(JC69Params::default()),
          alphabet: None,
        }),
      };
      let dates = TIP_DATES
        .iter()
        .chain(internal_dates)
        .map(|&(name, date)| (o!(name), Some(DateConstraint::exact(date))));
      let dates = dates_by_node(dates, &nwk_parsed.graph, &names);
      let input = TimetreeInput {
        graph: nwk_parsed.graph,
        names,
        sequences,
        dates: Some(dates),
        branch_lengths: nwk_parsed.branch_lengths,
      };
      let params = TimetreeParams {
        branch_length_mode: mode,
        dense,
        clock_rate: Some(CLOCK_RATE),
        collapse_short_branches,
        ..marginal_timetree_params()
      };
      let output = pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &NoopProgress)
        .map_err(OperationError::into_report)?;
      let no_lengths: EdgeMap<Option<f64>> = output.graph.get_edges().map(|edge| (edge.key(), None)).collect();
      nwk_write_str(
        &GraphSchedule::new(&output.graph)?.into_tree()?,
        &output.names,
        &no_lengths,
        &NwkWriteOptions::default(),
        &NwkNodeComments::new(),
      )
    }
  }
}
