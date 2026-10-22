#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{F81Params, GtrModel, JC69Params, Jtt92Params};
  use crate::progress::NoopProgress;
  use crate::test_utils::JC69_NUCLEOTIDE_SATURATED_DIFFERENCE;
  use crate::timetree::inference::input_branch_length::{InputBranchLikelihood, InputBranchModel, input_branch_model};
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_grid::MaxGridPoints;
  use treetime_utils::{pretty_assert_relative_eq, pretty_assert_ulps_eq};

  const SEQUENCE_LENGTH: usize = 1000;

  const BRANCH_RATE: f64 = 0.002;

  const AMINO_ACID_JC69_SATURATED_DIFFERENCE: f64 = 1.0 - 20.0 * 0.05 * 0.05;

  const V0_JTT92_SATURATED_DIFFERENCE: f64 = 0.9419118013829045;

  #[rustfmt::skip]
  #[rstest]
  #[case::at_zero(         0.0)]
  #[case::below_one_scale( 0.25)]
  #[case::inside_window(   1.0)]
  #[case::near_window_end( 3.0)]
  #[case::soft_tail(      10.0)]
  #[trace]
  fn test_input_branch_length_zero_length_decays_linearly_in_time(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::model(SEQUENCE_LENGTH)
      .likelihood(0.0)
      .distribution(BRANCH_RATE, MaxGridPoints::default())?;

    let expected = t * BRANCH_RATE * 1000.0;
    pretty_assert_relative_eq!(expected, distribution.eval(t)?, epsilon = 1e-9, max_relative = 1e-12);
    Ok(())
  }

  #[test]
  fn test_input_branch_length_gaussian_far_tail_follows_v0_linear_penalty() -> Result<(), Report> {
    let likelihood = helpers::model(7000).likelihood(1.5);
    let rate = 0.01;
    let distribution = likelihood.distribution(rate, MaxGridPoints::default())?;

    let (t_near, t_far) = (70.0, 10.0);
    let expected = likelihood.neg_log_lh(t_far * rate) - likelihood.neg_log_lh(t_near * rate);
    let actual = distribution.eval(t_far)? - distribution.eval(t_near)?;
    pretty_assert_relative_eq!(expected, actual, epsilon = 1e-9, max_relative = 1e-9);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::poisson_ten_mutations((  1000, 0.01), 0.01 - 1e-6)]
  #[case::poisson_long_genome(  ( 30000, 0.02), 0.02 - 1e-3 / 30000.0)]
  #[case::gaussian(             (  1000, 0.3 ), 0.3)]
  #[case::gaussian_narrow(      (1000000, 0.5), 0.5)]
  #[trace]
  fn test_input_branch_length_peak_lies_within_one_grid_step_of_the_mode(
    #[case] (sequence_length, mutation_length): (usize, f64),
    #[case] mode: f64,
  ) -> Result<(), Report> {
    let distribution = helpers::model(sequence_length)
      .likelihood(mutation_length)
      .distribution(BRANCH_RATE, MaxGridPoints::default())?;

    let peak = distribution.likely_time()?.expect("a branch distribution has a peak");
    let dx = helpers::function(&distribution).dx();
    let mode_time = mode / BRANCH_RATE;
    assert!(
      (peak - mode_time).abs() <= dx,
      "peak {peak} must lie within one grid step {dx} of the mode {mode_time}"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::poisson_ten_mutations((   1000, 0.01), 3.16e-3)]
  #[case::gaussian_narrow(      (1000000, 0.5 ), 9.37e-4)]
  #[trace]
  fn test_input_branch_length_grid_resolves_the_peak_width(
    #[case] (sequence_length, mutation_length): (usize, f64),
    #[case] width: f64,
  ) -> Result<(), Report> {
    let distribution = helpers::model(sequence_length)
      .likelihood(mutation_length)
      .distribution(BRANCH_RATE, MaxGridPoints::default())?;

    let dx = helpers::function(&distribution).dx() * BRANCH_RATE;
    assert!(dx <= width / 9.0, "grid step {dx} must resolve the peak width {width} with at least 9 steps");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_alignment(false)]
  #[case::with_alignment(   true)]
  #[trace]
  fn test_input_branch_model_takes_the_frequencies_of_the_model_parameters(#[case] with_alignment: bool) -> Result<(), Report> {
    let model = GtrModel::F81(F81Params {
      pi: Some(array![0.4, 0.3, 0.2, 0.1]),
      ..F81Params::default()
    });
    let (graph, branch_lengths, sources) = helpers::cherry(with_alignment)?;

    let actual = input_branch_model(&model, None, SEQUENCE_LENGTH, &graph, &branch_lengths, sources, &NoopProgress)?;

    pretty_assert_ulps_eq!(
      helpers::long_branch_variance(InputBranchModel::new(SEQUENCE_LENGTH, 0.7)),
      helpers::long_branch_variance(actual),
      max_ulps = 8
    );
    Ok(())
  }

  #[rstest]
  #[case::without_alignment(false)]
  #[case::with_alignment(true)]
  #[trace]
  fn test_input_branch_model_named_jc69_has_the_four_state_saturation(
    #[case] with_alignment: bool,
  ) -> Result<(), Report> {
    let (graph, branch_lengths, sources) = helpers::cherry(with_alignment)?;

    let actual = input_branch_model(
      &GtrModel::JC69(JC69Params::default()),
      None,
      SEQUENCE_LENGTH,
      &graph,
      &branch_lengths,
      sources,
      &NoopProgress,
    )?;

    assert_eq!(InputBranchModel::new(SEQUENCE_LENGTH, 0.75), actual);
    Ok(())
  }

  #[test]
  fn test_input_branch_model_without_alignment_infers_nothing_and_takes_jc69() -> Result<(), Report> {
    let (graph, branch_lengths, sources) = helpers::cherry(false)?;

    let actual = input_branch_model(
      &GtrModel::Infer,
      None,
      SEQUENCE_LENGTH,
      &graph,
      &branch_lengths,
      sources,
      &NoopProgress,
    )?;

    assert_eq!(InputBranchModel::new(SEQUENCE_LENGTH, 0.75), actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69_amino_acids(      GtrModel::JC69(JC69Params::default()),   Some(AlphabetName::AaNoStop), AMINO_ACID_JC69_SATURATED_DIFFERENCE)]
  #[case::infer_amino_acids(     GtrModel::Infer,                         Some(AlphabetName::AaNoStop), AMINO_ACID_JC69_SATURATED_DIFFERENCE)]
  #[case::jtt92_alphabet_given(  GtrModel::Jtt92(Jtt92Params::default()), Some(AlphabetName::AaNoStop), V0_JTT92_SATURATED_DIFFERENCE)]
  #[case::jtt92_alphabet_omitted(GtrModel::Jtt92(Jtt92Params::default()), None,                         V0_JTT92_SATURATED_DIFFERENCE)]
  #[trace]
  fn test_input_branch_model_without_alignment_builds_the_model_on_the_run_alphabet(
    #[case] model: GtrModel,
    #[case] alphabet: Option<AlphabetName>,
    #[case] saturated_difference: f64,
  ) -> Result<(), Report> {
    let (graph, branch_lengths, sources) = helpers::cherry(false)?;

    let actual = input_branch_model(&model, alphabet, SEQUENCE_LENGTH, &graph, &branch_lengths, sources, &NoopProgress)?;

    pretty_assert_ulps_eq!(
      helpers::long_branch_variance(InputBranchModel::new(SEQUENCE_LENGTH, saturated_difference)),
      helpers::long_branch_variance(actual),
      max_ulps = 8
    );
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::Alphabet;
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{leaf_seq_inputs, sequence_source};
    use std::collections::BTreeMap;
    use treetime_distribution::DistributionFunction;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    type Cherry = (Graph, BTreeMap<GraphEdgeKey, Option<f64>>, Vec<SequenceSource>);

    pub(super) fn cherry(with_alignment: bool) -> Result<Cherry, Report> {
      let nwk_parsed = nwk_read(b"(A:0.1,B:0.2)root;".as_slice())?;
      let sources = if with_alignment {
        let alphabet = Alphabet::new(AlphabetName::Nuc)?;
        let aln = fasta_read(b">A\nACGTACGTAC\n>B\nACGTTCGTAA\n".as_slice(), &alphabet)?
          .into_iter()
          .map(AlignmentRecord::from)
          .collect();
        let nodes = leaf_seq_inputs(&nwk_parsed.graph, &nwk_parsed.names, aln);
        vec![sequence_source(alphabet, nodes)?]
      } else {
        vec![]
      };
      Ok((nwk_parsed.graph, nwk_parsed.branch_lengths, sources))
    }

    pub(super) fn long_branch_variance(model: InputBranchModel) -> f64 {
      let InputBranchLikelihood::Gaussian { variance, .. } = model.likelihood(0.3) else {
        panic!("a branch of 0.3 substitutions per site has a Gaussian likelihood");
      };
      variance
    }

    pub(super) fn model(sequence_length: usize) -> InputBranchModel {
      InputBranchModel::new(sequence_length, JC69_NUCLEOTIDE_SATURATED_DIFFERENCE)
    }

    pub(super) fn function(distribution: &Distribution<NegLog>) -> &DistributionFunction<f64, NegLog> {
      let Distribution::Function(function) = distribution else {
        panic!("an input branch-length distribution must be a Function, got {distribution:?}");
      };
      function
    }
  }
}
