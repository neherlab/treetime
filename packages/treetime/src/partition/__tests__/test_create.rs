#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::{Representation, build_marginal_partition};
  use crate::partition::inferred_rate::InferredRate;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::progress::NoopProgress;
  use crate::test_utils::marginal_update_reconstruction;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::sparse_infer_unit(  Representation::Sparse, GtrModelName::Infer, InferredRate::Unit,   false, helpers::fitch_inferred_gtr_unit_rate)]
  #[case::sparse_infer_fitted(Representation::Sparse, GtrModelName::Infer, InferredRate::Fitted, false, helpers::fitch_inferred_gtr)]
  #[case::sparse_jc69(        Representation::Sparse, GtrModelName::JC69,  InferredRate::Unit,   false, helpers::jc69_gtr)]
  #[case::dense_infer_unit(   Representation::Dense,  GtrModelName::Infer, InferredRate::Unit,   true,  helpers::fitch_inferred_gtr_unit_rate)]
  #[case::dense_infer_fitted( Representation::Dense,  GtrModelName::Infer, InferredRate::Fitted, true,  helpers::fitch_inferred_gtr)]
  #[case::dense_jc69(         Representation::Dense,  GtrModelName::JC69,  InferredRate::Unit,   true,  helpers::jc69_gtr)]
  #[trace]
  fn test_build_marginal_partition_selects_representation_and_gtr(
    #[case] representation: Representation,
    #[case] model: GtrModelName,
    #[case] rate: InferredRate,
    #[case] expected_dense: bool,
    #[case]
    #[notrace]
    expected_gtr: helpers::GtrOracle,
  ) -> Result<(), Report> {
    let input = helpers::Input::star(&[("A", "ACGT"), ("B", "ACGT"), ("C", "ACGA")])?;

    let reconstruction = build_marginal_partition(
      representation,
      model,
      rate,
      &input.graph,
      Alphabet::default(),
      input.node_inputs.clone(),
      &input.branch_lengths,
      &NoopProgress,
    )?;

    assert_eq!(expected_dense, matches!(reconstruction, MarginalReconstruction::Dense(_)));
    assert_eq!(4, reconstruction.view().length());
    assert_eq!(
      format!("{:#?}", expected_gtr(&input)?),
      format!("{:#?}", reconstruction.gtr())
    );
    let (updated, _) = marginal_update_reconstruction(reconstruction, &input.graph, &input.branch_lengths)?;
    assert_eq!("ACGT", updated.view().root_sequence(&input.graph)?.to_string());
    Ok(())
  }

  #[test]
  fn test_build_marginal_partition_dense_named_model_reports_missing_leaf_sequence() -> Result<(), Report> {
    let input = helpers::Input::star(&[("A", "ACGT"), ("B", "ACGT")])?;

    let result = build_marginal_partition(
      Representation::Dense,
      GtrModelName::JC69,
      InferredRate::Unit,
      &input.graph,
      Alphabet::default(),
      input.node_inputs.clone(),
      &input.branch_lengths,
      &NoopProgress,
    );

    assert_error!(result, "Leaf sequence not found: 'C'");
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::gtr::gtr::GTR;
    use crate::partition::fitch::gtr_inference::infer_gtr_fitch;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::NodeSeqInput;
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};

    pub(super) type GtrOracle = fn(&Input) -> Result<GTR, Report>;

    pub(super) struct Input {
      pub(super) graph: Graph,
      pub(super) node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, f64>,
    }

    impl Input {
      pub(super) fn star(sequences: &[(&str, &str)]) -> Result<Self, Report> {
        let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1)root;".as_slice())?;
        let names = nwk_parsed.names();
        let aln = sequences
          .iter()
          .map(|(name, seq)| {
            Ok(AlignmentRecord {
              name: (*name).to_owned(),
              seq: Seq::try_from_str(seq)?,
            })
          })
          .collect::<Result<Vec<_>, Report>>()?;
        let node_inputs = leaf_seq_inputs(&nwk_parsed.graph, &names, aln);
        Ok(Self {
          branch_lengths: branch_lengths_or_zero(&nwk_parsed.branch_lengths),
          graph: nwk_parsed.graph,
          node_inputs,
        })
      }
    }

    pub(super) fn fitch_inferred_gtr(input: &Input) -> Result<GTR, Report> {
      let fitch = create_fitch_partition(&input.graph, Alphabet::default(), input.node_inputs.clone())?;
      infer_gtr_fitch(&fitch, &input.graph, &input.branch_lengths, &NoopProgress)
    }

    pub(super) fn fitch_inferred_gtr_unit_rate(input: &Input) -> Result<GTR, Report> {
      let mut gtr = fitch_inferred_gtr(input)?;
      assert!(
        (gtr.mu() - 1.0).abs() > 1e-3,
        "the fitted rate must differ from 1 to tell the rate rules apart"
      );
      gtr.set_mu(1.0)?;
      Ok(gtr)
    }

    pub(super) fn jc69_gtr(_input: &Input) -> Result<GTR, Report> {
      jc69(JC69Params::default())
    }
  }
}
