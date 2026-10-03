#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::params::{AncestralParams, MethodAncestral};
  use crate::ancestral::partition::AncestralPartition;
  use crate::ancestral::pipeline::run;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::marginal::sample::SampleMode;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::{AncestralInput, EdgeSeqInput, node_seq_inputs};
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::collections::container::get_exactly_one;

  const LONG_BRANCH_TREE: &str = "((A1:0.01,A2:0.01)X:1.0,(C1:0.01,C2:0.01)Y:0.001)root;";
  const LONG_BRANCH_ALIGNMENT: &str = indoc! {r#"
    >A1
    AAAAGGGGAA
    >A2
    AAAAGGGGAA
    >C1
    CAAAGGGGAT
    >C2
    CAAAGGGGAT
  "#};

  #[rstest]
  #[case::sparse(Some(false))]
  #[case::dense(Some(true))]
  fn test_root_sequence_marginal_root_follows_the_short_branch(#[case] dense: Option<bool>) -> Result<(), Report> {
    let (graph, partition) = helpers::reconstruct(
      LONG_BRANCH_TREE,
      LONG_BRANCH_ALIGNMENT,
      MethodAncestral::Marginal,
      dense,
    )?;
    assert_eq!("CAAAGGGGAT", partition.root_sequence(&graph)?.to_string());
    Ok(())
  }

  #[rstest]
  #[case::parsimony(MethodAncestral::Parsimony, None)]
  #[case::sparse(MethodAncestral::Marginal, Some(false))]
  #[case::dense(MethodAncestral::Marginal, Some(true))]
  fn test_root_sequence_equals_augur_root(
    #[case] method: MethodAncestral,
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let (graph, partition) = helpers::reconstruct(LONG_BRANCH_TREE, LONG_BRANCH_ALIGNMENT, method, dense)?;
    assert_eq!(partition.augur_root_sequence(&graph)?, partition.root_sequence(&graph)?);
    Ok(())
  }

  #[rstest]
  #[case::parsimony(MethodAncestral::Parsimony, None)]
  #[case::sparse(MethodAncestral::Marginal, Some(false))]
  #[case::dense(MethodAncestral::Marginal, Some(true))]
  fn test_root_sequence_plus_edge_subs_rebuilds_every_node(
    #[case] method: MethodAncestral,
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let (graph, partition) = helpers::reconstruct(LONG_BRANCH_TREE, LONG_BRANCH_ALIGNMENT, method, dense)?;
    let rebuilt = helpers::rebuild_from_root(&graph, &partition)?;
    let expected: BTreeMap<GraphNodeKey, Seq> = graph
      .get_nodes()
      .map(|node| Ok((node.key(), partition.augur_node_sequence(&graph, node.key())?)))
      .collect::<Result<_, Report>>()?;
    assert_eq!(expected, rebuilt);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn reconstruct(
      nwk: &str,
      fasta: &str,
      method: MethodAncestral,
      dense: Option<bool>,
    ) -> Result<(Graph, AncestralPartition), Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read_str(nwk)?;
      let names = parse.names();
      let sequences = read_many_fasta_str(fasta, &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect_vec();
      let mask = vec![false; sequences[0].seq.len()];
      let input = AncestralInput {
        nodes: node_seq_inputs(&parse.graph, &names, sequences),
        edges: parse
          .branch_lengths
          .into_iter()
          .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
          .collect(),
        graph: parse.graph,
      };
      let params = AncestralParams {
        method,
        model: GtrModelName::JC69,
        dense,
        include_leaves: true,
        impute_missing_data: false,
        gtr_iterations: 0,
        site_specific_gtr: false,
        seed: Some(0),
        sample_from_profile: SampleMode::Argmax,
      };
      let result = run(
        &params,
        &input,
        alphabet,
        mask,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let partition = result.partition.expect("ancestral reconstruction returns a partition");
      Ok((input.graph, partition))
    }

    pub(super) fn rebuild_from_root(
      graph: &Graph,
      partition: &AncestralPartition,
    ) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
      let mut rebuilt: BTreeMap<GraphNodeKey, Seq> = BTreeMap::new();
      graph.iter_depth_first_preorder_forward(|node| {
        let seq = if node.is_root {
          partition.root_sequence(graph)?
        } else {
          let (parent_key, edge_key) = get_exactly_one(&node.parent_keys)?;
          let mut seq: Seq = rebuilt[parent_key].clone();
          for sub in partition.edge_subs(graph, *edge_key)? {
            seq[sub.pos()] = sub.qry();
          }
          seq
        };
        rebuilt.insert(node.key, seq);
        Ok(())
      })?;
      Ok(rebuilt)
    }
  }
}
