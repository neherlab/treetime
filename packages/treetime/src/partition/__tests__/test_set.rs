#[cfg(test)]
mod tests {
  use crate::partition::create::Representation;
  use crate::partition::set::PartitionSet;
  use crate::test_utils::{TestDataset, split_at_gap_free_column};
  use eyre::Report;
  use helpers::{Resolution, clade_edge_log_lh, resolved_polytomy, updated_polytomy_set, updated_set};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_map_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Representation::Dense )]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_partition_set_of_two_column_blocks_sums_to_the_whole_alignment(
    #[case] representation: Representation,
  ) -> Result<(), Report> {
    let whole = TestDataset::Flu20.alignment()?;
    let [first, second] = split_at_gap_free_column(&whole)?;
    let (graph, whole) = updated_set(representation, vec![whole])?;
    let (_, split) = updated_set(representation, vec![first, second])?;

    pretty_assert_abs_diff_eq!(
      whole.graph_log_lh(&graph)?.value(),
      split.graph_log_lh(&graph)?.value(),
      epsilon = 1e-8
    );
    let sub_counts = |set: &PartitionSet| -> Result<Vec<Option<usize>>, Report> {
      graph.get_edges().map(|edge| set.edge_sub_count(&graph, edge.key())).collect()
    };
    assert_eq!(
      (whole.total_length(), sub_counts(&whole)?),
      (split.total_length(), sub_counts(&split)?)
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_merger(                        Resolution::MergeTwo)]
  #[case::merger_chain_removes_the_polytomy( Resolution::MergeAll)]
  #[trace]
  fn test_partition_set_sparse_log_lh_matches_dense_after_polytomy_resolution(
    #[case] resolution: Resolution,
  ) -> Result<(), Report> {
    let dense = resolved_polytomy(Representation::Dense, resolution)?;

    let sparse = resolved_polytomy(Representation::Sparse, resolution)?;

    pretty_assert_abs_diff_eq!(
      dense.set.graph_log_lh(&dense.graph)?.value(),
      sparse.set.graph_log_lh(&sparse.graph)?.value(),
      epsilon = 1e-12
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_merger(                        Resolution::MergeTwo, "(((A:0.01,B:0.02):0.005,C:0.015)P:0.01,D:0.03)R;")]
  #[case::merger_chain_removes_the_polytomy( Resolution::MergeAll, "(((A:0.01,B:0.02):0.005,C:0.015):0.015,D:0.03)R;")]
  #[trace]
  fn test_partition_set_sparse_edge_likelihoods_after_polytomy_resolution_match_a_sparse_build_on_the_resolved_tree(
    #[case] resolution: Resolution,
    #[case] resolved_tree: &str,
  ) -> Result<(), Report> {
    let built = updated_polytomy_set(Representation::Sparse, resolved_tree)?;

    let reconciled = resolved_polytomy(Representation::Sparse, resolution)?;

    pretty_assert_map_abs_diff_eq!(
      &clade_edge_log_lh(&built)?,
      &clade_edge_log_lh(&reconciled)?,
      epsilon = 1e-12
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::gather::EdgeContributions;
    use crate::optimize::indel::IndelTerm;
    use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
    use crate::partition::create::Representation;
    use crate::partition::set::PartitionSet;
    use crate::test_utils::{
      TestDataset, find_node_key_by_name, leaf_seq_inputs, parent_edge_key, sequence_source, whole_alignment_set,
    };
    use eyre::Report;
    use indoc::indoc;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::key_map::EdgeMap;

    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::reroot::remove_node_if_trivial;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    const POLYTOMY_TREE: &str = "((A:0.01,B:0.02,C:0.015)P:0.01,D:0.03)R;";

    const POLYTOMY_ALIGNMENT: &str = indoc! {"
      >A
      ACGTATGTACGTGCGTACGTCCGTACGTACGTACGTACGT
      >B
      ACGTATGTACGTGCGTACGTACGTACGTACATACGTACGT
      >C
      ACGTACGTCCGTACGTACGTACGTACGTACGTACGTACGT
      >D
      ACGTACGTACGTACGAACGTACGTAGGTACGTACGTACGT
    "};

    const MERGER_BRANCH_LENGTH: f64 = 0.005;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Resolution {
      MergeTwo,
      MergeAll,
    }

    pub(super) struct PolytomySet {
      pub(super) graph: Graph,
      pub(super) names: NodeNames,
      pub(super) branch_lengths: EdgeMap<f64>,
      pub(super) set: PartitionSet,
    }

    pub(super) fn updated_polytomy_set(representation: Representation, tree: &str) -> Result<PolytomySet, Report> {
      let parse = nwk_read(tree.as_bytes())?;
      let names = parse.names.clone();
      let graph = parse.graph;
      let records = fasta_read(POLYTOMY_ALIGNMENT.as_bytes(), &Alphabet::default())?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let source = sequence_source(Alphabet::default(), leaf_seq_inputs(&graph, &names, records))?;
      let branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let set = whole_alignment_set(
        representation,
        &GtrModel::JC69(JC69Params::default()),
        &graph,
        vec![source],
        &branch_lengths,
      )?;
      let (set, _) = set.marginal_update(&graph, &branch_lengths)?;
      Ok(PolytomySet {
        graph,
        names,
        branch_lengths,
        set,
      })
    }

    pub(super) fn resolved_polytomy(
      representation: Representation,
      resolution: Resolution,
    ) -> Result<PolytomySet, Report> {
      let PolytomySet {
        mut graph,
        names,
        mut branch_lengths,
        set,
      } = updated_polytomy_set(representation, POLYTOMY_TREE)?;
      let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("named node exists");
      let (p, a, b, c) = (key("P"), key("A"), key("B"), key("C"));
      let mut merge = |graph: &mut Graph, parent, children: [GraphNodeKey; 2]| -> Result<GraphNodeKey, Report> {
        let merger = graph.add_node();
        for child in children {
          graph.reparent_edge(parent_edge_key(graph, child), merger)?;
        }
        let edge = graph.add_edge(parent, merger)?;
        branch_lengths.insert(edge, MERGER_BRANCH_LENGTH);
        Ok(merger)
      };
      let merges = match resolution {
        Resolution::MergeTwo => {
          merge(&mut graph, p, [a, b])?;
          graph.build()?;
          vec![]
        },
        Resolution::MergeAll => {
          let ab = merge(&mut graph, p, [a, b])?;
          merge(&mut graph, p, [ab, c])?;
          graph.build()?;
          let parent_branch = Some(branch_lengths[&parent_edge_key(&graph, p)]);
          let merged =
            remove_node_if_trivial(&mut graph, p, parent_branch, Some(MERGER_BRANCH_LENGTH))?.expect("P has one child");
          branch_lengths.insert(
            merged.merged_edge_key,
            merged.merged_branch_length.expect("merged length"),
          );
          vec![merged]
        },
      };
      branch_lengths.retain(|key, _| graph.get_edge(key).is_some());
      let (set, _) = set
        .reconcile_topology(&graph, &merges)?
        .marginal_update(&graph, &branch_lengths)?;
      Ok(PolytomySet {
        graph,
        names,
        branch_lengths,
        set,
      })
    }

    pub(super) fn clade_edge_log_lh(polytomy: &PolytomySet) -> Result<BTreeMap<(Vec<&str>, u8), f64>, Report> {
      let PolytomySet { graph, names, set, .. } = polytomy;
      let mut ids = TipIds::default();
      let keys = graph_clade_keys(graph, names, &mut ids)?;
      let mut values = BTreeMap::new();
      for edge in graph.get_edges() {
        let clade = keys[&edge.target()]
          .iter()
          .map(|&id| ids.name(id))
          .sorted()
          .collect_vec();
        for (index, branch_length) in [(0_u8, 0.001), (1, 0.01), (2, 0.1)] {
          let log_lh = set.with_edge(edge.key(), |contributions| {
            evaluate_with_indels_log_lh_only(contributions, IndelTerm::new(0, 0.0), branch_length)
          })?;
          values.insert((clade.clone(), index), log_lh.value());
        }
      }
      Ok(values)
    }

    pub(super) fn updated_set(
      representation: Representation,
      files: Vec<Vec<AlignmentRecord>>,
    ) -> Result<(Graph, PartitionSet), Report> {
      let parse = nwk_read(TestDataset::Flu20.tree().as_bytes())?;
      let names = parse.names.clone();
      let sources = files
        .into_iter()
        .map(|records| sequence_source(Alphabet::default(), leaf_seq_inputs(&parse.graph, &names, records)))
        .collect::<Result<Vec<_>, Report>>()?;
      let branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let set = whole_alignment_set(
        representation,
        &GtrModel::JC69(JC69Params::default()),
        &parse.graph,
        sources,
        &branch_lengths,
      )?;
      let (set, _) = set.marginal_update(&parse.graph, &branch_lengths)?;
      Ok((parse.graph, set))
    }
  }
}
