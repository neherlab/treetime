#[cfg(test)]
mod tests {
  use crate::optimize::dispatch::optimize_edge_lengths;
  use crate::optimize::gather::gather_edge_contributions;
  use crate::optimize::indel::EdgeIndels;
  use crate::optimize::params::BranchOptMethod;
  use crate::pretty_assert_abs_diff_eq;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use helpers::{DenseSetup, optimize_once, root_edge_branch_lengths};
  use indoc::indoc;
  use itertools::Itertools;
  use ndarray::Array1;
  use rstest::rstest;

  #[rustfmt::skip]
  const BIFURCATING_TREE: &str = "((A:0.05,B:0.05)AB:0.15,(C:0.05,D:0.05)CD:0.05)root:0.0;";

  const ALIGNMENT: &str = indoc! {r#"
    >A
    ACGTACGTAAGGACGT
    >B
    ACGTACGTAAGGACGA
    >C
    TCGATCGATTCCACGT
    >D
    TCGATCGATTCCACGC
  "#};

  #[test]
  fn test_root_preservation_ratio_preserved_after_optimization() -> Result<(), Report> {
    let DenseSetup {
      graph,
      reconstruction,
      mut branch_lengths,
      ..
    } = DenseSetup::new(BIFURCATING_TREE, ALIGNMENT)?;

    let (bl0_before, bl1_before) = root_edge_branch_lengths(&graph, &branch_lengths);
    let total_before = bl0_before + bl1_before;
    let ratio_before = bl0_before / total_before;

    let total_length = reconstruction.total_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &EdgeIndels::Excluded,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;

    let (bl0_after, bl1_after) = root_edge_branch_lengths(&graph, &branch_lengths);
    let total_after = bl0_after + bl1_after;
    let ratio_after = bl0_after / total_after;

    assert_abs_diff_eq!(ratio_before, ratio_after, epsilon = 1e-12);
    assert!(total_after > 0.0, "optimized total should be positive");
    Ok(())
  }

  #[test]
  fn test_root_preservation_both_edges_zero_uses_equal_split() -> Result<(), Report> {
    let tree = "((A:0.1,B:0.2)AB:0.0,(C:0.15,D:0.12)CD:0.0)root:0.0;";
    let DenseSetup {
      graph,
      reconstruction,
      mut branch_lengths,
      ..
    } = DenseSetup::new(tree, ALIGNMENT)?;

    let total_length = reconstruction.total_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &EdgeIndels::Excluded,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;

    let (bl0, bl1) = root_edge_branch_lengths(&graph, &branch_lengths);
    let total = bl0 + bl1;
    assert!(
      total > 0.0,
      "optimizer should produce positive total from zero-start edges"
    );
    let ratio = bl0 / total;
    assert_abs_diff_eq!(ratio, 0.5, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_root_preservation_skips_trifurcating_root() -> Result<(), Report> {
    let tree = "((A:0.05,B:0.05)AB:0.1,(C:0.05,D:0.05)CD:0.1,(E:0.05,F:0.05)EF:0.1)root:0.0;";
    let fasta = indoc! {r#"
      >A
      AAAAACGTACGTACGT
      >B
      AAAAACGTACGTACGA
      >C
      CCCCACGTACGTACGG
      >D
      CCCCACGTACGTACGC
      >E
      GGGGACGTACGTACGT
      >F
      GGGGACGTACGTACGA
    "#};
    let DenseSetup {
      graph,
      reconstruction,
      mut branch_lengths,
      ..
    } = DenseSetup::new(tree, fasta)?;

    {
      let root = graph.get_exactly_one_root()?;
      let children = graph.children_of(root).collect::<Vec<_>>();
      assert_eq!(3, children.len());
    }

    let total_length = reconstruction.total_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &EdgeIndels::Excluded,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;
    Ok(())
  }

  #[test]
  fn test_root_preservation_asymmetric_ratio() -> Result<(), Report> {
    let tree = "((A:0.05,B:0.05)AB:0.9,(C:0.05,D:0.05)CD:0.01)root:0.0;";
    let DenseSetup {
      graph,
      reconstruction,
      mut branch_lengths,
      ..
    } = DenseSetup::new(tree, ALIGNMENT)?;

    let (bl0_before, bl1_before) = root_edge_branch_lengths(&graph, &branch_lengths);
    let total_before = bl0_before + bl1_before;
    let ratio_before = bl0_before / total_before;

    assert!(ratio_before > 0.9, "pre-condition: asymmetric ratio");

    let total_length = reconstruction.total_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &EdgeIndels::Excluded,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;

    let (bl0_after, bl1_after) = root_edge_branch_lengths(&graph, &branch_lengths);
    let total_after = bl0_after + bl1_after;
    let ratio_after = bl0_after / total_after;

    assert_abs_diff_eq!(ratio_before, ratio_after, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::bifurcating_root(
    "((A:0.05,B:0.05)AB:0.9,(C:0.05,D:0.05)CD:0.3)root:0.0;",
    &[("AB", 0), ("CD", 0)],
    "(A:0.05,B:0.05,(C:0.05,D:0.05)CD:1.2)AB:0.0;",
    0.0,
  )]
  #[case::single_child_node(
    "((A:0.05,B:0.05)AB:0.6,((C:0.05,D:0.05)CD:0.2)U:0.1)root:0.0;",
    &[("AB", 0), ("U", 0), ("CD", 0)],
    "(A:0.05,B:0.05,(C:0.05,D:0.05)CD:0.9)AB:0.0;",
    0.0,
  )]
  #[case::bifurcating_root_with_indels(
    "((A:0.05,B:0.05)AB:0.9,(C:0.05,D:0.05)CD:0.3)root:0.0;",
    &[("AB", 2), ("CD", 1)],
    "(A:0.05,B:0.05,(C:0.05,D:0.05)CD:1.2)AB:0.0;",
    3.0,
  )]
  #[case::single_child_node_with_indels(
    "((A:0.05,B:0.05)AB:0.6,((C:0.05,D:0.05)CD:0.2)U:0.1)root:0.0;",
    &[("AB", 1), ("U", 1), ("CD", 1)],
    "(A:0.05,B:0.05,(C:0.05,D:0.05)CD:0.9)AB:0.0;",
    3.0,
  )]
  #[trace]
  fn test_root_preservation_chain_sum_matches_single_edge_of_rerooted_tree(
    #[case] tree: &str,
    #[case] chain: &[(&str, usize)],
    #[case] rerooted_tree: &str,
    #[case] indel_rate: f64,
  ) -> Result<(), Report> {
    let chain_indel_count: usize = chain.iter().map(|(_, count)| count).sum();
    let mut rerooted = DenseSetup::new(rerooted_tree, ALIGNMENT)?;
    let rerooted_indel_counts = rerooted.indel_counts(&[("CD", chain_indel_count)]);
    optimize_once(&mut rerooted, &rerooted_indel_counts, indel_rate)?;
    let expected_sum = rerooted.branch_length("CD");

    let mut setup = DenseSetup::new(tree, ALIGNMENT)?;
    let indel_counts = setup.indel_counts(chain);
    optimize_once(&mut setup, &indel_counts, indel_rate)?;
    let actual_sum: f64 = chain.iter().map(|(name, _)| setup.branch_length(name)).sum();

    assert_abs_diff_eq!(expected_sum, actual_sum, epsilon = 1e-5);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_indels(0, [0.0,                 0.0                ])]
  #[case::with_indels(   1, [0.3 * (0.5 / 0.8),   0.3 * (0.3 / 0.8)  ])]
  #[trace]
  fn test_root_preservation_chain_with_zero_optimum_of_longest_edge(
    #[case] indel_count: usize,
    #[case] expected: [f64; 2],
  ) -> Result<(), Report> {
    let tree = "((A:0.05,B:0.05)AB:0.5,(C:0.05,D:0.05)CD:0.3)root:0.0;";
    let identical = indoc! {r#"
      >A
      ACGTACGTAAGGACGT
      >B
      ACGTACGTAAGGACGT
      >C
      ACGTACGTAAGGACGT
      >D
      ACGTACGTAAGGACGT
    "#};
    let mut setup = DenseSetup::new(tree, identical)?;
    let indel_counts = setup.indel_counts(&[("AB", indel_count)]);

    optimize_once(&mut setup, &indel_counts, 1.0)?;

    let actual = [setup.branch_length("AB"), setup.branch_length("CD")];
    pretty_assert_abs_diff_eq!(Array1::from_iter(expected), Array1::from_iter(actual), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_root_preservation_single_child_chain_keeps_ratios() -> Result<(), Report> {
    let tree = "((A:0.05,B:0.05)AB:0.15,((C:0.05,D:0.05)CD:0.03)U:0.02)root:0.0;";
    let chain = ["AB", "U", "CD"];
    let mut setup = DenseSetup::new(tree, ALIGNMENT)?;
    let indel_counts = setup.indel_counts(&[]);
    let ratios = |setup: &DenseSetup| {
      let lengths = chain.iter().map(|name| setup.branch_length(name)).collect_vec();
      let total: f64 = lengths.iter().sum();
      Array1::from_iter(lengths.iter().map(|bl| bl / total))
    };
    let expected = ratios(&setup);

    optimize_once(&mut setup, &indel_counts, 0.0)?;

    pretty_assert_abs_diff_eq!(expected, ratios(&setup), epsilon = 1e-12);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::dispatch::optimize_edge_lengths;
    use crate::optimize::gather::gather_edge_contributions;
    use crate::optimize::indel::EdgeIndels;
    use crate::optimize::params::BranchOptMethod;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::reconstruction::DenseReconstruction;
    use crate::partition::set::PartitionSet;
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) struct DenseSetup {
      pub(super) graph: Graph,
      pub(super) reconstruction: PartitionSet,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      pub(super) edges_by_child_name: BTreeMap<String, GraphEdgeKey>,
    }

    impl DenseSetup {
      pub(super) fn new(newick: &str, fasta: &str) -> Result<Self, Report> {
        let alphabet = Alphabet::new(AlphabetName::Nuc)?;
        let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &alphabet)?
          .into_iter()
          .map(AlignmentRecord::from)
          .collect();
        let nwk_parsed = nwk_read(newick.as_bytes())?;
        let names = nwk_parsed.names.clone();
        let graph = nwk_parsed.graph;
        let branch_lengths = nwk_parsed.branch_lengths;
        let gtr = jc69(JC69Params::default())?;
        let partition = PartitionMarginalDense::new(alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?;
        let reconstruction = PartitionSet::Dense(vec![DenseReconstruction::seeded(partition, gtr)]);
        let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
        let edges_by_child_name = graph
          .get_edges()
          .map(|edge| (names[edge.target()].to_owned(), edge.key()))
          .collect();
        Ok(Self {
          graph,
          reconstruction,
          branch_lengths,
          edges_by_child_name,
        })
      }

      pub(super) fn branch_length(&self, child_name: &str) -> f64 {
        self.branch_lengths[&self.edges_by_child_name[child_name]].unwrap()
      }

      pub(super) fn indel_counts(&self, counts_by_child_name: &[(&str, usize)]) -> BTreeMap<GraphEdgeKey, usize> {
        let mut indel_counts: BTreeMap<GraphEdgeKey, usize> =
          self.graph.get_edges().map(|edge| (edge.key(), 0)).collect();
        for (name, count) in counts_by_child_name {
          indel_counts.insert(self.edges_by_child_name[*name], *count);
        }
        indel_counts
      }
    }

    pub(super) fn optimize_once(
      setup: &mut DenseSetup,
      indel_counts: &BTreeMap<GraphEdgeKey, usize>,
      indel_rate: f64,
    ) -> Result<(), Report> {
      let contributions = gather_edge_contributions(&setup.graph, &setup.reconstruction)?;
      optimize_edge_lengths(
        &setup.graph,
        setup.reconstruction.total_length(),
        &contributions,
        &EdgeIndels::Included {
          counts: indel_counts.clone(),
          partition_rates: vec![indel_rate],
        },
        BranchOptMethod::Brent,
        &mut setup.branch_lengths,
      )?;
      Ok(())
    }

    pub(super) fn root_edge_branch_lengths(
      graph: &Graph,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    ) -> (f64, f64) {
      let root = graph.get_exactly_one_root().unwrap();
      let children = graph.children_of(root).collect::<Vec<_>>();
      assert_eq!(2, children.len());
      let bl0 = branch_lengths[&children[0].1.key()].unwrap_or(0.0);
      let bl1 = branch_lengths[&children[1].1.key()].unwrap_or(0.0);
      (bl0, bl1)
    }
  }
}
