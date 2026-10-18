#[cfg(test)]
mod tests {
  use crate::optimize::indel::IndelPolicy;
  use crate::test_utils::TestDataset;
  use eyre::Report;
  use helpers::{
    Dataset, SmallSet, blocks_from_strings, build_set, column_blocks, internal_edge_above, optimize, reversed,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_map_abs_diff_eq;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_dense(      TestDataset::Zika20,  true )]
  #[case::zika_sparse(     TestDataset::Zika20,  false)]
  #[case::flu_dense(       TestDataset::Flu20,   true )]
  #[case::flu_sparse(      TestDataset::Flu20,   false)]
  #[case::ebola_dense(     TestDataset::Ebola20, true )]
  #[case::ebola_sparse(    TestDataset::Ebola20, false)]
  #[trace]
  fn test_optimize_two_column_blocks_give_the_branch_lengths_of_one_alignment(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let dataset = Dataset::block_of(dataset)?;
    let whole = optimize(&dataset, &[dataset.whole()], dense)?;
    let split = optimize(&dataset, &column_blocks(&dataset)?, dense)?;
    assert_eq!(whole.edges, split.edges);
    pretty_assert_map_abs_diff_eq!(whole.branch_lengths, split.branch_lengths, epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_dense(   TestDataset::Flu20,  true )]
  #[case::ebola_sparse(TestDataset::Ebola20,false)]
  #[trace]
  fn test_optimize_results_do_not_depend_on_the_order_of_the_alignments(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let dataset = Dataset::block_of(dataset)?;
    let blocks = column_blocks(&dataset)?;
    let forward = optimize(&dataset, &blocks, dense)?;
    let backward = optimize(&dataset, &reversed(&blocks), dense)?;
    assert_eq!(forward.edges, backward.edges);
    pretty_assert_map_abs_diff_eq!(forward.branch_lengths, backward.branch_lengths, epsilon = 1e-10);
    assert_eq!(
      forward.paths.iter().rev().collect::<Vec<_>>(),
      backward.paths.iter().collect::<Vec<_>>()
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( true )]
  #[case::sparse(false)]
  #[trace]
  fn test_optimize_deletion_across_the_cut_counts_once_per_alignment_with_the_summed_rate(
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let blocks = blocks_from_strings(&[
      &[("A", "AC-"), ("B", "AC-"), ("C", "ACG"), ("D", "ACG")],
      &[("A", "-TT"), ("B", "-TT"), ("C", "ATT"), ("D", "ATT")],
    ])?;
    let SmallSet {
      graph,
      names,
      branch_lengths,
      partitions,
    } = build_set(&blocks, dense)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &partitions, &branch_lengths);
    let term = indels.term(internal_edge_above(&graph, &names, "AB"));
    assert_eq!(
      (6, 2, vec![Some(2.0), Some(2.0)]),
      (partitions.total_length(), term.count, indels.partition_rates(partitions.len()))
    );
    pretty_assert_ulps_eq!(4.0, term.rate, max_ulps = 4);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
    use crate::partition::create::Representation;
    use crate::partition::set::{PartitionSet, build_partition_set};
    use crate::progress::NoopProgress;
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{
      TestDataset, find_node_key_by_name, leaf_seq_inputs, sequence_source, split_at_gap_free_column,
    };
    use eyre::Report;
    use std::collections::{BTreeMap, BTreeSet};
    use std::path::PathBuf;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::{NwkParse, nwk_read};
    use treetime_primitives::{AlignmentRecord, Seq};
    use treetime_utils::o;

    const BLOCK_COLUMNS: usize = 400;

    pub(super) struct Dataset {
      dataset: TestDataset,
      records: Vec<AlignmentRecord>,
    }

    impl Dataset {
      pub(super) fn block_of(dataset: TestDataset) -> Result<Self, Report> {
        Ok(Self {
          dataset,
          records: dataset.leading_columns(BLOCK_COLUMNS)?,
        })
      }

      pub(super) fn whole(&self) -> Block {
        Block {
          path: o!("whole.fasta"),
          records: self.records.clone(),
        }
      }

      fn tree(&self) -> Result<NwkParse, Report> {
        nwk_read(self.dataset.tree().as_bytes())
      }
    }

    #[derive(Clone)]
    pub(super) struct Block {
      path: String,
      records: Vec<AlignmentRecord>,
    }

    pub(super) fn column_blocks(dataset: &Dataset) -> Result<Vec<Block>, Report> {
      let [first, second] = split_at_gap_free_column(&dataset.records)?;
      Ok(vec![
        Block {
          path: o!("first.fasta"),
          records: first,
        },
        Block {
          path: o!("second.fasta"),
          records: second,
        },
      ])
    }

    pub(super) fn reversed(blocks: &[Block]) -> Vec<Block> {
      blocks.iter().rev().cloned().collect()
    }

    pub(super) fn blocks_from_strings(files: &[&[(&str, &str)]]) -> Result<Vec<Block>, Report> {
      files
        .iter()
        .enumerate()
        .map(|(index, file)| {
          let records = file
            .iter()
            .map(|(name, seq)| {
              Ok(AlignmentRecord {
                name: (*name).to_owned(),
                seq: Seq::try_from_str(seq)?,
              })
            })
            .collect::<Result<Vec<_>, Report>>()?;
          Ok(Block {
            path: format!("{}.fasta", index + 1),
            records,
          })
        })
        .collect()
    }

    pub(super) struct Optimized {
      pub(super) edges: BTreeSet<(Option<String>, Option<String>)>,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, f64>,
      pub(super) paths: Vec<PathBuf>,
    }

    pub(super) fn optimize(dataset: &Dataset, blocks: &[Block], dense: bool) -> Result<Optimized, Report> {
      let parse = dataset.tree()?;
      let names = parse.names();
      let sources = sources(&parse.graph, &names, blocks)?;
      let params = OptimizeParams {
        model: GtrModelName::JC69,
        dense: Some(dense),
        initial_guess: InitialGuessMode::Auto,
        reroot_spec: None,
        loop_settings: OptimizeLoopSettings {
          max_iter: 3,
          dp: 0.1,
          damping: 0.0,
          opt_method: BranchOptMethod::Brent,
          indel_policy: IndelPolicy::Include,
          topology_ops: TopologyOps::default(),
        },
      };
      let input = OptimizeInput {
        graph: parse.graph,
        sources,
        branch_lengths: parse.branch_lengths,
      };
      let output = run(&params, input, &names, &NoopCancel, &NoopProgress, &NoopProgress)?;
      let edges = output
        .graph
        .get_edges()
        .map(|edge| {
          (
            output.names[&edge.source()].clone(),
            output.names[&edge.target()].clone(),
          )
        })
        .collect();
      Ok(Optimized {
        edges,
        branch_lengths: branch_lengths_or_zero(&output.branch_lengths),
        paths: output.results.into_iter().map(|result| result.path).collect(),
      })
    }

    pub(super) struct SmallSet {
      pub(super) graph: Graph,
      pub(super) names: BTreeMap<GraphNodeKey, Option<String>>,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      pub(super) partitions: PartitionSet,
    }

    pub(super) fn build_set(blocks: &[Block], dense: bool) -> Result<SmallSet, Report> {
      let parse = nwk_read(b"((A:0.1,B:0.1)AB:0.1,C:0.1,D:0.1)root;".as_slice())?;
      let names = parse.names();
      let sources = sources(&parse.graph, &names, blocks)?;
      let representation = if dense {
        Representation::Dense
      } else {
        Representation::Sparse
      };
      let branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let (_, partitions) = build_partition_set(
        representation,
        GtrModelName::JC69,
        &parse.graph,
        sources,
        &branch_lengths,
        &NoopProgress,
      )?;
      let (partitions, _) = partitions.marginal_update(&parse.graph, &branch_lengths)?;
      Ok(SmallSet {
        graph: parse.graph,
        names,
        branch_lengths: parse.branch_lengths,
        partitions,
      })
    }

    pub(super) fn internal_edge_above(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      name: &str,
    ) -> GraphEdgeKey {
      let node_key = find_node_key_by_name(graph, names, name).expect("the tree has the named node");
      graph.get_node(node_key).expect("node exists").inbound()[0]
    }

    fn sources(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      blocks: &[Block],
    ) -> Result<Vec<SequenceSource>, Report> {
      blocks
        .iter()
        .map(|block| {
          let source = sequence_source(
            Alphabet::default(),
            leaf_seq_inputs(graph, names, block.records.clone()),
          )?;
          Ok(SequenceSource {
            path: PathBuf::from(&block.path),
            ..source
          })
        })
        .collect()
    }
  }
}
