#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::GtrModel;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::partition_rates::{fix_equal_rates, normalize_partition_rates, update_partition_rates};
  use crate::optimize::run_loop::record_if_best;
  use crate::progress::NoopProgress;
  use crate::test_utils::{Columns, TestDataset};
  use eyre::Report;
  use helpers::{
    Dataset, SmallSet, blocks_from_strings, build_set, build_set_with_model, column_blocks, internal_edge_above,
    objective, optimize, optimize_codon_positions, rate_blocks, rate_step_from_rates, reversed,
  };
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::LogLh;
  use treetime_utils::pretty_assert_abs_diff_eq;
  use treetime_utils::pretty_assert_map_abs_diff_eq;
  use treetime_utils::pretty_assert_ulps_eq;

  const BLOCK: Columns = Columns::Leading(400);

  const CODON_BLOCK: Columns = Columns::Leading(399);

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_dense(  TestDataset::Zika20,  true,  BLOCK       )]
  #[case::zika_sparse( TestDataset::Zika20,  false, Columns::All)]
  #[case::flu_dense(   TestDataset::Flu20,   true,  BLOCK       )]
  #[case::flu_sparse(  TestDataset::Flu20,   false, Columns::All)]
  #[case::ebola_dense( TestDataset::Ebola20, true,  BLOCK       )]
  #[case::ebola_sparse(TestDataset::Ebola20, false, Columns::All)]
  #[trace]
  fn test_optimize_two_column_blocks_give_the_branch_lengths_of_one_alignment(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
    #[case] selection: Columns,
  ) -> Result<(), Report> {
    let dataset = Dataset::new(dataset, selection)?;
    let whole = optimize(&dataset, &[dataset.whole()], dense, PartitionRates::Equal)?;
    let split = optimize(&dataset, &column_blocks(&dataset)?, dense, PartitionRates::Equal)?;
    assert_eq!(whole.edges, split.edges);
    pretty_assert_map_abs_diff_eq!(whole.branch_lengths, split.branch_lengths, epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_dense(  TestDataset::Zika20,  true,  CODON_BLOCK          )]
  #[case::flu_dense(   TestDataset::Flu20,   true,  CODON_BLOCK          )]
  #[case::ebola_dense( TestDataset::Ebola20, true,  CODON_BLOCK          )]
  #[case::zika_sparse( TestDataset::Zika20,  false, Columns::WholeCodons )]
  #[case::flu_sparse(  TestDataset::Flu20,   false, Columns::WholeCodons )]
  #[case::ebola_sparse(TestDataset::Ebola20, false, Columns::WholeCodons )]
  #[trace]
  fn test_optimize_codon_positions_give_the_branch_lengths_and_mutations_of_one_alignment(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
    #[case] selection: Columns,
  ) -> Result<(), Report> {
    let dataset = Dataset::new(dataset, selection)?;
    let whole = optimize(&dataset, &[dataset.whole()], dense, PartitionRates::Equal)?;
    let split = optimize_codon_positions(&dataset, dense)?;
    assert_eq!(
      (&whole.edges, &whole.sequences, &whole.indel_rates),
      (&split.edges, &split.sequences, &split.indel_rates)
    );
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
    let dataset = Dataset::new(dataset, BLOCK)?;
    let blocks = column_blocks(&dataset)?;
    let forward = optimize(&dataset, &blocks, dense, PartitionRates::Estimate)?;
    let backward = optimize(&dataset, &reversed(&blocks), dense, PartitionRates::Estimate)?;
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
    let indels = IndelPolicy::Include.edge_indels(&graph, &partitions, &branch_lengths)?;
    let term = indels.term(internal_edge_above(&graph, &names, "AB"));
    assert_eq!(
      (6, 2, btreemap! {1 => 2.0, 2 => 2.0}),
      (partitions.total_length(), term.count, indels.alignment_rates())
    );
    pretty_assert_ulps_eq!(4.0, term.rate, max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( true )]
  #[case::sparse(false)]
  #[trace]
  fn test_optimize_rate_normalization_keeps_the_objective_on_branches_above_the_floor(
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let SmallSet {
      graph,
      branch_lengths,
      mut partitions,
      ..
    } = build_set(&rate_blocks()?, dense)?;
    partitions.set_rates(&[0.5, 3.0])?;
    let mut indels = IndelPolicy::Include.edge_indels(&graph, &partitions, &branch_lengths)?;
    let indel_rate = indels.rate();
    let before = objective(&graph, &partitions, &branch_lengths, &indels)?;

    let mut normalized_branch_lengths = branch_lengths;
    normalize_partition_rates(&mut partitions, &mut normalized_branch_lengths, Some(&mut indels))?;

    let mean = (8.0 * 0.5 + 6.0 * 3.0) / 14.0;
    let rates = partitions.rates();
    pretty_assert_ulps_eq!(0.5 / mean, rates[0], max_ulps = 4);
    pretty_assert_ulps_eq!(3.0 / mean, rates[1], max_ulps = 4);
    pretty_assert_ulps_eq!(indel_rate / mean, indels.rate(), max_ulps = 4);
    let after = objective(&graph, &partitions, &normalized_branch_lengths, &indels)?;
    pretty_assert_abs_diff_eq!(before, after, epsilon = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_sparse_double(TestDataset::Zika20, false, 2.0)]
  #[case::zika_sparse_half(  TestDataset::Zika20, false, 0.5)]
  #[case::flu_dense_double(  TestDataset::Flu20,  true,  2.0)]
  #[case::flu_dense_half(    TestDataset::Flu20,  true,  0.5)]
  #[trace]
  fn test_optimize_rate_step_does_not_depend_on_the_scale_of_the_starting_rates(
    #[case] dataset: TestDataset,
    #[case] dense: bool,
    #[case] factor: f64,
  ) -> Result<(), Report> {
    let dataset = Dataset::new(dataset, BLOCK)?;
    let blocks = column_blocks(&dataset)?;
    let reference = rate_step_from_rates(&dataset, &blocks, dense, 1.0)?;
    let scaled = rate_step_from_rates(&dataset, &blocks, dense, factor)?;
    pretty_assert_abs_diff_eq!(reference.rates[0], scaled.rates[0], epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(reference.rates[1], scaled.rates[1], epsilon = 1e-7);
    pretty_assert_map_abs_diff_eq!(reference.branch_lengths, scaled.branch_lengths, epsilon = 1e-8);
    Ok(())
  }

  #[test]
  fn test_optimize_restore_of_the_best_iteration_returns_its_rates_branch_lengths_and_indel_rates() -> Result<(), Report>
  {
    let SmallSet {
      graph,
      branch_lengths,
      partitions,
      ..
    } = build_set(&rate_blocks()?, false)?;
    let mut indels = IndelPolicy::Include.edge_indels(&graph, &partitions, &branch_lengths)?;
    let log_lh = LogLh::new(objective(&graph, &partitions, &branch_lengths, &indels)?);
    let mut best = None;
    record_if_best(&mut best, 0, log_lh, &branch_lengths, &partitions, &indels);
    let recorded = (partitions.rates(), branch_lengths.clone(), indels.alignment_rates());

    let mut current_branch_lengths = branch_lengths;
    let partitions = update_partition_rates(
      &graph,
      partitions,
      PartitionRates::Estimate,
      &mut current_branch_lengths,
      Some(&mut indels),
      &NoopProgress,
    )?;
    assert_ne!(recorded.0, partitions.rates());

    let best = best.expect("the first iteration is recorded");
    let partitions = best.restore(&graph, partitions, &mut current_branch_lengths, &mut indels)?;
    assert_eq!(
      recorded,
      (
        partitions.rates(),
        current_branch_lengths.clone(),
        indels.alignment_rates()
      )
    );
    let restored_log_lh = partitions.graph_log_lh(&graph)?
      + indels
        .recount(&graph, &partitions)?
        .total_log_lh(&graph, &current_branch_lengths)?;
    pretty_assert_abs_diff_eq!(log_lh.value(), restored_log_lh.value(), epsilon = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( true )]
  #[case::sparse(false)]
  #[trace]
  fn test_optimize_equal_rates_after_a_fitted_model_keep_the_messages_of_the_current_rates(
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let SmallSet {
      graph,
      mut branch_lengths,
      mut partitions,
      ..
    } = build_set_with_model(&rate_blocks()?, dense, &GtrModel::Infer)?;
    normalize_partition_rates(&mut partitions, &mut branch_lengths, None)?;
    assert_ne!(vec![1.0, 1.0], partitions.rates());

    let partitions = fix_equal_rates(&graph, partitions, PartitionRates::Equal, &branch_lengths)?;

    let (fresh, fresh_log_lh) = partitions
      .clone()
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    assert_eq!(vec![1.0, 1.0], fresh.rates());
    pretty_assert_abs_diff_eq!(
      fresh_log_lh.value(),
      partitions.graph_log_lh(&graph)?.value(),
      epsilon = 1e-10
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::indel::{EdgeIndels, IndelPolicy};
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::partition_rates::update_partition_rates;
    use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
    use crate::partition::create::Representation;
    use crate::partition::set::PartitionSet;
    use crate::partition::split::SequenceInput;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::SequenceSource;
    use crate::seq::mutation::Mutation;
    use crate::test_utils::{
      Columns, TestDataset, codon_positions_input, find_node_key_by_name, leaf_seq_inputs, sequence_source,
      split_at_gap_free_column, whole_alignment_set,
    };
    use eyre::Report;
    use std::collections::BTreeSet;
    use std::path::PathBuf;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;

    use treetime_io::nwk::{NwkParse, nwk_read};
    use treetime_primitives::{AlignmentRecord, Seq};
    use treetime_utils::o;

    pub(super) struct Dataset {
      dataset: TestDataset,
      records: Vec<AlignmentRecord>,
    }

    impl Dataset {
      pub(super) fn new(dataset: TestDataset, selection: Columns) -> Result<Self, Report> {
        Ok(Self {
          dataset,
          records: dataset.alignment_columns(selection)?,
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
      pub(super) edges: BTreeSet<(String, String)>,
      pub(super) branch_lengths: EdgeMap<f64>,
      pub(super) sequences: Vec<(Seq, EdgeMap<Vec<Mutation>>)>,
      pub(super) indel_rates: Vec<Option<f64>>,
      pub(super) paths: Vec<PathBuf>,
    }

    pub(super) fn optimize(
      dataset: &Dataset,
      blocks: &[Block],
      dense: bool,
      partition_rates: PartitionRates,
    ) -> Result<Optimized, Report> {
      run_optimize(dataset, dense, partition_rates, |graph, names| {
        SequenceInput::whole_alignments(sources(graph, names, blocks)?, &jc69())
      })
    }

    pub(super) fn optimize_codon_positions(dataset: &Dataset, dense: bool) -> Result<Optimized, Report> {
      run_optimize(dataset, dense, PartitionRates::Equal, |graph, names| {
        codon_positions_input(sources(graph, names, &[dataset.whole()])?, &jc69())
      })
    }

    fn run_optimize(
      dataset: &Dataset,
      dense: bool,
      partition_rates: PartitionRates,
      sequences: impl FnOnce(&Graph, &NodeNames) -> Result<SequenceInput, Report>,
    ) -> Result<Optimized, Report> {
      let parse = dataset.tree()?;
      let names = parse.names.clone();
      let sequences = sequences(&parse.graph, &names)?;
      let params = OptimizeParams {
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
          partition_rates,
        },
      };
      let input = OptimizeInput {
        graph: parse.graph,
        names,
        sequences,
        branch_lengths: parse.branch_lengths,
      };
      let output = run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)?;
      let edges = output
        .graph
        .get_edges()
        .map(|edge| {
          (
            output.names[edge.source()].to_owned(),
            output.names[edge.target()].to_owned(),
          )
        })
        .collect();
      let sequences = output
        .alignments
        .iter()
        .map(|alignment| {
          let sequences = alignment.sequences.as_ref().expect("optimize returns sequences");
          (sequences.root_sequence.clone(), sequences.edge_mutations.clone())
        })
        .collect();
      Ok(Optimized {
        edges,
        branch_lengths: branch_lengths_or_zero(&output.branch_lengths),
        sequences,
        indel_rates: output.alignments.iter().map(|alignment| alignment.indel_rate).collect(),
        paths: output.alignments.into_iter().map(|result| result.path).collect(),
      })
    }

    fn jc69() -> GtrModel {
      GtrModel::JC69(JC69Params::default())
    }

    pub(super) struct SmallSet {
      pub(super) graph: Graph,
      pub(super) names: NodeNames,
      pub(super) branch_lengths: EdgeMap<Option<f64>>,
      pub(super) partitions: PartitionSet,
    }

    pub(super) fn build_set(blocks: &[Block], dense: bool) -> Result<SmallSet, Report> {
      build_set_with_model(blocks, dense, &GtrModel::JC69(JC69Params::default()))
    }

    pub(super) fn build_set_with_model(blocks: &[Block], dense: bool, model: &GtrModel) -> Result<SmallSet, Report> {
      let parse = nwk_read(b"((A:0.1,B:0.1)AB:0.1,C:0.1,D:0.1)root;".as_slice())?;
      let names = parse.names.clone();
      let sources = sources(&parse.graph, &names, blocks)?;
      let representation = if dense {
        Representation::Dense
      } else {
        Representation::Sparse
      };
      let branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let partitions = whole_alignment_set(representation, model, &parse.graph, sources, &branch_lengths)?;
      let (partitions, _) = partitions.marginal_update(&parse.graph, &branch_lengths)?;
      Ok(SmallSet {
        graph: parse.graph,
        names,
        branch_lengths: parse.branch_lengths,
        partitions,
      })
    }

    pub(super) fn rate_blocks() -> Result<Vec<Block>, Report> {
      blocks_from_strings(&[
        &[
          ("A", "ACGTA-GT"),
          ("B", "ACGTA-GA"),
          ("C", "ACCTACGT"),
          ("D", "TCCTACGT"),
        ],
        &[("A", "GGT--A"), ("B", "GGTTCA"), ("C", "GATTCA"), ("D", "CATTCA")],
      ])
    }

    pub(super) fn objective(
      graph: &Graph,
      partitions: &PartitionSet,
      branch_lengths: &EdgeMap<Option<f64>>,
      indels: &EdgeIndels,
    ) -> Result<f64, Report> {
      let (partitions, log_lh) = partitions
        .clone()
        .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
      let indel_log_lh = indels
        .recount(graph, &partitions)?
        .total_log_lh(graph, branch_lengths)?;
      Ok((log_lh + indel_log_lh).value())
    }

    pub(super) struct RateStep {
      pub(super) rates: Vec<f64>,
      pub(super) branch_lengths: EdgeMap<f64>,
    }

    pub(super) fn rate_step_from_rates(
      dataset: &Dataset,
      blocks: &[Block],
      dense: bool,
      start_rate: f64,
    ) -> Result<RateStep, Report> {
      let parse = dataset.tree()?;
      let names = parse.names.clone();
      let sources = sources(&parse.graph, &names, blocks)?;
      let input_branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let mut partitions = whole_alignment_set(
        Representation::resolve(Some(dense)),
        &GtrModel::JC69(JC69Params::default()),
        &parse.graph,
        sources,
        &input_branch_lengths,
      )?;
      partitions.set_rates(&[start_rate, start_rate])?;
      let (partitions, _) = partitions.marginal_update(&parse.graph, &input_branch_lengths)?;
      let mut branch_lengths = parse.branch_lengths;
      let partitions = update_partition_rates(
        &parse.graph,
        partitions,
        PartitionRates::Estimate,
        &mut branch_lengths,
        None,
        &NoopProgress,
      )?;
      Ok(RateStep {
        rates: partitions.rates(),
        branch_lengths: branch_lengths_or_zero(&branch_lengths),
      })
    }

    pub(super) fn internal_edge_above(graph: &Graph, names: &NodeNames, name: &str) -> GraphEdgeKey {
      let node_key = find_node_key_by_name(graph, names, name).expect("the tree has the named node");
      graph.get_node(node_key).expect("node exists").inbound()[0]
    }

    fn sources(graph: &Graph, names: &NodeNames, blocks: &[Block]) -> Result<Vec<SequenceSource>, Report> {
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
