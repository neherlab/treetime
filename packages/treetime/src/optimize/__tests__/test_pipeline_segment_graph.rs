#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
  use crate::partition::split::SequenceInput;
  use crate::progress::NoopProgress;
  use crate::test_utils::{NUC_ALPHABET, complete_sequence_source};
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use regex::Regex;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_graph::merge_by_name::{MergedGraph, NamedTree, merge_trees_by_name};
  use treetime_io::fasta::fasta_read_file;
  use treetime_io::nwk::TREE_DIALECT_DEFAULT;
  use treetime_io::tree::tree_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::{make_report, pretty_assert_abs_diff_eq};

  const SEGMENTS: [&str; 3] = ["seg1", "seg2", "seg3"];

  #[rstest]
  #[case::dense(true)]
  #[case::sparse(false)]
  #[trace]
  fn test_optimize_segment_graph_where_every_segment_carries_every_edge_equals_the_multi_partition_tree(
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let tree = toy_tree("seg1")?;
    let trees = [tree.clone(), tree.clone(), tree];
    let alignments = SEGMENTS
      .iter()
      .map(|segment| alignment(segment))
      .collect::<Result<Vec<_>, _>>()?;
    let model = GtrModel::Infer;
    let merged = optimize(
      &trees,
      alignments.clone(),
      &model,
      dense,
      PartitionRates::Estimate,
      true,
    )?;
    let single = optimize(&trees, alignments, &model, dense, PartitionRates::Estimate, false)?;
    assert_eq!(single.len(), merged.len());
    for (edge, length) in &single {
      assert_eq!(length.to_bits(), merged[edge].to_bits(), "edge {edge:?}");
    }
    Ok(())
  }

  #[rstest]
  #[case::dense(true)]
  #[case::sparse(false)]
  #[trace]
  fn test_optimize_segment_graph_without_shared_edges_gives_each_segment_the_lengths_of_its_own_tree(
    #[case] dense: bool,
  ) -> Result<(), Report> {
    let trees = SEGMENTS
      .iter()
      .map(|segment| own_internal_names(segment))
      .collect::<Result<Vec<_>, _>>()?;
    let alignments = SEGMENTS
      .iter()
      .map(|segment| alignment(segment))
      .collect::<Result<Vec<_>, _>>()?;
    let model = GtrModel::JC69(JC69Params::default());
    let merged = optimize(&trees, alignments.clone(), &model, dense, PartitionRates::Equal, true)?;
    let mut compared = 0;
    for (tree, records) in trees.iter().zip(alignments) {
      let separate = optimize(
        std::slice::from_ref(tree),
        vec![records],
        &model,
        dense,
        PartitionRates::Equal,
        false,
      )?;
      for (edge, length) in &separate {
        pretty_assert_abs_diff_eq!(*length, merged[edge], epsilon = 1e-6);
        compared += 1;
      }
    }
    assert_eq!(merged.len(), compared);
    Ok(())
  }

  fn toy_path(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .join("../../data/toy/reassortment")
      .join(file)
  }

  fn toy_tree(segment: &str) -> Result<String, Report> {
    Ok(fs::read_to_string(toy_path(&format!("tree_{segment}.nwk")))?)
  }

  fn own_internal_names(segment: &str) -> Result<String, Report> {
    let internal = Regex::new(r"\)([A-Z0-9]+)")?;
    Ok(
      internal
        .replace_all(&toy_tree(segment)?, format!(")${{1}}_{segment}"))
        .into_owned(),
    )
  }

  fn alignment(segment: &str) -> Result<Vec<AlignmentRecord>, Report> {
    let records = fasta_read_file(toy_path(&format!("aln_{segment}.fasta")), &*NUC_ALPHABET)?;
    Ok(records.into_iter().map(AlignmentRecord::from).collect())
  }

  fn optimize(
    trees: &[String],
    alignments: Vec<Vec<AlignmentRecord>>,
    model: &GtrModel,
    dense: bool,
    partition_rates: PartitionRates,
    as_segments: bool,
  ) -> Result<BTreeMap<(String, String), f64>, Report> {
    let parsed = trees
      .iter()
      .map(|tree| tree_read(tree.as_bytes(), TREE_DIALECT_DEFAULT))
      .collect::<Result<Vec<_>, _>>()?;
    let labels = (1..=parsed.len()).map(|index| format!("tree {index}")).collect_vec();
    let named = parsed
      .iter()
      .zip(&labels)
      .map(|(tree, label)| NamedTree {
        label,
        graph: &tree.graph,
        names: &tree.names,
        unnamed_nodes: &tree.unnamed_nodes,
        branch_lengths: &tree.branch_lengths,
      })
      .collect_vec();
    let MergedGraph {
      graph,
      names,
      branch_lengths,
      scopes,
      ..
    } = merge_trees_by_name(&named)?;
    let sources = alignments
      .into_iter()
      .enumerate()
      .map(|(index, records)| {
        complete_sequence_source(&graph, &names, records, &NUC_ALPHABET, &format!("aln_{index}.fasta"))
      })
      .collect::<Result<Vec<_>, _>>()?;
    let params = OptimizeParams {
      dense: Some(dense),
      initial_guess: InitialGuessMode::Auto,
      reroot_spec: None,
      loop_settings: OptimizeLoopSettings {
        max_iter: 50,
        dp: 1e-10,
        damping: 0.0,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::NONE,
        partition_rates,
      },
    };
    let input = OptimizeInput {
      graph,
      names,
      sequences: SequenceInput::whole_alignments(sources, model)?,
      branch_lengths,
      scopes: as_segments.then(|| Arc::new(scopes)),
    };
    let output = run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress).map_err(|error| error.into_report())?;
    let lengths = branch_lengths_or_zero(&output.branch_lengths);
    output
      .graph
      .get_edges()
      .map(|edge| {
        let length = lengths
          .get(edge.key())
          .copied()
          .ok_or_else(|| make_report!("Edge {} has no branch length", edge.key()))?;
        Ok((
          (
            output.names[edge.source()].to_owned(),
            output.names[edge.target()].to_owned(),
          ),
          length,
        ))
      })
      .collect()
  }
}
