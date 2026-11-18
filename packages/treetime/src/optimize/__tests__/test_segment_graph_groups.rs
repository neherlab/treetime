#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::optimize::gather::{EdgeContributions, series_groups};
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::{BranchOptMethod, OptimizeLoopSettings, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::partition::create::{Representation, segment_inheritances};
  use crate::partition::set::{PartitionSet, build_inherited_partition_set};
  use crate::partition::split::SequenceInput;
  use crate::progress::NoopProgress;
  use crate::test_utils::{NUC_ALPHABET, complete_sequence_source};
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use std::fs;
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::merge_by_name::{MergedGraph, NamedTree, merge_trees_by_name};
  use treetime_io::fasta::fasta_read_file;
  use treetime_io::nwk::TREE_DIALECT_DEFAULT;
  use treetime_io::tree::tree_read;
  use treetime_primitives::{AlignmentRecord, Seq};

  const RECOMBINANT: [&str; 2] = [
    "(((A:0.1,B:0.1)AB:0.1,(C:0.1)R:0.1)X:0.1,D:0.2)root1;",
    "(((C:0.1)R:0.05,D:0.1)Y:0.1,(A:0.1,B:0.1)AB:0.2)root2;",
  ];

  fn merge(trees: &[&str]) -> Result<MergedGraph, Report> {
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
    merge_trees_by_name(&named)
  }

  fn partitions(
    merged: &MergedGraph,
    alignments: Vec<Vec<AlignmentRecord>>,
    segments: bool,
  ) -> Result<PartitionSet, Report> {
    let sources = alignments
      .into_iter()
      .enumerate()
      .map(|(index, records)| {
        complete_sequence_source(
          &merged.graph,
          &merged.names,
          records,
          &NUC_ALPHABET,
          &format!("aln_{index}.fasta"),
        )
      })
      .collect::<Result<Vec<_>, _>>()?;
    let split = SequenceInput::whole_alignments(sources, &GtrModel::JC69(JC69Params::default()))?.split()?;
    let scopes = segments.then(|| Arc::new(merged.scopes.clone()));
    let inheritances = segment_inheritances(
      scopes.as_ref(),
      split.alignments.len(),
      split.partitions.iter().map(|partition| &partition.data.layout),
    )?;
    let lengths = merged
      .branch_lengths
      .iter()
      .map(|(edge, length)| (edge, length.unwrap_or(0.0)))
      .collect();
    build_inherited_partition_set(
      Representation::Dense,
      &merged.graph,
      split.partitions,
      inheritances,
      split.models,
      &lengths,
      &NoopProgress,
    )
  }

  fn constant_alignment() -> Result<Vec<AlignmentRecord>, Report> {
    ["A", "B", "C", "D"]
      .iter()
      .map(|name| {
        Ok(AlignmentRecord {
          name: (*name).to_owned(),
          seq: Seq::try_from_str("ACGTACGTAC")?,
        })
      })
      .collect()
  }

  fn edge(parent: &str, child: &str) -> String {
    format!("{parent}->{child}")
  }

  fn edge_name(merged: &MergedGraph, key: GraphEdgeKey) -> Result<String, Report> {
    let (parent, child) = merged.graph.edge_endpoints(key)?;
    Ok(format!("{}->{}", &merged.names[parent], &merged.names[child]))
  }

  #[test]
  fn test_segment_graph_reassortment_node_has_a_timeless_child_edge() -> Result<(), Report> {
    let merged = merge(&RECOMBINANT)?;
    let set = partitions(&merged, vec![constant_alignment()?, constant_alignment()?], true)?;
    let timeless = set
      .timeless_nodes(&merged.graph)?
      .into_iter()
      .map(|node| {
        Ok((
          edge_name(&merged, node.held)?,
          node
            .others
            .into_iter()
            .map(|other| edge_name(&merged, other))
            .collect::<Result<BTreeSet<_>, Report>>()?,
        ))
      })
      .collect::<Result<Vec<_>, Report>>()?;
    assert_eq!(
      vec![(edge("R", "C"), BTreeSet::from([edge("X", "R"), edge("Y", "R")]))],
      timeless
    );
    Ok(())
  }

  #[test]
  fn test_segment_graph_series_groups_separate_edges_in_series_through_a_pass_through_node() -> Result<(), Report> {
    let merged = merge(&RECOMBINANT)?;
    let set = partitions(&merged, vec![constant_alignment()?, constant_alignment()?], true)?;
    let groups = series_groups(&merged.graph, &set)?
      .into_iter()
      .map(|group| {
        group
          .into_iter()
          .map(|key| edge_name(&merged, key))
          .collect::<Result<BTreeSet<_>, Report>>()
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let group_of = |name: &str| groups.iter().position(|group| group.contains(name));
    let apart = |a: &str, b: &str| group_of(a).is_some() && group_of(a) != group_of(b);
    assert_eq!(
      (merged.graph.get_edges().count(), true, true),
      (
        groups.iter().map(BTreeSet::len).sum::<usize>(),
        apart("X->R", "R->C") && apart("Y->R", "R->C"),
        groups.iter().all(|group| !group.is_empty()),
      )
    );
    Ok(())
  }

  #[test]
  fn test_segment_graph_series_groups_separate_the_two_edges_below_an_inner_segment_root() -> Result<(), Report> {
    let merged = merge(&["((A:0.1,B:0.1)S:0.1,(C:0.1,D:0.1)T:0.1)root;", "(A:0.1,B:0.1)S;"])?;
    let set = partitions(&merged, vec![constant_alignment()?, constant_alignment()?], true)?;
    let groups = series_groups(&merged.graph, &set)?
      .into_iter()
      .map(|group| {
        group
          .into_iter()
          .map(|key| edge_name(&merged, key))
          .collect::<Result<BTreeSet<_>, Report>>()
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let group_of = |name: &str| groups.iter().position(|group| group.contains(name));
    assert_ne!(group_of("S->A"), group_of("S->B"));
    Ok(())
  }

  #[test]
  fn test_tree_series_groups_are_one_group_of_every_edge() -> Result<(), Report> {
    let merged = merge(&[RECOMBINANT[0]])?;
    let set = partitions(&merged, vec![constant_alignment()?], false)?;
    let groups = series_groups(&merged.graph, &set)?;
    assert_eq!(
      vec![merged.graph.get_edges().count()],
      groups.iter().map(BTreeSet::len).collect_vec()
    );
    Ok(())
  }

  fn toy(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .join("../../data/toy/reassortment")
      .join(file)
  }

  #[test]
  fn test_segment_graph_optimize_loop_does_not_lower_the_likelihood() -> Result<(), Report> {
    let segments = ["seg1", "seg2", "seg3"];
    let trees = segments
      .iter()
      .map(|segment| fs::read_to_string(toy(&format!("tree_{segment}.nwk"))))
      .collect::<Result<Vec<_>, _>>()?;
    let merged = merge(&trees.iter().map(String::as_str).collect_vec())?;
    let alignments = segments
      .iter()
      .map(|segment| {
        Ok(
          fasta_read_file(toy(&format!("aln_{segment}.fasta")), &*NUC_ALPHABET)?
            .into_iter()
            .map(AlignmentRecord::from)
            .collect(),
        )
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let set = partitions(&merged, alignments, true)?;
    let (set, _) = set.marginal_update(&merged.graph, &branch_lengths_or_zero(&merged.branch_lengths))?;
    let MergedGraph {
      mut graph,
      names,
      branch_lengths,
      ..
    } = merged;
    let result = run_optimize_loop(
      &mut graph,
      set,
      &OptimizeLoopSettings {
        max_iter: 20,
        dp: 1e-6,
        damping: 0.0,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::NONE,
        partition_rates: PartitionRates::Equal,
      },
      branch_lengths,
      names,
      &NoopProgress,
    )?;
    let steps = result
      .lh_history
      .iter()
      .tuple_windows()
      .map(|(before, after)| after.value() - before.value())
      .collect_vec();
    assert!(steps.len() > 1 && steps.iter().all(|&step| step > -1e-6), "{steps:?}");
    Ok(())
  }
}
