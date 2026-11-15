#[cfg(test)]
mod tests {
  use crate::ancestral::params::{AncestralParams, MethodAncestral};
  use crate::ancestral::pipeline;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::partition::marginal::sample::SampleMode;
  use crate::partition::split::SequenceInput;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::{AncestralInput, EdgeSeqInput};
  use crate::seq::mutation::Mutation;
  use crate::seq::sink::{SeqItem, SeqSink, SeqStream, SeqTrack};
  use crate::test_utils::{NUC_ALPHABET, complete_sequence_source};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_graph::graph::Graph;
  use treetime_graph::key_map::EdgeMap;
  use treetime_graph::node_names::NodeNames;
  use treetime_graph::scope::EdgeScopes;
  use treetime_io::fasta::fasta_read_file;
  use treetime_io::nwk::TREE_DIALECT_DEFAULT;
  use treetime_io::tree::{segment_trees_read_files, tree_read_file};
  use treetime_primitives::{AlignmentRecord, LogLh, Seq};
  use treetime_utils::{make_report, pretty_assert_abs_diff_eq};

  const SEGMENTS: [&str; 3] = ["seg1", "seg2", "seg3"];

  #[derive(Clone, Copy, Debug)]
  enum Case {
    Parsimony,
    MarginalDenseJc69,
    MarginalDenseInfer,
  }

  impl Case {
    fn method(self) -> (MethodAncestral, Option<bool>) {
      match self {
        Self::Parsimony => (MethodAncestral::Parsimony, None),
        Self::MarginalDenseJc69 | Self::MarginalDenseInfer => (MethodAncestral::Marginal, Some(true)),
      }
    }

    fn model(self) -> GtrModel {
      match self {
        Self::Parsimony | Self::MarginalDenseJc69 => GtrModel::JC69(JC69Params::default()),
        Self::MarginalDenseInfer => GtrModel::Infer,
      }
    }
  }

  #[rstest]
  #[case::dense_jc69(Case::MarginalDenseJc69)]
  #[case::dense_infer(Case::MarginalDenseInfer)]
  #[trace]
  fn test_ancestral_marginal_on_segment_graph_gives_each_segment_the_result_of_its_own_tree(
    #[case] case: Case,
  ) -> Result<(), Report> {
    let merged = run_segment_graph(case)?;
    for (index, segment) in SEGMENTS.iter().enumerate() {
      let separate = run_segment_tree(segment, case)?;
      assert_eq!(
        (&separate.sequences, &separate.mutations, &separate.root_sequence),
        (
          &merged[index].sequences,
          &merged[index].mutations,
          &merged[index].root_sequence
        ),
        "segment {segment}"
      );
      pretty_assert_abs_diff_eq!(separate.log_lh, merged[index].log_lh, epsilon = 1e-9);
    }
    Ok(())
  }

  #[test]
  fn test_ancestral_parsimony_on_segment_graph_gives_each_segment_the_result_of_its_own_tree() -> Result<(), Report> {
    let merged = run_segment_graph(Case::Parsimony)?;
    assert_eq!(
      vec![(20, true), (20, true), (21, true)],
      merged
        .iter()
        .map(|result| (result.sequences.len(), result.mutations.values().any(|m| !m.is_empty())))
        .collect::<Vec<_>>()
    );
    for (index, segment) in SEGMENTS.iter().enumerate() {
      let separate = run_segment_tree(segment, Case::Parsimony)?;
      assert_eq!(
        (&separate.sequences, &separate.mutations, &separate.root_sequence),
        (
          &merged[index].sequences,
          &merged[index].mutations,
          &merged[index].root_sequence
        ),
        "segment {segment}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_ancestral_parsimony_on_segment_graph_reconstructs_only_the_nodes_of_each_segment() -> Result<(), Report> {
    let merged = run_segment_graph(Case::Parsimony)?;
    let internal = |index: usize| -> BTreeSet<&str> {
      merged[index]
        .sequences
        .keys()
        .map(String::as_str)
        .filter(|name| name.len() > 1 || "RKMNPQUVXYZ".contains(*name))
        .collect()
    };
    assert_eq!(
      [
        BTreeSet::from(["K", "M", "P", "Q", "R", "U", "V", "X", "Y", "Z"]),
        BTreeSet::from(["K", "M", "P", "Q", "R", "U", "V", "X", "Y", "Z"]),
        BTreeSet::from(["K", "M", "N", "P", "Q", "R3", "U", "V", "X", "Y", "Z"]),
      ],
      [internal(0), internal(1), internal(2)]
    );
    Ok(())
  }

  #[test]
  fn test_ancestral_parsimony_on_segment_graph_without_scopes_is_error() -> Result<(), Report> {
    let paths = SEGMENTS.map(|segment| toy_path(&format!("tree_{segment}.nwk")));
    let merged = segment_trees_read_files(&paths, TREE_DIALECT_DEFAULT)?;
    let alignments = SEGMENTS
      .iter()
      .map(|segment| alignment(segment))
      .collect::<Result<_, _>>()?;
    let result = run(
      merged.graph,
      &merged.names,
      merged.branch_lengths,
      None,
      alignments,
      Case::Parsimony,
    );
    let message = result.err().map(|error| error.to_string()).unwrap_or_default();
    assert!(
      message.contains("parents that pass on all sites, but every site has one parent"),
      "{message}"
    );
    Ok(())
  }

  struct SegmentResult {
    sequences: BTreeMap<String, Seq>,
    mutations: BTreeMap<(String, String), Vec<Mutation>>,
    root_sequence: Seq,
    log_lh: f64,
  }

  fn toy_path(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .join("../../data/toy/reassortment")
      .join(file)
  }

  fn alignment(segment: &str) -> Result<Vec<AlignmentRecord>, Report> {
    let records = fasta_read_file(toy_path(&format!("aln_{segment}.fasta")), &*NUC_ALPHABET)?;
    Ok(records.into_iter().map(AlignmentRecord::from).collect())
  }

  fn run_segment_graph(case: Case) -> Result<Vec<SegmentResult>, Report> {
    let paths = SEGMENTS.map(|segment| toy_path(&format!("tree_{segment}.nwk")));
    let merged = segment_trees_read_files(&paths, TREE_DIALECT_DEFAULT)?;
    let alignments = SEGMENTS
      .iter()
      .map(|segment| alignment(segment))
      .collect::<Result<_, _>>()?;
    run(
      merged.graph,
      &merged.names,
      merged.branch_lengths,
      Some(Arc::new(merged.scopes)),
      alignments,
      case,
    )
  }

  fn run_segment_tree(segment: &str, case: Case) -> Result<SegmentResult, Report> {
    let tree = tree_read_file(toy_path(&format!("tree_{segment}.nwk")), TREE_DIALECT_DEFAULT)?;
    let mut results = run(
      tree.graph,
      &tree.names,
      tree.branch_lengths,
      None,
      vec![alignment(segment)?],
      case,
    )?;
    results
      .pop()
      .ok_or_else(|| make_report!("The run of segment {segment} returned no alignment"))
  }

  fn run(
    graph: Graph,
    names: &NodeNames,
    branch_lengths: EdgeMap<Option<f64>>,
    scopes: Option<Arc<EdgeScopes>>,
    alignments: Vec<Vec<AlignmentRecord>>,
    case: Case,
  ) -> Result<Vec<SegmentResult>, Report> {
    let sources = alignments
      .into_iter()
      .enumerate()
      .map(|(index, records)| {
        complete_sequence_source(&graph, names, records, &NUC_ALPHABET, &format!("aln_{index}.fasta"))
      })
      .collect::<Result<Vec<_>, _>>()?;
    let input = AncestralInput {
      edges: branch_lengths
        .into_iter()
        .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
        .collect(),
      graph,
      sequences: SequenceInput::whole_alignments(sources, &case.model())?,
      scopes,
    };
    let (method, dense) = case.method();
    let params = AncestralParams {
      method,
      dense,
      include_leaves: true,
      report_ambiguous: false,
      impute_missing_data: false,
      seed: 0,
      sample_from_profile: SampleMode::Argmax,
    };
    let mut sink = TrackSink::default();
    let output = pipeline::run(
      &params,
      input,
      Some(SeqStream::new(&mut sink, names)),
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )
    .map_err(|error| error.into_report())?;
    let log_lh = output
      .partitions
      .iter()
      .map(|partition| partition.log_lh.map_or(0.0, LogLh::value))
      .collect::<Vec<_>>();
    output
      .alignments
      .into_iter()
      .zip(log_lh)
      .map(|(result, log_lh)| {
        let sequences = result
          .sequences
          .ok_or_else(|| make_report!("Alignment {} has no sequences", result.number))?;
        let mutations = sequences
          .edge_mutations
          .into_iter()
          .map(|(edge_key, mutations)| {
            let (parent, child) = output.graph.edge_endpoints(edge_key)?;
            Ok(((names[parent].to_owned(), names[child].to_owned()), mutations))
          })
          .collect::<Result<_, Report>>()?;
        Ok(SegmentResult {
          sequences: sink.by_track.remove(&result.number).unwrap_or_default(),
          mutations,
          root_sequence: sequences.root_sequence,
          log_lh,
        })
      })
      .collect()
  }

  #[derive(Default)]
  struct TrackSink {
    by_track: BTreeMap<usize, BTreeMap<String, Seq>>,
  }

  impl SeqSink for TrackSink {
    fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
      if let SeqTrack::Nuc(alignment) = item.track {
        self
          .by_track
          .entry(alignment)
          .or_default()
          .insert(item.name.to_owned(), item.seq.clone());
      }
      Ok(())
    }
  }
}
