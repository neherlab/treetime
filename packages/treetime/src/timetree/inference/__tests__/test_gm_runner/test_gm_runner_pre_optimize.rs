#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::cancel::NoopCancel;
  use crate::progress::NoopProgress;
  use crate::test_utils::{dates_by_node, leaf_seq_inputs, marginal_timetree_params};
  use crate::timetree::params::TimetreeParams;
  use crate::timetree::pipeline::{self, TimetreeInput};
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  fn test_gm_runner_pre_optimize_changes_branch_lengths_and_dates_every_node(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let nwk_parsed = nwk_read(case.rerooted_tree_nwk().as_bytes())?;
    let names = nwk_parsed.names();
    let dates = dates_by_node(load_dates_for_dataset(dataset)?, &nwk_parsed.graph, &names);
    let input_branch_lengths = nwk_parsed.branch_lengths.clone();
    let aln: Vec<AlignmentRecord> = load_alignment_for_dataset(dataset)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let aln_nodes = leaf_seq_inputs(&nwk_parsed.graph, &names, aln);
    let input = TimetreeInput {
      graph: nwk_parsed.graph,
      names,
      alphabet: ALPHABET.clone(),
      sequences: Some(aln_nodes),
      dates: Some(dates),
      branch_lengths: nwk_parsed.branch_lengths,
    };
    let params = TimetreeParams {
      clock_rate: Some(case.clock_rate()),
      ..marginal_timetree_params()
    };

    let output = pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let n_changed = input_branch_lengths
      .iter()
      .filter(|(key, before)| (before.unwrap_or(0.0) - output.branch_lengths[key].unwrap_or(0.0)).abs() > 1e-10)
      .count();
    assert!(
      n_changed > 0,
      "Expected at least one branch length to change after ML optimization, but none did"
    );

    let times: BTreeMap<GraphNodeKey, f64> = output
      .node_dates
      .iter()
      .filter_map(|(&key, &time)| time.map(|time| (key, time)))
      .collect();
    assert_eq!(output.graph.num_nodes(), times.len(), "every node must be dated");
    let non_finite: BTreeMap<&GraphNodeKey, &f64> = times.iter().filter(|(_, time)| !time.is_finite()).collect();
    assert_eq!(BTreeMap::new(), non_finite, "every node time must be finite");
    Ok(())
  }
}
