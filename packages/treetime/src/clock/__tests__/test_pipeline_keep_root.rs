#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::pipeline::{self, ClockInput, ClockParams};
  use crate::progress::NoopProgress;
  use crate::test_utils::dates_by_node;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::date::DateConstraint;
  use treetime_utils::o;

  #[test]
  fn test_pipeline_keep_root_keeps_an_undated_stem_root() -> Result<(), Report> {
    let parsed = nwk_read(b"(((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)R:0.001)STEM;".as_slice())?;
    let names = parsed.names.clone();
    let dates = [
      (o!("A"), 2013.0),
      (o!("B"), 2022.0),
      (o!("C"), 2017.0),
      (o!("D"), 2005.0),
    ]
    .into_iter()
    .map(|(name, date)| (name, Some(DateConstraint::exact(date))));
    let input = ClockInput {
      dates: dates_by_node(dates, &parsed.graph, &names),
      graph: parsed.graph,
      names,
      branch_lengths: parsed.branch_lengths,
    };
    let params = ClockParams {
      clock_filter: 0.0,
      allow_negative_rate: true,
      reroot_spec: None,
      prune_filtered_tips: None,
    };

    let output = pipeline::run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let root_key = output.graph.get_exactly_one_root()?.key();
    assert_eq!(("STEM", 8), (&output.names[root_key], output.graph.get_nodes().count()));
    Ok(())
  }
}
