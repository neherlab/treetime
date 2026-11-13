#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{HKY85_PARAMS, hky85_model, project_root};
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw, homoplasy_ancestral_params};
  use crate::commands::shared::sequence_inputs::{TreeAndSources, ancestral_input, read_tree_and_sources};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use tempfile::tempdir;
  use treetime::alphabet::alphabet::AlphabetName;
  use treetime::ancestral;
  use treetime::cancel::NoopCancel;
  use treetime::partition::split::SequenceInput;
  use treetime::progress::NoopProgress;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_homoplasy_model_params_reach_the_model() -> Result<(), Report> {
    let out = tempdir()?;
    let data = project_root().join("data/flu/h3n2/20");
    let tree = format!("--tree={}", data.join("tree.nwk").display());
    let alignment = format!("--alignment={}", data.join("aln.fasta.xz").display());
    let output = format!("--output-homoplasy-stats={}", out.path().join("stats.json").display());
    let mut argv = vec![
      "homoplasy",
      &tree,
      &alignment,
      "--method-anc=marginal",
      "--model=hky85",
      "--model-params",
    ];
    argv.extend(HKY85_PARAMS);
    argv.push(&output);
    let args = TreetimeHomoplasyArgs::try_from(TreetimeHomoplasyArgsRaw::try_parse_from(argv)?)?;
    let TreeAndSources {
      tree,
      sources,
      partitions,
    } = read_tree_and_sources(
      &args.tree,
      args.newick_dialect,
      &args.sequences,
      &args.model,
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;
    let input = ancestral_input(
      tree.graph,
      tree.branch_lengths,
      SequenceInput {
        alignments: sources.into_iter().map(|alignment| alignment.source).collect(),
        partitions,
      },
    );

    let output = ancestral::pipeline::run(
      &homoplasy_ancestral_params(&args, 0),
      input,
      None,
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )
    .map_err(|err| err.into_report())?;

    let [partition] = output.partitions.as_slice() else {
      panic!("one alignment gives one partition");
    };
    let gtr = partition.gtr.as_ref().expect("the partition has a substitution model");
    let model = hky85_model();
    let expected = model.build(AlphabetName::Nuc)?;
    assert_eq!(model.name(), partition.model_name);
    pretty_assert_ulps_eq!(expected.pi().clone(), gtr.pi().clone(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected.W().clone(), gtr.W().clone(), max_ulps = 4);
    Ok(())
  }
}
