#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
  use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};

  use crate::cancel::NoopCancel;
  use crate::progress::NoopProgress;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use std::path::Path;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::AlignmentRecord;

  #[test]
  fn test_optimize_pipeline_infer_gtr_mu_normalized() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(Path::parent)
      .expect("workspace root");

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");

    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let sequences: Vec<AlignmentRecord> = read_many_fasta_path(&[aln_path.to_str().expect("utf-8 path")], &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let params = OptimizeParams {
      model: GtrModelName::Infer,
      dense: Some(false),
      max_iter: 2,
      dp: 0.1,
      damping: 0.0,
      opt_method: BranchOptMethod::default(),
      initial_guess: InitialGuessMode::default(),
      no_indels: false,
      reroot_spec: None,
      topology_ops: TopologyOps::default(),
    };
    let input = OptimizeInput {
      graph,
      alphabet,
      sequences,
      branch_lengths,
    };

    let output = run(&params, input, &names, &NoopCancel, &NoopProgress)?;

    assert_ulps_eq!(output.gtr.mu, 1.0, max_ulps = 4);

    let partition_mu = output.sparse_partitions[0].gtr.mu;
    assert_ulps_eq!(output.gtr.mu, partition_mu, max_ulps = 4);

    Ok(())
  }
}
