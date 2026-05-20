use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr, write_gtr_json};
use crate::optimize::run_loop::{apply_initial_guess_mode, normalize_partition_rates};
use crate::optimize::run_loop::{collect_optimize_partitions, run_optimize_loop};
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::traits::MutationCommentProvider;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::write_graph_files_with;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::make_error;

#[derive(Clone, Debug, Default)]
pub struct TreetimeOptimizeParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

pub fn run_optimize(args: &TreetimeOptimizeArgs) -> Result<(), Report> {
  let TreetimeOptimizeArgs {
    input_fastas,
    tree,
    alphabet,
    model_name,
    dense,
    outdir,
    max_iter,
    dp,
    damping,
    branch_length_initial_guess,
    opt_method,
    no_indels,
  } = args;

  if !(0.0..1.0).contains(damping) {
    return make_error!("--damping must be in [0.0, 1.0), got {damping}");
  }

  let dense = dense.unwrap_or_else(infer_dense);
  let alphabet = Alphabet::new(alphabet.unwrap_or_default())?;
  let aln = read_many_fasta(input_fastas, &alphabet)?;
  let mut graph: GraphAncestral = nwk_read_file(tree)?;

  let sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if !dense {
    let fitch = create_fitch_partition(&graph, 0, alphabet.clone(), &aln)?;
    let gtr = match *model_name {
      GtrModelName::Infer => infer_gtr_fitch(&fitch, &graph)?,
      _ => get_gtr_by_name(*model_name)?,
    };
    log_gtr(&gtr, *model_name);
    let partition = fitch.into_marginal_sparse(gtr, &graph)?;
    vec![Arc::new(RwLock::new(partition))]
  } else {
    vec![]
  };

  let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = if dense {
    if *model_name == GtrModelName::Infer {
      let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
      let gtr = infer_gtr_fitch(&fitch, &graph)?;
      log_gtr(&gtr, *model_name);
      let partition = fitch.into_marginal_dense(gtr);
      vec![Arc::new(RwLock::new(partition))]
    } else {
      let length = get_common_length(&aln)?;
      let gtr = get_gtr_by_name(*model_name)?;
      log_gtr(&gtr, *model_name);
      let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);
      vec![Arc::new(RwLock::new(partition))]
    }
  } else {
    vec![]
  };

  // Initialize marginal data for whichever partitions are active
  update_marginal(&graph, &sparse_partitions)?;
  if !dense_partitions.is_empty() {
    initialize_marginal(&graph, &dense_partitions, &aln)?;
    update_marginal(&graph, &dense_partitions)?;
  }

  let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

  if *model_name == GtrModelName::Infer {
    normalize_partition_rates(&graph, &sparse_partitions);
    normalize_partition_rates(&graph, &dense_partitions);
  }

  for partition in &sparse_partitions {
    write_gtr_json(&partition.read_arc().gtr, *model_name, outdir, None)?;
  }
  for partition in &dense_partitions {
    write_gtr_json(&partition.read_arc().data.gtr, *model_name, outdir, None)?;
  }

  apply_initial_guess_mode(&graph, &mixed_partitions, *branch_length_initial_guess, *no_indels)?;

  run_optimize_loop(
    &mut graph,
    &sparse_partitions,
    &dense_partitions,
    &mixed_partitions,
    *max_iter,
    *dp,
    *damping,
    *opt_method,
    *no_indels,
  )?;

  // Re-run marginal to populate subs_ml after the loop (the last iteration
  // may have changed topology, clearing ML subs on affected edges).
  update_marginal(&graph, &sparse_partitions)?;
  update_marginal(&graph, &dense_partitions)?;

  if dense {
    let guard = dense_partitions[0].read_arc();
    let provider = MutationCommentProvider::new(&*guard, &graph);
    let providers = CommentProviders::new().with(&provider);
    write_graph_files_with(outdir, "annotated_tree", &graph, &providers)?;
  } else {
    let guard = sparse_partitions[0].read_arc();
    let provider = MutationCommentProvider::new(&*guard, &graph);
    let providers = CommentProviders::new().with(&provider);
    write_graph_files_with(outdir, "annotated_tree", &graph, &providers)?;
  }

  Ok(())
}
