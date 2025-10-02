#![allow(clippy::todo)]
use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::timetree::clock_model::{infer_clock_model, update_clock_model};
use crate::commands::timetree::date_constraints::load_date_constraints;
use crate::commands::timetree::timetree_args::{BranchLengthMode, TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::timetree_unified::run_timetree;
use crate::gtr::get_gtr::{JC69Params, get_gtr, jc69};
use crate::io::fasta::read_many_fasta;
use crate::io::fs::ensure_dir;
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition_timetree::PartitionTimetreeOps;
use crate::representation::partition_timetree_dense::PartitionTimetreeDense;
use crate::representation::partition_timetree_sparse::PartitionTimetreeSparse;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;

enum PartitionVariant {
  Sparse(Vec<Arc<RwLock<PartitionMarginalSparse>>>),
  Dense(Vec<Arc<RwLock<PartitionMarginalDense>>>),
}

pub fn run_timetree_estimation(args: &TreetimeTimetreeArgs) -> Result<(), Report> {
  // Validate and prepare output directory
  ensure_dir(&args.outdir).wrap_err("Failed to create output directory")?;

  // ============================================================================
  // LOAD INPUT DATA
  // ============================================================================

  // Load phylogenetic tree (using GraphAncestral for compatibility with ancestral reconstruction)
  let graph: GraphAncestral = if let Some(tree_path) = &args.tree {
    nwk_read_file(tree_path).wrap_err("Failed to load tree from file")?
  } else {
    todo!("Tree inference from alignment using FastTree/IQ-TREE")
  };

  // Load alignment if provided (needed for branch length optimization)
  let aln = if !args.input_fastas.is_empty() {
    let alphabet_name = args.alphabet.unwrap_or_default();
    let dense = args.dense.unwrap_or_else(infer_dense);
    let treat_gap_as_unknown = dense;
    let alphabet = Alphabet::new(alphabet_name, treat_gap_as_unknown)?;

    Some((read_many_fasta(&args.input_fastas, &alphabet)?, alphabet))
  } else if args.input_fastas.is_empty() && args.branch_length_mode != BranchLengthMode::Input {
    return Err(eyre::eyre!(
      "Alignment required when branch_length_mode is not 'input'. \
       Provide FASTA files or use --branch-length-mode=input"
    ));
  } else {
    None
  };

  // Parse date constraints from metadata files
  let constraints = load_date_constraints(args, &graph).wrap_err("Failed to load date constraints")?;

  // ============================================================================
  // INITIALIZATION PHASE
  // ============================================================================

  // Estimate initial molecular clock rate via root-to-tip regression
  let mut clock_model = infer_clock_model(args, &graph, &constraints).wrap_err("Failed to infer clock model")?;

  // Create marginal partitions for sequence data if alignment is provided
  // These are used for ancestral reconstruction to get branch length distributions

  let partitions_marginal = if let Some((aln_data, alphabet)) = &aln {
    let dense = args.dense.unwrap_or_else(infer_dense);
    let length = get_common_length(aln_data)?;

    if !dense {
      // Sparse representation: store only mutations relative to reference
      #[allow(clippy::iter_on_single_items)]
      let partitions_sparse = [PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: alphabet.clone(),
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
      }]
      .into_iter()
      .map(|p| Arc::new(RwLock::new(p)))
      .collect_vec();

      if !partitions_sparse.is_empty() {
        compress_sequences(&graph, &partitions_sparse, aln_data)?;

        for partition in &partitions_sparse {
          let gtr = get_gtr(&args.gtr, partition, &graph)?;
          partition.write_arc().gtr = gtr;
        }
      }

      Some(PartitionVariant::Sparse(partitions_sparse))
    } else {
      // Dense representation: store full probability distributions
      #[allow(clippy::iter_on_single_items)]
      let partitions_dense = [PartitionMarginalDense {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: alphabet.clone(),
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
      }]
      .into_iter()
      .map(|p| Arc::new(RwLock::new(p)))
      .collect_vec();

      if !partitions_dense.is_empty() {
        for partition in &partitions_dense {
          let gtr = get_gtr(&args.gtr, partition, &graph)?;
          partition.write_arc().gtr = gtr;
        }
      }

      Some(PartitionVariant::Dense(partitions_dense))
    }
  } else {
    None
  };

  // Initial ancestral sequence reconstruction if needed
  // This establishes branch length distributions for timetree inference
  if let Some(partitions) = &partitions_marginal {
    match args.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Using input branch lengths for timetree inference");
      },
      BranchLengthMode::Marginal | BranchLengthMode::Auto => {
        info!("Running marginal ancestral reconstruction for branch length distributions");
        match partitions {
          PartitionVariant::Sparse(parts) => {
            run_marginal(&graph, parts, None)?;
          },
          PartitionVariant::Dense(parts) => {
            run_marginal(&graph, parts, aln.as_ref().map(|(a, _)| a.as_slice()))?;
          },
        }
      },
      BranchLengthMode::Joint => {
        info!("Running joint ancestral reconstruction for branch length distributions");
        todo!("Joint ancestral reconstruction not yet implemented")
      },
    }
  }

  // Create timetree partitions for temporal inference
  // These hold date constraints, time distributions, and clock lengths
  // Using trait objects to allow mixing sparse and dense partitions
  let partitions_timetree: Vec<Arc<RwLock<dyn PartitionTimetreeOps>>> = {
    let dense = args.dense.unwrap_or_else(infer_dense);
    let sequence_length = if let Some((aln_data, _)) = &aln {
      Some(get_common_length(aln_data)?)
    } else {
      args.sequence_length
    };

    #[allow(clippy::iter_on_single_items, trivial_casts)]
    if !dense {
      [Arc::new(RwLock::new(PartitionTimetreeSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: aln.as_ref().map_or_else(
          || Alphabet::new(args.alphabet.unwrap_or_default(), false).unwrap(),
          |(_, a)| a.clone(),
        ),
        length: sequence_length.unwrap_or(0),
        sequence_length,
        nodes: btreemap! {},
        edges: btreemap! {},
      })) as Arc<RwLock<dyn PartitionTimetreeOps>>]
    } else {
      [Arc::new(RwLock::new(PartitionTimetreeDense {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: aln.as_ref().map_or_else(
          || Alphabet::new(args.alphabet.unwrap_or_default(), false).unwrap(),
          |(_, a)| a.clone(),
        ),
        length: sequence_length.unwrap_or(0),
        sequence_length,
        nodes: btreemap! {},
        edges: btreemap! {},
      })) as Arc<RwLock<dyn PartitionTimetreeOps>>]
    }
    .into_iter()
    .map(|partition| -> Result<_, Report> {
      partition
        .write_arc()
        .attach_date_constraints(&graph, &constraints.per_node)?;
      Ok(partition)
    })
    .try_collect()?
  };

  // Optional rerooting based on temporal signal
  // Python v0 reroots BEFORE initial make_time_tree() if root is set or n_iqd (clock filter)
  if !args.keep_root {
    {
      info!("Rerooting tree to optimize temporal signal");
      todo!("Rerooting algorithms not yet implemented")
    }
  }

  // Second ancestral reconstruction pass after rerooting (if rerooted)
  // TODO: Re-run ancestral reconstruction if tree topology changed

  // ============================================================================
  // INITIAL TIMETREE ESTIMATION (before iteration loop)
  // ============================================================================
  info!("### TreeTime: INITIAL ROUND");

  // Run timetree inference using trait-based unified runner (analogous to run_marginal)
  // This replaces the old manual backward/forward pass calls
  run_timetree(&graph, &partitions_timetree)?;

  // TODO: Build initial branch length distributions from input tree or ancestral reconstruction
  // let branch_mode = match args.branch_length_mode {
  //   BranchLengthMode::Auto => BranchDistributionContext::InputPrior,
  //   BranchLengthMode::Input => BranchDistributionContext::InputPrior,
  //   BranchLengthMode::Joint => BranchDistributionContext::Joint,
  //   BranchLengthMode::Marginal => BranchDistributionContext::Marginal,
  // };
  // let mut branch_distributions =
  //   build_branch_distributions(&graph, branch_mode).wrap_err("Failed to build initial branch distributions")?;

  // TODO: Old manual backward/forward pass calls - replaced by run_timetree above
  // let backward_state =
  //   run_backward_pass(&graph, &constraints, &mut branch_distributions).wrap_err("Initial backward pass failed")?;
  // let _forward_state = run_forward_pass(&graph, &constraints, &mut branch_distributions, &backward_state)
  //   .wrap_err("Initial forward pass failed")?;

  // TODO: Extract node time estimates (numdate) from posteriors and store in graph

  // TODO: If root was specified AND max_iter > 0:
  //   - Reroot again using least-squares on inferred node times
  //   - Re-run ancestral reconstruction if branch_length_mode != input
  //   - Re-run make_time_tree (backward + forward pass)
  // Python does this "second rerooting" after initial timetree at line 276

  // ============================================================================
  // ITERATION LOOP: Refine timetree, resolve polytomies, update clock model
  // ============================================================================
  let max_iter = args.max_iter.unwrap_or(0);
  for iteration in 1..=max_iter {
    info!("### Timetree iteration {iteration}/{max_iter}");

    // TODO: Add coalescent model if requested (args.coalescent)
    // Python adds merger model as additional factor in node probability calculation
    // For 'skyline', use 'const' until last iteration

    // TODO: Apply relaxed clock if requested (args.relax)
    // Estimates branch-specific rate variation

    // TODO: Resolve polytomies using temporal information if requested
    // Python's resolve_polytomies() uses time constraints to break ties
    let n_resolved = 0; // Placeholder

    // If polytomies were resolved, need to re-optimize branch lengths
    let mut need_new_time_tree = n_resolved > 0;
    if n_resolved > 0 {
      // TODO: prepare_tree() - reset internal state after tree structure changes
      // TODO: optimize_tree(max_iter=0) - re-run ancestral reconstruction
      need_new_time_tree = true;
    }

    // TODO: assign_gamma() if custom branch length scaling function provided
    // Modifies branch lengths by local clock rate multipliers

    // Conditional reconstruction order (Python lines 328-337):
    // If tree was modified (polytomies/coalescent/relaxed clock):
    //   1. make_time_tree() - update node times based on tree changes
    //   2. infer_ancestral_sequences() - re-infer sequences given new times
    // Otherwise (just iterating to convergence):
    //   1. infer_ancestral_sequences() - update sequences first
    //   2. make_time_tree() - update times based on new branch lengths

    let ndiff = if need_new_time_tree {
      // Tree structure or rates changed - update times first
      run_timetree(&graph, &partitions_timetree)
        .wrap_err_with(|| format!("Timetree inference failed (iteration {iteration})"))?;

      // TODO: Re-infer ancestral sequences and return number of changes
      if let Some(partitions) = &partitions_marginal {
        match partitions {
          PartitionVariant::Sparse(parts) => {
            run_marginal(&graph, parts, None)?;
          },
          PartitionVariant::Dense(parts) => {
            run_marginal(&graph, parts, aln.as_ref().map(|(a, _)| a.as_slice()))?;
          },
        }
      }

      0 // Placeholder
    } else {
      // No tree changes - update sequences first
      let sequence_changes = if let Some(partitions) = &partitions_marginal {
        match partitions {
          PartitionVariant::Sparse(parts) => {
            run_marginal(&graph, parts, None)?;
          },
          PartitionVariant::Dense(parts) => {
            run_marginal(&graph, parts, aln.as_ref().map(|(a, _)| a.as_slice()))?;
          },
        }
        0 // TODO: track actual sequence changes
      } else {
        0
      };

      // Then update time tree
      run_timetree(&graph, &partitions_timetree)
        .wrap_err_with(|| format!("Timetree inference failed (iteration {iteration})"))?;

      sequence_changes
    };

    // TODO: Update tracelog with iteration stats (niter, ndiff, n_resolved, LH values)

    // Update clock model based on inferred node times
    let new_clock_model = update_clock_model(&graph, &constraints, &clock_model)
      .wrap_err_with(|| format!("Failed to update clock model (iteration {iteration})"))?;

    // Check convergence: no sequence changes AND no polytomies resolved
    // Python line 354: if ndiff==0 and n_resolved==0 and Tc!='skyline': break
    let converged = ndiff == 0 && n_resolved == 0;

    clock_model = new_clock_model;

    if converged {
      info!("### TreeTime: CONVERGED at iteration {iteration}");
      break;
    }
  }

  // ============================================================================
  // POST-ITERATION: Rate variation and final confidence estimation
  // ============================================================================

  // TODO: Rate variation analysis if args.vary_rate is set
  // Python's calc_rate_susceptibility() re-runs timetree with rate ± std deviation
  // Useful for understanding sensitivity to clock rate uncertainty

  // TODO: Rate variation analysis if args.vary_rate is set
  // Python's calc_rate_susceptibility() re-runs timetree with rate ± std deviation
  // Useful for understanding sensitivity to clock rate uncertainty

  // FINAL CONFIDENCE ESTIMATION: Run timetree with marginal mode for confidence intervals
  if args.time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: timetree with marginal reconstruction for confidence intervals");
    run_timetree(&graph, &partitions_timetree).wrap_err("Final timetree inference failed")?;

    // TODO: Extract confidence intervals (HPD or quantiles) from marginal posteriors in partitions
    // TODO: Update tracelog with final marginal reconstruction stats
  }

  // TODO: Mark and report bad branches (outliers that violate molecular clock)
  // Python identifies tips where inferred date differs significantly from input constraint
  // These are marked with bad_branch=True and logged as warnings

  // TODO: Write posterior time estimates to node comments
  // TODO: Extract and write confidence intervals (HPD or quantiles)

  // Write output files
  let out_base = args.outdir.join("timetree");
  nwk_write_file(out_base.with_extension("nwk"), &graph, &NwkWriteOptions::default())
    .wrap_err("Failed to write Newick output")?;
  nex_write_file(out_base.with_extension("nexus"), &graph, &NexWriteOptions::default())
    .wrap_err("Failed to write Nexus output")?;

  // TODO: Write temporal data (node times, distributions, etc.) from partitions_timetree
  // TODO: Write confidence intervals to CSV
  // TODO: Write clock model parameters to JSON
  // TODO: Plot root-to-tip regression if requested (args.plot_rtt)
  // TODO: Plot time tree if requested (args.plot_tree)

  Ok(())
}
