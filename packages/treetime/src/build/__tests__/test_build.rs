#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::build::evaluate::{TreeWithCounts, compare_splits, fitch_edge_counts};
  use crate::build::mat::Mat;
  use crate::build::pipeline::{BuildParams, BuiltMat, SampleOrder, SearchMode, build_mat, mat_to_tree};
  use crate::build::query::Query;
  use crate::build::rebuild::RebuildParams;
  use crate::build::references::{ReferenceSelection, eccentricities, local_branching_index};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{CaseSpec, ModelFamily};
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rand::RngExt;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_simulation::coalescent::SampleBatch;
  use treetime_simulation::observe::{TipObservation, observe_tips};
  use treetime_simulation::tc_schedule::TcSchedule;
  use treetime_utils::pretty_assert_abs_diff_eq;
  use treetime_utils::sync::random::get_random_number_generator;

  const A: AsciiChar = AsciiChar::from_byte_unchecked(b'A');
  const C: AsciiChar = AsciiChar::from_byte_unchecked(b'C');

  struct PerfectPhylogeny {
    graph: Graph,
    names: NodeNames,
    sample_names: Vec<String>,
    seqs: Vec<Seq>,
    variable_mutations: usize,
  }

  fn perfect_phylogeny(n_tips: usize, seed: u64) -> Result<PerfectPhylogeny, Report> {
    let mut rng = get_random_number_generator(seed);
    let mut graph = Graph::new();
    let root = graph.add_node();
    let first = graph.add_node();
    let mut edges = vec![graph.add_edge(root, first)?];
    let mut tips = vec![first];
    while tips.len() < n_tips {
      let edge = edges.swap_remove(rng.random_range(0..edges.len()));
      let (parent, child) = graph.edge_endpoints(edge)?;
      graph.remove_edge(edge)?;
      let node = graph.add_node();
      let tip = graph.add_node();
      edges.extend([
        graph.add_edge(parent, node)?,
        graph.add_edge(node, child)?,
        graph.add_edge(node, tip)?,
      ]);
      tips.push(tip);
    }
    graph.build()?;

    let mut positions: BTreeMap<GraphNodeKey, Vec<usize>> = BTreeMap::new();
    let mut length = 0;
    let mut variable_mutations = 0;
    for edge in graph.get_edges() {
      let count = rng.random_range(1..=3);
      positions.insert(edge.target(), (length..length + count).collect());
      length += count;
      if edge.source() != root {
        variable_mutations += count;
      }
    }
    let seqs = tips
      .iter()
      .map(|&tip| {
        let mut seq: Seq = std::iter::repeat_n(A, length).collect();
        let mut node = tip;
        while let Some((parent, _)) = graph.node_parent(node)? {
          for &pos in &positions[&node] {
            seq[pos] = C;
          }
          node = parent;
        }
        Ok(seq)
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let sample_names = (0..tips.len()).map(|i| format!("t{i}")).collect_vec();
    let input_names = izip_names(&tips, &sample_names);
    let names = NodeNames::from_input(input_names, &graph)?.names;
    Ok(PerfectPhylogeny {
      graph,
      names,
      sample_names,
      seqs,
      variable_mutations,
    })
  }

  fn izip_names(tips: &[GraphNodeKey], names: &[String]) -> BTreeMap<GraphNodeKey, Option<String>> {
    tips.iter().copied().zip(names.iter().cloned().map(Some)).collect()
  }

  fn lbi_references() -> ReferenceSelection {
    ReferenceSelection::Lbi {
      max_references: 4,
      tau: None,
      pseudo_length: 0.5,
    }
  }

  fn walk_params(batch_size: usize, slack: usize, exact: bool) -> BuildParams {
    BuildParams {
      search: SearchMode::Walk {
        references: lbi_references(),
        slack,
        exact,
      },
      batch_size,
      order: SampleOrder::Input,
      spr_rounds: 3,
      spr_slack: Some(1),
      fitch_refresh: false,
      compare_exhaustive: true,
      rebuild: None,
    }
  }

  fn exhaustive_params(batch_size: usize) -> BuildParams {
    BuildParams {
      search: SearchMode::Exhaustive,
      ..walk_params(batch_size, 0, false)
    }
  }

  fn reference_params(batch_size: usize, references: ReferenceSelection) -> BuildParams {
    BuildParams {
      search: SearchMode::FromReferences { references },
      ..walk_params(batch_size, 0, false)
    }
  }

  fn rebuilding_params(batch_size: usize) -> BuildParams {
    BuildParams {
      rebuild: Some(RebuildParams {
        radius: 3,
        max_leaves: 16,
        orders: 3,
      }),
      ..refreshing_params(batch_size)
    }
  }

  fn refreshing_params(batch_size: usize) -> BuildParams {
    BuildParams {
      spr_slack: None,
      fitch_refresh: true,
      ..exhaustive_params(batch_size)
    }
  }

  fn assert_samples_reproduced(built: &BuiltMat, seqs: &[Seq]) {
    let alphabet = Alphabet::default();
    let mat = &built.mat;
    assert_eq!(mat.samples().len(), seqs.len());
    for (&leaf, &sample) in mat.samples() {
      let query = Query::from_seq(&seqs[sample], mat.reference(), &alphabet);
      assert_eq!(
        query.cost_against(&mat.node_state(leaf), mat.reference()),
        0,
        "sample {sample}"
      );
    }
    let recount: usize = mat
      .graph()
      .get_edges()
      .map(|edge| mat.edge_subs(edge.key()).len())
      .sum();
    assert_eq!(recount, mat.parsimony());
  }

  #[rstest]
  #[case::exhaustive_one_by_one(exhaustive_params(1))]
  #[case::exhaustive_batches(exhaustive_params(8))]
  #[case::walk_one_by_one(walk_params(1, 0, false))]
  #[case::walk_batches(walk_params(8, 0, false))]
  #[case::seeded(walk_params(8, 0, true))]
  #[case::lbi_references(reference_params(8, lbi_references()))]
  #[case::long_branch_references(reference_params(8, ReferenceSelection::LongBranches { cuts: 4 }))]
  #[case::refreshing(refreshing_params(8))]
  #[case::rebuilding(rebuilding_params(8))]
  fn test_build_recovers_perfect_phylogeny(#[case] params: BuildParams) -> Result<(), Report> {
    let alphabet = Alphabet::default();
    for seed in 0..5 {
      let truth = perfect_phylogeny(60, seed)?;
      let built = build_mat(&truth.seqs, &alphabet, &params)?;
      assert_samples_reproduced(&built, &truth.seqs);
      assert_eq!(built.stats.placement_parsimony, truth.variable_mutations, "seed {seed}");
      assert_eq!(built.stats.final_parsimony, truth.variable_mutations, "seed {seed}");
      let comparison = built.stats.comparison.unwrap_or_default();
      assert_eq!(comparison.worse, 0, "seed {seed}");

      let tree = mat_to_tree(built.mat, &truth.sample_names)?;
      let seqs = izip_seqs(&truth.sample_names, &truth.seqs);
      let built_counts = fitch_edge_counts(&tree.graph, &tree.names, &seqs, &alphabet)?;
      let true_counts = fitch_edge_counts(&truth.graph, &truth.names, &seqs, &alphabet)?;
      let splits = compare_splits(
        &TreeWithCounts {
          graph: &tree.graph,
          names: &tree.names,
          edge_counts: &built_counts,
        },
        &TreeWithCounts {
          graph: &truth.graph,
          names: &truth.names,
          edge_counts: &true_counts,
        },
      )?;
      assert_eq!((splits.only_first, splits.only_second), (0, 0), "seed {seed}");
      assert_eq!(splits.shared, 60 - 3, "seed {seed}");
    }
    Ok(())
  }

  fn izip_seqs<'a>(names: &'a [String], seqs: &'a [Seq]) -> BTreeMap<&'a str, &'a Seq> {
    names.iter().map(String::as_str).zip(seqs).collect()
  }

  #[rstest]
  #[case::exhaustive(exhaustive_params(4))]
  #[case::walk(walk_params(4, 1, false))]
  #[case::seeded(walk_params(4, 0, true))]
  #[case::lbi_references(reference_params(4, lbi_references()))]
  #[case::long_branch_references(reference_params(4, ReferenceSelection::LongBranches { cuts: 4 }))]
  #[case::refreshing(refreshing_params(4))]
  #[case::rebuilding(rebuilding_params(4))]
  fn test_build_keeps_mat_consistent_with_homoplasy_and_missing_data(
    #[case] params: BuildParams,
  ) -> Result<(), Report> {
    let alphabet = Alphabet::default();
    for seed in 0..3 {
      let spec = CaseSpec {
        seed,
        sampling: (0..8)
          .map(|i| SampleBatch {
            tau: f64::from(i) * 0.25,
            count: 12,
          })
          .collect(),
        schedule: TcSchedule::constant(2.0)?,
        model: ModelFamily::Jc69,
        mu: 2e-2,
        length: 800,
      };
      let case = spec.simulate()?;
      let mut rng = get_random_number_generator(seed);
      let observation = TipObservation {
        ambiguity_rate: 0.01,
        unknown_run_rate: 0.3,
        unknown_run_lengths: 10..=80,
      };
      let observed = observe_tips(&case.tree.graph, &case.sequences.sequences, &observation, &mut rng)?;
      let (sample_names, seqs): (Vec<String>, Vec<Seq>) = case
        .tree
        .graph
        .get_leaves()
        .map(|leaf| (case.tree.names[leaf.key()].to_owned(), observed[&leaf.key()].clone()))
        .unzip();

      let built = build_mat(&seqs, &alphabet, &params)?;
      assert_samples_reproduced(&built, &seqs);
      assert!(built.stats.final_parsimony <= built.stats.placement_parsimony);
      let exact = matches!(
        params.search,
        SearchMode::Exhaustive | SearchMode::Walk { exact: true, .. } | SearchMode::FromReferences { .. }
      );
      if exact {
        assert_eq!(built.stats.comparison.unwrap_or_default().worse, 0, "seed {seed}");
      }

      let final_parsimony = built.stats.final_parsimony;
      let tree = mat_to_tree(built.mat, &sample_names)?;
      let named = izip_seqs(&sample_names, &seqs);
      let fitch: usize = fitch_edge_counts(&tree.graph, &tree.names, &named, &alphabet)?
        .values()
        .sum();
      assert!(
        fitch <= final_parsimony,
        "seed {seed}: Fitch {fitch} > MAT {final_parsimony}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_lbi_of_a_cherry_matches_the_integral_over_its_branches() -> Result<(), Report> {
    let reference: Seq = std::iter::repeat_n(A, 4).collect();
    let mut mat = Mat::new(reference);
    let root = mat.root();
    let long = mat.add_sample(root, vec![Sub::new(A, 0_usize, C)?, Sub::new(A, 1_usize, C)?], 0)?;
    let short = mat.add_sample(root, vec![Sub::new(A, 2_usize, C)?], 1)?;
    let tau = 1.5;
    let lbi = local_branching_index(&mat, tau, 0.0);
    let integral = |length: f64| tau * (1.0 - (-length / tau).exp());
    pretty_assert_abs_diff_eq!(lbi[&root], integral(2.0) + integral(1.0), epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(
      lbi[&long],
      integral(2.0) + (-2.0 / tau).exp() * integral(1.0),
      epsilon = 1e-12
    );
    pretty_assert_abs_diff_eq!(
      lbi[&short],
      integral(1.0) + (-1.0 / tau).exp() * integral(2.0),
      epsilon = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_eccentricities_of_a_cherry_are_its_longest_paths() -> Result<(), Report> {
    let reference: Seq = std::iter::repeat_n(A, 4).collect();
    let mut mat = Mat::new(reference);
    let root = mat.root();
    let long = mat.add_sample(root, vec![Sub::new(A, 0_usize, C)?, Sub::new(A, 1_usize, C)?], 0)?;
    let short = mat.add_sample(root, vec![Sub::new(A, 2_usize, C)?], 1)?;
    let eccentricity = eccentricities(&mat, &BTreeSet::new());
    assert_eq!(
      (eccentricity[&root], eccentricity[&long], eccentricity[&short]),
      (2, 3, 3)
    );
    Ok(())
  }
}
