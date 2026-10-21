use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::{GTR, Mu};
use crate::optimize::pipeline::OptimizeInput;
use crate::seq::alignment::{AncestralInput, EdgeSeqInput, SequenceSource};
use crate::test_utils::{leaf_seq_inputs, marginal_timetree_params, sequence_source};
use crate::timetree::params::{TimeMarginalMode, TimetreeParams};
use crate::timetree::pipeline::TimetreeInput;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use ndarray::aview1;
use proptest::arbitrary::any;
use proptest::collection::vec;
use proptest::sample::select;
use proptest::strategy::{Just, Strategy};
use rayon::ThreadPoolBuilder;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::distribution::{ContinuousCDF, StudentsT};
use statrs::stats_tests::ttest_onesample::ttest_onesample;
use statrs::stats_tests::{Alternative, NaNPolicy};
use std::collections::BTreeMap;
use treetime_graph::clades::{TipIds, first_node_by_key, graph_clade_keys};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::date::DateConstraint;
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_simulation::branch_lengths::realized;
use treetime_simulation::coalescent::{SampleBatch, SimulatedTree, simulate_coalescent};
use treetime_simulation::evolve::{SimulatedSequences, evolve_sequences};
use treetime_simulation::substitution_model::{NUCLEOTIDES, SubstitutionModel};
use treetime_simulation::tc_schedule::TcSchedule;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating, usize_to_u64};
use treetime_utils::sync::random::get_random_number_generator;
use treetime_utils::{make_error, make_internal_error, make_report};

pub(crate) const TIMETREE_CELL_TC: f64 = 5.0;
pub(crate) const TIMETREE_CELL_MU: f64 = 3e-3;
pub(crate) const TIMETREE_CELL_LENGTH: usize = 1000;
const TIMETREE_CELL_BATCH_TAUS: [f64; 5] = [0.0, 0.5, 1.0, 1.5, 2.0];
const TIMETREE_CELL_TIPS_PER_BATCH: usize = 4;
const CASE_SAMPLING_WINDOW: f64 = 2.0;
const CASE_MUS: [f64; 2] = [3e-3, 1e-2];
const CASE_MIN_LENGTH: usize = 300;
const CASE_MAX_LENGTH: usize = 1500;
const REPLICATE_THREADS: usize = 8;

pub(crate) fn timetree_cell(seed: u64) -> Result<CaseSpec, Report> {
  Ok(CaseSpec {
    seed,
    sampling: TIMETREE_CELL_BATCH_TAUS
      .iter()
      .map(|&tau| SampleBatch {
        tau,
        count: TIMETREE_CELL_TIPS_PER_BATCH,
      })
      .collect(),
    schedule: TcSchedule::constant(TIMETREE_CELL_TC)?,
    model: ModelFamily::RandomGtr,
    mu: TIMETREE_CELL_MU,
    length: TIMETREE_CELL_LENGTH,
  })
}

pub(crate) fn timetree_cell_params() -> TimetreeParams {
  TimetreeParams {
    model: GtrModelName::Infer,
    keep_root: true,
    clock_filter: 0.0,
    max_iter: 2,
    time_marginal: TimeMarginalMode::OnlyFinal,
    ..marginal_timetree_params()
  }
}

pub(crate) fn arb_simulated_case() -> impl Strategy<Value = CaseSpec> {
  let schedules = vec![
    TcSchedule::constant(5.0).expect("constant Tc schedule is valid"),
    TcSchedule::exponential(5.0, 2.5).expect("exponential Tc schedule is valid"),
    TcSchedule::piecewise(vec![2.5, 3.0], vec![5.0, 0.5, 1.5]).expect("piecewise Tc schedule is valid"),
  ];
  let mu_and_length = select(CASE_MUS.to_vec()).prop_flat_map(|mu| (Just(mu), case_min_length(mu)..=CASE_MAX_LENGTH));
  (any::<u64>(), vec(2..=4_usize, 3..=5), select(schedules), mu_and_length).prop_map(
    |(seed, tip_counts, schedule, (mu, length))| CaseSpec {
      seed,
      sampling: spread_over_window(&tip_counts),
      schedule,
      model: ModelFamily::RandomGtr,
      mu,
      length,
    },
  )
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct CaseSpec {
  pub seed: u64,
  pub sampling: Vec<SampleBatch>,
  pub schedule: TcSchedule,
  pub model: ModelFamily,
  pub mu: f64,
  pub length: usize,
}

impl CaseSpec {
  pub(crate) fn simulate(&self) -> Result<SimulatedCase, Report> {
    let mut rng = get_random_number_generator(self.seed);
    let tree = simulate_coalescent(&self.schedule, &self.sampling, &mut rng)?;
    let model = match self.model {
      ModelFamily::RandomGtr => SubstitutionModel::random(self.mu, &mut rng)?,
      ModelFamily::Jc69 => SubstitutionModel::jc69(self.mu)?,
    };
    let sequences = evolve_sequences(&tree.graph, &tree.edge_durations(), &model, self.length, &mut rng)?;
    Ok(SimulatedCase { tree, model, sequences })
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ModelFamily {
  RandomGtr,
  Jc69,
}

#[derive(Debug)]
pub(crate) struct SimulatedCase {
  pub tree: SimulatedTree,
  pub model: SubstitutionModel,
  pub sequences: SimulatedSequences,
}

impl SimulatedCase {
  pub(crate) fn timetree_cell_input(&self, columns: usize) -> Result<TimetreeInput, Report> {
    let sequences = first_columns(&self.sequences.sequences, columns)?;
    let branch_lengths = realized(&self.tree.graph, &sequences)?;
    timetree_input(&self.tree, &sequences, &branch_lengths)
  }
}

pub(crate) fn gtr_from_model(model: &SubstitutionModel, mu: f64) -> Result<GTR, Report> {
  GTR::builder()
    .n_states(NUCLEOTIDES.len())
    .mu(Mu::Normalized(mu))
    .W(model.w().clone())
    .pi(model.pi().clone())
    .build()
}

pub(crate) fn first_columns(
  sequences: &BTreeMap<GraphNodeKey, Seq>,
  count: usize,
) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
  sequences
    .iter()
    .map(|(&node, seq)| {
      if count == 0 || count > seq.len() {
        return make_error!(
          "Cannot take the first {count} columns of the sequence of node {node}, which has {} columns",
          seq.len()
        );
      }
      Ok((node, Seq::from(&seq.as_slice()[..count])))
    })
    .collect()
}

pub(crate) fn simulated_sequence_source(
  graph: &Graph,
  names: &NodeNames,
  sequences: &BTreeMap<GraphNodeKey, Seq>,
) -> Result<SequenceSource, Report> {
  let records = graph
    .get_leaves()
    .map(|leaf| {
      let Some(seq) = sequences.get(&leaf.key()) else {
        return make_error!("Leaf {} has no simulated sequence", &names[leaf.key()]);
      };
      Ok(AlignmentRecord {
        name: names[leaf.key()].to_owned(),
        seq: seq.clone(),
      })
    })
    .try_collect()?;
  sequence_source(Alphabet::default(), leaf_seq_inputs(graph, names, records))
}

pub(crate) fn exact_tip_dates(tree: &SimulatedTree) -> BTreeMap<GraphNodeKey, DateConstraint> {
  tree
    .graph
    .get_leaves()
    .map(|leaf| (leaf.key(), DateConstraint::exact(tree.date(leaf.key()))))
    .collect()
}

pub(crate) fn timetree_input(
  tree: &SimulatedTree,
  sequences: &BTreeMap<GraphNodeKey, Seq>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<TimetreeInput, Report> {
  Ok(TimetreeInput {
    graph: copy_graph(&tree.graph)?,
    names: tree.names.clone(),
    sources: vec![simulated_sequence_source(&tree.graph, &tree.names, sequences)?],
    dates: Some(exact_tip_dates(tree)),
    branch_lengths: edge_lengths(&tree.graph, branch_lengths)?,
  })
}

pub(crate) fn optimize_input(
  tree: &SimulatedTree,
  sequences: &BTreeMap<GraphNodeKey, Seq>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<OptimizeInput, Report> {
  Ok(OptimizeInput {
    graph: copy_graph(&tree.graph)?,
    names: tree.names.clone(),
    sources: vec![simulated_sequence_source(&tree.graph, &tree.names, sequences)?],
    branch_lengths: edge_lengths(&tree.graph, branch_lengths)?,
  })
}

pub(crate) fn ancestral_input(
  tree: &SimulatedTree,
  sequences: &BTreeMap<GraphNodeKey, Seq>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<AncestralInput, Report> {
  Ok(AncestralInput {
    graph: copy_graph(&tree.graph)?,
    edges: edge_lengths(&tree.graph, branch_lengths)?
      .into_iter()
      .map(|(edge, branch_length)| (edge, EdgeSeqInput { branch_length }))
      .collect(),
    sources: vec![simulated_sequence_source(&tree.graph, &tree.names, sequences)?],
  })
}

pub(crate) fn match_internal_nodes(
  truth: &SimulatedTree,
  graph: &Graph,
  names: &NodeNames,
) -> Result<BTreeMap<GraphNodeKey, GraphNodeKey>, Report> {
  let mut ids = TipIds::default();
  let true_keys = graph_clade_keys(&truth.graph, &truth.names, &mut ids)?;
  let keys = graph_clade_keys(graph, names, &mut ids)?;
  let nodes_by_key = first_node_by_key(keys.iter().map(|(&node, key)| (node, key.as_slice())));
  truth
    .graph
    .get_internal_nodes()
    .map(|node| {
      let true_node = node.key();
      nodes_by_key
        .get(true_keys[&true_node].as_slice())
        .map(|&inferred| (true_node, inferred))
        .ok_or_else(|| {
          make_report!(
            "Node {} of the true tree has no node with the same clade in the inferred tree",
            &truth.names[true_node]
          )
        })
    })
    .collect()
}

pub(crate) fn run_replicates<T: Send>(
  base_seed: u64,
  count: usize,
  replicate: impl Fn(u64) -> Result<T, Report> + Sync,
) -> Result<Vec<T>, Report> {
  let results = ThreadPoolBuilder::new()
    .num_threads(REPLICATE_THREADS)
    .build()?
    .install(|| {
      (0..count)
        .into_par_iter()
        .map(|index| {
          let seed = base_seed
            .checked_add(usize_to_u64(index))
            .ok_or_else(|| make_report!("Seed of replicate {index} overflows: base seed {base_seed}"))?;
          replicate(seed).wrap_err_with(|| format!("Replicate {index} with seed {seed} failed"))
        })
        .collect::<Vec<_>>()
    });
  results.into_iter().try_collect()
}

pub(crate) fn recovery_check(name: &str, errors: &[f64], margin: f64, alpha: f64) -> Result<MeanCheck, Report> {
  validate_margin(margin)?;
  mean_check(name, errors, alpha, Some(-margin), Some(margin))
}

pub(crate) fn upper_recovery_check(name: &str, errors: &[f64], margin: f64, alpha: f64) -> Result<MeanCheck, Report> {
  validate_margin(margin)?;
  mean_check(name, errors, alpha, None, Some(margin))
}

pub(crate) fn decrease_check(name: &str, log_ratios: &[f64], alpha: f64) -> Result<MeanCheck, Report> {
  mean_check(name, log_ratios, alpha, None, Some(0.0))
}

pub(crate) fn require_passed(checks: &[MeanCheck]) -> Result<(), Report> {
  if checks.iter().all(MeanCheck::passed) {
    return Ok(());
  }
  let lines = checks
    .iter()
    .map(|check| {
      let result = if check.passed() { "passed" } else { "FAILED" };
      format!("{result} {}", check.summary())
    })
    .join("\n");
  make_error!("Statistical checks failed:\n{lines}")
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct MeanCheck {
  pub name: String,
  pub n: usize,
  pub mean: f64,
  pub std_dev: f64,
  pub alpha: f64,
  pub interval: MeanInterval,
  pub tests: Vec<BoundTest>,
}

impl MeanCheck {
  pub(crate) fn passed(&self) -> bool {
    self.tests.iter().all(|test| test.p_value < self.alpha)
  }

  pub(crate) fn summary(&self) -> String {
    let Self {
      name,
      n,
      mean,
      std_dev,
      alpha,
      interval,
      tests,
    } = self;
    let tests = tests
      .iter()
      .map(|test| {
        let relation = match test.side {
          BoundSide::Above => ">",
          BoundSide::Below => "<",
        };
        format!("p(mean {relation} {}) = {:.3e}", test.bound, test.p_value)
      })
      .join(", ");
    format!(
      "{name}: n = {n}, mean = {mean:.6}, standard deviation = {std_dev:.6}, interval [{:.6}, {:.6}], {tests}, alpha = {alpha:.3e}",
      interval.lower, interval.upper
    )
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MeanInterval {
  pub lower: f64,
  pub upper: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BoundTest {
  pub side: BoundSide,
  pub bound: f64,
  pub p_value: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BoundSide {
  Above,
  Below,
}

fn case_min_length(mu: f64) -> usize {
  f64_to_usize_saturating((1.0 / mu).ceil()).max(CASE_MIN_LENGTH)
}

fn spread_over_window(tip_counts: &[usize]) -> Vec<SampleBatch> {
  let gaps = tip_counts.len().saturating_sub(1).max(1).to_f64_exact();
  tip_counts
    .iter()
    .enumerate()
    .map(|(index, &count)| SampleBatch {
      tau: CASE_SAMPLING_WINDOW * index.to_f64_exact() / gaps,
      count,
    })
    .collect()
}

fn copy_graph(graph: &Graph) -> Result<Graph, Report> {
  let mut copy = Graph::new();
  for node in graph.get_nodes() {
    let key = copy.add_node();
    if key != node.key() {
      return make_internal_error!("Copy of node {} got key {key}; the graph has removed nodes", node.key());
    }
  }
  for edge in graph.get_edges() {
    let key = copy.add_edge(edge.source(), edge.target())?;
    if key != edge.key() {
      return make_internal_error!("Copy of edge {} got key {key}; the graph has removed edges", edge.key());
    }
  }
  copy.build()?;
  Ok(copy)
}

fn edge_lengths(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<BTreeMap<GraphEdgeKey, Option<f64>>, Report> {
  graph
    .get_edges()
    .map(|edge| {
      let Some(&length) = branch_lengths.get(&edge.key()) else {
        return make_error!("Edge {} has no branch length", edge.key());
      };
      Ok((edge.key(), Some(length)))
    })
    .collect()
}

fn validate_margin(margin: f64) -> Result<(), Report> {
  if !(margin.is_finite() && margin > 0.0) {
    return make_error!("Margin must be finite and positive, got {margin}");
  }
  Ok(())
}

fn mean_check(
  name: &str,
  values: &[f64],
  alpha: f64,
  lower_bound: Option<f64>,
  upper_bound: Option<f64>,
) -> Result<MeanCheck, Report> {
  if !(alpha > 0.0 && alpha < 0.5) {
    return make_error!("{name}: significance level must be in (0, 0.5), got {alpha}");
  }
  let n = values.len();
  if n < 2 {
    return make_error!("{name}: a check of the mean needs at least 2 values, got {n}");
  }
  if let Some(value) = values.iter().find(|value| !value.is_finite()) {
    return make_error!("{name}: values must be finite, got {value}");
  }
  let view = aview1(values);
  let mean = view.mean().ok_or_else(|| make_report!("{name}: no values"))?;
  let std_dev = view.std(1.0);
  let degrees_of_freedom = (n - 1).to_f64_exact();
  let quantile = StudentsT::new(0.0, 1.0, degrees_of_freedom)
    .map_err(|error| {
      make_report!("{name}: Student's t distribution with {degrees_of_freedom} degrees of freedom: {error}")
    })?
    .inverse_cdf(1.0 - alpha);
  let half_width = quantile * std_dev / n.to_f64_exact().sqrt();
  let interval = MeanInterval {
    lower: lower_bound.map_or(f64::NEG_INFINITY, |_| mean - half_width),
    upper: upper_bound.map_or(f64::INFINITY, |_| mean + half_width),
  };
  let tests = [(BoundSide::Above, lower_bound), (BoundSide::Below, upper_bound)]
    .into_iter()
    .filter_map(|(side, bound)| Some((side, bound?)))
    .map(|(side, bound)| {
      let alternative = match side {
        BoundSide::Above => Alternative::Greater,
        BoundSide::Below => Alternative::Less,
      };
      let (_, p_value) = ttest_onesample(values.to_vec(), bound, alternative, NaNPolicy::Error)
        .map_err(|error| make_report!("{name}: one-sample t-test against {bound}: {error}"))?;
      Ok(BoundTest { side, bound, p_value })
    })
    .collect::<Result<Vec<_>, Report>>()?;
  Ok(MeanCheck {
    name: name.to_owned(),
    n,
    mean,
    std_dev,
    alpha,
    interval,
    tests,
  })
}

#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::progress::NoopProgress;
  use crate::test_utils::{
    CaseSpec, ModelFamily, SimulatedCase, TIMETREE_CELL_LENGTH, TIMETREE_CELL_MU, TIMETREE_CELL_TC, ancestral_input,
    arb_simulated_case, decrease_check, exact_tip_dates, first_columns, gtr_from_model, match_internal_nodes,
    optimize_input, recovery_check, require_passed, run_replicates, timetree_cell, timetree_cell_params,
    timetree_input,
  };
  use crate::timetree::params::TimetreeParams;
  use crate::timetree::pipeline;
  use eyre::Report;
  use helpers::{Margins, leaves_with_sequence, simulated_tree_from_newick};
  use itertools::{Itertools, izip};
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::{assert_eq, assert_ne};
  use proptest::prelude::{ProptestConfig, prop_assert, prop_assert_eq, proptest};
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::iter::repeat_n;
  use treetime_graph::node::GraphNodeKey;
  use treetime_primitives::Seq;
  use treetime_primitives::date::DateConstraint;
  use treetime_simulation::branch_lengths::{exact, realized};
  use treetime_simulation::substitution_model::SubstitutionModel;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{assert_error, make_error, o, pretty_assert_ulps_eq};

  #[test]
  fn test_simulation_gtr_from_model_keeps_the_rates_of_the_model() -> Result<(), Report> {
    let model = SubstitutionModel::random(0.01, &mut get_random_number_generator(7))?;

    let gtr = gtr_from_model(&model, 1.0)?;

    pretty_assert_ulps_eq!(1.0, gtr.mu(), max_ulps = 0);
    pretty_assert_ulps_eq!(model.pi().clone(), gtr.pi().clone(), max_ulps = 2);
    pretty_assert_ulps_eq!(model.w().clone(), gtr.W().clone(), max_ulps = 2);
    Ok(())
  }

  #[test]
  fn test_simulation_exact_tip_dates_date_every_tip_of_the_cell_at_its_batch() -> Result<(), Report> {
    let case = timetree_cell(1)?.simulate()?;

    let dates = exact_tip_dates(&case.tree);

    let dates_by_name = dates
      .into_iter()
      .map(|(node, date)| (case.tree.names[node].to_owned(), date))
      .collect::<BTreeMap<_, _>>();
    let batch_dates = [2025.0, 2024.5, 2024.0, 2023.5, 2023.0];
    let expected = izip!(1..=20, batch_dates.iter().flat_map(|&date| repeat_n(date, 4)))
      .map(|(tip, date)| (format!("tip_{tip:04}"), DateConstraint::exact(date)))
      .collect::<BTreeMap<_, _>>();
    assert_eq!(expected, dates_by_name);
    Ok(())
  }

  #[test]
  fn test_simulation_timetree_cell_input_uses_realized_lengths_of_the_first_columns() -> Result<(), Report> {
    let case = timetree_cell(2)?.simulate()?;

    let input = case.timetree_cell_input(500)?;

    let first = first_columns(&case.sequences.sequences, 500)?;
    let expected = realized(&case.tree.graph, &first)?
      .into_iter()
      .map(|(edge, length)| (edge, Some(length)))
      .collect::<BTreeMap<_, _>>();
    assert_eq!(expected, input.branch_lengths);
    assert_eq!(Some(exact_tip_dates(&case.tree)), input.dates);
    assert_eq!(
      (500, 20),
      (input.sources[0].length, leaves_with_sequence(&input.sources[0]))
    );
    Ok(())
  }

  #[test]
  fn test_simulation_inputs_carry_the_given_lengths_and_every_leaf_sequence() -> Result<(), Report> {
    let case: SimulatedCase = CaseSpec {
      model: ModelFamily::Jc69,
      ..timetree_cell(4)?
    }
    .simulate()?;
    let sequences = &case.sequences.sequences;
    let lengths = exact(&case.tree.edge_durations(), &case.model);

    let timetree = timetree_input(&case.tree, sequences, &lengths)?;
    let optimize = optimize_input(&case.tree, sequences, &lengths)?;
    let ancestral = ancestral_input(&case.tree, sequences, &lengths)?;

    let expected = lengths
      .into_iter()
      .map(|(edge, length)| (edge, Some(length)))
      .collect::<BTreeMap<_, _>>();
    let leaf_counts =
      [&timetree.sources, &optimize.sources, &ancestral.sources].map(|sources| leaves_with_sequence(&sources[0]));
    assert_eq!(expected, timetree.branch_lengths);
    assert_eq!(expected, optimize.branch_lengths);
    assert_eq!(expected, ancestral.branch_lengths());
    assert_eq!([20; 3], leaf_counts);
    Ok(())
  }

  #[test]
  fn test_simulation_model_family_changes_the_model_but_not_the_tree() -> Result<(), Report> {
    let gtr_case = timetree_cell(5)?.simulate()?;
    let jc69_case = CaseSpec {
      model: ModelFamily::Jc69,
      ..timetree_cell(5)?
    }
    .simulate()?;

    assert_eq!(gtr_case.tree.tau, jc69_case.tree.tau);
    assert_eq!(&array![0.25, 0.25, 0.25, 0.25], jc69_case.model.pi());
    assert_ne!(gtr_case.model.pi(), jc69_case.model.pi());
    Ok(())
  }

  #[test]
  fn test_simulation_timetree_cell_runs_through_the_timetree_pipeline_and_matches_every_clade() -> Result<(), Report> {
    let case = timetree_cell(3)?.simulate()?;
    let params = TimetreeParams {
      clock_rate: Some(TIMETREE_CELL_MU),
      coalescent: Some(TIMETREE_CELL_TC),
      ..timetree_cell_params()
    };

    let output = pipeline::run(
      &params,
      case.timetree_cell_input(TIMETREE_CELL_LENGTH)?,
      None,
      None,
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    let matched = match_internal_nodes(&case.tree, &output.graph, &output.names)?;
    assert_eq!(19, matched.len());
    Ok(())
  }

  #[test]
  fn test_simulation_first_columns_keeps_the_leading_columns() -> Result<(), Report> {
    let sequences = btreemap! {
      GraphNodeKey(0) => Seq::try_from_str("ACGT")?,
      GraphNodeKey(1) => Seq::try_from_str("TTGA")?,
    };

    let actual = first_columns(&sequences, 2)?;

    let expected = btreemap! {
      GraphNodeKey(0) => Seq::try_from_str("AC")?,
      GraphNodeKey(1) => Seq::try_from_str("TT")?,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_columns(0, "Cannot take the first 0 columns of the sequence of node 0, which has 4 columns")]
  #[case::too_many(    5, "Cannot take the first 5 columns of the sequence of node 0, which has 4 columns")]
  #[trace]
  fn test_simulation_first_columns_rejects_a_count_outside_the_sequence(
    #[case] count: usize,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let sequences = btreemap! { GraphNodeKey(0) => Seq::try_from_str("ACGT")? };

    assert_error!(first_columns(&sequences, count), expected);
    Ok(())
  }

  #[test]
  fn test_simulation_match_internal_nodes_pairs_clades_in_any_child_order() -> Result<(), Report> {
    let truth = simulated_tree_from_newick("((A,B)ab,(C,D)cd)root;")?;
    let inferred = simulated_tree_from_newick("((D,C)x,(B,A)y)z;")?;

    let matched = match_internal_nodes(&truth, &inferred.graph, &inferred.names)?;

    let matched_names = matched
      .into_iter()
      .map(|(true_node, node)| (truth.names[true_node].to_owned(), inferred.names[node].to_owned()))
      .collect::<BTreeMap<_, _>>();
    let expected = btreemap! {
      o!("ab") => o!("y"),
      o!("cd") => o!("x"),
      o!("root") => o!("z"),
    };
    assert_eq!(expected, matched_names);
    Ok(())
  }

  #[test]
  fn test_simulation_match_internal_nodes_rejects_a_missing_clade() -> Result<(), Report> {
    let truth = simulated_tree_from_newick("((A,B)ab,(C,D)cd)root;")?;
    let inferred = simulated_tree_from_newick("((A,C)x,(B,D)y)z;")?;

    assert_error!(
      match_internal_nodes(&truth, &inferred.graph, &inferred.names),
      "Node ab of the true tree has no node with the same clade in the inferred tree"
    );
    Ok(())
  }

  #[test]
  fn test_simulation_run_replicates_returns_results_in_replicate_order() -> Result<(), Report> {
    let actual = run_replicates(100, 12, |seed| Ok(seed * 2))?;

    let expected = (100..112).map(|seed| seed * 2).collect_vec();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_simulation_run_replicates_reports_the_first_failed_replicate_with_its_seed() {
    let result = run_replicates(100, 12, |seed| {
      if seed >= 103 {
        return make_error!("no estimate");
      }
      Ok(seed)
    });

    assert_error!(result, "Replicate 3 with seed 103 failed: no estimate");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::two_sided_within_margin(    Margins::TwoSided, [ -1.0,  1.0], 2.0, true )]
  #[case::two_sided_wider_than_margin(Margins::TwoSided, [ -1.0,  1.0], 0.5, false)]
  #[case::two_sided_below_lower_margin(Margins::TwoSided, [-11.0, -9.0], 2.0, false)]
  #[case::two_sided_above_upper_margin(Margins::TwoSided, [  9.0, 11.0], 2.0, false)]
  #[case::upper_below_lower_margin(   Margins::Upper,    [-11.0, -9.0], 2.0, true )]
  #[case::upper_above_upper_margin(   Margins::Upper,    [  9.0, 11.0], 2.0, false)]
  #[trace]
  fn test_simulation_recovery_checks_pass_only_when_the_interval_is_inside_the_margin(
    #[case] margins: Margins,
    #[case] errors: [f64; 2],
    #[case] margin: f64,
    #[case] expected: bool,
  ) -> Result<(), Report> {
    let actual = margins.check("error", &errors, margin, 0.25)?.passed();

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::clear_decrease([-11.0, -9.0], true )]
  #[case::no_change(     [ -1.0,  1.0], false)]
  #[case::increase(      [  9.0, 11.0], false)]
  #[trace]
  fn test_simulation_decrease_check_passes_only_for_a_mean_below_zero(
    #[case] log_ratios: [f64; 2],
    #[case] expected: bool,
  ) -> Result<(), Report> {
    let actual = decrease_check("log ratio", &log_ratios, 0.25)?.passed();

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_simulation_require_passed_states_every_check_when_one_fails() -> Result<(), Report> {
    let checks = [
      recovery_check("rate 1", &[-1.0, 1.0], 2.0, 0.25)?,
      recovery_check("rate 2", &[-1.0, 1.0], 0.5, 0.25)?,
    ];

    assert_error!(
      require_passed(&checks),
      "Statistical checks failed:\n\
       passed rate 1: n = 2, mean = 0.000000, standard deviation = 1.414214, interval [-1.000000, 1.000000], p(mean > -2) = 1.476e-1, p(mean < 2) = 1.476e-1, alpha = 2.500e-1\n\
       FAILED rate 2: n = 2, mean = 0.000000, standard deviation = 1.414214, interval [-1.000000, 1.000000], p(mean > -0.5) = 3.524e-1, p(mean < 0.5) = 3.524e-1, alpha = 2.500e-1"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_value(     &[0.1],           0.1, "error: a check of the mean needs at least 2 values, got 1")]
  #[case::not_finite(    &[0.1, f64::NAN], 0.1, "error: values must be finite, got NaN")]
  #[case::alpha_zero(    &[0.1, 0.2],      0.0, "error: significance level must be in (0, 0.5), got 0")]
  #[case::alpha_too_high(&[0.1, 0.2],      0.5, "error: significance level must be in (0, 0.5), got 0.5")]
  #[trace]
  fn test_simulation_recovery_check_rejects_invalid_input(
    #[case] errors: &[f64],
    #[case] alpha: f64,
    #[case] expected: &str,
  ) {
    assert_error!(recovery_check("error", errors, 0.03, alpha), expected);
  }

  #[test]
  fn test_simulation_recovery_check_rejects_a_non_positive_margin() {
    assert_error!(
      recovery_check("error", &[0.1, 0.2], 0.0, 0.1),
      "Margin must be finite and positive, got 0"
    );
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    #[test]
    fn test_prop_simulation_arb_simulated_case_follows_the_case_ranges(spec in arb_simulated_case()) {
      let taus = spec.sampling.iter().map(|batch| batch.tau).collect_vec();
      let tip_count: usize = spec.sampling.iter().map(|batch| batch.count).sum();
      prop_assert!((3..=5).contains(&spec.sampling.len()), "batch count {}", spec.sampling.len());
      prop_assert!(spec.sampling.iter().all(|batch| (2..=4).contains(&batch.count)), "tip counts {:?}", spec.sampling);
      prop_assert_eq!(Some(&0.0), taus.first());
      prop_assert_eq!(Some(&2.0), taus.last());
      prop_assert!([3e-3, 1e-2].contains(&spec.mu), "mu {}", spec.mu);
      prop_assert!((300..=1500).contains(&spec.length), "length {}", spec.length);
      prop_assert!(spec.length.to_f64_exact() * spec.mu >= 1.0, "length {} with mu {}", spec.length, spec.mu);

      let case = spec.simulate().expect("a case of the strategy simulates");

      prop_assert_eq!(tip_count, case.tree.graph.num_leaves());
      let lengths = case.sequences.sequences.values().map(Seq::len).unique().collect_vec();
      prop_assert_eq!(vec![spec.length], lengths);
    }
  }

  mod helpers {
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{MeanCheck, recovery_check, upper_recovery_check};
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_io::nwk::nwk_read;
    use treetime_simulation::coalescent::SimulatedTree;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Margins {
      TwoSided,
      Upper,
    }

    impl Margins {
      pub(super) fn check(self, name: &str, errors: &[f64], margin: f64, alpha: f64) -> Result<MeanCheck, Report> {
        match self {
          Self::TwoSided => recovery_check(name, errors, margin, alpha),
          Self::Upper => upper_recovery_check(name, errors, margin, alpha),
        }
      }
    }

    pub(super) fn leaves_with_sequence(source: &SequenceSource) -> usize {
      source.nodes.values().filter(|node| node.seq.is_some()).count()
    }

    pub(super) fn simulated_tree_from_newick(newick: &str) -> Result<SimulatedTree, Report> {
      let parse = nwk_read(newick.as_bytes())?;
      Ok(SimulatedTree {
        graph: parse.graph,
        names: parse.names,
        tau: BTreeMap::new(),
      })
    }
  }
}
