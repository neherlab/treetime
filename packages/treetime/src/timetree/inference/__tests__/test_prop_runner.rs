#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::dates_by_node;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::optimization::relaxed_clock::unit_gammas;
  use eyre::Report;
  use generators::{TimetreeCase, gen_timetree_case};
  use helpers::run_once;
  use itertools::Itertools;
  use proptest::prelude::*;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_grid::MaxGridPoints;
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  const CLOCK_RATE: f64 = 0.001;

  const MIN_DATED_LEAVES: usize = 3;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn test_prop_runner_leaf_bad_branch_flags_and_exact_dates_pass_through(case in gen_timetree_case()) {
      let leaves = run_once(&case).unwrap();

      for (name, leaf) in &leaves {
        prop_assert_eq!(leaf.expected_bad, leaf.bad, "bad-branch flag of leaf {}", name);
        if let Some(date) = leaf.date {
          prop_assert_eq!(Some(date), leaf.time, "posterior time of exactly dated leaf {}", name);
        }
      }
    }
  }

  mod generators {
    use super::*;

    #[derive(Debug, Clone)]
    pub(super) struct TimetreeCase {
      pub newick: String,
      pub fasta: String,
      pub dates: Vec<(String, Option<f64>)>,
      pub outliers: Vec<String>,
    }

    pub(super) fn gen_timetree_case() -> impl Strategy<Value = TimetreeCase> {
      (4_usize..=7).prop_flat_map(|n_leaves| {
        (
          prop::collection::vec(
            (any::<prop::sample::Index>(), any::<prop::sample::Index>()),
            n_leaves - 1,
          ),
          prop::collection::vec(0.001_f64..0.02, 2 * (n_leaves - 1)),
          prop::collection::vec(2000.0_f64..2020.0, n_leaves),
          prop::collection::vec(any::<bool>(), n_leaves),
          prop::collection::vec(any::<bool>(), n_leaves),
          prop::collection::vec("[ACGT]{24}", n_leaves),
        )
          .prop_map(move |(picks, lengths, dates, undated, outlier, sequences)| {
            let names = (0..n_leaves).map(|index| format!("L{index}")).collect::<Vec<_>>();
            let always_good = |index: usize| index < MIN_DATED_LEAVES;
            TimetreeCase {
              newick: random_newick(&names, &picks, &lengths),
              fasta: names
                .iter()
                .zip(&sequences)
                .map(|(name, sequence)| format!(">{name}\n{sequence}\n"))
                .join(""),
              dates: names
                .iter()
                .enumerate()
                .map(|(index, name)| {
                  let date = (always_good(index) || !undated[index]).then_some(dates[index]);
                  (name.clone(), date)
                })
                .collect(),
              outliers: names
                .iter()
                .enumerate()
                .filter(|(index, _)| !always_good(*index) && outlier[*index])
                .map(|(_, name)| name.clone())
                .collect(),
            }
          })
      })
    }

    fn random_newick(
      names: &[String],
      picks: &[(prop::sample::Index, prop::sample::Index)],
      lengths: &[f64],
    ) -> String {
      let mut subtrees = names.to_vec();
      for (merge, (left, right)) in picks.iter().enumerate() {
        let first = subtrees.remove(left.index(subtrees.len()));
        let second = subtrees.remove(right.index(subtrees.len()));
        let (first_length, second_length) = (lengths[2 * merge], lengths[2 * merge + 1]);
        subtrees.push(format!("({first}:{first_length},{second}:{second_length})"));
      }
      format!("{}root;", subtrees.concat())
    }
  }

  mod helpers {
    use super::*;

    pub(super) struct LeafOutcome {
      pub expected_bad: bool,
      pub bad: bool,
      pub date: Option<f64>,
      pub time: Option<f64>,
    }

    pub(super) fn run_once(case: &TimetreeCase) -> Result<BTreeMap<String, LeafOutcome>, Report> {
      let nwk_parsed = nwk_read(case.newick.as_bytes())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(case.fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let partition = MarginalReconstruction::Dense(DenseReconstruction {
        partition: PartitionMarginalDense::new(0, alphabet, &graph, &node_seq_inputs(&graph, &names, aln))?,
        gtr: jc69(JC69Params::default())?,
        node_states: BTreeMap::new(),
        edges: MarginalEdges::default(),
      });
      let branch_model = BranchModel::Marginal(
        partition
          .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?
          .0,
      );

      let dates: BTreeMap<String, Option<DateConstraint>> = case
        .dates
        .iter()
        .map(|(name, date)| (name.clone(), date.map(DateConstraint::exact)))
        .collect();
      let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;
      let outliers: BTreeSet<_> = case
        .outliers
        .iter()
        .map(|name| find_node_key_by_name(&graph, &names, name).expect("generated leaf must exist"))
        .collect();
      let leaf_bad_branches = bad_leaves(&graph, &constraints, &outliers);
      let clock_model = ClockModel::for_testing(CLOCK_RATE, 0.0);
      let gammas = unit_gammas(&graph);

      let inference = run_timetree(
        &TimeInferenceInputs {
          graph: &graph,
          date_constraints: &constraints,
          leaf_bad_branches: &leaf_bad_branches,
          gammas: &gammas,
          branch_model: &branch_model,
          branch_lengths: &branch_lengths,
          names: &names,
          clock_model: &clock_model,
          clock_rate_fixed: false,
          no_indels: false,
          max_grid_points: MaxGridPoints::default(),
        },
        None,
        &NoopProgress,
      )?;

      let given_dates: BTreeMap<&String, Option<f64>> = case.dates.iter().map(|(name, date)| (name, *date)).collect();
      Ok(
        graph
          .get_leaves()
          .map(|leaf| {
            let key = leaf.key();
            let name = names[&key].clone().expect("generated leaves are named");
            let outcome = LeafOutcome {
              expected_bad: given_dates[&name].is_none() || case.outliers.contains(&name),
              bad: inference.bad_branches[&key],
              date: given_dates[&name],
              time: inference.posterior[&key].time,
            };
            (name, outcome)
          })
          .collect(),
      )
    }
  }
}
