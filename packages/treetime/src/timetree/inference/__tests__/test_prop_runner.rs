#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::undated_leaves;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::inference::time_inference::{TimeInference, unit_gammas};
  use eyre::Report;
  use generators::{TimetreeCase, gen_timetree_case};
  use helpers::run_twice_around_other_gammas;
  use itertools::Itertools;
  use proptest::prelude::*;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  const CLOCK_RATE: f64 = 0.001;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn test_prop_runner_run_timetree_idempotent(case in gen_timetree_case()) {
      let runs = run_twice_around_other_gammas(&case).unwrap();
      prop_assert_eq!(runs.first, runs.after_other_gammas);
    }
  }

  mod generators {
    use super::*;

    #[derive(Debug, Clone)]
    pub(super) struct TimetreeCase {
      pub newick: String,
      pub fasta: String,
      pub dates: Vec<(String, Option<f64>)>,
      pub other_gamma: f64,
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
          prop::collection::vec("[ACGT]{24}", n_leaves),
          0.5_f64..2.0,
        )
          .prop_map(move |(picks, lengths, dates, undated, sequences, other_gamma)| {
            let names = (0..n_leaves).map(|index| format!("L{index}")).collect::<Vec<_>>();
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
                  let date = (index < 3 || !undated[index]).then_some(dates[index]);
                  (name.clone(), date)
                })
                .collect(),
              other_gamma,
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

    type RunResult = Result<TimeInference, String>;

    pub(super) struct TwoRuns {
      pub first: RunResult,
      pub after_other_gammas: RunResult,
    }

    pub(super) fn run_twice_around_other_gammas(case: &TimetreeCase) -> Result<TwoRuns, Report> {
      let nwk_parsed = nwk_read_str(&case.newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = read_many_fasta_str(&case.fasta, &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let partition = PartitionTimetree::Dense(DenseReconstruction {
        partition: PartitionMarginalDense::new(0, alphabet, &graph, &node_seq_inputs(&graph, &names, aln))?,
        gtr: jc69(JC69Params::default())?,
        node_states: BTreeMap::new(),
        edges: MarginalEdges::default(),
      });
      let branch_model =
        BranchModel::Marginal(partition.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?);

      let dates: DatesMap = case
        .dates
        .iter()
        .map(|(name, date)| (name.clone(), date.map(DateConstraint::exact)))
        .collect();
      let constraints = load_date_constraints(&dates, &graph, &names, &NoopProgress)?;
      let leaf_bad_branches = undated_leaves(&graph, &constraints);
      let clock_model = ClockModel::for_testing(CLOCK_RATE, 0.0);

      let unit = unit_gammas(&graph);
      let other: BTreeMap<GraphEdgeKey, f64> = unit.keys().map(|key| (*key, case.other_gamma)).collect();
      let run = |gammas: &BTreeMap<GraphEdgeKey, f64>| -> RunResult {
        run_timetree(
          &graph,
          &constraints,
          &leaf_bad_branches,
          gammas,
          &branch_model,
          &branch_lengths,
          &names,
          &clock_model,
          None,
          false,
          &NoopProgress,
        )
        .map_err(|report| format!("{report:?}"))
      };

      let first = run(&unit);
      let _other_gammas_result = run(&other);
      let after_other_gammas = run(&unit);
      Ok(TwoRuns {
        first,
        after_other_gammas,
      })
    }
  }
}
