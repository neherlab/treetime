#[cfg(test)]
mod tests {
  use crate::partition::create::Representation;
  use eyre::Report;
  use helpers::{LARGE_STAR, LONG_BRANCHES, Placement, PolytomyTree, SHORT_BRANCHES, Star};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::long_branches_at_root(     LONG_BRANCHES,  Placement::AtRoot)]
  #[case::short_branches_at_root(    SHORT_BRANCHES, Placement::AtRoot)]
  #[case::short_branches_below_root( SHORT_BRANCHES, Placement::BelowRoot)]
  #[trace]
  fn test_root_polytomy_sparse_matches_dense(#[case] star: Star, #[case] placement: Placement) -> Result<(), Report> {
    let dense = PolytomyTree::new(star, placement, Representation::Dense)?;
    let sparse = PolytomyTree::new(star, placement, Representation::Sparse)?;

    let expected = dense.summary()?;
    let actual = sparse.summary()?;

    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_root_polytomy_child_contributions_evaluate_the_sites_of_the_same_polytomy_below_the_root()
  -> Result<(), Report> {
    let below_root = PolytomyTree::new(LARGE_STAR, Placement::BelowRoot, Representation::Sparse)?;
    let at_root = PolytomyTree::new(LARGE_STAR, Placement::AtRoot, Representation::Sparse)?;

    let expected = below_root.child_contribution_sites()?;
    let actual = at_root.child_contribution_sites()?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::short_branches_at_root(    SHORT_BRANCHES, Placement::AtRoot)]
  #[case::short_branches_below_root( SHORT_BRANCHES, Placement::BelowRoot)]
  #[case::large_star_at_root(        LARGE_STAR,     Placement::AtRoot)]
  #[trace]
  fn test_root_polytomy_child_contributions_hold_each_distinct_site_once_and_count_every_site_with_data(
    #[case] star: Star,
    #[case] placement: Placement,
  ) -> Result<(), Report> {
    let tree = PolytomyTree::new(star, placement, Representation::Sparse)?;

    let summaries = tree.child_contribution_summaries()?;
    let expected = summaries
      .iter()
      .map(|(child, summary)| (child.clone(), (summary.distinct_coefficients, summary.sites_with_data)))
      .collect::<BTreeMap<_, _>>();
    let actual = summaries
      .iter()
      .map(|(child, summary)| (child.clone(), (summary.entries, summary.multiplicity)))
      .collect::<BTreeMap<_, _>>();

    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::make_report;
    use crate::optimize::dispatch::optimize_edge_lengths;
    use crate::optimize::indel::EdgeIndels;
    use crate::optimize::likelihood::evaluate_partitions_log_lh_only;
    use crate::optimize::params::BranchOptMethod;
    use crate::optimize::solver_outcome::UnconvergedEdges;
    use crate::partition::create::Representation;
    use crate::partition::marginal::reconstruction::MarginalReconstruction;
    use crate::partition::optimize::contribution::OptimizationContribution;
    use crate::test_utils::{NUC_ALPHABET, run_dense_marginal, run_sparse_marginal};
    use eyre::Report;
    use itertools::Itertools;
    use pretty_assertions::assert_eq;
    use rand::RngExt;
    use std::collections::{BTreeMap, BTreeSet};
    use std::iter::repeat_with;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::{NwkParse, nwk_read};
    use treetime_primitives::AlignmentRecord;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::sync::random::get_random_number_generator;

    const NUCS: [u8; 4] = *b"ACGT";

    const TRIAL_LENGTHS: [f64; 3] = [1e-3, 0.05, 0.5];

    #[rustfmt::skip]
    pub(super) const LONG_BRANCHES:  Star = Star { children: 8,   length: 60,  branch_length: 0.8,  n_runs: 1, run_lengths: (4, 12), seed: 11 };
    #[rustfmt::skip]
    pub(super) const SHORT_BRANCHES: Star = Star { children: 12,  length: 60,  branch_length: 0.01, n_runs: 2, run_lengths: (4, 12), seed: 12 };
    #[rustfmt::skip]
    pub(super) const LARGE_STAR:     Star = Star { children: 150, length: 120, branch_length: 1e-4, n_runs: 2, run_lengths: (4, 12), seed: 13 };

    #[derive(Clone, Copy, Debug)]
    pub(super) struct Star {
      children: usize,
      length: usize,
      branch_length: f64,
      n_runs: usize,
      run_lengths: (usize, usize),
      seed: u64,
    }

    impl Star {
      fn newick(&self, placement: Placement) -> String {
        let children = (0..self.children)
          .map(|i| format!("c{i}:{}", self.branch_length))
          .join(",");
        match placement {
          Placement::AtRoot => format!("({children})root;"),
          Placement::BelowRoot => format!("(({children})P:0.02,X:0.02)root;"),
        }
      }

      fn fasta(&self) -> String {
        let mut rng = get_random_number_generator(self.seed);
        let root = repeat_with(|| rng.random_range(0..NUCS.len()))
          .take(self.length)
          .collect_vec();
        let mut child = |n_runs: usize| {
          let mut states = root.clone();
          let pos = rng.random_range(0..self.length);
          states[pos] = (states[pos] + 1) % NUCS.len();
          let mut seq = states.into_iter().map(|state| NUCS[state]).collect_vec();
          for _ in 0..n_runs {
            let width = rng.random_range(self.run_lengths.0..self.run_lengths.1);
            let start = rng.random_range(0..self.length - width);
            seq[start..start + width].fill(b'N');
          }
          String::from_utf8_lossy(&seq).into_owned()
        };
        let children = (0..self.children)
          .map(|i| {
            let n_runs = if i == 0 { 0 } else { self.n_runs };
            format!(">c{i}\n{}\n", child(n_runs))
          })
          .join("");
        format!("{children}>X\n{}\n", child(0))
      }
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Placement {
      AtRoot,
      BelowRoot,
    }

    pub(super) struct PolytomyTree {
      graph: Graph,
      names: NodeNames,
      branch_lengths: EdgeMap<Option<f64>>,
      length: usize,
      sites_with_data: BTreeMap<String, usize>,
      reconstruction: MarginalReconstruction,
      log_lh: f64,
    }

    impl PolytomyTree {
      pub(super) fn new(star: Star, placement: Placement, representation: Representation) -> Result<Self, Report> {
        let NwkParse {
          graph,
          branch_lengths,
          names,
          ..
        } = nwk_read(star.newick(placement).as_bytes())?;
        let aln: Vec<AlignmentRecord> = fasta_read(star.fasta().as_bytes(), &*NUC_ALPHABET)?
          .into_iter()
          .map(AlignmentRecord::from)
          .collect();
        let unknown = NUC_ALPHABET.unknown();
        let sites_with_data = aln
          .iter()
          .map(|record| {
            let with_data = record.seq.iter().filter(|&&state| state != unknown).count();
            (record.name.clone(), with_data)
          })
          .collect();
        let gtr = jc69(JC69Params::default())?;
        let (log_lh, reconstruction) = match representation {
          Representation::Dense => {
            let (log_lh, dense) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr)?;
            (log_lh, MarginalReconstruction::Dense(dense))
          },
          Representation::Sparse => {
            let (log_lh, sparse) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr)?;
            (log_lh, MarginalReconstruction::Sparse(sparse))
          },
        };
        Ok(Self {
          graph,
          names,
          branch_lengths,
          length: star.length,
          sites_with_data,
          reconstruction,
          log_lh,
        })
      }

      pub(super) fn summary(&self) -> Result<BTreeMap<String, f64>, Report> {
        let root = [(String::from("root log_lh"), self.log_lh)];
        let objectives = self
          .graph
          .get_edges()
          .map(|edge| edge.key())
          .cartesian_product(TRIAL_LENGTHS)
          .map(|(edge_key, length)| -> Result<_, Report> {
            let child = &self.names[self.graph.edge_endpoints(edge_key)?.1];
            let current = self.branch_lengths[&edge_key].unwrap_or(0.0);
            let contribution = [self.reconstruction.view().create_edge_contribution(edge_key)?];
            let at_length = evaluate_partitions_log_lh_only(&contribution, length)?;
            let at_current = evaluate_partitions_log_lh_only(&contribution, current)?;
            Ok((
              format!("edge {child} t={length} objective change"),
              at_length.value() - at_current.value(),
            ))
          })
          .collect::<Result<Vec<_>, Report>>()?;
        let optima = self
          .optimal_branch_lengths()?
          .into_iter()
          .map(|(child, length)| (format!("edge {child} optimum"), length));
        let posteriors = self.internal_posteriors()?;
        Ok(
          root
            .into_iter()
            .chain(objectives)
            .chain(optima)
            .chain(posteriors)
            .collect(),
        )
      }

      pub(super) fn child_contribution_sites(&self) -> Result<BTreeMap<String, usize>, Report> {
        self
          .graph
          .get_edges()
          .map(|edge| -> Result<_, Report> {
            let child = self.names[self.graph.edge_endpoints(edge.key())?.1].to_owned();
            Ok((child, edge.key()))
          })
          .filter_ok(|(child, _)| child.starts_with('c'))
          .map(|item| -> Result<_, Report> {
            let (child, edge_key) = item?;
            match self.reconstruction.view().create_edge_contribution(edge_key)? {
              OptimizationContribution::Sparse(sparse) => Ok((child, sparse.site_contributions.len())),
              OptimizationContribution::Dense(_) => Err(make_report!("edge {child} has a dense contribution")),
            }
          })
          .collect()
      }

      pub(super) fn child_contribution_summaries(&self) -> Result<BTreeMap<String, ContributionSummary>, Report> {
        self
          .graph
          .get_edges()
          .map(|edge| -> Result<_, Report> {
            let child = self.names[self.graph.edge_endpoints(edge.key())?.1].to_owned();
            Ok((child, edge.key()))
          })
          .filter_ok(|(child, _)| child.starts_with('c'))
          .map(|item| -> Result<_, Report> {
            let (child, edge_key) = item?;
            let OptimizationContribution::Sparse(sparse) =
              self.reconstruction.view().create_edge_contribution(edge_key)?
            else {
              return Err(make_report!("edge {child} has a dense contribution"));
            };
            let sites = &sparse.site_contributions;
            let distinct: BTreeSet<Vec<u64>> = sites
              .iter()
              .map(|site| site.coefficients.iter().map(|value| value.to_bits()).collect())
              .collect();
            let summary = ContributionSummary {
              entries: sites.len(),
              distinct_coefficients: distinct.len(),
              multiplicity: sites.iter().map(|site| site.multiplicity).sum(),
              sites_with_data: self.sites_with_data[&child].to_f64_exact(),
            };
            Ok((child, summary))
          })
          .collect()
      }

      fn optimal_branch_lengths(&self) -> Result<BTreeMap<String, f64>, Report> {
        let contributions: EdgeMap<Vec<OptimizationContribution>> = self
          .graph
          .get_edges()
          .map(|edge| -> Result<_, Report> {
            let contribution = self.reconstruction.view().create_edge_contribution(edge.key())?;
            Ok((edge.key(), vec![contribution]))
          })
          .try_collect()?;
        let mut branch_lengths = self.branch_lengths.clone();
        let unconverged = optimize_edge_lengths(
          &self.graph,
          self.length,
          &contributions,
          &EdgeIndels::Excluded,
          BranchOptMethod::Brent,
          &mut branch_lengths,
        )?;
        assert_eq!(UnconvergedEdges::default(), unconverged);
        branch_lengths
          .iter()
          .map(|(edge_key, length)| -> Result<_, Report> {
            let child = &self.names[self.graph.edge_endpoints(edge_key)?.1];
            let length = length.ok_or_else(|| make_report!("edge {child} has no optimized length"))?;
            Ok((child.to_owned(), length))
          })
          .collect()
      }

      fn internal_posteriors(&self) -> Result<Vec<(String, f64)>, Report> {
        self
          .graph
          .get_internal_nodes()
          .map(|node| node.key())
          .cartesian_product(0..self.length)
          .map(|(key, pos)| -> Result<_, Report> {
            let name = &self.names[key];
            let posterior = match &self.reconstruction {
              MarginalReconstruction::Dense(dense) => dense.node_states[&key].profile.dis.row(pos).to_vec(),
              MarginalReconstruction::Sparse(sparse) => {
                let node = &sparse.node_states[&key];
                node.profile.variable.get(pos).map_or_else(
                  || node.profile.fixed[node.sequence[pos]].dis.to_vec(),
                  |var| var.dis.to_vec(),
                )
              },
            };
            Ok(
              posterior
                .into_iter()
                .enumerate()
                .map(|(state, p)| (format!("node {name} position {pos} state {state}"), p))
                .collect_vec(),
            )
          })
          .flatten_ok()
          .collect()
      }
    }

    pub(super) struct ContributionSummary {
      pub(super) entries: usize,
      pub(super) distinct_coefficients: usize,
      pub(super) multiplicity: f64,
      pub(super) sites_with_data: f64,
    }
  }
}
