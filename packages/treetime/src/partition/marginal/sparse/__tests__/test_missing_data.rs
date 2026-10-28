#[cfg(test)]
mod tests {
  use crate::partition::create::Representation;
  use eyre::Report;
  use helpers::{
    ALL_N, CLADE_GAP_RUN, CLADE_N_RUN, CLEAN, GAP_RUN, MissingDataTree, N_AT_SIBLING_SUB, N_RUN, SINGLE_N,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::clean(            CLEAN)]
  #[case::leaf_n_run(       N_RUN)]
  #[case::leaf_gap_run(     GAP_RUN)]
  #[case::leaf_all_n(       ALL_N)]
  #[case::clade_n_run(      CLADE_N_RUN)]
  #[case::clade_gap_run(    CLADE_GAP_RUN)]
  #[case::single_n(         SINGLE_N)]
  #[case::n_at_sibling_sub( N_AT_SIBLING_SUB)]
  #[trace]
  fn test_missing_data_sparse_matches_dense(#[case] seqs: [&str; 5]) -> Result<(), Report> {
    let dense = MissingDataTree::new(seqs, Representation::Dense)?;
    let sparse = MissingDataTree::new(seqs, Representation::Sparse)?;

    let expected = dense.summary()?;
    let actual = sparse.summary()?;

    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::clean(            CLEAN)]
  #[case::leaf_n_run(       N_RUN)]
  #[case::leaf_gap_run(     GAP_RUN)]
  #[case::leaf_all_n(       ALL_N)]
  #[case::clade_n_run(      CLADE_N_RUN)]
  #[case::clade_gap_run(    CLADE_GAP_RUN)]
  #[case::single_n(         SINGLE_N)]
  #[case::n_at_sibling_sub( N_AT_SIBLING_SUB)]
  #[trace]
  fn test_missing_data_node_log_lh_equals_root_log_lh(
    #[case] seqs: [&str; 5],
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
  ) -> Result<(), Report> {
    let tree = MissingDataTree::new(seqs, representation)?;

    let root_log_lh = tree.log_lh;
    let expected: BTreeMap<String, f64> = tree.node_names().map(|name| (name, root_log_lh)).collect();
    let actual = tree.node_log_lh();

    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::clean(            CLEAN)]
  #[case::leaf_n_run(       N_RUN)]
  #[case::leaf_gap_run(     GAP_RUN)]
  #[case::leaf_all_n(       ALL_N)]
  #[case::clade_n_run(      CLADE_N_RUN)]
  #[case::clade_gap_run(    CLADE_GAP_RUN)]
  #[case::single_n(         SINGLE_N)]
  #[case::n_at_sibling_sub( N_AT_SIBLING_SUB)]
  #[trace]
  fn test_missing_data_dense_internal_posteriors_stay_below_folding_threshold(
    #[case] seqs: [&str; 5],
  ) -> Result<(), Report> {
    let dense = MissingDataTree::new(seqs, Representation::Dense)?;

    let actual = dense.internal_sites_resolved_above(1.0 - 1e-4)?;

    assert_eq!(Vec::<(String, usize)>::new(), actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_n_run(       N_RUN)]
  #[case::leaf_all_n(       ALL_N)]
  #[case::single_n(         SINGLE_N)]
  #[case::n_at_sibling_sub( N_AT_SIBLING_SUB)]
  #[trace]
  fn test_missing_data_sparse_leaf_posteriors_at_n_sites_match_dense(#[case] seqs: [&str; 5]) -> Result<(), Report> {
    let dense = MissingDataTree::new(seqs, Representation::Dense)?;
    let sparse = MissingDataTree::new(seqs, Representation::Sparse)?;

    let expected = dense.leaf_posteriors_at_n_sites()?;
    let actual = sparse.leaf_posteriors_at_n_sites()?;

    pretty_assert_map_abs_diff_eq!(expected, actual, epsilon = 1e-10);
    Ok(())
  }

  mod helpers {
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::gtr::gtr::GTR;
    use crate::make_report;
    use crate::optimize::likelihood::evaluate_partitions_log_lh_only;
    use crate::partition::create::Representation;
    use crate::partition::marginal::reconstruction::MarginalReconstruction;
    use crate::partition::marginal::sparse::reconstruct::leaf_posterior;
    use crate::test_utils::{
      NUC_ALPHABET, dense_reconstruction, marginal_update_reconstruction, run_dense_marginal, run_sparse_marginal,
    };
    use eyre::Report;
    use itertools::{Itertools, izip};
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::{NwkParse, nwk_read};
    use treetime_primitives::AlignmentRecord;

    pub(super) const TREE: &str = "((A:0.3,B:0.4)AB:0.2,(C:0.35,(D:0.25,E:0.3)DE:0.2)CDE:0.25)root;";

    pub(super) const LEAF_NAMES: [&str; 5] = ["A", "B", "C", "D", "E"];

    pub(super) const BRANCH_LENGTHS: [f64; 3] = [0.01, 0.1, 1.0];

    #[rustfmt::skip]
    pub(super) const CLEAN:            [&str; 5] = ["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGAACGTACCTTCGTACGT", "ACGAACGTACGTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const N_RUN:            [&str; 5] = ["ACGTNNNNNNGTACGTACGT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGAACGTACCTTCGTACGT", "ACGAACGTACGTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const GAP_RUN:          [&str; 5] = ["ACGT------GTACGTACGT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGAACGTACCTTCGTACGT", "ACGAACGTACGTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const ALL_N:            [&str; 5] = ["NNNNNNNNNNNNNNNNNNNN", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGAACGTACCTTCGTACGT", "ACGAACGTACGTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const CLADE_N_RUN:      [&str; 5] = ["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGANNNNNNCTTCGTACGT", "ACGANNNNNNGTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const CLADE_GAP_RUN:    [&str; 5] = ["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGA------CTTCGTACGT", "ACGA------GTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const SINGLE_N:         [&str; 5] = ["ACGTACGTANGTACGTACGT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGAACGTACCTTCGTACGT", "ACGAACGTACGTTCGTACGC"];
    #[rustfmt::skip]
    pub(super) const N_AT_SIBLING_SUB: [&str; 5] = ["ACGTACGTACGTACGTACNT", "ACGTACGTACGTACGTACAA", "ACGAACGTACGTTCGTACGT", "ACGAACGTACCTTCGTACGT", "ACGAACGTACGTTCGTACGC"];

    fn gtr() -> Result<GTR, Report> {
      jc69(JC69Params::default())
    }

    pub(super) struct MissingDataTree {
      graph: Graph,
      names: NodeNames,
      branch_lengths: EdgeMap<f64>,
      reconstruction: MarginalReconstruction,
      pub(super) log_lh: f64,
    }

    impl MissingDataTree {
      pub(super) fn new(seqs: [&str; 5], representation: Representation) -> Result<Self, Report> {
        let NwkParse {
          graph,
          branch_lengths,
          names,
          ..
        } = nwk_read(TREE.as_bytes())?;
        let fasta = izip!(LEAF_NAMES, seqs)
          .map(|(name, seq)| format!(">{name}\n{seq}\n"))
          .join("");
        let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &*NUC_ALPHABET)?
          .into_iter()
          .map(AlignmentRecord::from)
          .collect();
        let (log_lh, reconstruction) = match representation {
          Representation::Dense => {
            let (log_lh, dense) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr()?)?;
            (log_lh, MarginalReconstruction::Dense(dense))
          },
          Representation::Sparse => {
            let (log_lh, sparse) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr()?)?;
            (log_lh, MarginalReconstruction::Sparse(sparse))
          },
        };
        let branch_lengths = branch_lengths
          .iter()
          .map(|(key, length)| (key, length.unwrap_or(0.0)))
          .collect();
        Ok(Self {
          graph,
          names,
          branch_lengths,
          reconstruction,
          log_lh,
        })
      }

      pub(super) fn summary(&self) -> Result<BTreeMap<String, f64>, Report> {
        let root = [(String::from("root log_lh"), self.log_lh)];
        let changes = self
          .graph
          .get_edges()
          .map(|edge| edge.key())
          .cartesian_product(BRANCH_LENGTHS)
          .map(|(edge_key, length)| -> Result<_, Report> {
            let child = &self.names[self.graph.edge_endpoints(edge_key)?.1];
            Ok([
              (
                format!("edge {child} t={length} objective change"),
                self.objective_change(edge_key, length)?,
              ),
              (
                format!("edge {child} t={length} tree change"),
                self.tree_change(edge_key, length)?,
              ),
            ])
          })
          .flatten_ok()
          .collect::<Result<Vec<_>, Report>>()?;
        let nodes = self
          .node_log_lh()
          .into_iter()
          .map(|(name, log_lh)| (format!("node {name} log_lh"), log_lh));
        Ok(root.into_iter().chain(changes).chain(nodes).collect())
      }

      pub(super) fn node_names(&self) -> impl Iterator<Item = String> + '_ {
        self.graph.get_nodes().map(|node| self.names[node.key()].to_owned())
      }

      pub(super) fn node_log_lh(&self) -> BTreeMap<String, f64> {
        match &self.reconstruction {
          MarginalReconstruction::Dense(dense) => dense
            .node_states
            .iter()
            .map(|(key, node)| (self.names[key].to_owned(), node.profile.log_lh.value()))
            .collect(),
          MarginalReconstruction::Sparse(sparse) => sparse
            .node_states
            .iter()
            .map(|(key, node)| (self.names[key].to_owned(), node.profile.log_lh.value()))
            .collect(),
        }
      }

      pub(super) fn internal_sites_resolved_above(&self, threshold: f64) -> Result<Vec<(String, usize)>, Report> {
        let dense = dense_reconstruction(&self.reconstruction);
        Ok(
          self
            .graph
            .get_internal_nodes()
            .flat_map(|node| {
              let name = &self.names[node.key()];
              dense.node_states[&node.key()]
                .profile
                .dis
                .rows()
                .into_iter()
                .enumerate()
                .filter(|(_, row)| row.iter().any(|&p| p >= threshold))
                .map(|(pos, _)| (name.to_owned(), pos))
                .collect_vec()
            })
            .collect(),
        )
      }

      pub(super) fn leaf_posteriors_at_n_sites(&self) -> Result<BTreeMap<String, Array1<f64>>, Report> {
        let unknown = NUC_ALPHABET.unknown();
        self
          .graph
          .get_leaves()
          .map(|leaf| leaf.key())
          .map(|key| -> Result<_, Report> {
            let name = &self.names[key];
            let posteriors = match &self.reconstruction {
              MarginalReconstruction::Dense(dense) => {
                let node = &dense.node_states[&key];
                node
                  .seq
                  .sequence
                  .iter()
                  .positions(|&c| c == unknown)
                  .map(|pos| (format!("{name} {pos}"), node.profile.dis.row(pos).to_owned()))
                  .collect_vec()
              },
              MarginalReconstruction::Sparse(sparse) => {
                let (parent_key, edge_key) = self
                  .graph
                  .node_parent(key)?
                  .ok_or_else(|| make_report!("leaf {name} has no parent"))?;
                let down = &sparse.edges.forward[&edge_key].msg_from_parent;
                let parent = &sparse.node_states[&parent_key];
                sparse.node_states[&key]
                  .sequence
                  .iter()
                  .positions(|&c| c == unknown)
                  .map(|pos| -> Result<_, Report> {
                    let posterior = leaf_posterior(unknown, pos, down, parent, &NUC_ALPHABET)
                      .ok_or_else(|| make_report!("leaf {name} has no posterior at site {pos}"))?;
                    Ok((format!("{name} {pos}"), posterior))
                  })
                  .collect::<Result<Vec<_>, Report>>()?
              },
            };
            Ok(posteriors)
          })
          .flatten_ok()
          .collect()
      }

      fn objective_change(&self, edge_key: GraphEdgeKey, length: f64) -> Result<f64, Report> {
        let contribution = [self.reconstruction.view().create_edge_contribution(edge_key)?];
        let at_length = evaluate_partitions_log_lh_only(&contribution, length)?;
        let at_current = evaluate_partitions_log_lh_only(&contribution, self.branch_lengths[&edge_key])?;
        Ok(at_length.value() - at_current.value())
      }

      fn tree_change(&self, edge_key: GraphEdgeKey, length: f64) -> Result<f64, Report> {
        let mut branch_lengths = self.branch_lengths.clone();
        branch_lengths.insert(edge_key, length);
        let (_, log_lh) = marginal_update_reconstruction(self.reconstruction.clone(), &self.graph, &branch_lengths)?;
        Ok(log_lh.value() - self.log_lh)
      }
    }
  }
}
