#[cfg(test)]
pub(crate) mod tests {
  pub(crate) mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
    use crate::partition::storage::sparse::SparseNodeState;
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::{EdgeMap, NodeMap};
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};

    pub(crate) const TREE: &str = "((A:0.1,B:0.25)AB:0.05,(C:0.15,(D:0.3,E:0.12)DE:0.08)CDE:0.2)root;";

    pub(crate) const GAP_FREE: [(&str, &str); 5] = [
      ("A", "ACGTACGTTGCAACGTAAGGCCTT"),
      ("B", "ACGTACGATGCAACTTAAGGCCTA"),
      ("C", "ACCTACGTTGGAACGTAAGCCCTT"),
      ("D", "TCCTACGTTGGAACGAAAGCCGTT"),
      ("E", "TCCTAGGTTGGAACGAAAGCCGTC"),
    ];

    pub(crate) const LONG_BRANCHES: f64 = 6.0;

    pub(crate) const MISSING_DATA: [(&str, &str); 5] = [
      ("A", "ACGTACGTTGCAACGTAAGGCCTT"),
      ("B", "ACGTNNNATGCAACTTAA--CCTA"),
      ("C", "ACCTACGTTGRAACGTAAGCCCTT"),
      ("D", "TCCTAC---GGAACGAAAGCCGTT"),
      ("E", "TCCTAGGTTGGAAYGAAAGCCGNN"),
    ];

    pub(crate) struct Fixture {
      pub(crate) graph: Graph,
      pub(crate) branch_lengths: EdgeMap<f64>,
      pub(crate) dense: PartitionMarginalDense,
      pub(crate) sparse: PartitionMarginalSparse,
      pub(crate) sparse_nodes: NodeMap<SparseNodeState>,
    }

    pub(crate) fn fixture(alignment: &[(&str, &str)]) -> Result<Fixture, Report> {
      scaled_fixture(alignment, 1.0)
    }

    pub(crate) fn scaled_fixture(alignment: &[(&str, &str)], scale: f64) -> Result<Fixture, Report> {
      let parsed = nwk_read(TREE.as_bytes())?;
      let graph = parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&parsed.branch_lengths)
        .into_iter()
        .map(|(key, length)| (key, length * scale))
        .collect();
      let records = alignment
        .iter()
        .map(|(name, seq)| {
          Ok(AlignmentRecord {
            name: (*name).to_owned(),
            seq: Seq::try_from_str(seq)?,
          })
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let inputs = leaf_seq_inputs(&graph, &parsed.names, records);
      let dense = PartitionMarginalDense::new(Alphabet::default(), &graph, &inputs)?;
      let fitch = create_fitch_partition(&graph, Alphabet::default(), inputs)?;
      let (sparse, sparse_nodes) = fitch.into_marginal_sparse(&graph)?;
      Ok(Fixture {
        graph,
        branch_lengths,
        dense,
        sparse: sparse.with_site_scales(true),
        sparse_nodes,
      })
    }
  }
}
