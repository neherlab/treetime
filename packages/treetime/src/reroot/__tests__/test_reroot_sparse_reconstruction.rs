#[cfg(test)]
mod tests {
  use crate::test_utils::find_edge_key;
  use eyre::Report;
  use helpers::{Case, edge_subs, node_sequences};
  use itertools::Itertools;
  use ordered_float::OrderedFloat;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::{o, vec_of_owned};

  const FOUR_TAXA_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
  const FOUR_TAXA_ALIGNMENT: &[(&str, &str)] = &[
    ("A", "ACATCGCCNNA--GAC"),
    ("B", "GCATCCCTGTA-NG--"),
    ("C", "CCGGCGATGTRTTG--"),
    ("D", "TCGGCCGTGTRTTG--"),
  ];

  const SIX_TAXA_ALIGNMENT: &[(&str, &str)] = &[
    ("A", "TC--ACGTAN"),
    ("B", "TC--ACGTAC"),
    ("C", "GCGTACGTAC"),
    ("D", "GCGTACGAAC"),
    ("E", "ACGT--GTAC"),
    ("F", "ACGT--GTCC"),
  ];

  #[rustfmt::skip]
  #[rstest]
  #[case::path_of_several_edges(("AB", "A"),    0.5)]
  #[case::path_of_one_edge(     ("root", "AB"), 1.0)]
  #[trace]
  fn test_reroot_sparse_reconstruction_inverts_the_old_root_edge_that_the_merge_consumes(
    #[case] (source, target): (&str, &str),
    #[case] split: f64,
  ) -> Result<(), Report> {
    let mut case = Case::new(FOUR_TAXA_TREE, FOUR_TAXA_ALIGNMENT)?;

    let changes = case.reroot(source, target, split)?;

    assert!(changes.edge_merge.is_some(), "the old root is merged away");
    let partition = &case.sparse().partition;
    let merged_edge = &partition.obs_edges[&find_edge_key(&case.graph, &case.names, "AB", "CD").expect("merged edge")];
    assert_eq!(
      (
        o!("ACATCCCTGTA--G--"),
        vec_of_owned!["A1C", "A3G", "C7A", "T4G"],
        vec_of_owned!["11--13: -- -> TT"],
      ),
      (
        partition.root_sequence.as_str().to_owned(),
        edge_subs(&case, "AB", "CD").into_iter().sorted().collect_vec(),
        merged_edge.indels.iter().map(ToString::to_string).collect_vec(),
      )
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::bifurcating_old_root(  "(((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)W:0.1)Z:0.1,(E:0.1,F:0.1)Y:0.1)root;")]
  #[case::multifurcating_old_root("(((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)W:0.1)Z:0.1,E:0.1,F:0.1)root;")]
  #[trace]
  fn test_reroot_sparse_reconstruction_keeps_the_sequence_of_every_surviving_node(
    #[case] newick: &str,
  ) -> Result<(), Report> {
    let before = Case::new(newick, SIX_TAXA_ALIGNMENT)?;
    assert_eq!(
      (vec_of_owned!["A1G"], vec_of_owned!["G1T"]),
      (edge_subs(&before, "root", "Z"), edge_subs(&before, "Z", "X")),
      "the first site mutates on two edges of the paths through Z"
    );
    let sequences_before = node_sequences(&before)?;
    let targets = before
      .graph
      .get_edges()
      .flat_map(|edge| {
        let (source, target) = (&before.names[edge.source()], &before.names[edge.target()]);
        let splits = if before.graph.is_leaf(edge.target()) { vec![0.5] } else { vec![0.5, 1.0] };
        splits.into_iter().map(move |split| (source.to_owned(), target.to_owned(), split))
      })
      .collect_vec();

    let old_root = before.names[before.graph.root_key()?].to_owned();
    let mut expected = BTreeMap::new();
    let mut actual = BTreeMap::new();
    for (source, target, split) in targets {
      let mut case = Case::new(newick, SIX_TAXA_ALIGNMENT)?;
      let changes = case.reroot(&source, &target, split)?;
      let mut sequences = sequences_before.clone();
      if changes.edge_merge.is_some() {
        sequences.remove(&old_root);
      }
      if let Some(info) = &changes.edge_split {
        sequences.insert(case.names[info.new_node_key].to_owned(), sequences_before[&source].clone());
      }
      let key = (source, target, OrderedFloat(split));
      expected.insert(key.clone(), sequences);
      actual.insert(key, node_sequences(&case)?);
    }
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::reconstruction::SparseReconstruction;
    use crate::partition::marginal::sparse::reconstruct::parsimony_seq;
    use crate::partition::set::PartitionSet;
    use crate::reroot::find_root::EdgePoint;
    use crate::reroot::orchestrate::{RerootTopologyParams, apply_root_at_edge};
    use crate::test_utils::{find_edge_key, leaf_seq_inputs};
    use eyre::Report;
    use std::collections::BTreeMap;
    use std::mem;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::reroot::RerootResult;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_graph::split_fraction::SplitFraction;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};

    pub(super) struct Case {
      pub(super) graph: Graph,
      pub(super) names: NodeNames,
      branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      partitions: PartitionSet,
    }

    impl Case {
      pub(super) fn new(newick: &str, records: &[(&str, &str)]) -> Result<Self, Report> {
        let parse = nwk_read(newick.as_bytes())?;
        let records = records
          .iter()
          .map(|&(name, seq)| {
            Ok(AlignmentRecord {
              name: name.to_owned(),
              seq: Seq::try_from_str(seq)?,
            })
          })
          .collect::<Result<Vec<_>, Report>>()?;
        let fitch = create_fitch_partition(
          &parse.graph,
          Alphabet::new(AlphabetName::Nuc)?,
          leaf_seq_inputs(&parse.graph, &parse.names, records),
        )?;
        let (partition, node_states) = fitch.into_marginal_sparse(&parse.graph)?;
        let reconstruction = SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);
        Ok(Self {
          graph: parse.graph,
          names: parse.names,
          branch_lengths: parse.branch_lengths,
          partitions: PartitionSet::Sparse(vec![reconstruction]),
        })
      }

      pub(super) fn reroot(&mut self, source: &str, target: &str, split: f64) -> Result<RerootResult, Report> {
        let edge_key = find_edge_key(&self.graph, &self.names, source, target).expect("the tree has the edge");
        let changes = apply_root_at_edge(
          &mut self.graph,
          Some(EdgePoint {
            edge: edge_key,
            split: SplitFraction::new(split)?,
          }),
          RerootTopologyParams::default(),
          &mut self.branch_lengths,
          None,
          |_| true,
          &mut self.names,
        )?;
        let partitions = mem::replace(&mut self.partitions, PartitionSet::Sparse(vec![]));
        self.partitions = partitions.apply_reroot(&changes)?;
        Ok(changes)
      }

      pub(super) fn sparse(&self) -> &SparseReconstruction {
        let PartitionSet::Sparse(partitions) = &self.partitions else {
          panic!("the case builds a sparse partition");
        };
        &partitions[0]
      }
    }

    pub(super) fn edge_subs(case: &Case, source: &str, target: &str) -> Vec<String> {
      let edge_key = find_edge_key(&case.graph, &case.names, source, target).expect("the tree has the edge");
      case.sparse().partition.obs_edges[&edge_key]
        .fitch_subs()
        .iter()
        .map(ToString::to_string)
        .collect()
    }

    pub(super) fn node_sequences(case: &Case) -> Result<BTreeMap<String, Seq>, Report> {
      let partition = &case.sparse().partition;
      let schedule = GraphSchedule::new(&case.graph)?;
      let mut sequences: BTreeMap<_, Seq> = BTreeMap::new();
      for &key in schedule.preorder() {
        let sequence = match schedule.parents(key) {
          [] => partition.root_sequence.clone(),
          [(parent_key, edge_key)] => parsimony_seq(
            &sequences[&case.names[*parent_key]],
            &partition.obs_edges[edge_key],
            &partition.obs_nodes[&key],
            &partition.alphabet,
          ),
          _ => panic!("the case is a tree"),
        };
        sequences.insert(case.names[key].to_owned(), sequence);
      }
      Ok(sequences)
    }
  }
}
