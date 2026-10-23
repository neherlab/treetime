#[cfg(test)]
mod tests {
  use crate::gtr::jc_distance::jukes_cantor_distance;
  use crate::optimize::topology::collapse::collapse_edge;
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::partition::marginal::sparse::alignment_indels::sparse_alignment_groups;
  use crate::test_utils::{deletion, find_edge_key, insertion};
  use eyre::Report;
  use helpers::{SplitCase, alignment_indels};
  use pretty_assertions::assert_eq;
  use treetime_primitives::Seq;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_sibling_merge_of_a_codon_split_alignment_gives_the_branch_lengths_of_the_unsplit_alignment()
  -> Result<(), Report> {
    let records = [
      ("A", "TCG---GTA"),
      ("B", "TCGTACGTA"),
      ("C", "ACGTACGTA"),
      ("D", "ACGTACGTA"),
    ];
    let mut case = SplitCase::new("(A:0.1,B:0.1,C:0.1,D:0.1)root;", &records)?;
    for side in case.sides_mut() {
      merge_shared_mutation_branches(
        &mut side.graph,
        &mut side.names,
        side.partitions.sparse_mut(),
        &mut side.branch_lengths,
      )?;
    }

    let one_event = jukes_cantor_distance(1.0 / 9_usize.to_f64_exact(), 4);
    assert_eq!(case.whole.named_branch_lengths(), case.split.named_branch_lengths());
    let merged = case.whole.parent_name("A");
    assert_eq!(
      (Some(one_event), Some(one_event), Some(0.0), merged.clone()),
      (
        case.whole.branch_length("root", &merged),
        case.whole.branch_length(&merged, "A"),
        case.whole.branch_length(&merged, "B"),
        case.whole.parent_name("B"),
      )
    );
    let deletion_on_a = vec![deletion((3, 6), Seq::try_from_str("TAC")?)];
    assert_eq!(
      (deletion_on_a.clone(), deletion_on_a),
      (
        alignment_indels(&case.whole, &merged, "A")?,
        alignment_indels(&case.split, &merged, "A")?
      )
    );
    Ok(())
  }

  #[test]
  fn test_collapse_in_a_codon_split_alignment_keeps_adjacent_insertions_of_two_edges_as_two_events()
  -> Result<(), Report> {
    let records = [
      ("A", "ACTGGTACG"),
      ("B", "ACT-GTACG"),
      ("C", "ACT-GTACG"),
      ("D", "ACT-GTACG"),
    ];
    let mut case = SplitCase::new("((A:0.1,B:0.1)X:0.1,C:0.1,D:0.1)root;", &records)?;
    let parent = insertion((2, 3), Seq::try_from_str("T")?);
    let child = insertion((3, 4), Seq::try_from_str("G")?);
    case.whole.set_pieces("root", "X", &[(0, vec![parent.clone()])]);
    case.whole.set_pieces("X", "A", &[(0, vec![child.clone()])]);
    case
      .split
      .set_pieces("root", "X", &[(2, vec![insertion((0, 1), Seq::try_from_str("T")?)])]);
    case
      .split
      .set_pieces("X", "A", &[(0, vec![insertion((1, 2), Seq::try_from_str("G")?)])]);
    for side in case.sides_mut() {
      let edge_key = find_edge_key(&side.graph, &side.names, "root", "X").expect("the tree has the edge");
      let groups = sparse_alignment_groups(side.partitions.sparse_mut())?;
      collapse_edge(
        &mut side.graph,
        side.partitions.sparse_mut(),
        &groups,
        edge_key,
        &mut side.branch_lengths,
      )?;
    }

    let expected = vec![parent, child];
    assert_eq!(
      (expected.clone(), expected),
      (
        alignment_indels(&case.whole, "root", "A")?,
        alignment_indels(&case.split, "root", "A")?
      )
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::partition::create::Representation;
    use crate::partition::marginal::sparse::alignment_indels::{assembled_edge_indels, sparse_alignment_groups};
    use crate::partition::set::PartitionSet;
    use crate::seq::alignment::SequenceSource;
    use crate::seq::indel::InDel;
    use crate::test_utils::{
      codon_partition_set, find_edge_key, find_node_key_by_name, leaf_seq_inputs, sequence_source, whole_alignment_set,
    };
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};

    pub(super) struct SplitCase {
      pub(super) whole: Side,
      pub(super) split: Side,
    }

    impl SplitCase {
      pub(super) fn new(newick: &str, records: &[(&str, &str)]) -> Result<Self, Report> {
        Ok(Self {
          whole: Side::new(newick, records, whole_alignment_set)?,
          split: Side::new(newick, records, codon_partition_set)?,
        })
      }

      pub(super) fn sides_mut(&mut self) -> [&mut Side; 2] {
        [&mut self.whole, &mut self.split]
      }
    }

    pub(super) struct Side {
      pub(super) graph: Graph,
      pub(super) names: NodeNames,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      pub(super) partitions: PartitionSet,
    }

    type SetBuilder = fn(
      Representation,
      &GtrModel,
      &Graph,
      Vec<SequenceSource>,
      &BTreeMap<GraphEdgeKey, f64>,
    ) -> Result<PartitionSet, Report>;

    impl Side {
      fn new(newick: &str, records: &[(&str, &str)], build: SetBuilder) -> Result<Self, Report> {
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
        let source = sequence_source(
          Alphabet::default(),
          leaf_seq_inputs(&parse.graph, &parse.names, records),
        )?;
        let partitions = build(
          Representation::Sparse,
          &GtrModel::JC69(JC69Params::default()),
          &parse.graph,
          vec![source],
          &branch_lengths_or_zero(&parse.branch_lengths),
        )?;
        Ok(Self {
          graph: parse.graph,
          names: parse.names,
          branch_lengths: parse.branch_lengths,
          partitions,
        })
      }

      pub(super) fn set_pieces(&mut self, source: &str, target: &str, pieces: &[(usize, Vec<InDel>)]) {
        let edge_key = find_edge_key(&self.graph, &self.names, source, target).expect("the tree has the edge");
        let PartitionSet::Sparse(partitions) = &mut self.partitions else {
          panic!("the case builds sparse partitions");
        };
        for (partition, indels) in pieces {
          partitions[*partition]
            .partition
            .obs_edges
            .entry(edge_key)
            .or_default()
            .indels
            .clone_from(indels);
        }
      }

      pub(super) fn parent_name(&self, name: &str) -> String {
        let key = find_node_key_by_name(&self.graph, &self.names, name).expect("the tree has the node");
        let parent = self.graph.get_node(key).expect("the node exists").inbound()[0];
        self.names[self.graph.get_source_node_key(parent).expect("the edge has a source")].to_owned()
      }

      pub(super) fn branch_length(&self, source: &str, target: &str) -> Option<f64> {
        self.branch_lengths[&find_edge_key(&self.graph, &self.names, source, target).expect("the tree has the edge")]
      }

      pub(super) fn named_branch_lengths(&self) -> Vec<(String, String, Option<f64>)> {
        self
          .graph
          .get_edges()
          .map(|edge| {
            (
              self.names[edge.source()].to_owned(),
              self.names[edge.target()].to_owned(),
              self.branch_lengths[&edge.key()],
            )
          })
          .sorted_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)))
          .collect()
      }
    }

    pub(super) fn alignment_indels(side: &Side, source: &str, target: &str) -> Result<Vec<InDel>, Report> {
      let PartitionSet::Sparse(partitions) = &side.partitions else {
        panic!("the case builds sparse partitions");
      };
      let [group] = <[_; 1]>::try_from(sparse_alignment_groups(partitions)?).expect("the case has one alignment");
      let edge_key = find_edge_key(&side.graph, &side.names, source, target).expect("the tree has the edge");
      assembled_edge_indels(&group, partitions, edge_key)
    }
  }
}
