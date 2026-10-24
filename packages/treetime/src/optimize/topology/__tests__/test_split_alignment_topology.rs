#[cfg(test)]
mod tests {
  use crate::gtr::jc_distance::jukes_cantor_distance;
  use crate::optimize::topology::collapse::collapse_edge;
  use crate::optimize::topology::hoist_reversions::hoist_reverting_child;
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::partition::marginal::sparse::alignment_indels::sparse_alignment_groups;
  use crate::reroot::find_root::EdgePoint;
  use crate::reroot::orchestrate::{RerootTopologyParams, apply_root_at_edge};
  use crate::reroot::placement::Placement;
  use crate::seq::indel::InDel;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;
  use helpers::{SplitCase, alignment_indels, del, ins, inverted};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::reroot::{record_merge, remove_node_if_trivial};
  use treetime_graph::split_fraction::SplitFraction;
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
      side.finish()?;
    }

    let one_event = jukes_cantor_distance(1.0 / 9_usize.to_f64_exact(), 4);
    assert_eq!(case.whole.snapshot()?, case.split.snapshot()?);
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
    let deletion_on_a = vec![del(3, 6, "TAC")];
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
    let parent = ins(2, 3, "T");
    let child = ins(3, 4, "G");
    case.whole.set_pieces("root", "X", &[(0, vec![parent.clone()])]);
    case.whole.set_pieces("X", "A", &[(0, vec![child.clone()])]);
    case.split.set_pieces("root", "X", &[(2, vec![ins(0, 1, "T")])]);
    case.split.set_pieces("X", "A", &[(0, vec![ins(1, 2, "G")])]);
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
      side.finish()?;
    }

    assert_eq!(case.whole.snapshot()?, case.split.snapshot()?);
    assert_eq!(vec![parent, child], alignment_indels(&case.whole, "root", "A")?);
    Ok(())
  }

  #[test]
  fn test_reroot_of_a_codon_split_alignment_inverts_and_merges_the_events_of_the_unsplit_alignment()
  -> Result<(), Report> {
    let records = [
      ("A", "-CG--CGTA"),
      ("B", "-CG--CGTA"),
      ("C", "ACGT-CGTA"),
      ("D", "ACGTACG--"),
      ("E", "ACGTACG--"),
    ];
    let mut case = SplitCase::new("(((A:0.1,B:0.1)X:0.1,C:0.1)Z:0.1,(D:0.1,E:0.1)Y:0.1)root;", &records)?;
    let deletions_on_x = vec![del(0, 1, "A"), del(3, 4, "T")];
    let deletion_on_z = del(4, 5, "A");
    let deletion_on_y = del(7, 9, "TA");
    assert_eq!(
      (
        vec![del(0, 2, "AT")],
        deletions_on_x.clone(),
        vec![del(1, 2, "A")],
        vec![deletion_on_z.clone()],
        vec![deletion_on_y.clone()]
      ),
      (
        case.split.pieces("Z", "X", 0),
        alignment_indels(&case.split, "Z", "X")?,
        case.split.pieces("root", "Z", 1),
        alignment_indels(&case.split, "root", "Z")?,
        alignment_indels(&case.split, "root", "Y")?
      )
    );
    for side in case.sides_mut() {
      let edge_key = find_edge_key(&side.graph, &side.names, "X", "A").expect("the tree has the edge");
      let changes = apply_root_at_edge(
        &mut side.graph,
        Some(EdgePoint {
          edge: edge_key,
          split: SplitFraction::MIDDLE,
        }),
        Placement::EdgeEndBand,
        RerootTopologyParams::default(),
        &mut side.branch_lengths,
        None,
        |_| true,
        &mut side.names,
      )?;
      assert!(changes.edge_split.is_some() && changes.edge_merge.is_some() && !changes.inverted_edge_keys.is_empty());
      side.apply(|_, partitions| partitions.apply_reroot(&changes))?;
    }

    assert_eq!(case.whole.snapshot()?, case.split.snapshot()?);
    let insertions_on_z = deletions_on_x.into_iter().map(inverted).collect::<Vec<_>>();
    assert_eq!(
      (insertions_on_z, vec![inverted(deletion_on_z), deletion_on_y]),
      (
        alignment_indels(&case.whole, "X", "Z")?,
        alignment_indels(&case.whole, "Z", "Y")?
      )
    );
    Ok(())
  }

  #[test]
  fn test_trivial_root_merge_of_a_codon_split_alignment_composes_the_events_of_the_unsplit_alignment()
  -> Result<(), Report> {
    let records = [
      ("A", "-CGTA-GTA"),
      ("B", "-CGTA-GTA"),
      ("C", "ACG-ACGTA"),
      ("D", "ACG-ACGTA"),
    ];
    let mut case = SplitCase::new("((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)Y:0.1)root;", &records)?;
    let deletion_on_y = del(3, 4, "T");
    assert_eq!(
      (
        vec![del(0, 1, "A"), del(5, 6, "C")],
        vec![deletion_on_y.clone()],
        vec![del(0, 1, "A")],
        vec![del(1, 2, "C")],
        vec![del(1, 2, "T")]
      ),
      (
        alignment_indels(&case.split, "root", "X")?,
        alignment_indels(&case.split, "root", "Y")?,
        case.split.pieces("root", "X", 0),
        case.split.pieces("root", "X", 2),
        case.split.pieces("root", "Y", 0)
      )
    );
    let before = [case.whole.parsimony_sequences()?, case.split.parsimony_sequences()?];
    for side in case.sides_mut() {
      let edge_key = find_edge_key(&side.graph, &side.names, "root", "X").expect("the tree has the edge");
      let changes = apply_root_at_edge(
        &mut side.graph,
        Some(EdgePoint {
          edge: edge_key,
          split: SplitFraction::CHILD_END,
        }),
        Placement::EdgeEndBand,
        RerootTopologyParams::default(),
        &mut side.branch_lengths,
        None,
        |_| true,
        &mut side.names,
      )?;
      assert!(changes.edge_merge.is_some() && changes.inverted_edge_keys.len() == 1);
      side.apply(|_, partitions| partitions.apply_reroot(&changes))?;
    }

    assert_eq!(case.whole.snapshot()?, case.split.snapshot()?);
    let [mut whole_before, mut split_before] = before;
    whole_before.remove("root");
    split_before.remove("root");
    assert_eq!(
      (whole_before, split_before),
      (case.whole.parsimony_sequences()?, case.split.parsimony_sequences()?)
    );
    let inverted_x_then_y = vec![inverted(del(0, 1, "A")), deletion_on_y, inverted(del(5, 6, "C"))];
    assert_eq!(inverted_x_then_y, alignment_indels(&case.whole, "X", "Y")?);
    Ok(())
  }

  #[test]
  fn test_edge_merge_of_a_codon_split_alignment_composes_the_events_of_the_unsplit_alignment() -> Result<(), Report> {
    let records = [("A", "-CG-A-GTA"), ("B", "ACGTACGTA"), ("C", "ACGTACGTA")];
    let mut case = SplitCase::new("((A:0.1)U:0.1,B:0.1,C:0.1)root;", &records)?;
    let deletions_on_a = vec![del(0, 1, "A"), del(3, 4, "T"), del(5, 6, "C")];
    assert_eq!(deletions_on_a, alignment_indels(&case.split, "root", "U")?);
    case
      .whole
      .set_pieces("root", "U", &[(0, vec![del(0, 1, "A"), del(5, 6, "C")])]);
    case.whole.set_pieces("U", "A", &[(0, vec![del(3, 4, "T")])]);
    case
      .split
      .set_pieces("root", "U", &[(0, vec![del(0, 1, "A")]), (2, vec![del(1, 2, "C")])]);
    case.split.set_pieces("U", "A", &[(0, vec![del(1, 2, "T")])]);
    let before = [case.whole.observe("U")?, case.split.observe("U")?];
    for side in case.sides_mut() {
      let node_key = find_node_key_by_name(&side.graph, &side.names, "U").expect("the tree has the node");
      let merge = remove_node_if_trivial(&mut side.graph, node_key, Some(0.1), Some(0.1))?.expect("U has one child");
      record_merge(&mut side.branch_lengths, &merge);
      side.apply(|graph, partitions| partitions.reconcile_topology(graph, &[merge]))?;
    }

    assert_eq!(case.whole.snapshot()?, case.split.snapshot()?);
    let [mut whole_before, mut split_before] = before;
    whole_before.remove("U");
    split_before.remove("U");
    assert_eq!(
      (whole_before, split_before),
      (case.whole.parsimony_sequences()?, case.split.parsimony_sequences()?)
    );
    assert_eq!(deletions_on_a, alignment_indels(&case.whole, "root", "A")?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::deletions_adjacent_only_in_a_partition(
    (vec![del(0, 1, "A")],  vec![(0, vec![del(0, 1, "A")])]),
    (vec![del(3, 4, "T")],  vec![(0, vec![del(1, 2, "T")])]),
    (vec![del(0, 1, "A")], vec![], vec![del(3, 4, "T")]),
  )]
  #[case::adjacent_insertions_of_two_partitions(
    (vec![ins(2, 3, "G")], vec![(2, vec![ins(0, 1, "G")])]),
    (vec![ins(3, 4, "T")], vec![(0, vec![ins(1, 2, "T")])]),
    (vec![], vec![ins(2, 3, "G")], vec![ins(2, 3, "G"), ins(3, 4, "T")]),
  )]
  #[trace]
  fn test_reversion_hoist_of_a_codon_split_alignment_tests_interaction_on_alignment_positions(
    #[case] parent: (Vec<InDel>, Vec<(usize, Vec<InDel>)>),
    #[case] child: (Vec<InDel>, Vec<(usize, Vec<InDel>)>),
    #[case] expected: (Vec<InDel>, Vec<InDel>, Vec<InDel>),
  ) -> Result<(), Report> {
    let records = [
      ("A", "ACGTACGTA"),
      ("B", "AGGTACGTA"),
      ("E", "AGGTACGTA"),
      ("C", "ACGTACGTA"),
      ("D", "ACGTACGTA"),
    ];
    let mut case = SplitCase::new("((A:0.1,B:0.1,E:0.1)V:0.2,C:0.1,D:0.1)U;", &records)?;
    case.whole.set_pieces("U", "V", &[(0, parent.0)]);
    case.whole.set_pieces("V", "A", &[(0, child.0)]);
    case.split.set_pieces("U", "V", &parent.1);
    case.split.set_pieces("V", "A", &child.1);
    for side in case.sides_mut() {
      let edge = |source, target| find_edge_key(&side.graph, &side.names, source, target).expect("the tree has the edge");
      let (parent_edge_key, child_edge_key) = (edge("U", "V"), edge("V", "A"));
      let groups = sparse_alignment_groups(side.partitions.sparse_mut())?;
      hoist_reverting_child(
        &mut side.graph,
        &mut side.names,
        side.partitions.sparse_mut(),
        &groups,
        parent_edge_key,
        child_edge_key,
        &mut side.branch_lengths,
      )?;
      side.finish()?;
    }

    assert_eq!(case.whole.snapshot()?, case.split.snapshot()?);
    let hoisted = case.whole.parent_name("V");
    assert_eq!(
      expected,
      (
        alignment_indels(&case.whole, "U", &hoisted)?,
        alignment_indels(&case.whole, &hoisted, "V")?,
        alignment_indels(&case.whole, &hoisted, "A")?,
      )
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::gather::gather_edge_indel_counts;
    use crate::partition::assembly::AlignmentGroup;
    use crate::partition::create::Representation;
    use crate::partition::marginal::reconstruction::SparseReconstruction;
    use crate::partition::marginal::sparse::alignment_indels::{assembled_edge_indels, sparse_alignment_groups};
    use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
    use crate::partition::marginal::sparse::reconstruct::parsimony_seq;
    use crate::partition::set::PartitionSet;
    use crate::partition::storage::sparse::SparseNodeObs;
    use crate::seq::alignment::SequenceSource;
    use crate::seq::indel::InDel;
    use crate::seq::mutation::Sub;
    use crate::test_utils::{
      codon_partition_set, deletion, find_edge_key, find_node_key_by_name, insertion, leaf_seq_inputs, sequence_source,
      whole_alignment_set,
    };
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use std::mem;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
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
        for (partition, indels) in pieces {
          self.sparse()[*partition]
            .partition
            .obs_edges
            .entry(edge_key)
            .or_default()
            .indels
            .clone_from(indels);
        }
      }

      pub(super) fn pieces(&self, source: &str, target: &str, partition: usize) -> Vec<InDel> {
        let edge_key = find_edge_key(&self.graph, &self.names, source, target).expect("the tree has the edge");
        self.sparse_ref()[partition].partition.obs_edges[&edge_key]
          .indels
          .clone()
      }

      pub(super) fn parsimony_sequences(&self) -> Result<BTreeMap<String, Seq>, Report> {
        let partitions = self.sparse_ref();
        let group = only_group(partitions)?;
        let members = group.members(partitions);
        let per_partition: Vec<_> = members
          .iter()
          .map(|member| partition_parsimony_sequences(&self.graph, &member.partition))
          .try_collect()?;
        self
          .graph
          .get_nodes()
          .map(|node| {
            let sequences = per_partition
              .iter()
              .map(|sequences| sequences[&node.key()].clone())
              .collect_vec();
            Ok((
              self.names[node.key()].to_owned(),
              group.layout.assemble_sequence(&sequences)?,
            ))
          })
          .collect()
      }

      pub(super) fn observe(&mut self, name: &str) -> Result<BTreeMap<String, Seq>, Report> {
        let key = find_node_key_by_name(&self.graph, &self.names, name).expect("the tree has the node");
        for reconstruction in self.partitions.sparse_mut() {
          let partition = &mut reconstruction.partition;
          let sequence = partition_parsimony_sequences(&self.graph, partition)?
            .remove(&key)
            .expect("the node exists");
          partition
            .obs_nodes
            .insert(key, SparseNodeObs::new(&sequence, &partition.alphabet));
        }
        self.parsimony_sequences()
      }

      pub(super) fn parent_name(&self, name: &str) -> String {
        let key = find_node_key_by_name(&self.graph, &self.names, name).expect("the tree has the node");
        let parent = self.graph.get_node(key).expect("the node exists").inbound()[0];
        self.names[self.graph.get_source_node_key(parent).expect("the edge has a source")].to_owned()
      }

      pub(super) fn branch_length(&self, source: &str, target: &str) -> Option<f64> {
        self.branch_lengths[&find_edge_key(&self.graph, &self.names, source, target).expect("the tree has the edge")]
      }

      pub(super) fn snapshot(&self) -> Result<Snapshot, Report> {
        let partitions = self.sparse_ref();
        let group = only_group(partitions)?;
        let members = group.members(partitions);
        let root_sequence = group.layout.assemble_sequence(
          &members
            .iter()
            .map(|member| member.partition.root_sequence.clone())
            .collect_vec(),
        )?;
        let counts = gather_edge_indel_counts(&self.graph, &self.partitions)?;
        let edges = self
          .graph
          .get_edges()
          .map(|edge| {
            let key = edge.key();
            let substitutions = members
              .iter()
              .map(|member| {
                member
                  .partition
                  .obs_edges
                  .get(&key)
                  .map_or_else(Vec::new, |edge| edge.fitch_subs().to_vec())
              })
              .collect_vec();
            Ok(EdgeSnapshot {
              source: self.names[edge.source()].to_owned(),
              target: self.names[edge.target()].to_owned(),
              branch_length: self.branch_lengths[&key],
              substitutions: group.layout.assemble_substitutions(&substitutions)?,
              indels: assembled_edge_indels(&group, partitions, key)?,
              indel_count: counts[&key],
            })
          })
          .collect::<Result<Vec<_>, Report>>()?
          .into_iter()
          .sorted_by(|a, b| (&a.source, &a.target).cmp(&(&b.source, &b.target)))
          .collect();
        let (updated, _) = self
          .partitions
          .clone()
          .marginal_update(&self.graph, &branch_lengths_or_zero(&self.branch_lengths))?;
        let reconstructions = updated.iter().collect_vec();
        let sequences = self
          .graph
          .get_nodes()
          .map(|node| {
            let sequences: Vec<Seq> = group
              .members(&reconstructions)
              .iter()
              .map(|reconstruction| reconstruction.node_sequence(&self.graph, false, node.key()))
              .try_collect()?;
            Ok((
              self.names[node.key()].to_owned(),
              group.layout.assemble_sequence(&sequences)?,
            ))
          })
          .collect::<Result<BTreeMap<_, _>, Report>>()?;
        Ok(Snapshot {
          root_sequence,
          edges,
          sequences,
        })
      }

      pub(super) fn apply(
        &mut self,
        operation: impl FnOnce(&Graph, PartitionSet) -> Result<PartitionSet, Report>,
      ) -> Result<(), Report> {
        let partitions = mem::replace(&mut self.partitions, PartitionSet::Sparse(vec![]));
        self.partitions = operation(&self.graph, partitions)?;
        Ok(())
      }

      pub(super) fn finish(&mut self) -> Result<(), Report> {
        self.graph.build()?;
        self.apply(|graph, partitions| Ok(partitions.reconcile_node_states(graph)))
      }

      fn sparse(&mut self) -> &mut [SparseReconstruction] {
        self.partitions.sparse_mut()
      }

      fn sparse_ref(&self) -> &[SparseReconstruction] {
        let PartitionSet::Sparse(partitions) = &self.partitions else {
          panic!("the case builds sparse partitions");
        };
        partitions
      }
    }

    #[derive(Debug, PartialEq)]
    pub(super) struct Snapshot {
      root_sequence: Seq,
      edges: Vec<EdgeSnapshot>,
      sequences: BTreeMap<String, Seq>,
    }

    #[derive(Debug, PartialEq)]
    struct EdgeSnapshot {
      source: String,
      target: String,
      branch_length: Option<f64>,
      substitutions: Vec<Sub>,
      indels: Vec<InDel>,
      indel_count: usize,
    }

    pub(super) fn alignment_indels(side: &Side, source: &str, target: &str) -> Result<Vec<InDel>, Report> {
      let partitions = side.sparse_ref();
      let group = only_group(partitions)?;
      let edge_key = find_edge_key(&side.graph, &side.names, source, target).expect("the tree has the edge");
      assembled_edge_indels(&group, partitions, edge_key)
    }

    pub(super) fn del(lo: usize, hi: usize, chars: &str) -> InDel {
      deletion((lo, hi), Seq::try_from_str(chars).expect("the fixture is a sequence"))
    }

    pub(super) fn ins(lo: usize, hi: usize, chars: &str) -> InDel {
      insertion((lo, hi), Seq::try_from_str(chars).expect("the fixture is a sequence"))
    }

    pub(super) fn inverted(mut indel: InDel) -> InDel {
      indel.invert();
      indel
    }

    fn partition_parsimony_sequences(
      graph: &Graph,
      partition: &PartitionMarginalSparse,
    ) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
      let schedule = GraphSchedule::new(graph)?;
      let mut sequences = BTreeMap::new();
      for &key in schedule.preorder() {
        let sequence = match schedule.parents(key) {
          [] => partition.root_sequence.clone(),
          [(parent_key, edge_key)] => parsimony_seq(
            &sequences[parent_key],
            &partition.obs_edges[edge_key],
            &partition.obs_nodes[&key],
            &partition.alphabet,
          ),
          _ => panic!("the case is a tree"),
        };
        sequences.insert(key, sequence);
      }
      Ok(sequences)
    }

    fn only_group(partitions: &[SparseReconstruction]) -> Result<AlignmentGroup, Report> {
      let [group] = <[_; 1]>::try_from(sparse_alignment_groups(partitions)?).expect("the case has one alignment");
      Ok(group)
    }
  }
}
