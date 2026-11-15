use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR};
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::fitch::sub::{
  discover_fixed_disagreements_backward, finalize_sequence_forward, resolve_nonroot_substitutions_forward,
  resolve_root_forward, resolve_variable_positions_backward,
};
use crate::partition::storage::sparse::{FitchNodeData, FitchSeqDistribution, FitchSeqInfo, SparseEdgeObs};
use crate::seq::alignment::NodeSeqInput;
use crate::seq::alignment::get_common_length_of_node_inputs;
use crate::seq::composition::Composition;
use crate::seq::indel::{compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use crate::{make_internal_report, make_report};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::keymap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::{NodeEdgeMaps, PerEdge};
use treetime_graph::scope::Inheritance;
use treetime_primitives::seq;
use treetime_utils::interval::range_union::range_union;

pub fn create_fitch_partition(
  graph: &Graph,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  create_inherited_fitch_partition(graph, Inheritance::ALL, alphabet, node_inputs)
}

pub fn create_inherited_fitch_partition(
  graph: &Graph,
  inheritance: Inheritance,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let length = get_common_length_of_node_inputs(&node_inputs)?;
  let partition = PartitionFitch {
    alphabet,
    length,
    nodes: keymap! {},
    edges: keymap! {},
    inheritance,
  };
  compress_sequences(graph, partition, node_inputs)
}

pub(crate) fn compress_sequences(
  graph: &Graph,
  partition: PartitionFitch,
  node_inputs: NodeMap<NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let partition = attach_seqs_to_graph(graph, partition, node_inputs)?;
  let partition = fitch_backward(graph, partition)?;
  let partition = fitch_forward(graph, partition)?;
  Ok(fitch_cleanup(graph, partition))
}

pub(crate) fn attach_seqs_to_graph(
  graph: &Graph,
  mut partition: PartitionFitch,
  mut node_inputs: NodeMap<NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let inheritance = &partition.inheritance;
  let leaf_records = schedule
    .leaves()
    .iter()
    .filter(|&&leaf_key| inheritance.holds(&schedule, leaf_key))
    .map(|&leaf_key| -> Result<_, Report> {
      let node = node_inputs
        .get_mut(leaf_key)
        .ok_or_else(|| make_internal_report!("Leaf {leaf_key} has no sequence input"))?;
      let seq = node
        .seq
        .take()
        .ok_or_else(|| make_report!("Leaf sequence not found after alignment completion: '{}'", node.name))?;
      Ok((leaf_key, seq))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let alphabet = partition.alphabet.clone();
  let nodes = leaf_records
    .into_par_iter()
    .map(|(leaf_key, seq)| (leaf_key, FitchNodeData::new(seq, &alphabet)))
    .collect::<BTreeMap<_, _>>();
  partition.nodes.extend(nodes);

  let edges = graph
    .get_edges()
    .map(|edge| edge.key())
    .filter(|&edge_key| partition.inheritance.carries(edge_key))
    .collect_vec();
  for edge_key in edges {
    partition.edges.insert(edge_key, SparseEdgeObs::default());
  }

  Ok(partition)
}

pub(crate) fn fitch_backward(graph: &Graph, partition: PartitionFitch) -> Result<PartitionFitch, Report> {
  let PartitionFitch {
    alphabet,
    length,
    mut nodes,
    edges,
    inheritance,
  } = partition;
  let schedule = GraphSchedule::new(graph)?;
  inheritance.require_one_parent(&schedule)?;
  if let Some(key) = nodes.keys().find(|&key| !schedule.contains(key)) {
    return Err(make_internal_report!(
      "Partition node {key} is not in the tree before the Fitch backward pass"
    ));
  }
  let inputs = schedule.preorder().iter().map(|&key| (key, nodes.remove(key)));
  let pass = FitchBackwardPass {
    schedule: &schedule,
    alphabet: &alphabet,
    length,
    edges: &edges,
    inheritance: &inheritance,
  };
  let nodes = schedule
    .scan_up(inputs, |key, leaf, children| pass.node(key, leaf, &children))?
    .unwrap_present()
    .into_map()?;
  Ok(PartitionFitch {
    alphabet,
    length,
    nodes,
    edges,
    inheritance,
  })
}

struct FitchBackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  alphabet: &'a Alphabet,
  length: usize,
  edges: &'a EdgeMap<SparseEdgeObs>,
  inheritance: &'a Inheritance,
}

type FitchBackwardOutput = Option<Arc<FitchNodeData>>;

impl FitchBackwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    leaf: Option<FitchNodeData>,
    children: &[FitchBackwardOutput],
  ) -> Result<(FitchBackwardOutput, FitchBackwardOutput), Report> {
    if !self.inheritance.holds(self.schedule, key) {
      return Ok((None, None));
    }
    let children = self
      .inheritance
      .child_messages(self.schedule, key, children)?
      .into_iter()
      .map(|(_, edge_key, child)| (&child.seq, &self.edges[&edge_key]))
      .collect_vec();
    let node = if children.is_empty() {
      leaf.ok_or_else(|| make_internal_report!("Leaf {key} has no sequence before the Fitch backward pass"))?
    } else {
      fitch_backward_internal_node(self.alphabet, self.length, &children)
    };
    let node = Arc::new(node);
    Ok((Some(Arc::clone(&node)), Some(node)))
  }
}

fn fitch_backward_internal_node(
  alphabet: &Alphabet,
  length: usize,
  children: &[(&FitchSeqInfo, &SparseEdgeObs)],
) -> FitchNodeData {
  let child_non_chars: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.non_char).collect_vec();
  let child_gaps: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.gaps).collect_vec();

  let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
  let unknown = ranges.unknown;

  let child_unknown: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.unknown).collect_vec();
  let child_variable_indels: Vec<&_> = children.iter().map(|(c, _)| &c.fitch.variable_indel).collect_vec();

  let indels_bw = resolve_indels_backward(&child_gaps, &child_unknown, &child_variable_indels, length);

  let non_char = range_union(&[ranges.non_char, indels_bw.resolved_gaps.clone()]);

  let mut sequence = seq![FILL_CHAR; length];
  for r in &non_char {
    sequence[r.0..r.1].fill(NON_CHAR);
  }

  let discovered = discover_fixed_disagreements_backward(children, alphabet, &mut sequence);
  let variable = resolve_variable_positions_backward(children, &discovered, &non_char, &mut sequence);

  FitchNodeData {
    seq: FitchSeqInfo {
      gaps: indels_bw.resolved_gaps,
      unknown,
      non_char,
      fitch: FitchSeqDistribution {
        variable,
        variable_indel: indels_bw.variable_indel,
        chosen_state: btreemap! {},
      },
      sequence,
      composition: Composition::new(alphabet),
    },
  }
}

pub(crate) fn fitch_forward(graph: &Graph, partition: PartitionFitch) -> Result<PartitionFitch, Report> {
  let PartitionFitch {
    alphabet,
    length,
    nodes,
    edges,
    inheritance,
  } = partition;
  let schedule = GraphSchedule::new(graph)?;
  inheritance.require_one_parent(&schedule)?;
  let outputs = schedule.scan_down(
    inheritance.join_parent_edge(&schedule, nodes, edges),
    |key, (node, edge), parents: PerEdge<Option<FitchForwardOutput>>| {
      let Some(node) = node else {
        return Ok((None, None));
      };
      let parent = inheritance.parent_message(&schedule, key, &parents)?;
      let (output, message) = fitch_forward_node(&alphabet, key, node, edge, parent.map(|(_, _, parent)| &parent.0))?;
      Ok((Some(output), Some(message)))
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = inheritance.split_parent_edge(&schedule, outputs.unwrap_present().into_map()?)?;
  Ok(PartitionFitch {
    alphabet,
    length,
    nodes,
    edges,
    inheritance,
  })
}

fn fitch_forward_node(
  alphabet: &Alphabet,
  key: GraphNodeKey,
  mut node: FitchNodeData,
  edge: Option<SparseEdgeObs>,
  parent: Option<&FitchNodeData>,
) -> Result<(FitchForwardOutput, FitchForwardOutput), Report> {
  let edge = match (edge, parent) {
    (Some(mut edge), Some(parent)) => {
      let parent = &parent.seq;
      let seq = &mut node.seq;
      seq.composition = parent.composition.clone();

      for r in &seq.non_char {
        seq.sequence[r.0..r.1].clone_from_slice(&parent.sequence[r.0..r.1]);
      }

      let subs = resolve_nonroot_substitutions_forward(
        &mut seq.sequence,
        &seq.gaps,
        &mut seq.fitch.variable,
        &mut seq.fitch.chosen_state,
        &mut seq.composition,
        parent,
        alphabet,
      )?;

      let (indels, new_gaps) = resolve_indels_forward(
        &seq.fitch.variable_indel,
        &seq.gaps,
        &seq.non_char,
        &parent.gaps,
        &parent.sequence,
        &seq.sequence,
      );
      seq.gaps = new_gaps;
      seq.non_char = range_union(&[seq.non_char.clone(), seq.gaps.clone()]);
      for indel in &indels {
        seq.composition.add_indel(indel);
      }
      for r in &seq.unknown {
        for pos in r.0..r.1 {
          seq.composition.adjust_count(seq.sequence[pos], -1);
        }
        seq
          .composition
          .adjust_count(alphabet.unknown(), isize::try_from(r.1 - r.0)?);
      }

      edge.extend_fitch_subs(subs);
      edge.indels.extend(indels);
      Some(edge)
    },
    (None, None) => {
      let seq = &mut node.seq;
      resolve_root_forward(
        &mut seq.sequence,
        &seq.fitch.variable,
        &mut seq.fitch.chosen_state,
        alphabet,
      );
      None
    },
    (Some(_), None) | (None, Some(_)) => {
      return Err(make_internal_report!(
        "Fitch forward pass received an inconsistent parent edge and parent message at node {key}"
      ));
    },
  };

  let seq = &mut node.seq;
  finalize_sequence_forward(
    &mut seq.sequence,
    &seq.gaps,
    &seq.unknown,
    &mut seq.composition,
    alphabet,
    edge.is_none(),
  );

  let output = Arc::new((node, edge.into_iter().collect()));
  Ok((Arc::clone(&output), output))
}

type FitchForwardOutput = Arc<(FitchNodeData, PerEdge<SparseEdgeObs>)>;

fn fitch_cleanup(graph: &Graph, mut partition: PartitionFitch) -> PartitionFitch {
  for (key, node) in &mut partition.nodes {
    if !graph.is_leaf(key) {
      node.seq.fitch.variable = btreemap! {};
    }
  }
  partition
}
