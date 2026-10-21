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
use itertools::{Itertools, izip};
use maplit::btreemap;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::{NodeEdgeMaps, PerEdge};
use treetime_primitives::{AlphabetLike, seq};
use treetime_utils::interval::range_union::range_union;

pub fn create_fitch_partition(
  graph: &Graph,
  alphabet: Alphabet,
  node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let length = get_common_length_of_node_inputs(&node_inputs)?;
  let partition = PartitionFitch {
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  };
  compress_sequences(graph, partition, node_inputs)
}

pub(crate) fn compress_sequences(
  graph: &Graph,
  partition: PartitionFitch,
  node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let partition = attach_seqs_to_graph(graph, partition, node_inputs)?;
  let partition = fitch_backward(graph, partition)?;
  let partition = fitch_forward(graph, partition)?;
  Ok(fitch_cleanup(graph, partition))
}

pub(crate) fn attach_seqs_to_graph(
  graph: &Graph,
  mut partition: PartitionFitch,
  mut node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let leaf_records = graph
    .get_leaves()
    .map(|leaf| -> Result<_, Report> {
      let leaf_key = leaf.key();
      let node = node_inputs
        .get_mut(&leaf_key)
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

  for edge in graph.get_edges() {
    let edge_key = edge.key();
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
  } = partition;
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  if let Some(key) = nodes.keys().find(|&&key| !schedule.contains(key)) {
    return Err(make_internal_report!(
      "Partition node {key} is not in the tree before the Fitch backward pass"
    ));
  }
  let inputs = schedule.preorder().iter().map(|&key| (key, nodes.remove(&key)));
  let pass = FitchBackwardPass {
    schedule: &schedule,
    alphabet: &alphabet,
    length,
    edges: &edges,
  };
  let nodes = schedule
    .scan_up(inputs, |key, leaf, children| pass.node(key, leaf, &children))?
    .unwrap_shared()
    .into_map()?;
  Ok(PartitionFitch {
    alphabet,
    length,
    nodes,
    edges,
  })
}

struct FitchBackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  alphabet: &'a Alphabet,
  length: usize,
  edges: &'a BTreeMap<GraphEdgeKey, SparseEdgeObs>,
}

impl FitchBackwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    leaf: Option<FitchNodeData>,
    children: &[Arc<FitchNodeData>],
  ) -> Result<(Arc<FitchNodeData>, Arc<FitchNodeData>), Report> {
    let node = if self.schedule.children(key).is_empty() {
      leaf.ok_or_else(|| make_internal_report!("Leaf {key} has no sequence before the Fitch backward pass"))?
    } else {
      let children = izip!(self.schedule.children(key), children)
        .map(|(&(_, edge_key), child)| (&child.seq, &self.edges[&edge_key]))
        .collect_vec();
      fitch_backward_internal_node(self.alphabet, self.length, &children)
    };
    let node = Arc::new(node);
    Ok((Arc::clone(&node), node))
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
      composition: Composition::new(alphabet.chars(), alphabet.gap()),
    },
  }
}

pub(crate) fn fitch_forward(graph: &Graph, partition: PartitionFitch) -> Result<PartitionFitch, Report> {
  let PartitionFitch {
    alphabet,
    length,
    nodes,
    edges,
  } = partition;
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let inputs = schedule.join_parent_edges(nodes, edges)?;
  let outputs = schedule.scan_down(
    inputs,
    |key, (node, mut edges), parents: PerEdge<FitchForwardOutput>| {
      fitch_forward_node(
        &alphabet,
        key,
        node,
        edges.pop(),
        parents.first().map(|parent| &parent.0),
      )
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = outputs.unwrap_shared().split_parent_edges()?;
  Ok(PartitionFitch {
    alphabet,
    length,
    nodes,
    edges,
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
    if !graph.is_leaf(*key) {
      node.seq.fitch.variable = btreemap! {};
    }
  }
  partition
}
