use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::{GraphEdgeKey, Weighted};
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::PartitionWithGtrInference;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::graph_ancestral::{EdgeAncestral, NodeAncestral};
use crate::representation::graph_sparse::{
  MarginalSparseSeqDistribution, SparseEdgePartition, SparseNodePartition, VarPos,
};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_compressed::PartitionCompressed;
use crate::representation::partition_marginal::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::seq::Seq;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;
use std::iter::zip;
use treetime_utils::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionCompressed for PartitionMarginalSparse {
  fn index(&self) -> usize {
    self.index
  }

  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn length(&self) -> usize {
    self.length
  }

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition> {
    &self.nodes
  }

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &self.edges
  }

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition> {
    &mut self.nodes
  }

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &mut self.edges
  }
}

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalSparse {}

impl crate::commands::timetree::partition_ops::PartitionTimetreeOps for PartitionMarginalSparse {
  fn create_edge_contribution(
    &self,
    edge_key: GraphEdgeKey,
  ) -> Result<crate::commands::optimize::optimize_unified::OptimizationContribution, Report> {
    crate::commands::optimize::optimize_unified::OptimizationContribution::from_sparse(edge_key, self)
  }
}

impl PartitionMarginalOps<NodeAncestral, EdgeAncestral> for PartitionMarginalSparse {
  fn attach_sequences(&mut self, _graph: &GraphAncestral, _aln: &[FastaRecord]) -> Result<(), Report> {
    // Sparse partitions get sequences attached during compression phase
    Ok(())
  }

  fn process_node_backward(
    &mut self,
    node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    let length = self.length;
    let mut seq_info = self.nodes.remove(&node.key).unwrap();

    let msg_to_parent = if node.is_leaf {
      // this is mostly a copy (or ref here) of the fitch state.
      let alphabet = &self.alphabet;
      let fixed = alphabet
        .determined()
        .map(|state| (state, alphabet.get_profile(state).clone()))
        .collect();

      // convert the parsimony variable states to sparseSeqDis format
      let variable = seq_info
        .seq
        .fitch
        .variable
        .iter()
        .map(|(pos, p)| {
          let dis = alphabet.construct_profile(p.chars()).unwrap();
          let state = p.get_one();
          (*pos, VarPos { dis, state })
        })
        .collect();

      MarginalSparseSeqDistribution {
        fixed_counts: seq_info.seq.composition.clone(),
        variable,
        variable_indel: btreemap! {},
        fixed,
        log_lh: 0.0,
      }
    } else {
      // for internal nodes, combine the messages from the children
      let mut variable_pos = btreemap! {};
      let mut child_states = vec![];
      let mut child_messages: Vec<MarginalSparseSeqDistribution> = vec![];
      for (ci, (_child_key, edge_key)) in node.child_keys.iter().enumerate() {
        child_states.push(btreemap! {});
        let edge_data = &self.edges[edge_key];
        for m in &edge_data.subs {
          variable_pos.insert(m.pos(), m.reff());
          child_states[ci].insert(m.pos(), m.qry());
        }
        for (pos, p) in &edge_data.msg_from_child.variable {
          variable_pos.entry(*pos).or_insert(p.state);
        }
        child_messages.push(edge_data.msg_from_child.clone());
      }

      let alphabet = &self.alphabet;
      for (ci, (child_key, _)) in node.child_keys.iter().enumerate() {
        let states = &mut child_states[ci];
        let child_data = &self.nodes[child_key];
        for pos in variable_pos.keys() {
          if !states.contains_key(pos) && range_contains(&child_data.seq.non_char, *pos) {
            if range_contains(&child_data.seq.gaps, *pos) {
              states.insert(*pos, alphabet.gap());
            } else {
              states.insert(*pos, alphabet.unknown());
            }
          }
        }
      }

      combine_messages(
        &seq_info.seq.composition,
        &child_messages,
        &variable_pos,
        &child_states,
        &self.alphabet,
        node.is_root.then_some(&self.gtr.pi),
      )?
    };

    if node.is_root {
      // data from children * gtr.pi as calculated above is the root profile.
      seq_info.profile = msg_to_parent;
    } else {
      // what was calculated above is what is sent to the parent. we also calculate the propagated message to the parent (we need it in the forward pass).
      let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
      let branch_length = node.parent_edges[0].weight().unwrap_or(0.0);
      let branch_length = fix_branch_length(length, branch_length);
      let mut edge_data = self.edges.remove(edge_key).unwrap();
      edge_data.msg_from_child = propagate_raw(
        &self.gtr.expQt(branch_length).t().to_owned(),
        &msg_to_parent,
        edge_data.transmission.as_ref(),
      );
      edge_data.msg_to_parent = msg_to_parent;
      self.edges.insert(*edge_key, edge_data);
    }
    self.nodes.insert(node.key, seq_info);
    Ok(())
  }

  fn process_node_forward(
    &mut self,
    graph: &GraphAncestral,
    node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  ) -> Result<(), Report> {
    if !node.is_root {
      let mut seq_info = self.nodes.remove(&node.key).unwrap();
      let mut variable_pos = btreemap! {};
      let mut ref_states: Vec<BTreeMap<usize, crate::representation::seq_char::AsciiChar>> = vec![];
      let mut msgs_to_combine: Vec<MarginalSparseSeqDistribution> = vec![];
      let mut removed_edges = vec![];
      for (_, edge_key) in &node.parent_keys {
        let mut parent_state: BTreeMap<usize, crate::representation::seq_char::AsciiChar> = btreemap! {};
        let mut child_state: BTreeMap<usize, crate::representation::seq_char::AsciiChar> = btreemap! {};
        let edge_data = self.edges.remove(edge_key).unwrap();
        removed_edges.push((*edge_key, edge_data.clone()));

        for (pos, p) in &edge_data.msg_to_child.variable {
          if !range_contains(&seq_info.seq.gaps, *pos) {
            variable_pos.entry(*pos).or_insert(p.state);
            parent_state.insert(*pos, p.state);
          }
        }
        for m in &edge_data.subs {
          variable_pos.insert(m.pos(), m.qry());
          parent_state.entry(m.pos()).or_insert_with(|| m.reff());
          child_state.entry(m.pos()).or_insert_with(|| m.qry());
        }
        for (pos, p) in &edge_data.msg_to_parent.variable {
          variable_pos.entry(*pos).or_insert(p.state);
          child_state.entry(*pos).or_insert(p.state);
        }

        let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
        let branch_length = edge_payload.branch_length.unwrap_or(0.0);

        msgs_to_combine.push(propagate_raw(
          &self.gtr.expQt(branch_length),
          &edge_data.msg_to_child,
          edge_data.transmission.as_ref(),
        ));
        msgs_to_combine.push(edge_data.msg_to_parent.clone());

        ref_states.push(parent_state);
        ref_states.push(child_state);
      }

      for (edge_key, edge_data) in removed_edges {
        self.edges.insert(edge_key, edge_data);
      }

      seq_info.profile = combine_messages(
        &seq_info.seq.composition,
        &msgs_to_combine,
        &variable_pos,
        &ref_states,
        &self.alphabet,
        None,
      )?;

      self.nodes.insert(node.key, seq_info);
    }

    // precalculate messages to children that summarize info from their siblings and the parent
    for child_edge_key in &node.child_edge_keys {
      let child_edge_data = self.edges.remove(child_edge_key).unwrap();
      let seq_info = &self.nodes[&node.key];
      let mut seq_dis = MarginalSparseSeqDistribution {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        fixed: btreemap! {},
        fixed_counts: seq_info.seq.composition.clone(),
        log_lh: seq_info.profile.log_lh - child_edge_data.msg_from_child.log_lh,
      };

      let child_dis = child_edge_data.msg_from_child.clone();
      let mut parent_states: BTreeMap<usize, crate::representation::seq_char::AsciiChar> = btreemap! {};
      let mut child_states: BTreeMap<usize, crate::representation::seq_char::AsciiChar> = btreemap! {};
      for (pos, p) in &seq_info.profile.variable {
        child_states.insert(*pos, p.state);
        parent_states.insert(*pos, p.state);
      }
      for (pos, p) in &child_dis.variable {
        child_states.insert(*pos, p.state);
        parent_states.entry(*pos).or_insert(p.state);
      }
      for sub in &child_edge_data.subs {
        child_states.insert(sub.pos(), sub.qry());
        parent_states.insert(sub.pos(), sub.reff());
      }

      let mut delta_ll = 0.0;
      for (pos, pstate) in parent_states {
        let divisor = if let Some(dis) = child_dis.variable.get(&pos) {
          &dis.dis
        } else {
          &child_dis.fixed[&child_states[&pos]]
        };
        let numerator = if let Some(dis) = seq_info.profile.variable.get(&pos) {
          &dis.dis
        } else {
          &seq_info.profile.fixed[&pstate]
        };
        let dis = numerator / divisor;
        let norm = dis.sum();
        delta_ll += norm.ln();
        seq_dis.variable.insert(
          pos,
          VarPos {
            dis: dis / norm,
            state: pstate,
          },
        );
        seq_dis.fixed_counts.adjust_count(pstate, -1);
      }
      for (s, p) in &seq_info.profile.fixed {
        let dis = p / &child_dis.fixed[s];
        let norm = dis.sum();
        delta_ll += norm.ln() * (seq_dis.fixed_counts.get(*s).unwrap() as f64);
        seq_dis.fixed.insert(*s, dis / norm);
      }
      seq_dis.log_lh += delta_ll;
      let mut updated_edge_data = child_edge_data;
      updated_edge_data.msg_to_child = seq_dis;
      self.edges.insert(*child_edge_key, updated_edge_data);
    }
    Ok(())
  }

  fn extract_ancestral_sequence(&mut self, node_key: GraphNodeKey) -> Seq {
    if let Some(node_data) = self.nodes.get(&node_key) {
      if !node_data.seq.sequence.is_empty() {
        node_data.seq.sequence.clone()
      } else {
        // Return empty seq to be filled later by the reconstruction algorithm
        crate::seq![]
      }
    } else {
      crate::seq![]
    }
  }

  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
    include_leaves: bool,
  ) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let mut node_data = self.nodes.remove(&node.key)?;

    let mut seq = if node.is_root {
      node_data.seq.sequence.clone()
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      let parent_data = self.nodes.get(parent_key)?;
      let edge_data = self.edges.get(edge_key)?;

      let mut seq = parent_data.seq.sequence.clone();

      // Implant mutations
      for m in &edge_data.subs {
        seq[m.pos()] = m.qry();
      }

      // Implant indels
      for indel in &edge_data.indels {
        if indel.deletion {
          seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
        } else {
          seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        }
      }

      seq
    };

    // At the node itself, mask whatever is unknown in the node.
    let alphabet = &self.alphabet;
    for r in &node_data.seq.unknown {
      let ambig_char = alphabet.unknown();
      seq[r.0..r.1].fill(ambig_char);
    }

    // change variable sites to their most likely state
    for (pos, states) in &node_data.profile.variable {
      seq[*pos] = alphabet.char(states.dis.argmax().unwrap());
    }

    node_data.seq.sequence = seq.clone();
    self.nodes.insert(node.key, node_data);

    Some(seq)
  }

  fn get_sequence_length(&self) -> Option<usize> {
    Some(self.length)
  }
}

impl PartitionWithGtrInference for PartitionMarginalSparse {
  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn get_seq_composition(&self, node_key: GraphNodeKey) -> &Composition {
    &self.nodes[&node_key].seq.composition
  }

  fn get_edge_substitutions(&self, edge_key: GraphEdgeKey, _graph: &GraphAncestral) -> Vec<Sub> {
    self.edges[&edge_key].subs.clone()
  }
}

const EPS: f64 = 1e-4;

fn combine_messages(
  composition: &Composition,
  messages: &[MarginalSparseSeqDistribution],
  variable_pos: &BTreeMap<usize, crate::representation::seq_char::AsciiChar>,
  reference_states: &[BTreeMap<usize, crate::representation::seq_char::AsciiChar>],
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<MarginalSparseSeqDistribution, Report> {
  let mut seq_dis = MarginalSparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: btreemap! {},
    fixed: btreemap! {},
    fixed_counts: composition.clone(),
    log_lh: messages.iter().map(|m| m.log_lh).sum(),
  };

  let mut fixed_counts = composition
    .counts()
    .iter()
    .map(|(k, v)| (*k, *v as f64))
    .collect::<BTreeMap<_, _>>();

  for (&pos, &state) in variable_pos {
    let mut all_states_equal = true;
    let mut vec = if let Some(gtr_weight) = gtr_weight {
      gtr_weight.clone()
    } else {
      Array1::from_elem(alphabet.n_canonical(), 1.0)
    };

    for (msg, states) in zip(messages, reference_states) {
      if let Some(var) = msg.variable.get(&pos) {
        vec *= &var.dis;
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = states.get(&pos) {
        if alphabet.is_canonical(*ref_state) {
          vec *= &msg.fixed[ref_state];
        }
        if ref_state != &state {
          all_states_equal = false;
        }
      } else {
        vec *= &msg.fixed[&state];
      }
    }

    let vec_norm = vec.sum();
    seq_dis.log_lh += vec_norm.ln();
    if let Some(count) = fixed_counts.get_mut(&state) {
      *count -= 1.0;
    }

    if (*vec.max()? < (1.0 - EPS) * vec_norm) || !all_states_equal {
      seq_dis.fixed_counts.adjust_count(state, -1);
      let dis = vec / vec_norm;
      seq_dis.variable.insert(pos, VarPos { dis, state });
    }
  }

  for state in alphabet.canonical() {
    let mut vec = if let Some(gtr_weight) = gtr_weight {
      gtr_weight.clone()
    } else {
      Array1::from_elem(alphabet.n_canonical(), 1.0)
    };

    for msg in messages {
      vec *= &msg.fixed[&state];
    }
    let vec_norm = vec.sum();

    seq_dis.log_lh += fixed_counts[&state] * vec_norm.ln();
    seq_dis.fixed.insert(state, vec / vec_norm);
  }
  Ok(seq_dis)
}

fn propagate_raw(
  exp_qt: &Array2<f64>,
  seq_dis: &MarginalSparseSeqDistribution,
  transmission: Option<&Vec<(usize, usize)>>,
) -> MarginalSparseSeqDistribution {
  let mut message = MarginalSparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: btreemap! {},
    fixed: btreemap! {},
    fixed_counts: seq_dis.fixed_counts.clone(),
    log_lh: seq_dis.log_lh,
  };
  for (pos, state) in &seq_dis.variable {
    if let Some(transmission) = &transmission {
      if !range_contains(transmission, *pos) {
        continue;
      }
    }

    let dis = exp_qt.dot(&state.dis);
    let child_state = state.state;
    message.variable.insert(
      *pos,
      VarPos {
        dis,
        state: child_state,
      },
    );
  }

  for (&s, p) in &seq_dis.fixed {
    message.fixed.insert(s, exp_qt.dot(p));
  }

  message
}
