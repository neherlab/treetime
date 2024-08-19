use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::graph::node::Named;
use crate::gtr::gtr::GTR;
use crate::port::constants::GAP_CHAR;
use crate::port::seq_partitions::SeqPartition;
use crate::port::seq_sparse::{SparseGraph, SparseNode, SparseSeqDis, SparseSeqEdge, SparseSeqNode, VarPos};
use crate::seq::range::range_contains;
use crate::utils::ndarray::{product_axis, stack_owned};
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use ndarray::{array, Array1, Array2, Axis};
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;

const EPS: f64 = 1e-6;

fn ingroup_profiles_sparse(graph: &SparseGraph) {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;

  graph.par_iter_breadth_first_backward(|mut node| {
    for (si, seq_info) in node.payload.sparse_partitions.iter_mut().enumerate() {
      let &SeqPartition { gtr, length } = &sparse_partitions[si];

      if node.is_leaf {
        // this is mostly a copy (or ref here) of the fitch state.
        let fixed = gtr
          .alphabet()
          .determined()
          .map(|state| (state, gtr.alphabet().get_profile(state).clone()))
          .collect();

        let variable = seq_info.seq.fitch.variable.clone();

        seq_info.msg_to_parents = SparseSeqDis {
          fixed_counts: seq_info.seq.composition.clone(),
          variable,
          fixed,
          ..SparseSeqDis::default()
        };
      } else {
        // internal nodes
        let child_expQt = node
          .children
          .iter()
          .map(|(c, e)| gtr.expQt(e.read_arc().weight().unwrap_or(0.0)).t().to_owned())
          .collect_vec();

        let child_seqs = node
          .children
          .iter()
          .map(|(c, _)| c.read_arc().sparse_partitions[si].clone()) // FIXME: avoid cloning
          .collect_vec();

        let child_edges = node
          .children
          .iter()
          .map(|(_, e)| e.read_arc().sparse_partitions[si].clone()) // FIXME: avoid cloning
          .collect_vec();

        // get all variable positions, the reference state, and the child states at these positions
        let (variable_pos, child_states) = get_variable_states_children(&child_seqs, &child_edges);

        let mut seq_dis = SparseSeqDis {
          fixed_counts: seq_info.seq.composition.clone(),
          ..SparseSeqDis::default()
        };

        seq_info.msgs_from_children = btreemap! {};
        for (ci, (child, _)) in node.children.iter().enumerate() {
          let name = child.read_arc().name().unwrap().as_ref().to_owned();
          seq_dis.log_lh += child_seqs[ci].msg_to_parents.log_lh;
          let message = propagate(
            &child_expQt[ci],
            &child_seqs[ci].msg_to_parents,
            &variable_pos,
            &child_states[ci],
            &child_seqs[ci].seq.non_char,
            &None,
          );
          seq_info.msgs_from_children.insert(name, message);
        }

        combine_messages(
          &mut seq_dis,
          &seq_info.msgs_from_children.values().cloned().collect_vec(), // FIXME: avoid cloning
          &variable_pos,
          gtr,
          None,
        )
        .unwrap();

        seq_info.msg_to_parents = seq_dis;
      }
    }

    GraphTraversalContinuation::Continue
  });
}

fn propagate(
  expQt: &Array2<f64>,
  seq_dis: &SparseSeqDis,
  variable_pos: &BTreeMap<usize, char>,
  child_states: &BTreeMap<usize, char>,
  non_char: &[(usize, usize)],
  transmission: &Option<Vec<(usize, usize)>>,
) -> SparseSeqDis {
  let mut message = SparseSeqDis {
    fixed_counts: seq_dis.fixed_counts.clone(),
    ..SparseSeqDis::default()
  };
  for (pos, &state) in variable_pos {
    if let Some(transmission) = &transmission {
      if !range_contains(transmission, *pos) {
        continue; // transmission field is not currently used
      }
    }

    if range_contains(non_char, *pos) {
      continue;
    }

    let v = if let Some(variable) = seq_dis.variable.get(pos) {
      variable.dis.clone()
    } else if let Some(&state) = child_states.get(pos) {
      message.fixed_counts.adjust_count(state, -1);
      seq_dis.fixed[&state].clone()
    } else {
      message.fixed_counts.adjust_count(state, -1);
      seq_dis.fixed[&state].clone()
    };

    let dis = expQt.dot(&v);
    let state = Some(state);
    message.variable.insert(*pos, VarPos { dis, state });
  }

  for (&s, p) in &seq_dis.fixed {
    message.fixed.insert(s, expQt.dot(p));
  }

  message
}

fn combine_messages(
  seq_dis: &mut SparseSeqDis,
  messages: &[SparseSeqDis],
  variable_pos: &BTreeMap<usize, char>,
  gtr: &GTR,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<(), Report> {
  // go over all putatively variable positions
  for (&pos, &state) in variable_pos {
    // collect the profiles of children to multiply
    let mut msg_dis: Vec<Array1<f64>> = gtr_weight
      .map(|gtr_weight| vec![gtr_weight.clone()])
      .unwrap_or_default();

    for msg in messages {
      if let Some(var) = msg.variable.get(&pos) {
        msg_dis.push(var.dis.clone());
      }
    }

    // calculate new profile and likelihood contribution
    let vec = if !msg_dis.is_empty() {
      let msg_dis = stack_owned(Axis(0), &msg_dis)?;
      product_axis(&msg_dis, Axis(0))
    } else {
      array![1.0]
    };
    let vec_norm = vec.sum();

    // add position to variable states if the subleading states have a probability exceeding eps
    if *vec.max()? < (1.0 - EPS) * vec_norm {
      if vec.ndim() > 1 {
        return make_internal_error!("Unexpected dimensionality in probability vector: {}", vec.ndim());
      }

      seq_dis.log_lh += vec_norm.ln();
      seq_dis.fixed_counts.adjust_count(state, -1);

      let dis = vec / vec_norm;
      let state = Some(state);
      seq_dis.variable.insert(pos, VarPos { dis, state });
    }
  }

  // collect contribution from the fixed sites
  for state in gtr.alphabet().canonical() {
    // indeterminate parts in some children are not handled correctly here.
    // they should not contribute to the product. This will require some additional
    // handling or could be handled by treating these positions as variable
    let mut msg_dis: Vec<Array1<f64>> = gtr_weight
      .map(|gtr_weight| vec![gtr_weight.clone()])
      .unwrap_or_default();

    for msg in messages {
      msg_dis.push(msg.fixed[&state].clone());
    }

    let msg_dis = stack_owned(Axis(0), &msg_dis)?;
    let vec: Array1<f64> = product_axis(&msg_dis, Axis(0));
    let vec_norm = vec.sum();

    let fixed_count = seq_dis
      .fixed_counts
      .get(state)
      .ok_or_else(|| make_internal_report!("Unable to find character count for {state}"))?;

    seq_dis.log_lh += (fixed_count as f64) * vec_norm.ln();
    seq_dis.fixed.insert(state, vec / vec_norm);
  }

  Ok(())
}

fn get_variable_states_children(
  child_seqs: &[SparseSeqNode],
  child_edges: &[SparseSeqEdge],
) -> (BTreeMap<usize, char>, Vec<BTreeMap<usize, char>>) {
  let mut variable_pos = btreemap! {};
  let mut child_states: Vec<BTreeMap<usize, char>> = vec![btreemap! {}; child_edges.len()];
  for (ci, edge) in child_edges.iter().enumerate() {
    // go over all mutations and get reference state
    for m in &edge.muts {
      variable_pos.insert(m.pos, m.reff); // this might be set multiple times, but the reference state should always be the same
      child_states[ci].insert(m.pos, m.qry);
    }
  }
  for seq_info in child_seqs {
    // go over child variable position and get reference state
    for (pos, p) in &seq_info.msg_to_parents.variable {
      if let Some(p_state) = p.state {
        variable_pos.entry(*pos).or_insert(p_state);
      }
    }
  }
  (variable_pos, child_states)
}

fn outgroup_profiles_sparse(graph: &SparseGraph) {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;

  graph.par_iter_breadth_first_forward(|mut node| {
    if node.is_root {
      return GraphTraversalContinuation::Continue;
    }

    let name = node.payload.name().unwrap().as_ref().to_owned();

    for (si, seq_info) in node.payload.sparse_partitions.iter_mut().enumerate() {
      let &SeqPartition { gtr, length } = &sparse_partitions[si];

      let parent_nodes = node
        .parents
        .iter()
        .map(|(p, _)| p.read_arc().sparse_partitions[si].profile.clone()) // FIXME: avoid cloning
        .collect_vec();

      let parent_edges = node
        .parents
        .iter()
        .map(|(_, e)| e.read_arc().sparse_partitions[si].clone()) // FIXME: avoid cloning
        .collect_vec();

      let (variable_pos, parent_states) =
        get_variable_states_parents(&seq_info.msg_to_parents, &parent_nodes, &parent_edges);

      let mut msgs_from_parents = vec![];
      for (p, e) in &node.parents {
        let pseq_info = &p.read_arc().sparse_partitions[si];
        let msg = propagate(
          &gtr.expQt(e.read_arc().weight().unwrap_or(0.0)),
          &pseq_info.msgs_to_children[&name],
          &variable_pos,
          &parent_states[si],
          &pseq_info.seq.non_char,
          &None,
        );
        msgs_from_parents.push(msg);
      }

      let mut seq_dis = SparseSeqDis {
        fixed_counts: seq_info.seq.composition.clone(),
        ..SparseSeqDis::default()
      };

      combine_messages(&mut seq_dis, &msgs_from_parents, &variable_pos, gtr, None).unwrap();
      seq_info.profile = seq_dis;

      // precalculate messages to children that summarize info from their siblings and the parent
      seq_info.msgs_to_children = btreemap! {};
      for child_name in seq_info.msgs_from_children.keys() {
        let msgs = seq_info
          .msgs_from_children
          .iter()
          .filter(|&(k, _)| k != child_name)
          .map(|(_, m)| m.to_owned()) // FIXME: avoid cloning
          .collect_vec();

        let mut seq_dis = SparseSeqDis {
          fixed_counts: seq_info.msg_to_parents.fixed_counts.clone(),
          ..SparseSeqDis::default()
        };

        combine_messages(&mut seq_dis, &msgs, &variable_pos, gtr, None).unwrap();
        seq_info.msgs_to_children.insert(child_name.clone(), seq_dis);
      }
    }

    GraphTraversalContinuation::Continue
  });
}

fn get_variable_states_parents(
  node_seq_info: &SparseSeqDis,
  parents: &[SparseSeqDis],
  parent_edges: &[SparseSeqEdge],
) -> (BTreeMap<usize, char>, Vec<BTreeMap<usize, char>>) {
  let mut variable_pos: BTreeMap<usize, char> = node_seq_info
    .variable
    .iter()
    .filter_map(|(pos, p)| p.state.map(|p_state| (*pos, p_state)))
    .collect();

  let mut parent_states: Vec<BTreeMap<usize, char>> = vec![btreemap! {}; parents.len()];

  for (pi, pseq) in parents.iter().enumerate() {
    // go over all mutations and get reference state
    for m in &parent_edges[pi].muts {
      variable_pos.insert(m.pos, m.qry); // this might be set multiple times, but the reference state should always be the same
      parent_states[pi].insert(m.pos, m.reff);
    }
    // go over child variable position and get reference state
    for (pos, p) in &pseq.variable {
      if !variable_pos.contains_key(pos) {
        if let Some(p_state) = p.state {
          variable_pos.entry(*pos).or_insert(p_state);
        }
      }
    }
  }

  (variable_pos, parent_states)
}

fn calculate_root_state_sparse(graph: &SparseGraph) -> f64 {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;

  let mut log_lh = 0.0;
  for root_node in graph.get_roots() {
    let root_partitions = &mut root_node.write_arc().payload().write_arc().sparse_partitions;
    for (si, seq_info) in root_partitions.iter_mut().enumerate() {
      let &SeqPartition { gtr, length } = &sparse_partitions[si];

      let mut seq_profile = SparseSeqDis {
        fixed_counts: seq_info.msg_to_parents.fixed_counts.clone(),
        log_lh: seq_info.msg_to_parents.log_lh,
        ..SparseSeqDis::default()
      };

      // multiply the info from the tree with the GTR equilibrium probabilities (variable and fixed)
      for (pos, p) in &seq_info.msg_to_parents.variable {
        let vec = &p.dis * &gtr.pi;
        let vec_norm = vec.sum();
        seq_profile.log_lh += vec_norm.ln();
        let dis = vec / vec_norm;
        let state = p.state;
        seq_profile.variable.insert(*pos, VarPos { dis, state });
      }
      for (state, p) in &seq_info.msg_to_parents.fixed {
        let vec = p * &gtr.pi;
        let vec_norm = vec.sum();
        let count = seq_info.msg_to_parents.fixed_counts.get(*state).unwrap_or_default();
        seq_profile.log_lh += vec_norm.ln() * count as f64;
        seq_profile.fixed.insert(*state, vec / vec_norm);
      }

      log_lh += seq_profile.log_lh;
      seq_info.profile = seq_profile;

      // calculate messages to children
      let variable_pos: BTreeMap<usize, char> = seq_info
        .msg_to_parents
        .variable
        .iter()
        .filter_map(|(pos, p)| p.state.map(|p_state| (*pos, p_state)))
        .collect();

      seq_info.msgs_to_children = btreemap! {};
      for (child, _) in graph.children_of(&*root_node.read_arc()) {
        let name = child
          .read_arc()
          .payload()
          .read_arc()
          .name()
          .unwrap()
          .as_ref()
          .to_owned();

        let msgs = seq_info
          .msgs_from_children
          .iter()
          .filter(|&(k, _)| k != &name)
          .map(|(_, m)| m)
          .cloned() // FIXME: avoid cloning
          .collect_vec();

        let mut seq_dis = SparseSeqDis {
          fixed_counts: seq_info.msg_to_parents.fixed_counts.clone(),
          log_lh: 0.0,
          ..SparseSeqDis::default()
        };

        combine_messages(&mut seq_dis, &msgs, &variable_pos, gtr, None).unwrap();
        seq_info.msgs_to_children.insert(name, seq_dis);
      }
    }
  }
  log_lh
}

pub fn run_marginal_sparse(graph: &SparseGraph) -> Result<f64, Report> {
  ingroup_profiles_sparse(graph);
  let log_lh = calculate_root_state_sparse(graph);
  outgroup_profiles_sparse(graph);
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal_sparse(
  graph: &SparseGraph,
  include_leaves: bool,
  mut visitor: impl FnMut(&SparseNode, Vec<char>),
) -> Result<(), Report> {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;
  let n_partitions = sparse_partitions.len();

  graph.iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq = (0..n_partitions)
      .flat_map(|si| {
        let &SeqPartition { gtr, length } = &sparse_partitions[si];
        let node_seq = &node.payload.sparse_partitions[si].seq;

        let mut seq = if node.is_root {
          node_seq.sequence.clone()
        } else {
          let (parent, edge) = node.get_exactly_one_parent().unwrap();
          let parent = &parent.read_arc().sparse_partitions[si];
          let edge = &edge.read_arc().sparse_partitions[si];

          let mut seq = parent.seq.sequence.clone();

          // Implant mutations
          for m in &edge.muts {
            seq[m.pos] = m.qry;
          }

          // Implant most likely state of variable sites
          for (&pos, vec) in &node.payload.sparse_partitions[si].seq.fitch.variable {
            seq[pos] = gtr.alphabet.char(vec.dis.argmax().unwrap());
          }

          // Implant indels
          for indel in &edge.indels {
            if indel.deletion {
              seq[indel.range.0..indel.range.1].fill(GAP_CHAR);
            } else {
              seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
            }
          }

          seq
        };

        // At the node itself, mask whatever is unknown in the node.
        for r in &node_seq.unknown {
          let ambig_char = sparse_partitions[si].gtr.alphabet.unknown();
          seq[r.0..r.1].fill(ambig_char);
        }

        for (pos, p) in &node_seq.fitch.variable {
          seq[*pos] = sparse_partitions[si].code(&p.dis);
        }

        seq
      })
      .collect();

    visitor(&node.payload, seq);
  });

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::port::fitch::{compress_sequences, PartitionModel};
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_ancestral_reconstruction_marginal_sparse() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let inputs = read_many_fasta_str(indoc! {r#"
      >root
      ACAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#})?;

    let expected = read_many_fasta_str(indoc! {r#"
      >root
      ACAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >CD
      CCGGCCATGTATTG--
    "#})?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let mut graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let gtr = &jc69(JC69Params::default())?;
    let partitions = vec![PartitionModel { gtr, aln: inputs }];
    compress_sequences(&mut graph, &partitions)?;
    run_marginal_sparse(&graph)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal_sparse(&graph, false, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq));
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }
}
