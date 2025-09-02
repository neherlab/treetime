use crate::alphabet::alphabet::Alphabet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::make_internal_error;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::graph_sparse::{MarginalSparseSeqDistribution, VarPos};
use crate::representation::log_lh::graph_log_lh;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use crate::seq::composition;
use crate::utils::container::get_exactly_one;
use crate::utils::interval::range::range_contains;
use eyre::Report;
use log::debug;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use ndarray_stats::QuantileExt;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::iter::zip;
use std::sync::Arc;

const EPS: f64 = 1e-4;

fn marginal_sparse_backward(graph: &GraphAncestral, partitions: &[Arc<RwLock<PartitionMarginalSparse>>]) {
  graph.par_iter_breadth_first_backward(|node| {
    run_marginal_sparse_backward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
}

/// Backward pass calculates ingroup profiles
fn run_marginal_sparse_backward(
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  for partition in partitions {
    let mut partition = partition.write_arc();
    let length = partition.length;

    let mut seq_info = partition.nodes.remove(&node.key).unwrap();
    let msg_to_parent = if node.is_leaf {
      // this is mostly a copy (or ref here) of the fitch state.
      let alphabet = &partition.alphabet;
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
      // to do so, we need to loop over incoming edges, collect variable positions and the child states at them
      let mut variable_pos = btreemap! {};
      let mut child_states = vec![];
      let mut child_messages: Vec<MarginalSparseSeqDistribution> = vec![];
      for (ci, (_child_key, edge_key)) in node.child_keys.iter().enumerate() {
        // go over all mutations and get reference and child state
        child_states.push(btreemap! {});
        let edge_data = &partition.edges[edge_key];
        for m in &edge_data.subs {
          variable_pos.insert(m.pos(), m.reff()); // this might be set multiple times, but the reference state should always be the same
          child_states[ci].insert(m.pos(), m.qry());
        }
        // go over child variable position and get reference state
        for (pos, p) in &edge_data.msg_from_child.variable {
          variable_pos.entry(*pos).or_insert(p.state);
        }
        // FIXME: avoid cloning. could move this loop over child_edges into combine_messages
        child_messages.push(edge_data.msg_from_child.clone());
      }

      // now that all variable positions are determined, check whether they are characters in each child
      let alphabet = &partition.alphabet;
      for (ci, (child_key, _)) in node.child_keys.iter().enumerate() {
        let states = &mut child_states[ci];
        let child_data = &partition.nodes[child_key];
        for pos in variable_pos.keys() {
          // test whether pos in states, otherwise check whether it is in non-char
          if !states.contains_key(pos) && range_contains(&child_data.seq.non_char, *pos) {
            if range_contains(&child_data.seq.gaps, *pos) {
              states.insert(*pos, alphabet.gap());
            } else {
              states.insert(*pos, alphabet.unknown());
            }
          }
        }
      }

      // messages are combined, if this is the root node the gtr.pi is used to multiply the data from children
      combine_messages(
        &seq_info.seq.composition,
        &child_messages,
        &variable_pos,
        &child_states,
        &partition.alphabet,
        node.is_root.then_some(&partition.gtr.pi),
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
      let mut edge_data = partition.edges.remove(edge_key).unwrap();
      edge_data.msg_from_child = propagate_raw(
        &partition.gtr.expQt(branch_length).t().to_owned(),
        &msg_to_parent,
        edge_data.transmission.as_ref(),
      );
      edge_data.msg_to_parent = msg_to_parent;
      partition.edges.insert(*edge_key, edge_data);
    }
    partition.nodes.insert(node.key, seq_info);
  }
  Ok(())
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
        continue; // transmission field is not currently used
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

fn combine_messages(
  composition: &composition::Composition,
  messages: &[MarginalSparseSeqDistribution],
  variable_pos: &BTreeMap<usize, AsciiChar>,
  reference_states: &[BTreeMap<usize, AsciiChar>],
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<MarginalSparseSeqDistribution, Report> {
  let mut seq_dis = MarginalSparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: btreemap! {},
    fixed: btreemap! {},
    fixed_counts: composition.clone(), //this comes from fitch without any variable positions
    log_lh: messages.iter().map(|m| m.log_lh).sum(),
  };

  // create a copy of the composition to keep track of the fixed states. this needs to be separate from the fixed_counts field in seq_dis
  // pulling out the counts as BTmap with f64 is sufficient
  // copy composition.counts and cast the values of f64
  let mut fixed_counts = composition
    .counts()
    .iter()
    .map(|(k, v)| (*k, *v as f64))
    .collect::<BTreeMap<_, _>>();
  // go over all putatively variable positions
  for (&pos, &state) in variable_pos {
    // collect the profiles of children to multiply
    let mut all_states_equal = true;
    // start with a vector given by the gtr_weight if it is provided otherwise use a vector of ones
    let mut vec = if let Some(gtr_weight) = gtr_weight {
      gtr_weight.clone()
    } else {
      Array1::from_elem(alphabet.n_canonical(), 1.0)
    };

    // element wise multiplication of the profiles of the messages with the vec
    // zip messages and the dict of corresponding states to iterate over them simultaneously
    for (msg, states) in zip(messages, reference_states) {
      //FIXME: handle transmission and non-char
      if let Some(var) = msg.variable.get(&pos) {
        // position variable in child
        vec *= &var.dis;
        if var.state != state {
          all_states_equal = false;
        }
      } else if let Some(ref_state) = states.get(&pos) {
        // position fixed in origin of the message, but different from focal node
        // If the state of the message is a non-character, skip multiplication
        if alphabet.is_canonical(*ref_state) {
          vec *= &msg.fixed[ref_state];
        }
        if ref_state != &state {
          all_states_equal = false;
        }
      } else {
        // position fixed in child and same as parent
        vec *= &msg.fixed[&state];
      }
    }

    let vec_norm = vec.sum();
    seq_dis.log_lh += vec_norm.ln();
    if let Some(count) = fixed_counts.get_mut(&state) {
      *count -= 1.0;
    }
    // add position to variable states if the subleading states have a probability exceeding eps
    if (*vec.max()? < (1.0 - EPS) * vec_norm) || !all_states_equal {
      if vec.ndim() > 1 {
        return make_internal_error!("Unexpected dimensionality in probability vector: {}", vec.ndim());
      }

      // FIXME: check what is supposed to happen here when state is not diverse
      seq_dis.fixed_counts.adjust_count(state, -1);

      let dis = vec / vec_norm;
      seq_dis.variable.insert(pos, VarPos { dis, state });
    }
  }

  // collect contribution from the fixed sites
  for state in alphabet.canonical() {
    // indeterminate parts in some messages are not handled correctly here.
    // they should not contribute to the product. This will require some additional
    // handling or could be handled by treating these positions as variable

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

/// Forward pass calculates outgroup profiles
fn marginal_sparse_forward(graph: &GraphAncestral, partitions: &[Arc<RwLock<PartitionMarginalSparse>>]) {
  graph.par_iter_breadth_first_forward(|node| {
    run_marginal_sparse_forward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
}

fn run_marginal_sparse_forward(
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  for partition in partitions {
    let mut partition = partition.write_arc();
    let length = partition.length;

    let mut seq_info = partition.nodes.remove(&node.key).unwrap();
    if !node.is_root {
      // the root has no input from parents, profile is already calculated
      let mut variable_pos = btreemap! {};
      let mut ref_states: Vec<BTreeMap<usize, AsciiChar>> = vec![];
      let mut msgs_to_combine: Vec<MarginalSparseSeqDistribution> = vec![];
      let mut removed_edges = vec![];
      for (_, edge_key) in &node.parent_keys {
        // go over all mutations and get reference state
        let mut parent_state: BTreeMap<usize, AsciiChar> = btreemap! {};
        let mut child_state: BTreeMap<usize, AsciiChar> = btreemap! {};
        let edge_data = partition.edges.remove(edge_key).unwrap();
        removed_edges.push((*edge_key, edge_data.clone()));
        // go over parent variable position and get reference state
        for (pos, p) in &edge_data.msg_to_child.variable {
          if !range_contains(&seq_info.seq.gaps, *pos) {
            // no need to track this position if the position is non-char
            variable_pos.entry(*pos).or_insert(p.state);
            parent_state.insert(*pos, p.state);
          }
        }
        // record all states that involve a mutation
        for m in &edge_data.subs {
          variable_pos.insert(m.pos(), m.qry());
          parent_state.entry(m.pos()).or_insert_with(|| m.reff());
          child_state.entry(m.pos()).or_insert_with(|| m.qry());
        }
        // go over variable position in children (info pushed to parent) and get reference state
        for (pos, p) in &edge_data.msg_to_parent.variable {
          variable_pos.entry(*pos).or_insert(p.state);
          child_state.entry(*pos).or_insert(p.state);
        }

        let branch_length = node.parents[0].1.read_arc().branch_length.unwrap_or(0.0);
        let branch_length = fix_branch_length(length, branch_length);

        msgs_to_combine.push(propagate_raw(
          &partition.gtr.expQt(branch_length),
          &edge_data.msg_to_child,
          edge_data.transmission.as_ref(),
        ));
        // add combined message from children (which is sent to the parent).
        msgs_to_combine.push(edge_data.msg_to_parent.clone());

        // NOTE: this empty parent_state is necessary since msgs and states are iterated over in a zip
        ref_states.push(parent_state);
        ref_states.push(child_state);
      }

      // Put the edges back
      for (edge_key, edge_data) in removed_edges {
        partition.edges.insert(edge_key, edge_data);
      }

      seq_info.profile = combine_messages(
        &seq_info.seq.composition,
        &msgs_to_combine,
        &variable_pos,
        &ref_states,
        &partition.alphabet,
        None,
      )?;
    }

    // precalculate messages to children that summarize info from their siblings and the parent
    for child_edge_key in &node.child_edge_keys {
      let child_edge_data = partition.edges.remove(child_edge_key).unwrap();
      let mut seq_dis = MarginalSparseSeqDistribution {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        fixed: btreemap! {},
        fixed_counts: seq_info.seq.composition.clone(),
        log_lh: seq_info.profile.log_lh - child_edge_data.msg_from_child.log_lh,
      };

      let child_dis = child_edge_data.msg_from_child.clone();
      let mut parent_states: BTreeMap<usize, AsciiChar> = btreemap! {};
      let mut child_states: BTreeMap<usize, AsciiChar> = btreemap! {};
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
      // subtract the message from the current child from the profile
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
      partition.edges.insert(*child_edge_key, updated_edge_data);
    }
    partition.nodes.insert(node.key, seq_info);
  }
  Ok(())
}

pub fn run_marginal_sparse(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> Result<f64, Report> {
  marginal_sparse_backward(graph, partitions);
  let log_lh = graph_log_lh(graph, partitions)?;
  debug!("Log likelihood: {log_lh}");
  marginal_sparse_forward(graph, partitions);
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal_sparse(
  graph: &GraphAncestral,
  include_leaves: bool,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  mut visitor: impl FnMut(&NodeAncestral, Seq),
) -> Result<(), Report> {
  graph.iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq: Seq = partitions
      .iter()
      .flat_map(|partition| {
        let mut partition = partition.write_arc();
        let mut node_data = partition.nodes.remove(&node.key).unwrap();

        let mut seq = if node.is_root {
          node_data.seq.sequence.clone()
        } else {
          let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).unwrap();
          let parent_data = &partition.nodes[parent_key];
          let edge_data = &partition.edges[edge_key];

          let mut seq = parent_data.seq.sequence.clone();

          // Implant mutations
          for m in &edge_data.subs {
            seq[m.pos()] = m.qry();
          }

          // Implant indels
          for indel in &edge_data.indels {
            if indel.deletion {
              seq[indel.range.0..indel.range.1].fill(partition.alphabet.gap());
            } else {
              seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
            }
          }

          seq
        };

        // At the node itself, mask whatever is unknown in the node.
        let alphabet = &partition.alphabet;
        for r in &node_data.seq.unknown {
          let ambig_char = alphabet.unknown();
          seq[r.0..r.1].fill(ambig_char);
        }

        // change variable sites to their most likely state
        for (pos, states) in &node_data.profile.variable {
          seq[*pos] = alphabet.char(states.dis.argmax().unwrap());
        }

        node_data.seq.sequence = seq.clone();
        partition.nodes.insert(node.key, node_data);

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
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{JsonPretty, json_write_str};
  use crate::io::nwk::nwk_read_str;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::graph_ancestral::GraphAncestral;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_ancestral_reconstruction_marginal_sparse() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTG--
      >AB
      ACATCGCTGTA-TG--
      >CD
      TCGGCGGTGTATTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_marginal_sparse = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_marginal_sparse, &aln)?;

    let log_lh = run_marginal_sparse(&graph, &partitions_marginal_sparse)?;

    // generate ancestral reconstruction and test against expectation
    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal_sparse(&graph, false, &partitions_marginal_sparse, |node, seq| {
      actual.insert(node.name.clone(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    // test overall likelihood
    pretty_assert_ulps_eq!(-55.55428499726621, log_lh, epsilon = 1e-6);
    Ok(())
  }

  //   #[test]
  //   fn test_root_state() -> Result<(), Report> {
  //     rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //     let aln = read_many_fasta_str(
  //       indoc! {r#"
  //       >A
  //       ACATCGCCNNA--GAC
  //       >B
  //       GCATCCCTGTA-NG--
  //       >C
  //       CCGGCGATGTRTTG--
  //       >D
  //       TCGGCCGTGTRTTG--
  //     "#},
  //       &NUC_ALPHABET,
  //     )?;
  //     let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  //
  //     let alphabet = Alphabet::default();
  //
  //     let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
  //     let partitions = compress_sequences(&graph, partitions)?;
  //
  //     // use non-trivial GTR with non-uniform stationary distribution
  //     let mu = 1.0;
  //     let pi = array![0.2, 0.3, 0.15, 0.35];
  //     let gtr = GTR::new(GTRParams {
  //       alphabet: Alphabet::default(),
  //       W: None,
  //       pi,
  //       mu,
  //     })?;
  //
  //     let partitions = partitions
  //       .into_iter()
  //       .map(|p| PartitionLikelihood::from_parsimony(gtr.clone(), p))
  //       .collect_vec();
  //
  //     let log_lh = run_marginal_sparse(&graph, &partitions)?;
  //     // note that this LH is slightly different from dense or python treetime due to
  //     // different handling of ambiguous characters (value from test_scripts/ancestral_sparse.py)
  //     pretty_assert_ulps_eq!(-56.946298878390444, log_lh, epsilon = 1e-6);
  //
  //     // test variable position distribution at the root (pos 0)
  //     let root = &graph
  //       .get_exactly_one_root()?
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .sparse_partitions[0];
  //     let pos: usize = 0;
  //     let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
  //     pretty_assert_ulps_eq!(&root.profile.variable[&pos].dis, &pos_zero_root, epsilon = 1e-6);
  //
  //     // pull out internal node AB for further testing
  //     let node_ab = &graph
  //       .get_node(GraphNodeKey(1))
  //       .unwrap()
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .sparse_partitions[0];
  //
  //     // test variable position distribution at internal node (pos 0)
  //     let pos: usize = 0;
  //     let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
  //     pretty_assert_ulps_eq!(&node_ab.profile.variable[&pos].dis, &pos_zero_ab, epsilon = 1e-6);
  //
  //     // test variable position distribution at internal node (pos 3)
  //     let dis_ab = array![
  //       0.0013914677323952813,
  //       0.002087201598592933,
  //       0.042827146239885545,
  //       0.9536941844291262
  //     ];
  //     let pos: usize = 3;
  //     pretty_assert_ulps_eq!(&node_ab.profile.variable[&pos].dis, &dis_ab, epsilon = 1e-6);
  //
  //     Ok(())
  //   }
}
