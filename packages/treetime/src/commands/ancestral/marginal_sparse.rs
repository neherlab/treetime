use crate::alphabet::alphabet::Alphabet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::representation::graph_sparse::{SparseGraph, SparseNode, SparseSeqDis, VarPos};
use crate::representation::partitions_likelihood::PartitionLikelihood;
use crate::seq::composition;
use crate::utils::container::get_exactly_one_mut;
use crate::utils::interval::range::range_contains;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use log::debug;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;
use std::iter::zip;

const EPS: f64 = 1e-3;

fn ingroup_profiles_sparse(graph: &SparseGraph, partitions: &[PartitionLikelihood]) {
  graph.par_iter_breadth_first_backward(|mut node| {
    for (si, seq_info) in node.payload.sparse_partitions.iter_mut().enumerate() {
      let PartitionLikelihood { gtr, alphabet, length } = &partitions[si];
      let msg_to_parent = if node.is_leaf {
        // this is mostly a copy (or ref here) of the fitch state.
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
            (
              *pos,
              VarPos {
                dis: p.dis.clone(),
                state: p.state.unwrap(),
              },
            )
          })
          .collect();

        SparseSeqDis {
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
        let mut child_states: Vec<BTreeMap<usize, char>> = vec![];
        let mut child_messages: Vec<SparseSeqDis> = vec![];
        for (ci, (c, edge)) in node.children.iter().enumerate() {
          // go over all mutations and get reference and child state
          child_states.push(btreemap! {});
          for m in &edge.read_arc().sparse_partitions[si].subs {
            variable_pos.insert(m.pos, m.reff); // this might be set multiple times, but the reference state should always be the same
            child_states[ci].insert(m.pos, m.qry);
          }
          // go over child variable position and get reference state
          for (pos, p) in &edge.read_arc().sparse_partitions[si].msg_from_child.variable {
            variable_pos.entry(*pos).or_insert(p.state);
          }
          // FIXME: avoid cloning. could move this loop over child_edges into combine_messages
          child_messages.push(edge.read_arc().sparse_partitions[si].msg_from_child.clone());
        }

        // now that all variable positions are determined, check whether they are characters in each child
        for (ci, (c, edge)) in node.children.iter().enumerate() {
          let states = &mut child_states[ci];
          for (pos, state) in variable_pos.iter() {
            // test whether pos in states, otherwise check whether it is in non-char
            if !states.contains_key(pos) {
              if range_contains(&c.read_arc().sparse_partitions[si].seq.non_char, *pos) {
                dbg!("non-char");
                if range_contains(&c.read_arc().sparse_partitions[si].seq.gaps, *pos) {
                  states.insert(*pos, alphabet.gap());
                } else {
                  states.insert(*pos, alphabet.unknown());
                }
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
          alphabet,
          node.is_root.then_some(&gtr.pi),
        )
        .unwrap()
      };

      if node.is_root {
        // data from children * gtr.pi as calculated above is the root profile.
        seq_info.profile = msg_to_parent;
      } else {
        // what was calculated above is what is sent to the parent. we also calculate the propagated message to the parent (we need it in the forward pass).
        let edge_to_parent =
          get_exactly_one_mut(&mut node.parent_edges).expect("Only nodes with exactly one parent are supported"); // HACK
        let branch_length = edge_to_parent.weight().unwrap_or(0.0);
        let branch_length = fix_branch_length(*length, branch_length);
        let edge_data = &mut edge_to_parent.sparse_partitions[si];
        edge_data.msg_from_child = propagate_raw(
          &gtr.expQt(branch_length).t().to_owned(),
          &msg_to_parent,
          &edge_data.transmission,
        );
        edge_data.msg_to_parent = msg_to_parent;
      }
    }

    GraphTraversalContinuation::Continue
  });
}

fn propagate_raw(
  expQt: &Array2<f64>,
  seq_dis: &SparseSeqDis,
  transmission: &Option<Vec<(usize, usize)>>,
) -> SparseSeqDis {
  let mut message = SparseSeqDis {
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

    let dis = expQt.dot(&state.dis);
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
    message.fixed.insert(s, expQt.dot(p));
  }

  message
}

fn combine_messages(
  composition: &composition::Composition,
  messages: &[SparseSeqDis],
  variable_pos: &BTreeMap<usize, char>,
  reference_states: &[BTreeMap<usize, char>],
  alphabet: &Alphabet,
  gtr_weight: Option<&Array1<f64>>,
) -> Result<SparseSeqDis, Report> {
  let mut seq_dis = SparseSeqDis {
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
    .into_iter()
    .map(|(k, v)| (*k, *v as f64))
    .collect::<BTreeMap<_, _>>();
  // go over all putatively variable positions
  for (&pos, &state) in variable_pos {
    // collect the profiles of children to multiply

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
      } else if let Some(ref_state) = states.get(&pos) {
        // position fixed in origin of the message, but different from focal node
        // If the state of the message is a non-character, skip multiplication
        if alphabet.is_canonical(*ref_state) {
          vec *= &msg.fixed[ref_state];
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
    if *vec.max()? < (1.0 - EPS) * vec_norm {
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

    let fixed_count = seq_dis
      .fixed_counts
      .get(state)
      .ok_or_else(|| make_internal_report!("Unable to find character count for {state}"))?;

    seq_dis.log_lh += fixed_counts[&state] * vec_norm.ln();
    seq_dis.fixed.insert(state, vec / vec_norm);
  }
  Ok(seq_dis)
}

fn outgroup_profiles_sparse(graph: &SparseGraph, partitions: &[PartitionLikelihood]) {
  graph.par_iter_breadth_first_forward(|mut node| {
    for (si, seq_info) in node.payload.sparse_partitions.iter_mut().enumerate() {
      let PartitionLikelihood { gtr, alphabet, length } = &partitions[si];
      if !node.is_root {
        // the root has no input from parents, profile is already calculated
        let mut variable_pos = btreemap! {};
        let mut ref_states: Vec<BTreeMap<usize, char>> = vec![];
        let mut msgs_to_combine: Vec<SparseSeqDis> = vec![];
        for (pi, (p, edge)) in node.parents.iter().enumerate() {
          // go over all mutations and get reference state
          let mut parent_state: BTreeMap<usize, char> = btreemap! {};
          let mut child_state: BTreeMap<usize, char> = btreemap! {};
          // go over parent variable position and get reference state
          for (pos, p) in &edge.read_arc().sparse_partitions[si].msg_to_child.variable {
            if !range_contains(&seq_info.seq.gaps, *pos) {
              // no need to track this position if the position is non-char
              variable_pos.entry(*pos).or_insert(p.state);
              parent_state.insert(*pos, p.state);
            }
          }
          // record all states that involve a mutation
          for m in &edge.read_arc().sparse_partitions[si].subs {
            variable_pos.insert(m.pos, m.qry);
            parent_state.entry(m.pos).or_insert(m.reff);
            child_state.entry(m.pos).or_insert(m.qry);
          }
          // go over variable position in children (info pushed to parent) and get reference state
          for (pos, p) in &edge.read_arc().sparse_partitions[si].msg_to_parent.variable {
            variable_pos.entry(*pos).or_insert(p.state);
            child_state.entry(*pos).or_insert(p.state);
          }

          let branch_length = edge.read_arc().branch_length.unwrap_or(0.0);
          let branch_length = fix_branch_length(*length, branch_length);

          msgs_to_combine.push(propagate_raw(
            &gtr.expQt(branch_length),
            &edge.read_arc().sparse_partitions[si].msg_to_child,
            &edge.read_arc().sparse_partitions[si].transmission,
          ));
          // add combined message from children (which is sent to the parent).
          msgs_to_combine.push(edge.read_arc().sparse_partitions[si].msg_to_parent.clone());
          // NOTE: this empty parent_state is necessary since msgs and states are iterated over in a zip
          ref_states.push(parent_state);
          ref_states.push(child_state);
        }

        seq_info.profile = combine_messages(
          &seq_info.seq.composition,
          &msgs_to_combine,
          &variable_pos,
          &ref_states,
          alphabet,
          None,
        )
        .unwrap();
      }

      // precalculate messages to children that summarize info from their siblings and the parent
      for child_edge in &mut node.child_edges {
        let mut seq_dis = SparseSeqDis {
          variable: btreemap! {},
          variable_indel: btreemap! {},
          fixed: btreemap! {},
          fixed_counts: seq_info.seq.composition.clone(),
          log_lh: seq_info.profile.log_lh - child_edge.sparse_partitions[si].msg_from_child.log_lh,
        };

        let child_dis = &child_edge.sparse_partitions[si].msg_from_child;
        let mut parent_states: BTreeMap<usize, char> = btreemap! {};
        let mut child_states: BTreeMap<usize, char> = btreemap! {};
        for (pos, p) in &seq_info.profile.variable {
          child_states.insert(*pos, p.state);
          parent_states.insert(*pos, p.state);
        }
        for (pos, p) in &child_dis.variable {
          child_states.insert(*pos, p.state);
          parent_states.entry(*pos).or_insert(p.state);
        }
        for sub in &child_edge.sparse_partitions[si].subs {
          child_states.insert(sub.pos, sub.qry);
          parent_states.insert(sub.pos, sub.reff);
        }

        let mut delta_ll = 0.0;
        // subtract the message from the current child from the profile
        for (pos, pstate) in parent_states {
          let divisor = if let Some(dis) = child_dis.variable.get(&pos) {
            &dis.dis
          } else {
            child_dis.fixed.get(child_states.get(&pos).unwrap()).unwrap()
          };
          let numerator = if let Some(dis) = seq_info.profile.variable.get(&pos) {
            &dis.dis
          } else {
            seq_info.profile.fixed.get(&pstate).unwrap()
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
          let dis = &*p / &child_dis.fixed[&s];
          let norm = dis.sum();
          delta_ll += norm.ln() * (seq_dis.fixed_counts.get(*s).unwrap() as f64);
          seq_dis.fixed.insert(*s, dis / norm);
        }
        seq_dis.log_lh += delta_ll;
        child_edge.sparse_partitions[si].msg_to_child = seq_dis;
      }
    }
    GraphTraversalContinuation::Continue
  });
}

pub fn run_marginal_sparse(graph: &SparseGraph, partitions: &[PartitionLikelihood]) -> Result<f64, Report> {
  ingroup_profiles_sparse(graph, partitions);
  let log_lh = graph
    .get_exactly_one_root()?
    .read_arc()
    .payload()
    .read_arc()
    .sparse_partitions
    .iter()
    .map(|p| p.profile.log_lh)
    .sum();
  debug!("Log likelihood: {}", log_lh);
  outgroup_profiles_sparse(graph, partitions);
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal_sparse(
  graph: &SparseGraph,
  include_leaves: bool,
  partitions: &[PartitionLikelihood],
  mut visitor: impl FnMut(&SparseNode, Vec<char>),
) -> Result<(), Report> {
  let n_partitions = partitions.len();

  graph.iter_depth_first_preorder_forward(|mut node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq = (0..n_partitions)
      .flat_map(|si| {
        let PartitionLikelihood { alphabet, .. } = &partitions[si];
        let node_seq = &node.payload.sparse_partitions[si].seq;

        let mut seq = if node.is_root {
          node_seq.sequence.clone()
        } else {
          let (parent, edge) = node.get_exactly_one_parent().unwrap();
          let parent = &parent.read_arc().sparse_partitions[si];
          let edge = &edge.read_arc().sparse_partitions[si];

          let mut seq = parent.seq.sequence.clone();

          // Implant mutations
          for m in &edge.subs {
            seq[m.pos] = m.qry;
          }

          // Implant most likely state of variable sites
          for (&pos, vec) in &node.payload.sparse_partitions[si].seq.fitch.variable {
            seq[pos] = alphabet.char(vec.dis.argmax().unwrap());
          }

          // Implant indels
          for indel in &edge.indels {
            if indel.deletion {
              seq[indel.range.0..indel.range.1].fill(alphabet.gap());
            } else {
              seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
            }
          }

          seq
        };

        // At the node itself, mask whatever is unknown in the node.
        for r in &node_seq.unknown {
          let ambig_char = partitions[si].alphabet.unknown();
          seq[r.0..r.1].fill(ambig_char);
        }

        for (pos, p) in &node_seq.fitch.variable {
          seq[*pos] = alphabet.get_code(&p.dis);
        }

        node.payload.sparse_partitions[si].seq.sequence = seq.clone();

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
  use crate::commands::ancestral::fitch::compress_sequences;
  use crate::graph::node::GraphNodeKey;
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use crate::gtr::gtr::{GTRParams, GTR};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partitions_parsimony::PartitionParsimonyWithAln;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use lazy_static::lazy_static;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_ancestral_reconstruction_marginal_sparse() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
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
    "#},
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      ACAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >CD
      CCGGCCATGTATTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let gtr = jc69(JC69Params::default())?;
    let partitions = partitions
      .into_iter()
      .map(|p| PartitionLikelihood::from_parsimony(gtr.clone(), p))
      .collect_vec();
    let log_lh = run_marginal_sparse(&graph, &partitions)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal_sparse(&graph, false, &partitions, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq));
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    // test overall likelihood
    pretty_assert_ulps_eq!(-55.55428499726621, log_lh, epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_root_state() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
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
    "#},
      &NUC_ALPHABET,
    )?;
    let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let pi = array![0.2, 0.3, 0.15, 0.35];

    // symmetric rate matrix
    let mu = 1.0;

    let gtr = GTR::new(GTRParams {
      alphabet: Alphabet::default(),
      W: None,
      pi: pi.clone(),
      mu,
    })
    .unwrap();

    let partitions = partitions
      .into_iter()
      .map(|p| PartitionLikelihood::from_parsimony(gtr.clone(), p))
      .collect_vec();
    let log_lh = run_marginal_sparse(&graph, &partitions)?;

    pretty_assert_ulps_eq!(-56.946298878390444, log_lh, epsilon = 1e-6);

    // test variable position distribution at the root
    let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
    let root = &graph
      .get_exactly_one_root()
      .unwrap()
      .read_arc()
      .payload()
      .read_arc()
      .sparse_partitions[0];
    let pos: usize = 0;
    pretty_assert_ulps_eq!(&root.profile.variable[&pos].dis, &pos_zero_root, epsilon = 1e-6);

    // test variable position distribution at internal node
    let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
    let node_ab = &graph
      .get_node(GraphNodeKey(1))
      .unwrap()
      .read_arc()
      .payload()
      .read_arc()
      .sparse_partitions[0];
    let pos: usize = 0;
    pretty_assert_ulps_eq!(&node_ab.profile.variable[&pos].dis, &pos_zero_ab, epsilon = 1e-6);

    // test variable position distribution at internal node
    let dis_ab = array![
      0.0013914677323952813,
      0.002087201598592933,
      0.042827146239885545,
      0.9536941844291262
    ];
    let pos: usize = 3;
    pretty_assert_ulps_eq!(&node_ab.profile.variable[&pos].dis, &dis_ab, epsilon = 1e-6);

    Ok(())
  }
}
