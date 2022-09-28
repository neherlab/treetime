#![allow(non_snake_case)]

use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::anc_args::TreetimeAncestralArgs;
use crate::ancestral::anc_graph::AncestralGraph;
use crate::ancestral::run_ancestral_reconstruction::TreetimeAncestralParams;
use crate::constants::{MIN_BRANCH_LENGTH, TINY_NUMBER};
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::graph::node::NodeType;
use crate::gtr::gtr::GTR;
use crate::seq_utils::seq2prof::{prof2seq, seq2prof, Prof2SeqParams, Prof2SeqResult};
use crate::utils::ndarray::{
  argmax_axis, choose1, choose2, exp, log, max_axis, maximum_scalar, sanitize_in_place, zeros,
};
use eyre::Report;
use itertools::Itertools;
use ndarray::{s, stack, Array2, ArrayBase, Axis};
use rand::Rng;
use std::fmt::Display;

pub fn run_anc_method_ml_joint(
  sequence_data: &SequenceData,
  gtr: &GTR,
  graph: &mut AncestralGraph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  traverse_backward(sequence_data, gtr, graph, rng, ancestral_args, ancestral_params);
  compute_roots(sequence_data, gtr, graph, rng, ancestral_args, ancestral_params)?;
  traverse_forward(sequence_data, gtr, graph, rng, ancestral_args, ancestral_params);
  Ok(())
}

fn traverse_backward(
  sequence_data: &SequenceData,
  gtr: &GTR,
  graph: &mut AncestralGraph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) {
  let L = sequence_data.len_compressed();
  let n_states = gtr.alphabet().len();

  graph.par_iter_breadth_first_backward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: node,
       children,
     }| {
      if is_root {
        return GraphTraversalContinuation::Continue;
      }

      let branch_length = {
        // TODO: Where does `node.branch_length` comes from on the graph?
        let branch_length = 1.0_f64;
        f64::max(MIN_BRANCH_LENGTH * sequence_data.one_mutation(), branch_length)
      };

      // transition matrix from parent states to the current node states.
      // denoted as Pij(i), where j - parent state, i - node state
      let Qt = gtr.expQt(branch_length);
      let log_transitions = maximum_scalar(&Qt, TINY_NUMBER);

      let msg_from_children: Array2<f64> = match &node.node_type {
        NodeType::Leaf(name) => {
          let mut msg_from_children = match sequence_data.get_compressed(name) {
            Some(compressed_alignment) => {
              let prof = seq2prof(&compressed_alignment, gtr.alphabet()).unwrap();
              log(&maximum_scalar(&prof, TINY_NUMBER))
            }
            None => zeros((L, n_states)),
          };
          sanitize_in_place(&mut msg_from_children);
          msg_from_children
        }
        NodeType::Internal(weight) | NodeType::Root(weight) => {
          let child_joint_Lxs = children
            .into_iter()
            .map(|(child, edge)| {
              let child = child.read();
              // TODO(perf): avoid this copy. It seems we need a version of `stack()` function accepting an iterator, not ` &[ArrayView<A, D>]`.
              child.joint_Lx.clone()
            })
            .collect_vec();
          let child_joint_Lx_views = child_joint_Lxs.iter().map(ArrayBase::view).collect_vec();
          let joint_Lx = stack(Axis(0), &child_joint_Lx_views).unwrap();
          joint_Lx.sum_axis(Axis(0))
        }
      };

      // For every possible state of the parent node, get the best state of the current node
      // and compute the likelihood of this state.
      // Preallocate storage:
      node.joint_Lx = zeros((L, n_states));
      node.joint_Cx = zeros((L, n_states));

      for (char_i, char) in gtr.alphabet().iter().enumerate() {
        // Pij(i) * L_ch(i) for given parent state j
        // if the node has a mask, P_ij is uniformly 1 at masked positions as no info is propagated
        let msg_to_parent = match &node.mask {
          Some(mask) => {
            unimplemented!("node.mask is not yet implemented");
            // ((log_transitions[:,char_i]*np.repeat([node.mask], self.gtr.n_states, axis=0).T) + msg_from_children)
          }
          None => &log_transitions.slice(s![.., char_i]).t() + &msg_from_children.view(),
        };

        // For this parent state, choose the best state of the current node
        let joint_Cx_target = &mut node.joint_Cx.slice_mut(s![.., char_i]);
        let joint_Cx_values = argmax_axis(&msg_to_parent, Axis(1));
        joint_Cx_target.assign(&joint_Cx_values);

        // and compute the likelihood of the best state of the current node
        // given the state of the parent (char_i) -- at masked position, there is no contribution
        let joint_Lx_target = &mut node.joint_Lx.slice_mut(s![.., char_i]);
        let joint_Lx_values = max_axis(&msg_to_parent, Axis(1));
        joint_Lx_target.assign(&joint_Lx_values);

        if let Some(mask) = &node.mask {
          unimplemented!("node.mask is not yet implemented");
          // joint_Lx_target *= mask
        }
      }

      GraphTraversalContinuation::Continue
    },
  );
}

fn compute_roots(
  sequence_data: &SequenceData,
  model: &GTR,
  graph: &mut AncestralGraph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  // root node profile = likelihood of the total tree
  let msg_from_children = {
    let non_root_joint_Lxs = graph.filter_map(|node| {
      if node.is_root {
        None
      } else {
        // TODO(perf): avoid this copy. It seems we need a version of `stack()` function accepting an iterator, not ` &[ArrayView<A, D>]`.
        Some(node.payload.read().joint_Lx.clone())
      }
    });
    let non_root_joint_Lx_views = non_root_joint_Lxs.iter().map(ArrayBase::view).collect_vec();
    let joint_Lx = stack(Axis(0), &non_root_joint_Lx_views).unwrap();
    joint_Lx.sum_axis(Axis(0))
  };

  // TODO: consider case when there's multiple roots
  let mut roots = graph.get_roots();
  if roots.len() > 1 {
    // TODO: consider case when there's multiple roots
    unimplemented!("The case of multiple roots is not implemented");
  }

  let root = &mut roots[0];
  let root = root.write();
  let root = root.payload();
  let mut root = root.write();

  // Pi(i) * Prod_ch Lch(i)
  let root_joint_Lx = &msg_from_children + &log(&model.pi).t();
  root.joint_Lx = root_joint_Lx;

  let normalized_profile = (&root.joint_Lx.t() - &max_axis(&root.joint_Lx, Axis(1))).t().to_owned();

  // choose sequence characters from this profile.
  // treat root node differently to avoid piling up mutations on the longer branch
  let normalized_profile = exp(&normalized_profile);
  let Prof2SeqResult {
    prof_values,
    seq,
    seq_ii,
  } = prof2seq(
    &normalized_profile,
    model,
    rng,
    &Prof2SeqParams {
      should_sample_from_profile: ancestral_params.sample_from_profile,
      should_normalize_profile: false,
    },
  )?;

  // compute the likelihood of the most probable root sequence
  let root_joint_Lx_t = root.joint_Lx.t();
  let sequence_LH = choose2(&seq_ii, &root_joint_Lx_t);
  let sequence_joint_LH = (&sequence_LH * &sequence_data.multiplicity()).sum();
  root.seq = seq;
  root.seq_ii = seq_ii;

  Ok(())
}

fn traverse_forward(
  sequence_data: &SequenceData,
  model: &GTR,
  graph: &mut AncestralGraph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) {
  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       is_root,
       is_leaf,
       key,
       payload: node,
       parents,
     }| {
      // root node has no mutations, everything else has been already set
      if is_root {
        return GraphTraversalContinuation::Continue;
      }

      if is_leaf && !ancestral_args.reconstruct_tip_states {
        return GraphTraversalContinuation::Continue;
      }

      if parents.len() > 1 {
        unimplemented!("Multiple parent nodes not handled yet");
      }

      let (parent, _) = &parents[0];
      let parent = parent.read();

      // choose the value of the Cx(i), corresponding to the state of the
      // parent node i. This is the state of the current node
      node.seq_ii = choose2(&parent.seq_ii, &node.joint_Cx.t());
      node.seq = choose1(&node.seq_ii, &model.alphabet().alphabet);

      GraphTraversalContinuation::Continue
    },
  );
}
