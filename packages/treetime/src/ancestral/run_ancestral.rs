use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::run_ancestral::NodeType::{Internal, Leaf};
use crate::cli::treetime_cli::{MethodAncestral, TreetimeAncestralArgs};
use crate::constants::{MIN_BRANCH_LENGTH, TINY_NUMBER};
use crate::graph::graph::{Graph as GenericGraph, GraphNodeBackward, GraphNodeForward, NodeEdgePair, Weighted};
use crate::gtr::get_gtr::get_gtr;
use crate::gtr::gtr::GTR;
use crate::io::fasta::{FastaRecord, FastaWriter};
use crate::io::file::create_file;
use crate::io::nwk::read_nwk;
use crate::seq_utils::seq2prof::{prof2seq, seq2prof, Prof2SeqParams, Prof2SeqResult};
use crate::utils::ndarray::{
  argmax_axis, choose1, choose2, exp, log, max_axis, maximum_scalar, sanitize_in_place, zeros,
};
use crate::utils::random::get_random_number_generator;
use color_eyre::Section;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::{array, s, stack, Array1, Array2, ArrayBase, Axis};
use rand::Rng;
use std::fmt::{Display, Formatter};
use std::path::Path;

#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Debug, Default)]
pub struct TreetimeAncestralParams {
  sample_from_profile: bool,
  fixed_pi: bool,
}

pub fn run_ancestral(ancestral_args: &TreetimeAncestralArgs) -> Result<(), Report> {
  let TreetimeAncestralArgs {
    input_fastas,
    aln,
    vcf_reference,
    tree,
    gtr,
    gtr_params,
    aa,
    keep_overhangs,
    zero_based,
    reconstruct_tip_states,
    report_ambiguous,
    method_anc,
    outdir,
    seed,
  } = ancestral_args;

  let mut rng = get_random_number_generator(*seed);

  // TODO: alphabet is hardcoded. Make it dynamic.
  let alphabet = Alphabet::new("nuc")?;

  let sequence_data = SequenceData::new(input_fastas, alphabet.ambiguous())?;

  let model = get_gtr(gtr)?;

  let mut graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree)?,
  };

  let ancestral_params = TreetimeAncestralParams::default();

  ml_anc(
    &sequence_data,
    &alphabet,
    &model,
    &mut graph,
    &mut rng,
    ancestral_args,
    &ancestral_params,
  )?;

  let fasta_records = get_sequences(&sequence_data, &graph, ancestral_args, &ancestral_params);

  let fasta_file = create_file(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);

  for FastaRecord { seq, seq_name, .. } in fasta_records {
    fasta_writer.write(&seq_name, &seq)?;
  }

  Ok(())
}

pub fn get_sequences(
  sequence_data: &SequenceData,
  graph: &Graph,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Vec<FastaRecord> {
  graph.filter_map(|node| {
    let payload = &node.payload.read();

    let seq = if let (Leaf(name), true) = (&payload.node_type, !ancestral_args.reconstruct_tip_states) {
      sequence_data
        .get_full(name)
        .ok_or_else(|| format!("Sequence not found: '{name}'"))
        .unwrap()
    } else {
      sequence_data.decompress(&payload.seq)
    }
    .iter()
    .join("");

    Some(FastaRecord {
      index: node.key,
      seq_name: payload.name.clone(),
      seq,
    })
  })
}

#[derive(Clone, Debug, PartialEq)]
pub enum NodeType {
  Leaf(String),
  Internal(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,
  pub joint_Lx: Array2<f64>,   // likelihood
  pub joint_Cx: Array2<usize>, // max likelihood indices
  pub mask: Option<usize>,
  pub seq: Array1<char>,
  pub seq_ii: Array1<usize>,
}

impl Node {
  pub fn leaf(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: Leaf(name.to_owned()),
      joint_Lx: Array2::<f64>::default((0, 0)),
      joint_Cx: Array2::<usize>::default((0, 0)),
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
    }
  }

  pub fn internal(name: &str, weight: f64) -> Self {
    Self {
      name: name.to_owned(),
      node_type: Internal(weight),
      joint_Lx: array![[]],
      joint_Cx: array![[]],
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
    }
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match &self.node_type {
      Leaf(name) => write!(f, "{name}"),
      Internal(weight) => write!(f, "{weight:1.4}"),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
  weight: f64,
}

impl Weighted for Edge {}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

type Graph = GenericGraph<Node, Edge>;

fn infer_graph() -> Result<Graph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

fn create_graph(tree_path: impl AsRef<Path>) -> Result<Graph, Report> {
  let nwk_tree = read_nwk(tree_path)
    .wrap_err("When parsing input tree")
    .with_section(|| "Note: only Newick format is currently supported")?;

  let mut graph = Graph::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, usize>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  let mut node_counter: usize = 0;
  for (nwk_idx, nwk_node) in nwk_tree.g.raw_nodes().iter().enumerate() {
    // Attempt to parse weight as float. If not a float, then it's a named leaf node, otherwise - internal node.
    let inserted_node_idx = if let Ok(weight) = nwk_node.weight.parse::<f64>() {
      let node = Node::internal(&format!("NODE_{node_counter:07}"), weight);
      node_counter += 1;
      graph.add_node(node)
    } else {
      let name = nwk_node.weight.as_str();
      if name == "N/A" {
        let node = Node::internal(&format!("NODE_{node_counter:07}"), 0.0);
        node_counter += 1;
        graph.add_node(node)
      } else {
        graph.add_node(Node::leaf(name))
      }
    };

    index_map.insert(nwk_idx, inserted_node_idx);
  }

  // Insert edges
  for (nwk_idx, nwk_edge) in nwk_tree.g.raw_edges().iter().enumerate() {
    let weight: f64 = nwk_edge.weight as f64;
    let source: usize = nwk_edge.source().index();
    let target: usize = nwk_edge.target().index();

    let source = index_map
      .get(&source)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {source} not found."))?;

    let target = index_map
      .get(&target)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {target} not found."))?;

    graph.add_edge(*source, *target, Edge { weight });
  }

  graph.build()?;

  Ok(graph)
}

fn ml_anc(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut Graph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  match ancestral_args.method_anc {
    MethodAncestral::MaximumLikelihoodJoint => ml_anc_joint(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
    MethodAncestral::MaximumLikelihoodMarginal => ml_anc_marginal(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
    MethodAncestral::Parsimony => ml_anc_fitch(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
  }
}

fn ml_anc_joint(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut Graph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  let L = sequence_data.len_compressed();
  let n_states = alphabet.len();

  graph.par_iter_breadth_first_backward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: node,
       children,
     }| {
      if is_root {
        return;
      }

      let branch_length = {
        // TODO: Where does `node.branch_length` comes from on the graph?
        let branch_length = 1.0_f64;
        f64::max(MIN_BRANCH_LENGTH * sequence_data.one_mutation(), branch_length)
      };

      // transition matrix from parent states to the current node states.
      // denoted as Pij(i), where j - parent state, i - node state
      let Qt = model.expQt(branch_length);
      let log_transitions = maximum_scalar(&Qt, TINY_NUMBER);

      let msg_from_children: Array2<f64> = match &node.node_type {
        Leaf(name) => {
          let mut msg_from_children = match sequence_data.get_compressed(name) {
            Some(compressed_alignment) => {
              let prof = seq2prof(&compressed_alignment, &model.profile_map).unwrap();
              log(&maximum_scalar(&prof, TINY_NUMBER))
            }
            None => zeros((L, n_states)),
          };
          sanitize_in_place(&mut msg_from_children);
          msg_from_children
        }
        Internal(weight) => {
          let child_joint_Lxs = children
            .into_iter()
            .map(|child| {
              let node = child.node.read();
              // TODO(perf): avoid this copy. It seems we need a version of `stack()` function accepting an iterator, not ` &[ArrayView<A, D>]`.
              node.joint_Lx.clone()
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

      for (char_i, char) in alphabet.iter().enumerate() {
        // Pij(i) * L_ch(i) for given parent state j
        // if the node has a mask, P_ij is uniformly 1 at masked positions as no info is propagated
        let msg_to_parent = match node.mask {
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
    },
  );

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
  {
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
  }

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
        return;
      }

      if is_leaf && !ancestral_args.reconstruct_tip_states {
        return;
      }

      if parents.len() > 1 {
        unimplemented!("Multiple parent nodes not handled yet");
      }

      let NodeEdgePair { edge, node: parent } = &parents[0];
      let parent = parent.read();

      // choose the value of the Cx(i), corresponding to the state of the
      // parent node i. This is the state of the current node
      node.seq_ii = choose2(&parent.seq_ii, &node.joint_Cx.t());
      node.seq = choose1(&node.seq_ii, &alphabet.alphabet);
    },
  );

  Ok(())
}

fn ml_anc_marginal(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut Graph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  unimplemented!("ml_anc_marginal: not yet implemented");
}

fn ml_anc_fitch(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut Graph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  unimplemented!("ml_anc_fitch: not yet implemented");
}
