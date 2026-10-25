use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::mask::create_mask;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::build::evaluate::fitch_edge_subs;
use crate::build::pipeline::BuiltTree;
use crate::gtr::get_gtr::GtrModel;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::make_internal_report;
use crate::optimize::indel::IndelPolicy;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
use crate::optimize::partition_rates::PartitionRates;
use crate::optimize::pipeline::{OptimizeParams, optimized_partitions};
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, CandidateSites, SiteSlot};
use crate::partition::marginal::sparse::message::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::resolve::{PartitionEntry, resolve_partitions};
use crate::partition::set::PartitionSet;
use crate::partition::sites::SiteRange;
use crate::partition::split::SequenceInput;
use crate::partition::storage::sparse::SparseSeqDistribution;
use crate::progress::NoopProgress;
use crate::seq::alignment::{NodeSeqInput, SequenceSource};
use eyre::Report;
use itertools::{EitherOrBoth, Itertools, izip};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::{AsciiChar, LogLh, Seq};
use treetime_utils::interval::range::range_contains;

/// Likelihood SPR on the sparse marginal messages of a fixed tree, in rounds: optimize the branch lengths
/// and partition rates, search a new place for every subtree in parallel, apply the non-conflicting
/// improving moves, and start again from the new tree
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LikelihoodParams {
  pub rounds: usize,
  /// Number of branches the search for a new place walks away from the old one
  pub max_depth: usize,
  /// The search stops in a direction where a place scores this many log-likelihood units below the best
  pub drop: f64,
  pub min_gain: f64,
  /// Split the sites into partitions with their own rate, by the number of Fitch changes on the input tree
  pub rate_classes: bool,
  pub optimize_iterations: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct LikelihoodRound {
  pub log_lh: f64,
  pub tried: usize,
  pub improving: usize,
  pub moves: usize,
  pub predicted_gain: f64,
  /// Largest difference between the log-likelihood of a subtree put back at its old place and that of the tree
  pub max_reattach_error: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct LikelihoodStats {
  pub rounds: Vec<LikelihoodRound>,
  pub partition_sites: Vec<usize>,
  pub final_log_lh: f64,
}

pub fn refine_likelihood(
  tree: BuiltTree,
  seqs: &[Seq],
  sample_names: &[String],
  alphabet: &Alphabet,
  params: &LikelihoodParams,
) -> Result<(BuiltTree, LikelihoodStats), Report> {
  let BuiltTree {
    mut graph,
    names,
    mut branch_lengths,
  } = tree;
  let index: BTreeMap<&str, usize> = sample_names.iter().map(String::as_str).zip(0..).collect();
  let leaves: BTreeMap<GraphNodeKey, (String, usize)> = graph
    .get_leaves()
    .map(|leaf| {
      let name = names[leaf.key()].to_owned();
      let sample = index
        .get(name.as_str())
        .copied()
        .ok_or_else(|| make_internal_report!("Leaf '{name}' has no sequence"))?;
      Ok((leaf.key(), (name, sample)))
    })
    .collect::<Result<_, Report>>()?;
  let length = seqs.first().map_or(0, Seq::len);
  let classes = if params.rate_classes {
    let named: BTreeMap<&str, &Seq> = sample_names.iter().map(String::as_str).zip(seqs).collect();
    let mut changes = vec![0_usize; length];
    for sub in fitch_edge_subs(&graph, &names, &named, alphabet)?.values().flatten() {
      changes[sub.pos()] += 1;
    }
    changes.into_iter().map(rate_class).collect_vec()
  } else {
    vec![0; length]
  };
  let mut stats = LikelihoodStats {
    partition_sites: classes
      .iter()
      .copied()
      .counts()
      .into_iter()
      .sorted()
      .map(|(_, n)| n)
      .collect(),
    ..LikelihoodStats::default()
  };

  let mut best: Option<(f64, Graph, BTreeMap<GraphEdgeKey, Option<f64>>)> = None;
  for round in 0..=params.rounds {
    let names = leaf_names(&graph, &leaves)?;
    let input = sequence_input(&graph, &leaves, seqs, alphabet, &classes)?;
    let (partitions, lengths, log_lh) = optimized_partitions(
      &optimize_params(params),
      &mut graph,
      names,
      input,
      branch_lengths,
      &NoopProgress,
    )?;
    branch_lengths = lengths;
    let log_lh = log_lh.value();
    if let Some((best_lh, best_graph, best_lengths)) = &best
      && log_lh < *best_lh
    {
      graph = best_graph.clone();
      branch_lengths = best_lengths.clone();
      stats.final_log_lh = *best_lh;
      break;
    }
    stats.final_log_lh = log_lh;
    best = Some((log_lh, graph.clone(), branch_lengths.clone()));
    if round == params.rounds {
      break;
    }
    let PartitionSet::Sparse(mut parts) = partitions else {
      return Err(make_internal_report!("Likelihood SPR needs sparse partitions"));
    };
    let lengths = branch_lengths_or_zero(&branch_lengths);
    let states = parts
      .iter()
      .map(|part| FitchStates::new(part, &graph))
      .collect::<Result<Vec<_>, Report>>()?;
    let corrections = Search {
      parts: &parts,
      states: &states,
      graph: &graph,
      lengths: &lengths,
      params,
      offsets: Offsets::default(),
    }
    .calibrations()?;
    for (edge, shifts) in corrections {
      for (part, shift) in izip!(&mut parts, shifts) {
        if let Some(forward) = part.edges.forward.get_mut(&edge) {
          forward.msg_to_child.log_lh += LogLh::new(shift);
        }
      }
    }
    let mut search = Search {
      parts: &parts,
      states: &states,
      graph: &graph,
      lengths: &lengths,
      params,
      offsets: Offsets::default(),
    };
    search.offsets = search.offsets(log_lh)?;
    let search = search;
    let candidates = graph
      .get_nodes()
      .map(|node| node.key())
      .filter(|&node| !graph.is_root(node))
      .collect_vec()
      .into_par_iter()
      .map(|node| search.best_move(node))
      .collect::<Result<Vec<_>, Report>>()?;
    let mut summary = LikelihoodRound {
      log_lh,
      tried: candidates.len(),
      ..LikelihoodRound::default()
    };
    let mut moves = vec![];
    for found in candidates.into_iter().flatten() {
      summary.max_reattach_error = summary.max_reattach_error.max((found.original - log_lh).abs());
      if found.gain > params.min_gain {
        summary.improving += 1;
        moves.push(found);
      }
    }
    moves.sort_by(|a, b| b.gain.total_cmp(&a.gain));
    let mut touched = BTreeSet::new();
    for found in moves {
      if apply_move(&mut graph, &mut branch_lengths, &found, &mut touched)? {
        summary.moves += 1;
        summary.predicted_gain += found.gain;
      }
    }
    let moved = summary.moves;
    stats.rounds.push(summary);
    if moved == 0 {
      break;
    }
  }
  let names = leaf_names(&graph, &leaves)?;
  Ok((
    BuiltTree {
      graph,
      names,
      branch_lengths,
    },
    stats,
  ))
}

fn children(graph: &Graph, node: GraphNodeKey) -> Result<Vec<(GraphNodeKey, GraphEdgeKey)>, Report> {
  let node = graph
    .get_node(node)
    .ok_or_else(|| make_internal_report!("Node {node} not found"))?;
  Ok(graph.children_keys_of(node).collect())
}

/// Classes of 0, 1 to 2, 3 to 8, 9 to 26, and 27 or more changes
fn rate_class(changes: usize) -> usize {
  match changes {
    0 => 0,
    1..=2 => 1,
    3..=8 => 2,
    9..=26 => 3,
    _ => 4,
  }
}

fn optimize_params(params: &LikelihoodParams) -> OptimizeParams {
  OptimizeParams {
    dense: Some(false),
    initial_guess: InitialGuessMode::Never,
    reroot_spec: None,
    loop_settings: OptimizeLoopSettings {
      max_iter: params.optimize_iterations,
      dp: 0.1,
      damping: 0.0,
      opt_method: BranchOptMethod::Brent,
      indel_policy: IndelPolicy::Exclude,
      topology_ops: TopologyOps {
        collapse_short_branches: false,
        merge_siblings: false,
        flip_parent_child: false,
      },
      partition_rates: PartitionRates::Estimate,
    },
  }
}

fn leaf_names(graph: &Graph, leaves: &BTreeMap<GraphNodeKey, (String, usize)>) -> Result<NodeNames, Report> {
  let input = leaves
    .iter()
    .map(|(&leaf, (name, _))| (leaf, Some(name.clone())))
    .collect();
  Ok(NodeNames::from_input(input, graph)?.names)
}

fn sequence_input(
  graph: &Graph,
  leaves: &BTreeMap<GraphNodeKey, (String, usize)>,
  seqs: &[Seq],
  alphabet: &Alphabet,
  classes: &[usize],
) -> Result<SequenceInput, Report> {
  let nodes: BTreeMap<GraphNodeKey, NodeSeqInput> = graph
    .get_leaves()
    .map(|leaf| {
      let (name, sample) = &leaves[&leaf.key()];
      (
        leaf.key(),
        NodeSeqInput {
          name: name.clone(),
          seq: Some(seqs[*sample].clone()),
        },
      )
    })
    .collect();
  let length = classes.len();
  let mask = create_mask(nodes.values().filter_map(|node| node.seq.as_ref()), length, alphabet);
  let entries = classes
    .iter()
    .copied()
    .sorted()
    .dedup()
    .map(|class| {
      let sites = izip!(1.., classes)
        .filter(|&(_, &c)| c == class)
        .map(|(site, _)| site)
        .collect_vec();
      let ranges = sites
        .chunk_by(|a, b| a + 1 == *b)
        .map(|run| format!("{}-{}", run[0], run[run.len() - 1]).parse::<SiteRange>())
        .collect::<Result<Vec<_>, _>>()?;
      Ok(PartitionEntry {
        name: Some(format!("rate{class}")),
        alignment: None,
        sites: Some(ranges),
        model: GtrModel::Infer,
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let entries = if entries.len() > 1 { entries } else { vec![] };
  let plans = resolve_partitions(&entries, &[length], &GtrModel::Infer, Some(alphabet.name()))?.plans;
  Ok(SequenceInput {
    alignments: vec![SequenceSource {
      path: PathBuf::from("built.fasta"),
      alphabet: alphabet.clone(),
      length,
      nodes,
      mask,
    }],
    partitions: plans,
  })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Place {
  Node(GraphNodeKey),
  Edge(GraphEdgeKey),
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Move {
  node: GraphNodeKey,
  place: Place,
  gain: f64,
  /// Uncalibrated score of the old place, for the consistency check against the tree
  original: f64,
}

/// Sparse distributions of one message, one per partition, with fixed vectors keyed by the Fitch states of
/// `key`. A message from a subtree contributes nothing where the subtree has no data (`mask`), as in the
/// backward pass; a message that contains the root prior has `prior`.
#[derive(Clone, Debug)]
struct Msg {
  dis: Vec<SparseSeqDistribution>,
  key: GraphNodeKey,
  prior: bool,
  mask: bool,
}

/// Fitch states of every node of a partition: the root sequence and, per node, the sorted positions where
/// the node differs from it
struct FitchStates {
  root: Seq,
  diffs: BTreeMap<GraphNodeKey, Vec<(usize, AsciiChar)>>,
}

impl FitchStates {
  /// Applies the Fitch substitutions and the indels of each edge, as `parsimony_seq()` does
  fn new(part: &SparseReconstruction, graph: &Graph) -> Result<Self, Report> {
    let root = part.partition.root_sequence.clone();
    let gap = part.partition.alphabet.gap();
    let schedule = GraphSchedule::new(graph)?;
    let mut diffs: BTreeMap<GraphNodeKey, Vec<(usize, AsciiChar)>> = BTreeMap::new();
    for &node in schedule.preorder() {
      let own = match schedule.parents(node).first() {
        None => vec![],
        Some(&(parent, edge)) => {
          let edge = &part.partition.obs_edges[&edge];
          let mut states: BTreeMap<usize, AsciiChar> = diffs[&parent].iter().copied().collect();
          let changes = edge
            .fitch_subs()
            .iter()
            .map(|sub| (sub.pos(), sub.qry()))
            .chain(edge.indels.iter().flat_map(|indel| {
              (indel.range.0..indel.range.1).map(move |pos| {
                let state = if indel.is_deletion() {
                  gap
                } else {
                  indel.seq[pos - indel.range.0]
                };
                (pos, state)
              })
            }));
          for (pos, state) in changes {
            if state == root[pos] {
              states.remove(&pos);
            } else {
              states.insert(pos, state);
            }
          }
          states.into_iter().collect()
        },
      };
      diffs.insert(node, own);
    }
    Ok(Self { root, diffs })
  }

  /// Positions where the Fitch states of `a` and `b` are different nucleotides, with the state of `a` and of
  /// `b`. These are the positions of the Fitch substitutions between them.
  fn differences(&self, a: GraphNodeKey, b: GraphNodeKey, alphabet: &Alphabet) -> Vec<(usize, AsciiChar, AsciiChar)> {
    self.diffs[&a]
      .iter()
      .merge_join_by(&self.diffs[&b], |(pa, _), (pb, _)| pa.cmp(pb))
      .filter_map(|pair| match pair {
        EitherOrBoth::Both(&(pos, sa), &(_, sb)) => (sa != sb).then_some((pos, sa, sb)),
        EitherOrBoth::Left(&(pos, sa)) => Some((pos, sa, self.root[pos])),
        EitherOrBoth::Right(&(pos, sb)) => Some((pos, self.root[pos], sb)),
      })
      .filter(|&(_, sa, sb)| alphabet.is_canonical(sa) && alphabet.is_canonical(sb))
      .collect()
  }
}

struct Search<'a> {
  parts: &'a [SparseReconstruction],
  states: &'a [FitchStates],
  graph: &'a Graph,
  lengths: &'a BTreeMap<GraphEdgeKey, f64>,
  params: &'a LikelihoodParams,
  offsets: Offsets,
}

/// How far the combination at a node, or at the middle of an edge, of the stored messages of the tree is
/// from the log-likelihood of the tree. Subtracting them makes the scores of all places comparable.
#[derive(Default)]
struct Offsets {
  nodes: BTreeMap<GraphNodeKey, f64>,
  edges: BTreeMap<GraphEdgeKey, f64>,
}

impl Search<'_> {
  /// Walks from the old place of the subtree at `node` through the tree without it. The messages along the
  /// walk are combined from the stored messages of the other edges, so they never contain the subtree.
  fn best_move(&self, node: GraphNodeKey) -> Result<Option<Move>, Report> {
    let graph = self.graph;
    let Some((parent, edge)) = graph.node_parent(node)? else {
      return Ok(None);
    };
    let pendant = self.lengths[&edge];
    let subtree = Msg {
      dis: self.parts_map(|part| Ok(part.edges.backward[&edge].msg_to_parent.clone()))?,
      key: node,
      prior: false,
      mask: true,
    };
    let subtree = self.propagate(&subtree, pendant, true);
    let remaining = children(graph, parent)?.len() - 1 + usize::from(!graph.is_root(parent));
    if remaining < 2 && !graph.is_root(parent) {
      return Ok(None);
    }
    let suppressed = remaining <= 2;

    let mut original = None;
    let mut best: Option<(f64, Place)> = None;
    let mut best_lh = f64::NEG_INFINITY;
    let mut stack: Vec<(GraphNodeKey, Option<(GraphNodeKey, Msg)>, usize)> = vec![(parent, None, 0)];
    while let Some((at, from, depth)) = stack.pop() {
      let from_key = from.as_ref().map(|(key, _)| *key);
      let mut msgs: Vec<(GraphNodeKey, Option<GraphEdgeKey>, Msg)> =
        from.into_iter().map(|(k, m)| (k, None, m)).collect();
      if let Some((up, up_edge)) = graph.node_parent(at)?
        && Some(up) != from_key
      {
        msgs.push((up, Some(up_edge), self.parent_msg(up, up_edge)?));
      }
      for (child, child_edge) in children(graph, at)? {
        if child != node && Some(child) != from_key {
          msgs.push((child, Some(child_edge), self.child_msg(child, child_edge)?));
        }
      }
      let is_root = graph.is_root(at);
      let mut all = msgs.iter().map(|(_, _, msg)| msg).collect_vec();
      all.push(&subtree);
      let (_, raw_lh) = self.combine(at, &all, is_root)?;
      let node_lh = raw_lh - self.offsets.nodes.get(&at).copied().unwrap_or(0.0);
      if at == parent {
        original = Some((node_lh, raw_lh));
        best_lh = best_lh.max(node_lh);
      } else if node_lh > best_lh {
        best_lh = node_lh;
        best = Some((node_lh, Place::Node(at)));
      }
      for (i, (next, next_edge, _)) in msgs.iter().enumerate() {
        let Some(next_edge) = *next_edge else {
          continue;
        };
        let others = msgs
          .iter()
          .enumerate()
          .filter(|&(j, _)| j != i)
          .map(|(_, (_, _, msg))| msg)
          .collect_vec();
        let (behind, _) = self.combine(at, &others, is_root)?;
        let length = self.lengths[&next_edge];
        let downward = graph.edge_endpoints(next_edge)?.0 == at;
        let (mid_lh, arriving) = if downward {
          let below = Msg {
            dis: self.parts_map(|part| Ok(part.edges.backward[&next_edge].msg_to_parent.clone()))?,
            key: *next,
            prior: false,
            mask: true,
          };
          let (_, lh) = self.combine(
            at,
            &[
              &self.propagate(&behind, length / 2.0, false),
              &self.propagate(&below, length / 2.0, true),
              &subtree,
            ],
            false,
          )?;
          (lh, self.propagate(&behind, length, false))
        } else {
          let above = Msg {
            dis: self.parts_map(|part| Ok(part.edges.forward[&next_edge].msg_to_child.clone()))?,
            key: *next,
            prior: true,
            mask: false,
          };
          let (_, lh) = self.combine(
            *next,
            &[
              &self.propagate(&above, length / 2.0, false),
              &self.propagate(&behind, length / 2.0, true),
              &subtree,
            ],
            false,
          )?;
          (lh, self.propagate(&behind, length, true))
        };
        let mid_lh = mid_lh - self.offsets.edges.get(&next_edge).copied().unwrap_or(0.0);
        if !(suppressed && at == parent) && mid_lh > best_lh {
          best_lh = mid_lh;
          best = Some((mid_lh, Place::Edge(next_edge)));
        }
        if !graph.is_leaf(*next) && depth < self.params.max_depth && mid_lh > best_lh - self.params.drop {
          stack.push((*next, Some((at, arriving)), depth + 1));
        }
      }
    }
    let (original, raw) =
      original.ok_or_else(|| make_internal_report!("The search did not score the old place of {node}"))?;
    Ok(Some(match best {
      Some((lh, place)) => Move {
        node,
        place,
        gain: lh - original,
        original: raw,
      },
      None => Move {
        node,
        place: Place::Node(parent),
        gain: 0.0,
        original: raw,
      },
    }))
  }

  /// Shift of the log-likelihood of each edge's message from the rest of the tree, per partition, that makes
  /// it combine with the subtree below the edge to the log-likelihood of the tree. The forward pass keeps the
  /// distributions but not this bookkeeping consistent on alignments with gaps.
  fn calibrations(&self) -> Result<Vec<(GraphEdgeKey, Vec<f64>)>, Report> {
    let graph = self.graph;
    let root = graph.root_key()?;
    let totals = self
      .parts
      .iter()
      .map(|part| part.node_states[&root].profile.log_lh.value())
      .collect_vec();
    graph
      .get_edges()
      .map(|edge| edge.key())
      .collect_vec()
      .into_par_iter()
      .map(|edge| {
        let (parent, child) = graph.edge_endpoints(edge)?;
        let mut msgs = vec![self.parent_msg(parent, edge)?];
        if graph.is_leaf(child) {
          msgs.push(Msg {
            dis: self.parts_map(|part| Ok(part.edges.backward[&edge].msg_to_parent.clone()))?,
            key: child,
            prior: false,
            mask: true,
          });
        } else {
          for (grandchild, child_edge) in children(graph, child)? {
            msgs.push(self.child_msg(grandchild, child_edge)?);
          }
        }
        let (combined, _) = self.combine(child, &msgs.iter().collect_vec(), false)?;
        let shifts = izip!(&totals, &combined.dis)
          .map(|(total, dis)| total - dis.log_lh.value())
          .collect();
        Ok((edge, shifts))
      })
      .collect()
  }

  fn offsets(&self, log_lh: f64) -> Result<Offsets, Report> {
    let graph = self.graph;
    let nodes = graph
      .get_nodes()
      .map(|node| node.key())
      .filter(|&node| !graph.is_leaf(node))
      .collect_vec()
      .into_par_iter()
      .map(|node| {
        let mut msgs = vec![];
        if let Some((parent, edge)) = graph.node_parent(node)? {
          msgs.push(self.parent_msg(parent, edge)?);
        }
        for (child, edge) in children(graph, node)? {
          msgs.push(self.child_msg(child, edge)?);
        }
        let (_, lh) = self.combine(node, &msgs.iter().collect_vec(), graph.is_root(node))?;
        Ok((node, lh - log_lh))
      })
      .collect::<Result<_, Report>>()?;
    let edges = graph
      .get_edges()
      .map(|edge| edge.key())
      .collect_vec()
      .into_par_iter()
      .map(|edge| {
        let (parent, child) = graph.edge_endpoints(edge)?;
        let half = self.lengths[&edge] / 2.0;
        let above = Msg {
          dis: self.parts_map(|part| Ok(part.edges.forward[&edge].msg_to_child.clone()))?,
          key: parent,
          prior: true,
          mask: false,
        };
        let below = Msg {
          dis: self.parts_map(|part| Ok(part.edges.backward[&edge].msg_to_parent.clone()))?,
          key: child,
          prior: false,
          mask: true,
        };
        let (_, lh) = self.combine(
          parent,
          &[
            &self.propagate(&above, half, false),
            &self.propagate(&below, half, true),
          ],
          false,
        )?;
        Ok((edge, lh - log_lh))
      })
      .collect::<Result<_, Report>>()?;
    Ok(Offsets { nodes, edges })
  }

  fn parts_map(
    &self,
    f: impl Fn(&SparseReconstruction) -> Result<SparseSeqDistribution, Report>,
  ) -> Result<Vec<SparseSeqDistribution>, Report> {
    self.parts.iter().map(f).collect()
  }

  /// Message of the subtree below `child`, at the parent end of its edge
  fn child_msg(&self, child: GraphNodeKey, edge: GraphEdgeKey) -> Result<Msg, Report> {
    Ok(Msg {
      dis: self.parts_map(|part| Ok(part.edges.backward[&edge].msg_from_child.clone()))?,
      key: child,
      prior: false,
      mask: true,
    })
  }

  /// Message of the tree outside the child of `edge`, at the child end of the edge
  fn parent_msg(&self, parent: GraphNodeKey, edge: GraphEdgeKey) -> Result<Msg, Report> {
    let above = Msg {
      dis: self.parts_map(|part| Ok(part.edges.forward[&edge].msg_to_child.clone()))?,
      key: parent,
      prior: true,
      mask: false,
    };
    Ok(self.propagate(&above, self.lengths[&edge], false))
  }

  fn propagate(&self, msg: &Msg, length: f64, upward: bool) -> Msg {
    let dis = izip!(self.parts, &msg.dis)
      .map(|(part, dis)| {
        let t = fix_branch_length(part.partition.total_length, length);
        if let Some(rates) = part.gtr.site_rates() {
          propagate_raw_per_site(&part.gtr, rates, t, upward, dis, None)
        } else {
          let exp_qt = part.gtr.expQt(t);
          let exp_qt = if upward { exp_qt.t().to_owned() } else { exp_qt };
          propagate_raw(&exp_qt, dis, None)
        }
      })
      .collect();
    Msg {
      dis,
      key: msg.key,
      prior: msg.prior,
      mask: msg.mask,
    }
  }

  /// Candidate sites of partition `i` for combining `msgs` at `at`
  fn sites(&self, i: usize, at: GraphNodeKey, msgs: &[&Msg]) -> CandidateSites {
    let part = &self.parts[i];
    let states = &self.states[i];
    let alphabet = &part.partition.alphabet;
    let mut events = CandidateSiteEvents::with_capacity(msgs.len(), 64);
    for (j, msg) in msgs.iter().enumerate() {
      if msg.key != at {
        for (pos, state, reference) in states.differences(at, msg.key, alphabet) {
          events.insert(pos, SiteSlot::State, state);
          events.insert(pos, SiteSlot::Reference(j), reference);
        }
      }
      for (&pos, var) in &msg.dis[i].variable {
        events.or_insert(pos, SiteSlot::State, var.state);
        events.or_insert(pos, SiteSlot::Reference(j), var.state);
      }
    }
    let mut sites = events.finish();
    sites.fill_missing_references(|_, _, state| state);
    sites.map_references(|j, pos, reference| {
      let obs = &part.partition.obs_nodes[&msgs[j].key];
      if !msgs[j].mask || !range_contains(&obs.non_char, pos) {
        reference
      } else if range_contains(&obs.gaps, pos) {
        alphabet.gap()
      } else {
        alphabet.unknown()
      }
    });
    sites
  }

  /// Product of messages at node `at`, and its log-likelihood summed over partitions. As in the marginal
  /// passes, the candidate sites are the explicit positions of the messages and the positions where the
  /// Fitch states keying a message differ from those of `at`.
  fn combine(&self, at: GraphNodeKey, msgs: &[&Msg], root: bool) -> Result<(Msg, f64), Report> {
    let prior = msgs.iter().any(|msg| msg.prior);
    let use_pi = root && !prior;
    let mut dis = Vec::with_capacity(self.parts.len());
    let mut log_lh = 0.0;
    for (i, part) in self.parts.iter().enumerate() {
      let alphabet = &part.partition.alphabet;
      let sites = self.sites(i, at, msgs);
      let messages = msgs.iter().map(|msg| &msg.dis[i]).collect_vec();
      let combined = combine_messages(
        &part.partition.obs_nodes[&at].composition,
        &messages,
        &sites,
        alphabet,
        use_pi.then_some(part.gtr.pi()),
      )?;
      log_lh += combined.log_lh.value();
      dis.push(combined);
    }
    Ok((
      Msg {
        dis,
        key: at,
        prior: prior || use_pi,
        mask: false,
      },
      log_lh,
    ))
  }
}

/// Moves the subtree unless an earlier move of the round touched its nodes, its target, or their ancestors.
/// The old parent is removed when it is left with one child.
fn apply_move(
  graph: &mut Graph,
  lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  found: &Move,
  touched: &mut BTreeSet<GraphNodeKey>,
) -> Result<bool, Report> {
  let node = found.node;
  if touched.contains(&node) || graph.get_node(node).is_none() {
    return Ok(false);
  }
  let Some((parent, edge)) = graph.node_parent(node)? else {
    return Ok(false);
  };
  let target = match found.place {
    Place::Node(target) => {
      if graph.get_node(target).is_none() {
        return Ok(false);
      }
      vec![target]
    },
    Place::Edge(target) => {
      if graph.get_edge(target).is_none() {
        return Ok(false);
      }
      let (a, b) = graph.edge_endpoints(target)?;
      vec![a, b]
    },
  };
  let grandparent = graph.node_parent(parent)?.map(|(up, _)| up);
  if [node, parent]
    .iter()
    .chain(&target)
    .chain(&grandparent)
    .any(|key| touched.contains(key))
  {
    return Ok(false);
  }
  let pendant = lengths.remove(&edge).flatten().unwrap_or(0.0);
  graph.remove_edge(edge)?;
  match found.place {
    Place::Node(target) => {
      let new = graph.add_edge(target, node)?;
      lengths.insert(new, Some(pendant));
    },
    Place::Edge(target) => {
      let (a, b) = graph.edge_endpoints(target)?;
      let length = lengths.remove(&target).flatten().unwrap_or(0.0);
      graph.remove_edge(target)?;
      let middle = graph.add_node();
      lengths.insert(graph.add_edge(a, middle)?, Some(length / 2.0));
      lengths.insert(graph.add_edge(middle, b)?, Some(length / 2.0));
      lengths.insert(graph.add_edge(middle, node)?, Some(pendant));
      touched.insert(middle);
    },
  }
  let mut marked = vec![node, parent];
  marked.extend(&target);
  if let [(child, child_edge)] = children(graph, parent)?[..] {
    let child_length = lengths.remove(&child_edge).flatten().unwrap_or(0.0);
    match graph.node_parent(parent)? {
      Some((up, up_edge)) => {
        let up_length = lengths.remove(&up_edge).flatten().unwrap_or(0.0);
        graph.remove_node(parent)?;
        lengths.insert(graph.add_edge(up, child)?, Some(up_length + child_length));
        marked.push(up);
      },
      None => {
        graph.remove_node(parent)?;
      },
    }
    marked.push(child);
  }
  graph.build()?;
  for start in marked {
    let mut current = Some(start);
    while let Some(key) = current {
      if !touched.insert(key) && key != start {
        break;
      }
      current = if graph.get_node(key).is_some() {
        graph.node_parent(key)?.map(|(up, _)| up)
      } else {
        None
      };
    }
  }
  Ok(true)
}
