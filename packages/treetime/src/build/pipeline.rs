use crate::alphabet::alphabet::Alphabet;
use crate::build::insert::{Insertee, insert};
use crate::build::mat::Mat;
use crate::build::placement::{Placement, Target, place_exhaustive, place_from, place_greedy};
use crate::build::query::Query;
use crate::build::rebuild::{RebuildParams, RebuildStats, rebuild_round};
use crate::build::references::{Reaches, ReferenceSelection, select_references};
use crate::build::refresh::fitch_refresh;
use crate::build::reprune::{RepruneParams, RepruneStats, reprune};
use crate::build::spr::{SprStats, spr_round};
use crate::make_error;
use eyre::Report;
use itertools::{Itertools, izip};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::reroot::remove_stem_root;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::cast::ToF64Exact;

#[derive(Clone, Debug, PartialEq)]
pub enum SearchMode {
  Exhaustive,
  Walk {
    references: ReferenceSelection,
    slack: usize,
    exact: bool,
  },
  FromReferences {
    references: ReferenceSelection,
  },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleOrder {
  Input,
  FewestMissingFirst,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BuildParams {
  pub search: SearchMode,
  pub batch_size: usize,
  pub order: SampleOrder,
  pub spr_rounds: usize,
  pub spr_slack: Option<usize>,
  pub fitch_refresh: bool,
  pub compare_exhaustive: bool,
  pub rebuild: Option<RebuildParams>,
  pub reprune: Option<RepruneParams>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ExhaustiveComparison {
  pub placed: usize,
  pub worse: usize,
  pub excess: usize,
  pub max_excess: usize,
  pub visited_search: usize,
  pub visited_exhaustive: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BuildStats {
  pub placement_parsimony: usize,
  pub final_parsimony: usize,
  pub replaced_in_batch: usize,
  pub visited: usize,
  pub references: usize,
  pub comparison: Option<ExhaustiveComparison>,
  pub spr: Vec<SprStats>,
  pub rebuild: Vec<RebuildStats>,
  pub reprune: Option<RepruneStats>,
  pub collapsed_edges: usize,
}

#[derive(Debug)]
pub struct BuiltMat {
  pub mat: Mat,
  pub stats: BuildStats,
}

pub fn build_mat(seqs: &[Seq], alphabet: &Alphabet, params: &BuildParams) -> Result<BuiltMat, Report> {
  if params.batch_size == 0 {
    return make_error!("The placement batch size must be positive");
  }
  let reference = consensus(seqs, alphabet)?;
  let queries = encode_queries(seqs, &reference, alphabet);
  let order = sample_order(&queries, params.order);
  let mut mat = Mat::new(reference);
  let mut stats = BuildStats::default();
  place_samples(&mut mat, &queries, &order, params, &mut stats)?;
  stats.placement_parsimony = mat.parsimony();
  refine_mat(mat, seqs, alphabet, params, stats)
}

pub(crate) fn place_samples(
  mat: &mut Mat,
  queries: &[Query],
  order: &[usize],
  params: &BuildParams,
  stats: &mut BuildStats,
) -> Result<(), Report> {
  let Some((&first, rest)) = order.split_first() else {
    return make_error!("Cannot build a tree without sequences");
  };
  let root = mat.root();
  insert(mat, &queries[first], Target::Node(root), Insertee::Sample(first))?;

  for batch in rest.chunks(params.batch_size) {
    let search = BatchSearch::new(mat, &params.search);
    stats.references += search.n_references();
    let placements = batch
      .par_iter()
      .map(|&sample| search.place(mat, &queries[sample]))
      .collect::<Result<Vec<_>, Report>>()?;
    if params.compare_exhaustive {
      let exhaustive = batch
        .par_iter()
        .map(|&sample| place_exhaustive(mat, &queries[sample], None, None, false))
        .collect::<Result<Vec<_>, Report>>()?;
      let comparison = stats.comparison.get_or_insert_default();
      for (search, oracle) in izip!(&placements, &exhaustive) {
        let excess = search.cost - oracle.cost;
        comparison.placed += 1;
        comparison.worse += usize::from(excess > 0);
        comparison.excess += excess;
        comparison.max_excess = comparison.max_excess.max(excess);
        comparison.visited_search += search.visited;
        comparison.visited_exhaustive += oracle.visited;
      }
    }

    let mut touched: BTreeSet<GraphNodeKey> = BTreeSet::new();
    let mut split_into: BTreeMap<GraphEdgeKey, GraphNodeKey> = BTreeMap::new();
    for (&sample, placement) in izip!(batch, placements) {
      let query = &queries[sample];
      let placement = match stale_starts(mat, placement.target, &touched, &split_into) {
        None => placement,
        Some(local) => {
          stats.replaced_in_batch += 1;
          search.replace(mat, query, local)?
        },
      };
      stats.visited += placement.visited;
      let inserted = insert(mat, query, placement.target, Insertee::Sample(sample))?;
      touched.extend(inserted.touched);
      if let (Target::Edge(edge), Some(node)) = (placement.target, inserted.created) {
        split_into.insert(edge, node);
      }
    }
  }
  Ok(())
}

pub fn mat_from_tree(
  graph: Graph,
  samples: BTreeMap<GraphNodeKey, usize>,
  seqs: &[Seq],
  alphabet: &Alphabet,
) -> Result<Mat, Report> {
  let root = graph.root_key()?;
  let reference = consensus(seqs, alphabet)?;
  fitch_refresh(Mat::from_parts(graph, root, reference, vec![], samples), seqs, alphabet)
}

pub fn refine_mat(
  mut mat: Mat,
  seqs: &[Seq],
  alphabet: &Alphabet,
  params: &BuildParams,
  mut stats: BuildStats,
) -> Result<BuiltMat, Report> {
  let mut queries = encode_queries(seqs, mat.reference(), alphabet);
  if stats.placement_parsimony == 0 {
    stats.placement_parsimony = mat.parsimony();
  }
  for _ in 0..params.spr_rounds {
    if params.fitch_refresh {
      mat = fitch_refresh(mat, seqs, alphabet)?;
      queries = encode_queries(seqs, mat.reference(), alphabet);
    }
    let mut round = spr_round(&mut mat, &queries, params.spr_slack)?;
    if let Some(rebuild) = &params.rebuild {
      stats.rebuild.push(rebuild_round(&mut mat, seqs, alphabet, rebuild)?);
    }
    round.parsimony_before = stats
      .spr
      .last()
      .map_or(stats.placement_parsimony, |last| last.parsimony_after);
    round.parsimony_after = mat.parsimony();
    stats.spr.push(round);
    if round.parsimony_after >= round.parsimony_before {
      break;
    }
  }
  if params.fitch_refresh {
    mat = fitch_refresh(mat, seqs, alphabet)?;
  }
  if let Some(reprune_params) = &params.reprune {
    let (repruned, reprune_stats) = reprune(mat, seqs, alphabet, reprune_params)?;
    mat = repruned;
    stats.reprune = Some(reprune_stats);
  }
  stats.collapsed_edges = mat.collapse_empty_internal_edges()?;
  stats.final_parsimony = mat.parsimony();
  Ok(BuiltMat { mat, stats })
}

fn encode_queries(seqs: &[Seq], reference: &Seq, alphabet: &Alphabet) -> Vec<Query> {
  seqs
    .par_iter()
    .map(|seq| Query::from_seq(seq, reference, alphabet))
    .collect()
}

enum BatchSearch {
  Exhaustive,
  Walk {
    starts: Vec<GraphNodeKey>,
    slack: usize,
    exact: bool,
  },
  FromReferences {
    starts: Vec<(GraphNodeKey, Vec<(usize, AsciiChar)>)>,
    reaches: Reaches,
  },
}

impl BatchSearch {
  fn new(mat: &Mat, search: &SearchMode) -> Self {
    match search {
      SearchMode::Exhaustive => Self::Exhaustive,
      SearchMode::Walk {
        references,
        slack,
        exact,
      } => {
        let mut starts = select_references(mat, references);
        starts.push(mat.root());
        Self::Walk {
          starts,
          slack: *slack,
          exact: *exact,
        }
      },
      SearchMode::FromReferences { references } => {
        let starts = select_references(mat, references)
          .into_iter()
          .chain([mat.root()])
          .map(|node| (node, mat.node_state(node)))
          .collect();
        Self::FromReferences {
          starts,
          reaches: Reaches::new(mat, &BTreeSet::new()),
        }
      },
    }
  }

  fn n_references(&self) -> usize {
    match self {
      Self::Exhaustive => 0,
      Self::Walk { starts, .. } => starts.len(),
      Self::FromReferences { starts, .. } => starts.len(),
    }
  }

  fn place(&self, mat: &Mat, query: &Query) -> Result<Placement, Report> {
    match self {
      Self::Exhaustive => place_exhaustive(mat, query, None, None, true),
      Self::Walk { starts, slack, exact } => walk(mat, query, starts, *slack, *exact),
      Self::FromReferences { starts, reaches } => {
        let start = starts
          .iter()
          .min_by_key(|(node, state)| (query.cost_against(state, mat.reference()), *node))
          .map_or_else(|| mat.root(), |&(node, _)| node);
        place_from(mat, query, start, reaches)
      },
    }
  }

  fn replace(&self, mat: &Mat, query: &Query, local: Vec<GraphNodeKey>) -> Result<Placement, Report> {
    match self {
      Self::Walk { starts, slack, exact } => {
        let starts = local.into_iter().chain(starts.iter().copied()).collect_vec();
        walk(mat, query, &starts, *slack, *exact)
      },
      Self::Exhaustive | Self::FromReferences { .. } => place_exhaustive(mat, query, None, None, true),
    }
  }
}

fn walk(mat: &Mat, query: &Query, starts: &[GraphNodeKey], slack: usize, exact: bool) -> Result<Placement, Report> {
  let greedy = place_greedy(mat, query, starts, slack, None)?;
  if !exact {
    return Ok(greedy);
  }
  let exact = place_exhaustive(mat, query, None, Some(greedy), true)?;
  Ok(Placement {
    visited: greedy.visited + exact.visited,
    ..exact
  })
}

fn stale_starts(
  mat: &Mat,
  target: Target,
  touched: &BTreeSet<GraphNodeKey>,
  split_into: &BTreeMap<GraphEdgeKey, GraphNodeKey>,
) -> Option<Vec<GraphNodeKey>> {
  match target {
    Target::Node(node) => touched.contains(&node).then(|| vec![node]),
    Target::Edge(edge) if mat.has_edge(edge) => {
      let (parent, child) = mat.edge_endpoints(edge);
      (touched.contains(&parent) || touched.contains(&child)).then(|| vec![parent, child])
    },
    Target::Edge(edge) => Some(split_into.get(&edge).copied().into_iter().collect()),
  }
}

fn sample_order(queries: &[Query], order: SampleOrder) -> Vec<usize> {
  match order {
    SampleOrder::Input => (0..queries.len()).collect(),
    SampleOrder::FewestMissingFirst => (0..queries.len())
      .sorted_by_key(|&sample| (queries[sample].n_missing(), sample))
      .collect(),
  }
}

pub fn consensus(seqs: &[Seq], alphabet: &Alphabet) -> Result<Seq, Report> {
  let canonical = alphabet.canonical().collect_vec();
  let Some(length) = seqs.first().map(Seq::len) else {
    return make_error!("Cannot compute the consensus of no sequences");
  };
  if let Some(seq) = seqs.iter().find(|seq| seq.len() != length) {
    return make_error!(
      "All sequences must have the same length, found {} and {length}",
      seq.len()
    );
  }
  let mut counts = vec![0_usize; length * canonical.len()];
  for seq in seqs {
    for (pos, state) in seq.as_slice().iter().enumerate() {
      if let Some(index) = canonical.iter().position(|c| c == state) {
        counts[pos * canonical.len() + index] += 1;
      }
    }
  }
  Ok(
    counts
      .chunks(canonical.len())
      .map(|site| canonical[site.iter().position_max().unwrap_or(0)])
      .collect(),
  )
}

#[derive(Debug)]
pub struct BuiltTree {
  pub graph: Graph,
  pub names: NodeNames,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

pub fn mat_to_tree(mat: Mat, sample_names: &[String]) -> Result<BuiltTree, Report> {
  let length = mat.reference().len().to_f64_exact();
  let root = mat.root();
  let (mut graph, subs, samples) = mat.into_parts();
  graph.build()?;
  remove_stem_root(&mut graph, root)?;
  let input_names = samples
    .iter()
    .map(|(&node, &sample)| (node, sample_names.get(sample).cloned()))
    .collect();
  let names = NodeNames::from_input(input_names, &graph)?.names;
  let branch_lengths = graph
    .get_edges()
    .map(|edge| {
      let count = subs.get(edge.key().as_usize()).map_or(0, Vec::len);
      (edge.key(), Some(count.to_f64_exact() / length))
    })
    .collect();
  Ok(BuiltTree {
    graph,
    names,
    branch_lengths,
  })
}
