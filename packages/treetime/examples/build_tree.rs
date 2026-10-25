use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;
use treetime::alphabet::alphabet::Alphabet;
use treetime::build::evaluate::{
  MissedSplit, TreeWithCounts, compare_splits, fitch_edge_counts, fitch_edge_subs, missed_splits,
};
use treetime::build::pipeline::{
  BuildParams, BuildStats, InitialTree, SampleOrder, ScaffoldParams, SearchMode, build_mat, mat_from_tree, mat_to_tree,
  refine_mat,
};
use treetime::build::rebuild::RebuildParams;
use treetime::build::references::ReferenceSelection;
use treetime::build::reprune::RepruneParams;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::GraphSchedule;
use treetime_io::fasta::fasta_read_file;
use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::make_error;

struct Args {
  aln: String,
  sites_out: Option<String>,
  flags: BTreeSet<String>,
  order_file: Option<String>,
  tree: Option<String>,
  start: Option<String>,
  out: Option<String>,
  params: BuildParams,
}

fn parse_args() -> Result<Args, Report> {
  let options: BTreeMap<String, String> = std::env::args()
    .skip(1)
    .map(|arg| {
      let (key, value) = arg.split_once('=').unwrap_or((arg.as_str(), "true"));
      (key.trim_start_matches('-').to_owned(), value.to_owned())
    })
    .collect();
  let get = |key: &str, default: &str| options.get(key).cloned().unwrap_or_else(|| default.to_owned());
  let Some(aln) = options.get("aln").cloned() else {
    return make_error!(
      "usage: build_tree aln=<fasta> [tree=<reference nwk>] [out=<nwk>] [search=exhaustive|walk|references] [refs=lbi|branches] [max_refs=32] [tau=<radius/ln n>] [pseudo=0.5] [cuts=16] [slack=1] [exact] [batch=64] [order=input|missing|random|scaffold] [anchors=50] [seed=0] [order_file=<names, one per line, placed first in this order>] [spr=5] [spr_slack=2|all] [refresh] [compare] [start_from_reference] [start=<nwk to refine>] [sites] [sites_out=<tsv>] [locality] [rebuild=<radius>] [rebuild_leaves=24] [rebuild_orders=4] [rebuild_nj] [initial=placement|nj] [top_hits=<m, e.g. sqrt n>] [scaffold=<size>] [scaffold_restarts=5] [reprune=<batches>] [reprune_rounds=1]"
    );
  };
  let references = match get("refs", "lbi").as_str() {
    "lbi" => ReferenceSelection::Lbi {
      max_references: get("max_refs", "32").parse()?,
      tau: options.get("tau").map(|tau| tau.parse()).transpose()?,
      pseudo_length: get("pseudo", "0.5").parse()?,
    },
    "branches" => ReferenceSelection::LongBranches {
      cuts: get("cuts", "16").parse()?,
    },
    other => return make_error!("Unknown reference selection '{other}'"),
  };
  let search = match get("search", "exhaustive").as_str() {
    "exhaustive" => SearchMode::Exhaustive,
    "walk" => SearchMode::Walk {
      references,
      slack: get("slack", "1").parse()?,
      exact: options.contains_key("exact"),
    },
    "references" => SearchMode::FromReferences { references },
    other => return make_error!("Unknown search mode '{other}'"),
  };
  let order = match get("order", "random").as_str() {
    "input" => SampleOrder::Input,
    "missing" => SampleOrder::FewestMissingFirst,
    "random" => SampleOrder::Random {
      seed: get("seed", "0").parse()?,
    },
    "scaffold" => SampleOrder::Scaffold {
      anchors: get("anchors", "50").parse()?,
      seed: get("seed", "0").parse()?,
    },
    other => return make_error!("Unknown order '{other}'"),
  };
  Ok(Args {
    aln,
    tree: options.get("tree").cloned(),
    start: options.get("start").cloned(),
    order_file: options.get("order_file").cloned(),
    sites_out: options.get("sites_out").cloned(),
    flags: options
      .iter()
      .filter(|(_, value)| value.as_str() == "true")
      .map(|(key, _)| key.clone())
      .collect(),
    out: options.get("out").cloned(),
    params: BuildParams {
      search,
      batch_size: get("batch", "64").parse()?,
      order,
      spr_rounds: get("spr", "5").parse()?,
      spr_slack: match get("spr_slack", "2").as_str() {
        "all" => None,
        slack => Some(slack.parse()?),
      },
      fitch_refresh: options.contains_key("refresh"),
      compare_exhaustive: options.contains_key("compare"),
      initial: match get("initial", "placement").as_str() {
        "placement" => InitialTree::Placement,
        "nj" => InitialTree::NeighborJoining {
          top_hits: options.get("top_hits").map(|m| m.parse()).transpose()?,
        },
        other => return make_error!("Unknown initial tree '{other}'"),
      },
      scaffold: options
        .get("scaffold")
        .map(|size| -> Result<ScaffoldParams, Report> {
          Ok(ScaffoldParams {
            size: size.parse()?,
            restarts: get("scaffold_restarts", "5").parse()?,
            seed: get("seed", "0").parse()?,
          })
        })
        .transpose()?,
      reprune: options
        .get("reprune")
        .map(|batches| -> Result<RepruneParams, Report> {
          Ok(RepruneParams {
            batches: batches.parse()?,
            rounds: get("reprune_rounds", "1").parse()?,
            seed: 0,
          })
        })
        .transpose()?,
      rebuild: options
        .get("rebuild")
        .map(|radius| -> Result<RebuildParams, Report> {
          Ok(RebuildParams {
            radius: radius.parse()?,
            max_leaves: get("rebuild_leaves", "24").parse()?,
            orders: get("rebuild_orders", "4").parse()?,
            neighbor_joining: options.contains_key("rebuild_nj"),
          })
        })
        .transpose()?,
    },
  })
}

#[expect(clippy::print_stdout, reason = "the example reports its measurements on stdout")]
fn report_locality(label: &str, missed: &[MissedSplit]) {
  let conflicting = missed.iter().filter(|split| split.conflicting).collect_vec();
  let histogram = |values: &[usize], edges: &[usize]| {
    let mut counts = vec![0; edges.len() + 1];
    for &value in values {
      counts[edges.iter().position(|&edge| value <= edge).unwrap_or(edges.len())] += 1;
    }
    edges
      .iter()
      .map(|edge| format!("<={edge}"))
      .chain([format!(">{}", edges[edges.len() - 1])])
      .zip(counts)
      .map(|(label, count)| format!("{label}:{count}"))
      .join(",")
  };
  let moved = histogram(
    &conflicting.iter().map(|split| split.misplaced.len()).collect_vec(),
    &[1, 2, 5, 20],
  );
  let tips = conflicting.iter().flat_map(|split| &split.misplaced).collect_vec();
  let branches = histogram(&tips.iter().map(|tip| tip.branches).collect_vec(), &[1, 2, 3, 5, 10]);
  let substitutions = tips.iter().map(|tip| tip.substitutions).sorted().collect_vec();
  let quantile = |q: usize| {
    substitutions
      .get((substitutions.len() * q).div_euclid(100))
      .copied()
      .unwrap_or(0)
  };
  println!(
    "{label}: missed={} unresolved={} conflicting={} tips_to_move={} misplaced_tips={} branches_from_split={} substitutions_from_split_median={} p90={} max={}",
    missed.len(),
    missed.len() - conflicting.len(),
    conflicting.len(),
    moved,
    tips.len(),
    branches,
    quantile(50),
    quantile(90),
    substitutions.last().copied().unwrap_or(0)
  );
}

fn is_nucleotide(state: AsciiChar) -> bool {
  matches!(u8::from(state), b'A' | b'C' | b'G' | b'T')
}

fn options_flag(args: &Args, key: &str) -> bool {
  args.flags.contains(key)
}

fn site_counts(
  graph: &Graph,
  names: &NodeNames,
  named: &BTreeMap<&str, &Seq>,
  alphabet: &Alphabet,
  length: usize,
) -> Result<Vec<usize>, Report> {
  let mut counts = vec![0; length];
  for sub in fitch_edge_subs(graph, names, named, alphabet)?.values().flatten() {
    counts[sub.pos()] += 1;
  }
  Ok(counts)
}

#[expect(clippy::print_stdout, reason = "the example reports its measurements on stdout")]
fn report_sites(
  built: &TreeWithCounts<'_>,
  reference: &TreeWithCounts<'_>,
  named: &BTreeMap<&str, &Seq>,
  seqs: &[Seq],
  alphabet: &Alphabet,
  out: Option<&str>,
) -> Result<(), Report> {
  let length = seqs.first().map_or(0, Seq::len);
  let built_counts = site_counts(built.graph, built.names, named, alphabet, length)?;
  let reference_counts = site_counts(reference.graph, reference.names, named, alphabet, length)?;
  if let Some(out) = out {
    let table = std::iter::once("position\tbuilt\treference".to_owned())
      .chain((0..length).map(|pos| format!("{}\t{}\t{}", pos + 1, built_counts[pos], reference_counts[pos])))
      .join("\n");
    std::fs::write(out, format!("{table}\n")).wrap_err_with(|| format!("When writing '{out}'"))?;
  }
  let minor = |pos: usize| {
    let counts = seqs
      .iter()
      .map(|seq| seq[pos])
      .filter(|&state| is_nucleotide(state))
      .counts();
    let total: usize = counts.values().sum();
    total - counts.values().max().copied().unwrap_or(0)
  };
  let missing = |pos: usize| seqs.iter().filter(|seq| !is_nucleotide(seq[pos])).count();
  let diffs = (0..length)
    .map(|pos| {
      (
        built_counts[pos].cast_signed() - reference_counts[pos].cast_signed(),
        pos,
      )
    })
    .filter(|&(diff, _)| diff != 0)
    .collect_vec();
  let excess: isize = diffs.iter().filter(|(d, _)| *d > 0).map(|(d, _)| d).sum();
  let deficit: isize = diffs.iter().filter(|(d, _)| *d < 0).map(|(d, _)| -d).sum();
  println!(
    "sites_built_more={} excess={excess} sites_reference_more={} deficit={deficit}",
    diffs.iter().filter(|(d, _)| *d > 0).count(),
    diffs.iter().filter(|(d, _)| *d < 0).count()
  );
  let histogram = diffs.iter().map(|&(d, _)| d.clamp(-3, 3)).counts();
  println!(
    "site_difference_histogram={}",
    histogram.iter().sorted().map(|(d, n)| format!("{d}:{n}")).join(",")
  );
  for &(diff, pos) in diffs
    .iter()
    .sorted_by_key(|&&(d, p)| (std::cmp::Reverse(d), p))
    .take(12)
  {
    println!(
      "site={} built={} reference={} diff={diff} minor_allele_count={} missing={}",
      pos + 1,
      built_counts[pos],
      reference_counts[pos],
      minor(pos),
      missing(pos)
    );
  }
  let excess_minor: Vec<usize> = diffs.iter().filter(|(d, _)| *d > 0).map(|&(_, p)| minor(p)).collect();
  let all_minor: Vec<usize> = (0..length).map(minor).filter(|&m| m > 0).collect();
  let median = |mut values: Vec<usize>| {
    values.sort_unstable();
    values.get(values.len().div_euclid(2)).copied().unwrap_or(0)
  };
  println!(
    "median_minor_allele_count excess_sites={} all_variable_sites={}",
    median(excess_minor),
    median(all_minor)
  );
  Ok(())
}

fn options_order(args: &Args) -> Result<Option<Vec<String>>, Report> {
  args
    .order_file
    .as_ref()
    .map(|path| {
      Ok(
        std::fs::read_to_string(path)
          .wrap_err_with(|| format!("When reading the order file '{path}'"))?
          .lines()
          .map(str::to_owned)
          .collect(),
      )
    })
    .transpose()
}

#[expect(clippy::print_stdout, reason = "the example reports its measurements on stdout")]
fn main() -> Result<(), Report> {
  let args = parse_args()?;
  let alphabet = Alphabet::default();
  let records = fasta_read_file(&args.aln, &alphabet).wrap_err("When reading the alignment")?;
  let reference = args.tree.as_ref().map(nwk_read_file).transpose()?;
  let keep: Option<BTreeSet<String>> = reference.as_ref().map(|parse| {
    parse
      .graph
      .get_leaves()
      .map(|leaf| parse.names[leaf.key()].to_owned())
      .collect()
  });
  let rank: BTreeMap<String, usize> = match options_order(&args)? {
    Some(order) => order.into_iter().enumerate().map(|(rank, name)| (name, rank)).collect(),
    None => BTreeMap::new(),
  };
  let (names, seqs): (Vec<String>, Vec<Seq>) = records
    .into_iter()
    .filter(|record| keep.as_ref().is_none_or(|keep| keep.contains(&record.seq_name)))
    .enumerate()
    .sorted_by_key(|(index, record)| (rank.get(&record.seq_name).copied().unwrap_or(usize::MAX), *index))
    .map(|(_, record)| (record.seq_name, record.seq))
    .unzip();
  println!("samples={} length={}", seqs.len(), seqs.first().map_or(0, Seq::len));

  let start = Instant::now();
  let start_tree = match (&args.tree, options_flag(&args, "start_from_reference")) {
    (_, false) => args.start.as_ref(),
    (tree, true) => tree.as_ref(),
  };
  let built = match start_tree {
    Some(path) => {
      let parse = nwk_read_file(path)?;
      let index: BTreeMap<&str, usize> = names.iter().map(String::as_str).zip(0..).collect();
      let samples = parse
        .graph
        .get_leaves()
        .map(|leaf| (leaf.key(), index[&parse.names[leaf.key()]]))
        .collect();
      let mat = mat_from_tree(parse.graph, samples, &seqs, &alphabet)?;
      refine_mat(mat, &seqs, &alphabet, &args.params, BuildStats::default())?
    },
    None => build_mat(&seqs, &alphabet, &args.params)?,
  };
  let elapsed = start.elapsed().as_secs_f64();
  let stats = &built.stats;
  println!(
    "build_seconds={elapsed:.2} placement_parsimony={} final_parsimony={} replaced_in_batch={} visited={} references={} collapsed={}",
    stats.placement_parsimony,
    stats.final_parsimony,
    stats.replaced_in_batch,
    stats.visited,
    stats.references,
    stats.collapsed_edges
  );
  for (round, spr) in stats.spr.iter().enumerate() {
    println!(
      "spr_round={round} tried={} moves={} gain={} visited={} parsimony_before={} parsimony_after={}",
      spr.tried, spr.moves, spr.gain, spr.visited, spr.parsimony_before, spr.parsimony_after
    );
  }
  for (round, rebuild) in stats.rebuild.iter().enumerate() {
    println!(
      "rebuild_round={round} tried={} improved={} gain={} too_large={}",
      rebuild.tried, rebuild.improved, rebuild.gain, rebuild.too_large
    );
  }
  if let Some(reprune) = &stats.reprune {
    println!(
      "reprune batches={} accepted={} gain={}",
      reprune.batches, reprune.accepted, reprune.gain
    );
  }
  if let Some(comparison) = &stats.comparison {
    println!(
      "vs_exhaustive placed={} worse={} excess={} max_excess={} visited_search={} visited_exhaustive={}",
      comparison.placed,
      comparison.worse,
      comparison.excess,
      comparison.max_excess,
      comparison.visited_search,
      comparison.visited_exhaustive
    );
  }

  let tree = mat_to_tree(built.mat, &names)?;
  if let Some(out) = &args.out {
    nwk_write_file(
      out,
      &GraphSchedule::new(&tree.graph)?.into_tree()?,
      &tree.names,
      &tree.branch_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    )?;
  }

  let named: BTreeMap<&str, &Seq> = names.iter().map(String::as_str).zip(&seqs).collect();
  let built_counts = fitch_edge_counts(&tree.graph, &tree.names, &named, &alphabet)?;
  let built_fitch: usize = built_counts.values().sum();
  let built_informative = built_counts.values().filter(|&&count| count > 0).count();
  println!("built_fitch_parsimony={built_fitch} built_edges_with_mutations={built_informative}");

  if let Some(reference) = &reference {
    let reference_counts = fitch_edge_counts(&reference.graph, &reference.names, &named, &alphabet)?;
    let reference_fitch: usize = reference_counts.values().sum();
    let splits = compare_splits(
      &TreeWithCounts {
        graph: &tree.graph,
        names: &tree.names,
        edge_counts: &built_counts,
      },
      &TreeWithCounts {
        graph: &reference.graph,
        names: &reference.names,
        edge_counts: &reference_counts,
      },
    )?;
    println!(
      "reference_fitch_parsimony={reference_fitch} splits_shared={} only_built={} only_reference={}",
      splits.shared, splits.only_first, splits.only_second
    );
    if options_flag(&args, "locality") {
      let built_tree = TreeWithCounts {
        graph: &tree.graph,
        names: &tree.names,
        edge_counts: &built_counts,
      };
      let reference_tree = TreeWithCounts {
        graph: &reference.graph,
        names: &reference.names,
        edge_counts: &reference_counts,
      };
      report_locality(
        "reference_splits_missing_in_built",
        &missed_splits(&built_tree, &reference_tree)?,
      );
      report_locality(
        "built_splits_missing_in_reference",
        &missed_splits(&reference_tree, &built_tree)?,
      );
    }
    if options_flag(&args, "sites") {
      report_sites(
        &TreeWithCounts {
          graph: &tree.graph,
          names: &tree.names,
          edge_counts: &built_counts,
        },
        &TreeWithCounts {
          graph: &reference.graph,
          names: &reference.names,
          edge_counts: &reference_counts,
        },
        &named,
        &seqs,
        &alphabet,
        args.sites_out.as_deref(),
      )?;
    }
  }
  Ok(())
}
