use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;
use treetime::alphabet::alphabet::Alphabet;
use treetime::build::evaluate::{TreeWithCounts, compare_splits, fitch_edge_counts};
use treetime::build::pipeline::{BuildParams, SampleOrder, SearchMode, build_mat, mat_to_tree};
use treetime::build::references::ReferenceSelection;
use treetime_graph::schedule::GraphSchedule;
use treetime_io::fasta::fasta_read_file;
use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_primitives::Seq;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;

struct Args {
  aln: String,
  order_file: Option<String>,
  tree: Option<String>,
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
      "usage: build_tree aln=<fasta> [tree=<reference nwk>] [out=<nwk>] [search=exhaustive|walk|regions] [refs=lbi|branches] [max_refs=32] [tau=<radius/ln n>] [pseudo=0.5] [cuts=16] [slack=1] [exact] [batch=64] [order=input|missing] [order_file=<names, one per line, placed first in this order>] [spr=5] [spr_slack=2|all] [refresh] [compare]"
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
    "regions" => SearchMode::Regions { references },
    other => return make_error!("Unknown search mode '{other}'"),
  };
  let order = match get("order", "missing").as_str() {
    "input" => SampleOrder::Input,
    "missing" => SampleOrder::FewestMissingFirst,
    other => return make_error!("Unknown order '{other}'"),
  };
  Ok(Args {
    aln,
    tree: options.get("tree").cloned(),
    order_file: options.get("order_file").cloned(),
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
    },
  })
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
  let built = build_mat(&seqs, &alphabet, &args.params)?;
  let elapsed = start.elapsed().as_secs_f64();
  let stats = &built.stats;
  println!(
    "build_seconds={elapsed:.2} placement_parsimony={} final_parsimony={} replaced_in_batch={} visited={} references={} certified={} collapsed={}",
    stats.placement_parsimony,
    stats.final_parsimony,
    stats.replaced_in_batch,
    stats.visited,
    stats.references,
    stats.certified,
    stats.collapsed_edges
  );
  println!(
    "mean_global_diffs={:.1} mean_local_diffs={:.1}",
    stats.global_diffs.to_f64_exact() / seqs.len().max(1).to_f64_exact(),
    stats.local_diffs.to_f64_exact() / seqs.len().max(1).to_f64_exact()
  );
  for (round, spr) in stats.spr.iter().enumerate() {
    println!(
      "spr_round={round} tried={} moves={} gain={} visited={} parsimony_before={} parsimony_after={}",
      spr.tried, spr.moves, spr.gain, spr.visited, spr.parsimony_before, spr.parsimony_after
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
  }
  Ok(())
}
