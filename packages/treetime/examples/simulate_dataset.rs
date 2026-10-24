use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime_graph::schedule::GraphSchedule;
use treetime_io::fasta::FastaWriter;
use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_write_file};
use treetime_simulation::coalescent::{SampleBatch, simulate_coalescent};
use treetime_simulation::evolve::{evolve_sequences_with_rates, gamma_site_rates};
use treetime_simulation::observe::{TipObservation, observe_tips};
use treetime_simulation::substitution_model::SubstitutionModel;
use treetime_simulation::tc_schedule::TcSchedule;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;
use treetime_utils::sync::random::get_random_number_generator;

#[expect(clippy::print_stdout, reason = "the example reports what it wrote on stdout")]
fn main() -> Result<(), Report> {
  let options: BTreeMap<String, String> = std::env::args()
    .skip(1)
    .map(|arg| {
      let (key, value) = arg.split_once('=').unwrap_or((arg.as_str(), "true"));
      (key.to_owned(), value.to_owned())
    })
    .collect();
  let get = |key: &str, default: &str| options.get(key).cloned().unwrap_or_else(|| default.to_owned());
  let Some(out) = options.get("out") else {
    return make_error!(
      "usage: simulate_dataset out=<dir> [tips_per_time=10] [times=40] [interval=0.1] [tc=10] [mu=0.003] [length=10000] [alpha=<gamma shape, none for equal rates>] [messy] [seed=0]"
    );
  };
  let seed: u64 = get("seed", "0").parse()?;
  let tips_per_time: usize = get("tips_per_time", "10").parse()?;
  let times: usize = get("times", "40").parse()?;
  let interval: f64 = get("interval", "0.1").parse()?;
  let mu: f64 = get("mu", "0.003").parse()?;
  let length: usize = get("length", "10000").parse()?;

  let mut rng = get_random_number_generator(seed);
  let sampling = (0..times)
    .map(|time| SampleBatch {
      tau: time.to_f64_exact() * interval,
      count: tips_per_time,
    })
    .collect::<Vec<_>>();
  let tree = simulate_coalescent(&TcSchedule::constant(get("tc", "10").parse()?)?, &sampling, &mut rng)?;
  let model = SubstitutionModel::random(mu, &mut rng)?;
  let rates = match options.get("alpha") {
    Some(alpha) => gamma_site_rates(alpha.parse()?, length, &mut rng)?,
    None => vec![1.0; length],
  };
  let durations = tree.edge_durations();
  let simulated = evolve_sequences_with_rates(&tree.graph, &durations, &model, &rates, &mut rng)?;
  let tips = if options.contains_key("messy") {
    let observation = TipObservation {
      ambiguity_rate: 0.002,
      unknown_run_rate: 0.05,
      unknown_run_lengths: 5..=100,
    };
    observe_tips(&tree.graph, &simulated.sequences, &observation, &mut rng)?
  } else {
    tree
      .graph
      .get_leaves()
      .map(|leaf| (leaf.key(), simulated.sequences[&leaf.key()].clone()))
      .collect()
  };

  std::fs::create_dir_all(out).wrap_err_with(|| format!("When creating '{out}'"))?;
  let mut writer = FastaWriter::create(format!("{out}/tips.fasta"))?;
  for (&leaf, seq) in &tips {
    writer.write(&tree.names[leaf], None, seq)?;
  }
  writer.finish()?;
  let branch_lengths = durations
    .iter()
    .map(|(&edge, &duration)| (edge, Some(duration)))
    .collect();
  nwk_write_file(
    format!("{out}/tree.nwk"),
    &GraphSchedule::new(&tree.graph)?.into_tree()?,
    &tree.names,
    &branch_lengths,
    &NwkWriteOptions::default(),
    &NwkNodeComments::new(),
  )?;
  let substitutions: usize = simulated.substitutions.values().sum();
  println!(
    "tips={} length={length} substitutions_on_tree={substitutions}",
    tips.len()
  );
  Ok(())
}
