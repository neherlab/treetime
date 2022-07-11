use ctor::ctor;
use eyre::Report;
use log::{info, LevelFilter};
use ndarray::Array2;
use ndarray_rand::RandomExt;
use rand::SeedableRng;
use rand_isaac::Isaac64Rng;
use treetime::nuc_models::jc69::jc69;
use treetime::seq_utils::normalize_profile::normalize_profile;
use treetime::seq_utils::seq2prof::{prof2seq, seq2prof, Prof2SeqParams};
use treetime::utils::global_init::{global_init, setup_logger};
use treetime::utils::ndarray::random;

#[cfg(all(target_family = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
  setup_logger(LevelFilter::Info);
}

const RANDOM_SEED: u64 = 42;

fn main() -> Result<(), Report> {
  let gtr = jc69()?;
  info!("gtr:\n{gtr:#?}\n");

  let mut rng = Isaac64Rng::seed_from_u64(RANDOM_SEED);
  let dummy_prof: Array2<f64> = random((100, 5), &mut rng);
  info!("dummy_prof:\n{dummy_prof}\n");

  // used a lot (300us)
  let (norm_prof, _) = normalize_profile(&dummy_prof, false)?;
  info!("norm_prof:\n{norm_prof}\n");

  // used less but still a lot (50us)
  let evolved = gtr.evolve(&norm_prof, 0.1, false);
  info!("evolved:\n{evolved}\n");

  // used less but still a lot (50us)
  let propagated = gtr.propagate_profile(&norm_prof, 0.1, false);
  info!("propagated:\n{propagated}\n");

  // used only in final, sample_from_prof=False speeds it up (600us or 300us)
  let prof2seq_result = prof2seq(
    &norm_prof,
    &gtr,
    &mut rng,
    &Prof2SeqParams {
      should_sample_from_profile: true,
      should_normalize_profile: false,
    },
  )?;

  info!("prof2seq_result:\n{prof2seq_result:#?}\n");

  let prof = seq2prof(&prof2seq_result.seq, &gtr.profile_map)?;
  info!("prof:\n{prof}\n");

  Ok(())
}
