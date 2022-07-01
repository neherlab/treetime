use ctor::ctor;
use eyre::Report;
use ndarray::Array2;
use ndarray_rand::RandomExt;
use treetime::nuc_models::jc69::jc69;
use treetime::seq_utils::normalize_profile::normalize_profile;
use treetime::seq_utils::seq2prof::{prof2seq, seq2prof};
use treetime::utils::global_init::global_init;

#[cfg(all(target_family = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let gtr = jc69()?;

  let dummy_prof = Array2::<f32>::random((10000, 5));

  // used a lot (300us)
  let (norm_prof, _) = normalize_profile(&dummy_prof, false);

  // used less but still a lot (50us)
  gtr.evolve(&norm_prof, 0.1, false);

  // used less but still a lot (50us)
  gtr.propagate_profile(norm_prof, 0.1, false);

  // // used only in final, sample_from_prof=False speeds it up (600us or 300us)
  // let  (seq, p, seq_ii) = prof2seq(norm_prof, &gtr, true, false);
  //
  // // used only initially (slow, 5ms)
  // let tmp_prof = seq2prof(seq, gtr.profile_map);

  Ok(())
}
