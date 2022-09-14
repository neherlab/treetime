pub mod alphabet;
pub mod ancestral;
pub mod clock;
pub mod constants;
pub mod graph;
pub mod gtr;
pub mod homoplasy;
pub mod io;
pub mod mugration;
pub mod nuc_models;
pub mod seq_utils;
pub mod timetree;
pub mod utils;

#[cfg(test)]
mod tests {
  use crate::utils::global_init::global_init;
  use ctor::ctor;

  #[ctor]
  fn init() {
    global_init();
  }
}
