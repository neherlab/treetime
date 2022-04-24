pub mod ancestral;
pub mod cli;
pub mod clock;
pub mod homoplasy;
pub mod io;
pub mod mugration;
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
