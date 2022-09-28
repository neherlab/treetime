pub mod alphabet;
pub mod cli;
pub mod commands;
pub mod constants;
pub mod graph;
pub mod gtr;
pub mod io;
pub mod seq_utils;
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
