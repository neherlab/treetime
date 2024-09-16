pub mod alphabet;
pub mod cli;
pub mod commands;
pub mod constants;
pub mod graph;
pub mod gtr;
pub mod hacks;
pub mod io;
pub mod representation;
pub mod seq;
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
