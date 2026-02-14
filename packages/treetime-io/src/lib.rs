mod __tests__;
pub mod auspice;
pub mod auspice_types;
pub mod concat;
pub mod csv;
pub mod dates_csv;
pub mod discrete_states_csv;
pub mod fasta;
pub mod graphviz;
pub mod json;
pub mod nex;
pub mod nwk;
pub mod parse_delimited;
pub mod phyloxml;
pub mod usher_mat;
pub mod yaml;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
