#[cfg(test)]
mod __tests__;
pub mod auspice_types;
pub mod csv;
pub mod dates_csv;
pub mod discrete_states_csv;
pub mod fasta;
pub mod gff;
pub mod graphviz;
pub mod name_list;
pub mod nex;
pub mod nwk;
pub mod tree;
pub mod usher_mat;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
