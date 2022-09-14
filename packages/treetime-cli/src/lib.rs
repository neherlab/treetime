pub mod cli;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime::utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
