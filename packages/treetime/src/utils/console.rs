pub fn is_tty() -> bool {
  #[cfg(not(target_arch = "wasm32"))]
  {
    atty::is(atty::Stream::Stderr)
  }

  #[cfg(target_arch = "wasm32")]
  {
    false
  }
}
