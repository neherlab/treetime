use std::io::{IsTerminal, stderr};

pub fn is_tty() -> bool {
  #[cfg(not(target_arch = "wasm32"))]
  {
    stderr().is_terminal()
  }

  #[cfg(target_arch = "wasm32")]
  {
    false
  }
}
