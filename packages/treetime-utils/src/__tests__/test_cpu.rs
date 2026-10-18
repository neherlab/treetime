#[cfg(test)]
mod tests {
  use crate::init::cpu::cpu_has_build_features;

  #[test]
  fn test_cpu_build_machine_has_the_features_the_build_targets() {
    assert!(cpu_has_build_features());
  }
}
