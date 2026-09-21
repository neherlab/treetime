#[cfg(test)]
mod tests {
  use crate::dependency_queue::{run_dependency_queue, validate_dependency_graph};
  use eyre::Report;
  use std::sync::atomic::{AtomicUsize, Ordering};
  use treetime_utils::{assert_error, make_report};

  #[test]
  fn test_dependency_queue_validate_accepts_acyclic_graph() -> Result<(), Report> {
    let prerequisites = [2, 0, 0];
    let successors = [vec![], vec![0], vec![0]];
    validate_dependency_graph(&prerequisites, &successors)?;
    Ok(())
  }

  #[test]
  fn test_dependency_queue_validate_rejects_cycle() {
    let prerequisites = [1, 1];
    let successors = [vec![1], vec![0]];
    let result = validate_dependency_graph(&prerequisites, &successors);
    assert_error!(
      result,
      "Dependency graph is cyclic: visited 0 of 2 nodes. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_dependency_queue_failing_visit_stops_and_returns_error() -> Result<(), Report> {
    let prerequisites = [0, 1];
    let successors = [vec![1], vec![]];
    let visited = AtomicUsize::new(0);

    let result = run_dependency_queue(&prerequisites, &successors, |index| {
      visited.fetch_add(1, Ordering::AcqRel);
      if index == 0 {
        Err(make_report!("injected visit failure"))
      } else {
        Ok(())
      }
    });

    assert_error!(result, "injected visit failure");
    assert_eq!(1, visited.load(Ordering::Acquire));
    Ok(())
  }
}
