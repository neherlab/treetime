use crossbeam_channel::{Receiver, Sender, select, unbounded};
use eyre::Report;
use std::collections::VecDeque;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use treetime_utils::make_internal_report;

pub fn run_dependency_queue(
  prerequisites: &[usize],
  successors: &[Vec<usize>],
  visit: impl Fn(usize) -> Result<(), Report> + Sync + Send,
) -> Result<(), Report> {
  if prerequisites.is_empty() {
    return Ok(());
  }
  debug_assert_eq!(prerequisites.len(), successors.len());

  let remaining = prerequisites.iter().copied().map(AtomicUsize::new).collect::<Vec<_>>();
  let completed = AtomicUsize::new(0);
  // First-error reporting without a lock: one worker is elected atomically (see `run_worker`) and is
  // the sole writer of this single-publication cell, so competing errors never contend on it.
  let error = OnceLock::new();
  let failed = AtomicBool::new(false);
  let workers = rayon::current_num_threads();
  let (work_sender, work_receiver) = unbounded();
  let (stop_sender, stop_receiver) = unbounded();
  prerequisites
    .iter()
    .enumerate()
    .filter(|(_, count)| **count == 0)
    .try_for_each(|(index, _)| work_sender.send(index).map_err(Report::new))?;

  let workers = DependencyWorkers {
    node_count: prerequisites.len(),
    successors,
    remaining: &remaining,
    completed: &completed,
    error: &error,
    failed: &failed,
    work_sender: &work_sender,
    work_receiver: &work_receiver,
    stop_sender: &stop_sender,
    stop_receiver: &stop_receiver,
    worker_count: workers,
    visit: &visit,
  };
  workers.run();

  error.into_inner().map_or(Ok(()), Err)
}

pub fn validate_dependency_graph(prerequisites: &[usize], successors: &[Vec<usize>]) -> Result<(), Report> {
  if prerequisites.len() != successors.len() {
    return Err(make_internal_report!(
      "Dependency graph has {} nodes but {} successor lists",
      prerequisites.len(),
      successors.len()
    ));
  }
  let mut remaining = prerequisites.to_vec();
  let mut ready = remaining
    .iter()
    .enumerate()
    .filter_map(|(index, count)| (*count == 0).then_some(index))
    .collect::<VecDeque<_>>();
  let mut visited = 0;
  while let Some(index) = ready.pop_front() {
    visited += 1;
    for successor in &successors[index] {
      let count = remaining
        .get_mut(*successor)
        .ok_or_else(|| make_internal_report!("Dependency successor index {successor} is outside the node set"))?;
      if *count == 0 {
        return Err(make_internal_report!(
          "Dependency graph contains duplicate readiness for node {successor}"
        ));
      }
      *count -= 1;
      if *count == 0 {
        ready.push_back(*successor);
      }
    }
  }
  if visited != prerequisites.len() {
    return Err(make_internal_report!(
      "Dependency graph is cyclic: visited {visited} of {} nodes",
      prerequisites.len()
    ));
  }
  Ok(())
}

struct DependencyWorkers<'a, F> {
  node_count: usize,
  successors: &'a [Vec<usize>],
  remaining: &'a [AtomicUsize],
  completed: &'a AtomicUsize,
  error: &'a OnceLock<Report>,
  failed: &'a AtomicBool,
  work_sender: &'a Sender<usize>,
  work_receiver: &'a Receiver<usize>,
  stop_sender: &'a Sender<()>,
  stop_receiver: &'a Receiver<()>,
  worker_count: usize,
  visit: &'a F,
}

impl<F> DependencyWorkers<'_, F>
where
  F: Fn(usize) -> Result<(), Report> + Sync,
{
  fn run(&self) {
    rayon::scope(|scope| {
      for _ in 0..self.worker_count {
        scope.spawn(|_| self.run_worker());
      }
    });
  }

  fn run_worker(&self) {
    loop {
      select! {
        recv(self.work_receiver) -> index => {
          let Ok(index) = index else { return };
          if self.failed.load(Ordering::Acquire) {
            continue;
          }
          if let Err(report) = (self.visit)(index) {
            // Atomically elect a single error publisher: exactly one worker sees `false` here and
            // becomes the sole writer of the single-publication `error` cell, then cancels the rest.
            // Later failures observe `true` and drop their report, so no two workers ever race to
            // initialize the cell.
            if !self.failed.swap(true, Ordering::AcqRel) {
              assert!(
                self.error.set(report).is_ok(),
                "The elected error publisher must publish the first error exactly once"
              );
              self.stop();
            }
            continue;
          }

          for successor in &self.successors[index] {
            let previous = self.remaining[*successor].fetch_sub(1, Ordering::AcqRel);
            debug_assert!(previous > 0);
            if previous == 1 {
              self.work_sender.send(*successor).expect("Dependency work channel must remain connected");
            }
          }
          if self.completed.fetch_add(1, Ordering::AcqRel) + 1 == self.node_count {
            self.stop();
          }
        }
        recv(self.stop_receiver) -> _ => return,
      }
    }
  }

  fn stop(&self) {
    for _ in 0..self.worker_count {
      self
        .stop_sender
        .send(())
        .expect("Dependency stop channel must remain connected");
    }
  }
}

#[cfg(test)]
mod tests {
  use super::{run_dependency_queue, validate_dependency_graph};
  use eyre::Report;
  use std::sync::atomic::{AtomicUsize, Ordering};
  use treetime_utils::{assert_error, make_report};

  #[test]
  fn test_dependency_queue_validate_accepts_acyclic_graph() -> Result<(), Report> {
    // A fork: node 0 depends on 1 and 2; both are ready immediately.
    let prerequisites = [2, 0, 0];
    let successors = [vec![], vec![0], vec![0]];
    validate_dependency_graph(&prerequisites, &successors)?;
    Ok(())
  }

  #[test]
  fn test_dependency_queue_validate_rejects_cycle() {
    // A two-node cycle: each node waits on the other, so neither is ever ready.
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
    // Node 1 depends on node 0. Node 0 fails, so node 1 must never be scheduled.
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
