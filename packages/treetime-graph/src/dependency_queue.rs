use crossbeam_channel::{Receiver, Sender, select, unbounded};
use eyre::Report;
use parking_lot::Mutex;
use std::collections::VecDeque;
use std::sync::atomic::{AtomicUsize, Ordering};
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
  let error = Mutex::new(None);
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
  error: &'a Mutex<Option<Report>>,
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
          if self.error.lock().is_some() {
            continue;
          }
          if let Err(report) = (self.visit)(index) {
            let mut error = self.error.lock();
            if error.is_none() {
              *error = Some(report);
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
