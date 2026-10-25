use crate::key_map::{EdgeMap, NodeMap};
use crate::node::GraphNodeKey;
use crate::schedule::GraphSchedule;
use eyre::Report;
use itertools::{Itertools, izip};
use parking_lot::Mutex;
use rayon::{Scope, current_num_threads, scope};
use smallvec::{SmallVec, smallvec};
use std::iter::{repeat_n, repeat_with};
use std::marker::PhantomData;
use std::mem;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};
use treetime_utils::fmt::string::counted;
use treetime_utils::make_internal_report;

impl GraphSchedule<'_> {
  pub fn scan_down<N, M, Out>(
    &self,
    inputs: impl IntoIterator<Item = (GraphNodeKey, N)>,
    step: impl Fn(GraphNodeKey, N, PerEdge<M>) -> Result<(Out, M), Report> + Sync,
  ) -> Result<ScanOutputs<'_, Out>, Report>
  where
    N: Send,
    M: Clone + Send,
    Out: Send + Sync,
  {
    let scan = Scan::new(self, inputs, step)?;
    let inboxes = self.inboxes(|key| self.parents(key).len());
    let (scan_ref, inboxes) = (&scan, inboxes.as_slice());
    scope(|scope| {
      if let Some((&first, rest)) = self.roots().split_first() {
        for &root in rest {
          scope.spawn(move |scope| scan_ref.down(scope, inboxes, root, PerEdge::new()));
        }
        scan_ref.down(scope, inboxes, first, PerEdge::new());
      }
    });
    scan.finish()
  }

  pub fn scan_up<N, M, Out>(
    &self,
    inputs: impl IntoIterator<Item = (GraphNodeKey, N)>,
    step: impl Fn(GraphNodeKey, N, PerEdge<M>) -> Result<(Out, M), Report> + Sync,
  ) -> Result<ScanOutputs<'_, Out>, Report>
  where
    N: Send,
    M: Clone + Send,
    Out: Send + Sync,
  {
    let scan = Scan::new(self, inputs, step)?;
    let inboxes = self.inboxes(|key| self.children(key).len());
    let leaves = self.leaves();
    let next = AtomicUsize::new(0);
    scope(|scope| {
      for _ in 0..current_num_threads().min(leaves.len()) {
        scope.spawn(|_| {
          let mut runnable = Vec::new();
          while let Some(&leaf) = leaves.get(next.fetch_add(1, Ordering::Relaxed)) {
            runnable.push((leaf, PerEdge::new()));
            scan.up(&inboxes, &mut runnable);
          }
        });
      }
    });
    scan.finish()
  }

  pub fn join_parent_edges<N, E>(
    &self,
    nodes: NodeMap<N>,
    mut edges: EdgeMap<E>,
  ) -> Result<impl Iterator<Item = (GraphNodeKey, (N, PerEdge<E>))>, Report> {
    for &key in self.preorder() {
      if !nodes.contains_key(key) {
        return Err(make_internal_report!("Scan input is missing node {key}"));
      }
      if let Some((_, edge_key)) = self
        .parents(key)
        .iter()
        .find(|(_, edge_key)| !edges.contains_key(edge_key))
      {
        return Err(make_internal_report!(
          "Scan input is missing edge {edge_key} above node {key}"
        ));
      }
    }
    if let Some(key) = nodes.keys().find(|&key| !self.contains(key)) {
      return Err(make_internal_report!(
        "Scan input names node {key}, which is not in the graph"
      ));
    }
    if edges.len() > self.edge_count()
      && let Some(edge_key) = edges.keys().find(|&edge_key| !self.contains_edge(edge_key))
    {
      return Err(make_internal_report!(
        "Scan input names edge {edge_key}, which is not in the graph"
      ));
    }
    Ok(nodes.into_iter().map(move |(key, node)| {
      let values = self
        .parents(key)
        .iter()
        .filter_map(|(_, edge_key)| edges.remove(edge_key))
        .collect();
      (key, (node, values))
    }))
  }

  fn inboxes<M>(&self, slots: impl Fn(GraphNodeKey) -> usize) -> Vec<Mutex<Inbox<M>>> {
    (0..self.key_capacity())
      .map(|index| Mutex::new(Inbox::new(slots(GraphNodeKey(index)))))
      .collect_vec()
  }
}

struct Scan<'s, 'g, N, M, Out, F> {
  schedule: &'s GraphSchedule<'g>,
  inputs: Vec<Mutex<Option<N>>>,
  outputs: Vec<OnceLock<Out>>,
  error: OnceLock<Report>,
  step: F,
  messages: PhantomData<fn(PerEdge<M>) -> M>,
}

impl<'s, 'g, N, M, Out, F> Scan<'s, 'g, N, M, Out, F>
where
  F: Fn(GraphNodeKey, N, PerEdge<M>) -> Result<(Out, M), Report>,
{
  fn new(
    schedule: &'s GraphSchedule<'g>,
    inputs: impl IntoIterator<Item = (GraphNodeKey, N)>,
    step: F,
  ) -> Result<Self, Report> {
    let capacity = schedule.key_capacity();
    let mut slots = repeat_with(|| Mutex::new(None)).take(capacity).collect_vec();
    for (key, input) in inputs {
      if !schedule.contains(key) {
        return Err(make_internal_report!(
          "Scan input names node {key}, which is not in the graph"
        ));
      }
      if slots[key.as_usize()].get_mut().replace(input).is_some() {
        return Err(make_internal_report!("Scan input names node {key} twice"));
      }
    }
    if let Some(key) = schedule
      .preorder()
      .iter()
      .find(|key| slots[key.as_usize()].get_mut().is_none())
    {
      return Err(make_internal_report!("Scan input is missing node {key}"));
    }
    Ok(Self {
      schedule,
      inputs: slots,
      outputs: repeat_with(OnceLock::new).take(capacity).collect_vec(),
      error: OnceLock::new(),
      step,
      messages: PhantomData,
    })
  }

  fn visit(&self, key: GraphNodeKey, received: PerEdge<M>) -> Option<M> {
    if self.error.get().is_some() {
      return None;
    }
    let result = self.inputs[key.as_usize()]
      .lock()
      .take()
      .ok_or_else(|| make_internal_report!("Scan visited node {key} twice"))
      .and_then(|input| (self.step)(key, input, received));
    match result {
      Ok((output, message)) => {
        if self.outputs[key.as_usize()].set(output).is_err() {
          self.fail(make_internal_report!("Scan stored the output of node {key} twice"));
          return None;
        }
        Some(message)
      },
      Err(report) => {
        self.fail(report);
        None
      },
    }
  }

  fn deliver(
    &self,
    inboxes: &[Mutex<Inbox<M>>],
    node: GraphNodeKey,
    slots: usize,
    slot: usize,
    message: M,
  ) -> Option<PerEdge<M>> {
    if (slots, slot) == (1, 0) {
      return Some(smallvec![message]);
    }
    match inboxes[node.as_usize()].lock().deliver(node, slot, message) {
      Ok(ready) => ready,
      Err(report) => {
        self.fail(report);
        None
      },
    }
  }

  fn fail(&self, report: Report) {
    self.error.get_or_init(|| report);
  }

  fn finish(self) -> Result<ScanOutputs<'s, Out>, Report> {
    if let Some(report) = self.error.into_inner() {
      return Err(report);
    }
    let outputs = self.outputs.into_iter().map(OnceLock::into_inner).collect_vec();
    if let Some(key) = self
      .schedule
      .preorder()
      .iter()
      .find(|key| outputs[key.as_usize()].is_none())
    {
      return Err(make_internal_report!("Scan finished without visiting node {key}"));
    }
    Ok(ScanOutputs {
      schedule: self.schedule,
      outputs,
    })
  }
}

impl<N, M, Out, F> Scan<'_, '_, N, M, Out, F>
where
  N: Send,
  M: Clone + Send,
  Out: Send + Sync,
  F: Fn(GraphNodeKey, N, PerEdge<M>) -> Result<(Out, M), Report> + Sync,
{
  fn down<'s>(&'s self, scope: &Scope<'s>, inboxes: &'s [Mutex<Inbox<M>>], key: GraphNodeKey, received: PerEdge<M>) {
    let mut current = Some((key, received));
    while let Some((key, received)) = current.take() {
      let Some(message) = self.visit(key, received) else {
        return;
      };
      let children = self.schedule.children(key);
      let ready = izip!(children, repeat_n(message, children.len())).filter_map(|(&(child, edge), message)| {
        let slot = self.schedule.parent_slot(edge);
        self
          .deliver(inboxes, child, self.schedule.parents(child).len(), slot, message)
          .map(|ready| (child, ready))
      });
      for (child, received) in ready {
        if current.is_none() {
          current = Some((child, received));
        } else {
          scope.spawn(move |scope| self.down(scope, inboxes, child, received));
        }
      }
    }
  }

  fn up(&self, inboxes: &[Mutex<Inbox<M>>], runnable: &mut Vec<(GraphNodeKey, PerEdge<M>)>) {
    while let Some((key, received)) = runnable.pop() {
      let Some(message) = self.visit(key, received) else {
        runnable.clear();
        return;
      };
      let parents = self.schedule.parents(key);
      let ready = izip!(parents, repeat_n(message, parents.len())).filter_map(|(&(parent, edge), message)| {
        let slot = self.schedule.child_slot(edge);
        self
          .deliver(inboxes, parent, self.schedule.children(parent).len(), slot, message)
          .map(|ready| (parent, ready))
      });
      runnable.extend(ready);
    }
  }
}

struct Inbox<M> {
  missing: usize,
  messages: Vec<Option<M>>,
}

impl<M> Inbox<M> {
  fn new(slots: usize) -> Self {
    Self {
      missing: slots,
      messages: Vec::new(),
    }
  }

  fn deliver(&mut self, node: GraphNodeKey, slot: usize, message: M) -> Result<Option<PerEdge<M>>, Report> {
    if self.messages.is_empty() {
      self.messages.resize_with(self.missing, || None);
    }
    let entry = self
      .messages
      .get_mut(slot)
      .ok_or_else(|| make_internal_report!("Scan found no inbox slot {slot} at node {node}"))?;
    *entry = Some(message);
    self.missing -= 1;
    if self.missing > 0 {
      return Ok(None);
    }
    mem::take(&mut self.messages)
      .into_iter()
      .enumerate()
      .map(|(slot, message)| {
        message.ok_or_else(|| {
          make_internal_report!("Scan completed the inbox of node {node} without a message in slot {slot}")
        })
      })
      .try_collect()
      .map(Some)
  }
}

pub struct ScanOutputs<'s, Out> {
  schedule: &'s GraphSchedule<'s>,
  outputs: Vec<Option<Out>>,
}

impl<Out> ScanOutputs<'_, Out> {
  pub fn into_map(self) -> NodeMap<Out> {
    NodeMap::from_slots(self.outputs)
  }
}

impl<N, E> ScanOutputs<'_, (N, PerEdge<E>)> {
  pub fn split_parent_edges(self) -> Result<NodeEdgeMaps<N, E>, Report> {
    split_parent_edges(self.schedule, self.outputs.into_iter().map(Ok))
  }
}

impl<'s, T> ScanOutputs<'s, Arc<T>> {
  pub fn unwrap_shared(self) -> SharedOutputs<'s, T> {
    SharedOutputs {
      schedule: self.schedule,
      outputs: self.outputs,
    }
  }
}

pub struct SharedOutputs<'s, T> {
  schedule: &'s GraphSchedule<'s>,
  outputs: Vec<Option<Arc<T>>>,
}

impl<T> SharedOutputs<'_, T> {
  pub fn into_map(self) -> Result<NodeMap<T>, Report> {
    unwrapped_slots(self.outputs).try_collect().map(NodeMap::from_slots)
  }
}

impl<N, E> SharedOutputs<'_, (N, PerEdge<E>)> {
  pub fn split_parent_edges(self) -> Result<NodeEdgeMaps<N, E>, Report> {
    split_parent_edges(self.schedule, unwrapped_slots(self.outputs))
  }
}

fn unwrapped_slots<T>(outputs: Vec<Option<Arc<T>>>) -> impl Iterator<Item = Result<Option<T>, Report>> {
  outputs.into_iter().enumerate().map(|(index, output)| {
    output
      .map(|output| {
        Arc::into_inner(output).ok_or_else(|| {
          let key = GraphNodeKey(index);
          make_internal_report!("Scan output of node {key} is still shared after the scan")
        })
      })
      .transpose()
  })
}

fn split_parent_edges<N, E>(
  schedule: &GraphSchedule<'_>,
  slots: impl Iterator<Item = Result<Option<(N, PerEdge<E>)>, Report>>,
) -> Result<NodeEdgeMaps<N, E>, Report> {
  let mut nodes = Vec::with_capacity(schedule.key_capacity());
  let mut edges = EdgeMap::new();
  for (index, slot) in slots.enumerate() {
    let Some((node, values)) = slot? else {
      nodes.push(None);
      continue;
    };
    let key = GraphNodeKey(index);
    let parents = schedule.parents(key);
    if values.len() != parents.len() {
      return Err(make_internal_report!(
        "Scan output of node {key} holds {} for {}",
        counted(values.len(), "parent edge value"),
        counted(parents.len(), "parent edge")
      ));
    }
    edges.extend(izip!(parents, values).map(|(&(_, edge_key), value)| (edge_key, value)));
    nodes.push(Some(node));
  }
  Ok(NodeEdgeMaps {
    nodes: NodeMap::from_slots(nodes),
    edges,
  })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NodeEdgeMaps<N, E> {
  pub nodes: NodeMap<N>,
  pub edges: EdgeMap<E>,
}

pub type PerEdge<T> = SmallVec<[T; 1]>;

#[cfg(test)]
mod tests {
  use super::Inbox;
  use crate::node::GraphNodeKey;
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::second_message_for_one_slot(2, &[0, 0], "Scan completed the inbox of node 5 without a message in slot 1")]
  #[case::slot_beyond_the_inputs(     2, &[2],    "Scan found no inbox slot 2 at node 5")]
  #[case::message_after_the_last_slot(1, &[0, 0], "Scan found no inbox slot 0 at node 5")]
  #[case::message_to_a_node_without_inputs(0, &[0], "Scan found no inbox slot 0 at node 5")]
  #[trace]
  fn test_schedule_scan_inbox_rejects_messages_that_do_not_match_the_slots(
    #[case] slots: usize,
    #[case] deliveries: &[usize],
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let node = GraphNodeKey(5);
    let mut inbox = Inbox::new(slots);
    let (&last, earlier) = deliveries.split_last().expect("every case delivers at least one message");
    for &slot in earlier {
      inbox.deliver(node, slot, ())?;
    }

    let result = inbox.deliver(node, last, ());

    assert_error!(
      result,
      format!("{expected}. This is an internal error. Please report it to developers.")
    );
    Ok(())
  }
}
