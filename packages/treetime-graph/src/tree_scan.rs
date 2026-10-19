use crate::node::GraphNodeKey;
use crate::tree_view::TreeView;
use eyre::Report;
use itertools::Itertools;
use parking_lot::Mutex;
use rayon::{Scope, current_num_threads, scope};
use std::collections::BTreeMap;
use std::iter::repeat_with;
use std::marker::PhantomData;
use std::mem;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};
use treetime_utils::make_internal_report;

impl TreeView<'_> {
  pub fn scan_down<N, M, Out>(
    &self,
    inputs: impl IntoIterator<Item = (GraphNodeKey, N)>,
    step: impl Fn(GraphNodeKey, N, Option<M>) -> Result<(Out, M), Report> + Sync,
  ) -> Result<ScanOutputs<Out>, Report>
  where
    N: Send,
    M: Clone + Send,
    Out: Send + Sync,
  {
    let scan = Scan::new(self, inputs, step)?;
    scope(|scope| scan.down(scope, self.root(), None));
    scan.finish()
  }

  pub fn scan_up<N, M, Out>(
    &self,
    inputs: impl IntoIterator<Item = (GraphNodeKey, N)>,
    step: impl Fn(GraphNodeKey, N, Vec<M>) -> Result<(Out, M), Report> + Sync,
  ) -> Result<ScanOutputs<Out>, Report>
  where
    N: Send,
    M: Send,
    Out: Send + Sync,
  {
    let scan = Scan::new(self, inputs, step)?;
    let inboxes = (0..self.key_capacity())
      .map(|index| Mutex::new(Inbox::new(self.children(GraphNodeKey(index)).len())))
      .collect_vec();
    let leaves = self
      .preorder()
      .iter()
      .copied()
      .filter(|&key| self.children(key).is_empty())
      .collect_vec();
    let next = AtomicUsize::new(0);
    scope(|scope| {
      for _ in 0..current_num_threads().min(leaves.len()) {
        scope.spawn(|_| {
          while let Some(&leaf) = leaves.get(next.fetch_add(1, Ordering::Relaxed)) {
            scan.up(&inboxes, leaf);
          }
        });
      }
    });
    scan.finish()
  }
}

pub struct ScanOutputs<Out> {
  outputs: Vec<Option<Out>>,
}

impl<Out> ScanOutputs<Out> {
  pub fn into_map(self) -> BTreeMap<GraphNodeKey, Out> {
    bulk_map(self.outputs.len(), entries(self.outputs))
  }
}

fn entries<Out>(outputs: Vec<Option<Out>>) -> impl Iterator<Item = (GraphNodeKey, Out)> {
  outputs
    .into_iter()
    .enumerate()
    .filter_map(|(index, output)| output.map(|output| (GraphNodeKey(index), output)))
}

fn bulk_map<V>(len: usize, entries: impl Iterator<Item = (GraphNodeKey, V)>) -> BTreeMap<GraphNodeKey, V> {
  let mut sorted = Vec::with_capacity(len);
  sorted.extend(entries);
  BTreeMap::from_iter(sorted)
}

struct Scan<'t, 'v, N, P, M, Out, F> {
  tree: &'t TreeView<'v>,
  inputs: Vec<Mutex<Option<N>>>,
  outputs: Vec<OnceLock<Out>>,
  error: OnceLock<Report>,
  step: F,
  messages: PhantomData<fn(P) -> M>,
}

impl<'t, 'v, N, P, M, Out, F> Scan<'t, 'v, N, P, M, Out, F>
where
  F: Fn(GraphNodeKey, N, P) -> Result<(Out, M), Report>,
{
  fn new(tree: &'t TreeView<'v>, inputs: impl IntoIterator<Item = (GraphNodeKey, N)>, step: F) -> Result<Self, Report> {
    let capacity = tree.key_capacity();
    let mut slots = repeat_with(|| Mutex::new(None)).take(capacity).collect_vec();
    for (key, input) in inputs {
      let is_node = key.as_usize() < capacity && (key == tree.root() || tree.parent(key).is_some());
      if !is_node {
        return Err(make_internal_report!(
          "Scan input names node {key}, which is not in the tree"
        ));
      }
      if slots[key.as_usize()].get_mut().replace(input).is_some() {
        return Err(make_internal_report!("Scan input names node {key} twice"));
      }
    }
    if let Some(key) = tree
      .preorder()
      .iter()
      .find(|key| slots[key.as_usize()].get_mut().is_none())
    {
      return Err(make_internal_report!("Scan input is missing node {key}"));
    }
    Ok(Self {
      tree,
      inputs: slots,
      outputs: repeat_with(OnceLock::new).take(capacity).collect_vec(),
      error: OnceLock::new(),
      step,
      messages: PhantomData,
    })
  }

  fn visit(&self, key: GraphNodeKey, received: P) -> Option<M> {
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

  fn fail(&self, report: Report) {
    self.error.get_or_init(|| report);
  }

  fn finish(self) -> Result<ScanOutputs<Out>, Report> {
    if let Some(report) = self.error.into_inner() {
      return Err(report);
    }
    let outputs = self.outputs.into_iter().map(OnceLock::into_inner).collect_vec();
    if let Some(key) = self
      .tree
      .preorder()
      .iter()
      .find(|key| outputs[key.as_usize()].is_none())
    {
      return Err(make_internal_report!("Scan finished without visiting node {key}"));
    }
    Ok(ScanOutputs { outputs })
  }
}

impl<N, M, Out, F> Scan<'_, '_, N, Option<M>, M, Out, F>
where
  N: Send,
  M: Clone + Send,
  Out: Send + Sync,
  F: Fn(GraphNodeKey, N, Option<M>) -> Result<(Out, M), Report> + Sync,
{
  fn down<'s>(&'s self, scope: &Scope<'s>, mut key: GraphNodeKey, mut received: Option<M>) {
    loop {
      let Some(message) = self.visit(key, received) else {
        return;
      };
      let Some((&(first, _), rest)) = self.tree.children(key).split_first() else {
        return;
      };
      for &(child, _) in rest {
        let message = message.clone();
        scope.spawn(move |scope| self.down(scope, child, Some(message)));
      }
      (key, received) = (first, Some(message));
    }
  }
}

impl<N, M, Out, F> Scan<'_, '_, N, Vec<M>, M, Out, F>
where
  N: Send,
  M: Send,
  Out: Send + Sync,
  F: Fn(GraphNodeKey, N, Vec<M>) -> Result<(Out, M), Report> + Sync,
{
  fn up(&self, inboxes: &[Mutex<Inbox<M>>], leaf: GraphNodeKey) {
    let (mut key, mut received) = (leaf, Vec::new());
    loop {
      let Some(message) = self.visit(key, received) else {
        return;
      };
      let Some((parent, _)) = self.tree.parent(key) else {
        return;
      };
      let position = self.tree.position_in_parent(key);
      let delivery = inboxes[parent.as_usize()].lock().deliver(parent, position, message);
      match delivery {
        Ok(Some(siblings)) => (key, received) = (parent, siblings),
        Ok(None) => return,
        Err(report) => {
          self.fail(report);
          return;
        },
      }
    }
  }
}

struct Inbox<M> {
  missing: usize,
  messages: Vec<Option<M>>,
}

impl<M> Inbox<M> {
  fn new(children: usize) -> Self {
    Self {
      missing: children,
      messages: repeat_with(|| None).take(children).collect_vec(),
    }
  }

  fn deliver(&mut self, parent: GraphNodeKey, position: usize, message: M) -> Result<Option<Vec<M>>, Report> {
    let slot = self.messages.get_mut(position).ok_or_else(|| {
      make_internal_report!("Scan found no inbox slot for child position {position} of node {parent}")
    })?;
    *slot = Some(message);
    self.missing -= 1;
    if self.missing > 0 {
      return Ok(None);
    }
    mem::take(&mut self.messages)
      .into_iter()
      .enumerate()
      .map(|(child, message)| {
        message.ok_or_else(|| {
          make_internal_report!(
            "Scan completed the inbox of node {parent} without a message from child position {child}"
          )
        })
      })
      .try_collect()
      .map(Some)
  }
}

#[cfg(test)]
mod tests {
  use super::Inbox;
  use crate::node::GraphNodeKey;
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::second_message_for_one_slot(  2, &[0, 0], "Scan completed the inbox of node 5 without a message from child position 1")]
  #[case::position_beyond_the_children( 2, &[2],    "Scan found no inbox slot for child position 2 of node 5")]
  #[case::message_after_the_last_child( 1, &[0, 0], "Scan found no inbox slot for child position 0 of node 5")]
  #[trace]
  fn test_tree_scan_inbox_rejects_messages_that_do_not_match_the_children(
    #[case] children: usize,
    #[case] positions: &[usize],
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let parent = GraphNodeKey(5);
    let mut inbox = Inbox::new(children);
    let (&last, earlier) = positions.split_last().expect("every case delivers at least one message");
    for &position in earlier {
      inbox.deliver(parent, position, ())?;
    }

    let result = inbox.deliver(parent, last, ());

    assert_error!(
      result,
      format!("{expected}. This is an internal error. Please report it to developers.")
    );
    Ok(())
  }
}
