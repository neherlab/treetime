use crate::model::graph::NewickGraph;
use std::mem;

pub struct Preorder<'g> {
  graph: &'g NewickGraph,
  stack: Vec<usize>,
  visited: Vec<bool>,
}

impl<'g> Preorder<'g> {
  pub(crate) fn new(graph: &'g NewickGraph) -> Self {
    let node_count = graph.node_count();
    let stack = if graph.root() < node_count {
      vec![graph.root()]
    } else {
      Vec::new()
    };
    Self {
      graph,
      stack,
      visited: vec![false; node_count],
    }
  }
}

impl Iterator for Preorder<'_> {
  type Item = usize;

  fn next(&mut self) -> Option<usize> {
    while let Some(node) = self.stack.pop() {
      if mem::replace(&mut self.visited[node], true) {
        continue;
      }
      let children = self.graph.child_edges(node);
      self
        .stack
        .extend(children.iter().rev().map(|&edge| self.graph.edge(edge).child()));
      return Some(node);
    }
    None
  }
}

pub struct Postorder<'g> {
  graph: &'g NewickGraph,
  stack: Vec<Frame>,
  marks: Vec<Mark>,
}

impl<'g> Postorder<'g> {
  pub(crate) fn new(graph: &'g NewickGraph) -> Self {
    let node_count = graph.node_count();
    let stack = if graph.root() < node_count {
      vec![Frame {
        node: graph.root(),
        next_child: 0,
      }]
    } else {
      Vec::new()
    };
    Self {
      graph,
      stack,
      marks: vec![Mark::New; node_count],
    }
  }

  pub(crate) fn has_cycle(graph: &NewickGraph) -> bool {
    let mut postorder = Postorder::new(graph);
    while postorder.step().is_some() {}
    postorder.marks.contains(&Mark::Cycle)
  }

  fn step(&mut self) -> Option<usize> {
    while let Some(frame) = self.stack.last_mut() {
      let node = frame.node;
      if frame.next_child == 0 && self.marks[node] == Mark::New {
        self.marks[node] = Mark::Open;
      }
      if let Some(&edge) = self.graph.child_edges(node).get(frame.next_child) {
        frame.next_child += 1;
        let child = self.graph.edge(edge).child();
        match self.marks[child] {
          Mark::New => self.stack.push(Frame {
            node: child,
            next_child: 0,
          }),
          Mark::Open => self.marks[child] = Mark::Cycle,
          Mark::Done | Mark::Cycle => {},
        }
        continue;
      }
      self.stack.pop();
      if self.marks[node] == Mark::Open {
        self.marks[node] = Mark::Done;
      }
      return Some(node);
    }
    None
  }
}

impl Iterator for Postorder<'_> {
  type Item = usize;

  fn next(&mut self) -> Option<usize> {
    self.step()
  }
}

#[derive(Clone, Copy)]
struct Frame {
  node: usize,
  next_child: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mark {
  New,
  Open,
  Done,
  Cycle,
}
