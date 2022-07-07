/// The Traverse enum is used when exploring edges in a graph. It's
/// states signify if an edge should be included in the search, skipped
/// or if the search should stop because a terminating condition has
/// been met such as finding a sink node.
pub enum Traverse {
  Include,
  Skip,
  Finish,
}

pub const OPEN: bool = false;
pub const CLOSED: bool = true;

pub enum Continue<T> {
  Yes(T),
  No(T),
}
