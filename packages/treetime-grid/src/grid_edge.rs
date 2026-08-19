/// The live anchor an edge-relative boundary law reads on evaluation: a grid edge's coordinate and
/// its stored neg-log ordinate.
///
/// Both [`HardApproachLaw`](crate::HardApproachLaw) and [`SoftTailLaw`](crate::SoftTailLaw) store
/// only shape and read this anchor from the current grid edge when evaluated, so the law survives
/// re-windowing and normalization with no refit. Grouping the coordinate and the ordinate into one
/// named struct makes the two impossible to transpose at a call site, where both are plain `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridEdge {
  /// Coordinate of the grid edge on the time axis (`t_edge`).
  pub t: f64,
  /// Stored neg-log ordinate at the edge (`y_edge`).
  pub y: f64,
}
