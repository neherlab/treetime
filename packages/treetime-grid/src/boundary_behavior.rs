use crate::edge_cell_law::EdgeCellLaw;
use crate::hard_approach_law::HardApproachLaw;
use crate::soft_tail_law::SoftTailLaw;

pub const DEFAULT_TAIL_FIT_POINTS: usize = 5;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub enum BoundaryBehavior {
  #[default]
  Error,
  Hard,
  HardApproach(HardApproachLaw),
  HardEdgeCell(EdgeCellLaw),
  Linear(SoftTailLaw),
}

impl BoundaryBehavior {
  pub fn is_soft(self) -> bool {
    matches!(self, BoundaryBehavior::Linear(_))
  }
}
