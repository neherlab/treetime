use crate::hard_approach_law::HardApproachLaw;
use crate::soft_tail_law::SoftTailLaw;
use serde::{Deserialize, Serialize};

pub const DEFAULT_TAIL_FIT_POINTS: usize = 5;

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryBehavior {
  #[default]
  Error,
  Hard,
  HardApproach(HardApproachLaw),
  Linear(SoftTailLaw),
}

impl BoundaryBehavior {
  pub fn is_soft(self) -> bool {
    matches!(self, BoundaryBehavior::Linear(_))
  }

  pub fn soft_law(self) -> Option<SoftTailLaw> {
    match self {
      BoundaryBehavior::Linear(law) => Some(law),
      _ => None,
    }
  }
}
