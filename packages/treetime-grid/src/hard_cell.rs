use crate::{GridEdge, HardApproachLaw};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HardCell {
  pub law: HardApproachLaw,
  pub anchor: GridEdge,
}

impl HardCell {
  pub fn eval(&self, t: f64) -> f64 {
    self.law.eval(self.anchor, t)
  }
}
