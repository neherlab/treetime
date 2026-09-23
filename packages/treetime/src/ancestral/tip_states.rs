#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TipStates {
  pub include_leaves: bool,
  pub impute: bool,
}
