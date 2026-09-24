#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TipStates {
  pub(crate) include_leaves: bool,
  pub(crate) impute: bool,
}
