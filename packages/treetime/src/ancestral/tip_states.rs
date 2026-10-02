#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct TipStates {
  pub(crate) include_leaves: bool,
  pub(crate) impute: bool,
}
