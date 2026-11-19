use std::marker::PhantomData;
use treetime_grid::YAxisPolicy;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PolicyMarker<Y: YAxisPolicy>(PhantomData<Y>);

impl<Y: YAxisPolicy> PolicyMarker<Y> {
  pub(crate) fn new() -> Self {
    Self(PhantomData)
  }
}
