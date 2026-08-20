use crate::Distribution;
use crate::distribution_ops::multiply::{distribution_multiplication, multiply_functions};
use crate::policy::YAxisPolicy;
use eyre::Report;
use treetime_utils::make_internal_error;

/// Product of N distribution factors: the N-ary generalization of [`distribution_multiplication`].
///
/// A product of independent factors is pointwise -- a sum of neg-log ordinates -- so it is exact,
/// associative, and independent of factor order. The gridded `Function` factors are co-located on one
/// common grid and reduced once ([`multiply_functions`], the same primitive the pairwise
/// `multiply_function_function` uses), so every function is interpolated at most once regardless of
/// fan-out and no accumulator is re-resampled per factor. The remaining exact factors (`Point`,
/// `Range`, `Formula`) need no grid: their product is formed with the pairwise
/// [`distribution_multiplication`] and multiplied into the function product in one final step.
///
/// An `Empty` factor makes the product `Empty` (an empty operand is disjoint from every domain).
pub fn distribution_product<Y: YAxisPolicy>(factors: &[&Distribution<Y>]) -> Result<Distribution<Y>, Report> {
  let mut functions = Vec::new();
  let mut others = Vec::new();
  for &factor in factors {
    match factor {
      Distribution::Empty => return Ok(Distribution::empty()),
      Distribution::Function(function) => functions.push(function),
      Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_) => others.push(factor),
    }
  }

  let function_product = if functions.is_empty() {
    None
  } else {
    Some(multiply_functions(&functions)?)
  };

  let mut other_product: Option<Distribution<Y>> = None;
  for factor in others {
    other_product = Some(match other_product {
      None => factor.clone(),
      Some(current) => distribution_multiplication(&current, factor)?,
    });
  }

  match (function_product, other_product) {
    (Some(functions), Some(others)) => distribution_multiplication(&functions, &others),
    (Some(single), None) | (None, Some(single)) => Ok(single),
    (None, None) => make_internal_error!("distribution_product requires at least one factor"),
  }
}
