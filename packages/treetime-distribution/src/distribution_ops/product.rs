use crate::Distribution;
use crate::distribution_ops::multiply::{distribution_multiplication, multiply_functions};
use crate::policy::YAxisPolicy;
use eyre::Report;
use treetime_grid::MaxGridPoints;
use treetime_utils::make_internal_error;

pub fn distribution_product<Y: YAxisPolicy>(
  factors: &[&Distribution<Y>],
  max_points: MaxGridPoints,
) -> Result<Distribution<Y>, Report> {
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
    Some(multiply_functions(&functions, max_points)?)
  };

  let mut other_product: Option<Distribution<Y>> = None;
  for factor in others {
    other_product = Some(match other_product {
      None => factor.clone(),
      Some(current) => distribution_multiplication(&current, factor, max_points)?,
    });
  }

  match (function_product, other_product) {
    (Some(functions), Some(others)) => distribution_multiplication(&functions, &others, max_points),
    (Some(single), None) | (None, Some(single)) => Ok(single),
    (None, None) => make_internal_error!("distribution_product requires at least one factor"),
  }
}
