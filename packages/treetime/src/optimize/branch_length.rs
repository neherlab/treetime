use crate::branch_lengths::is_valid_branch_length_value;
use crate::make_report;
use eyre::Report;

pub(crate) fn validate_branch_length_value(branch_length: f64) -> Result<(), Report> {
  if is_valid_branch_length_value(branch_length) {
    Ok(())
  } else {
    Err(make_report!(
      "branch length must be finite and non-negative, got {branch_length}"
    ))
  }
}
