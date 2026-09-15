#[cfg(test)]
mod tests {
  use crate::commands::timetree::result::TimetreeEdgeOut;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_result_profile_branch_length_prefers_clock_length() {
    let edge = TimetreeEdgeOut {
      branch_length: Some(0.5),
      time_length: Some(2.0),
      clock_branch_length: Some(0.25),
      gamma: 1.0,
    };
    assert_eq!(Some(0.25), edge.profile_branch_length());
  }

  #[test]
  fn test_result_profile_branch_length_falls_back_to_branch_length() {
    let edge = TimetreeEdgeOut {
      branch_length: Some(0.5),
      time_length: Some(2.0),
      clock_branch_length: None,
      gamma: 1.0,
    };
    assert_eq!(Some(0.5), edge.profile_branch_length());
  }

  #[test]
  fn test_result_profile_branch_length_none_when_both_absent() {
    let edge = TimetreeEdgeOut {
      branch_length: None,
      time_length: None,
      clock_branch_length: None,
      gamma: 1.0,
    };
    assert_eq!(None, edge.profile_branch_length());
  }
}
