//! Proptest-compatible approximate equality assertions with pretty diffs.
//!
//! Use `return Err(TestCaseError::fail(...))` instead of panicking,
//! avoiding noisy backtraces during proptest shrinking.
//!
//! Dependencies resolved at expansion site: `approx`, `proptest`, `pretty_assertions`.

#[macro_export]
macro_rules! prop_assert_abs_diff_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::abs_diff_eq!(lhs, rhs $(, $opt = $val)*) {
      let lhs_fmt = format!("{lhs:?}");
      let rhs_fmt = format!("{rhs:?}");
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_abs_diff_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_ulps_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::ulps_eq!(lhs, rhs $(, $opt = $val)*) {
      let lhs_fmt = format!("{lhs:?}");
      let rhs_fmt = format!("{rhs:?}");
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_ulps_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_relative_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::relative_eq!(lhs, rhs $(, $opt = $val)*) {
      let lhs_fmt = format!("{lhs:?}");
      let rhs_fmt = format!("{rhs:?}");
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_relative_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_abs_diff_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::abs_diff_eq!(lhs, rhs $(, $opt = $val)*) {
      let lhs_fmt = $crate::testing::assert::format_array(format!("{lhs:#?}"));
      let rhs_fmt = $crate::testing::assert::format_array(format!("{rhs:#?}"));
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_array_abs_diff_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_ulps_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::ulps_eq!(lhs, rhs $(, $opt = $val)*) {
      let lhs_fmt = $crate::testing::assert::format_array(format!("{lhs:#?}"));
      let rhs_fmt = $crate::testing::assert::format_array(format!("{rhs:#?}"));
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_array_ulps_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_relative_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::relative_eq!(lhs, rhs $(, $opt = $val)*) {
      let lhs_fmt = $crate::testing::assert::format_array(format!("{lhs:#?}"));
      let rhs_fmt = $crate::testing::assert::format_array(format!("{rhs:#?}"));
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_array_relative_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_map_abs_diff_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !$crate::testing::assert::maps_eq(lhs, rhs, |a, b| approx::abs_diff_eq!(a, b $(, $opt = $val)*)) {
      let lhs_fmt = $crate::testing::assert::format_map(lhs);
      let rhs_fmt = $crate::testing::assert::format_map(rhs);
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_map_abs_diff_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_map_ulps_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !$crate::testing::assert::maps_eq(lhs, rhs, |a, b| approx::ulps_eq!(a, b $(, $opt = $val)*)) {
      let lhs_fmt = $crate::testing::assert::format_map(lhs);
      let rhs_fmt = $crate::testing::assert::format_map(rhs);
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_map_ulps_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_map_relative_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !$crate::testing::assert::maps_eq(lhs, rhs, |a, b| approx::relative_eq!(a, b $(, $opt = $val)*)) {
      let lhs_fmt = $crate::testing::assert::format_map(lhs);
      let rhs_fmt = $crate::testing::assert::format_map(rhs);
      return Err(proptest::test_runner::TestCaseError::fail(format!(
        "prop_assert_map_relative_eq!({}, {}{})\n\n{}",
        stringify!($lhs),
        stringify!($rhs),
        stringify!($(, $opt = $val)*),
        pretty_assertions::StrComparison::new(&lhs_fmt, &rhs_fmt),
      )));
    }
  }};
}
