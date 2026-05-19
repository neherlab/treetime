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
macro_rules! prop_assert_array_finite {
  ($arr:expr $(,)?) => {{
    let arr = &$arr;
    for (idx, val) in arr.indexed_iter() {
      if !val.is_finite() {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_finite!({}): element {idx:?} = {val}",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_nonneg {
  ($arr:expr $(,)?) => {{
    let arr = &$arr;
    for (idx, val) in arr.indexed_iter() {
      if *val < 0.0 {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_nonneg!({}): element {idx:?} = {val} < 0",
          stringify!($arr),
        )));
      }
    }
  }};
  ($arr:expr, epsilon = $eps:expr $(,)?) => {{
    let arr = &$arr;
    let eps = $eps;
    for (idx, val) in arr.indexed_iter() {
      if *val < -eps {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_nonneg!({}): element {idx:?} = {val} < -{eps}",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_positive {
  ($arr:expr $(,)?) => {{
    let arr = &$arr;
    for (idx, val) in arr.indexed_iter() {
      if *val <= 0.0 {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_positive!({}): element {idx:?} = {val} <= 0",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_upper_bounded {
  ($arr:expr, bound = $bound:expr, epsilon = $eps:expr $(,)?) => {{
    let arr = &$arr;
    let bound = $bound;
    let eps = $eps;
    for (idx, val) in arr.indexed_iter() {
      if *val > bound + eps {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_upper_bounded!({}): element {idx:?} = {val} > {bound} + {eps}",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_diag_abs {
  ($arr:expr, epsilon = $eps:expr $(,)?) => {{
    let arr = &$arr;
    let eps = $eps;
    for (i, val) in arr.diag().iter().enumerate() {
      if val.abs() >= eps {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_diag_abs!({}): diag[{i}] = {val}, |{val}| >= {eps}",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_offdiag_nonneg {
  ($arr:expr $(,)?) => {{
    let arr = &$arr;
    for ((i, j), val) in arr.indexed_iter() {
      if i != j && *val < 0.0 {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_offdiag_nonneg!({}): [{i},{j}] = {val} < 0",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_offdiag_positive {
  ($arr:expr $(,)?) => {{
    let arr = &$arr;
    for ((i, j), val) in arr.indexed_iter() {
      if i != j && *val <= 0.0 {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_offdiag_positive!({}): [{i},{j}] = {val} <= 0",
          stringify!($arr),
        )));
      }
    }
  }};
}

#[macro_export]
macro_rules! prop_assert_array_diag_nonpositive {
  ($arr:expr $(,)?) => {{
    let arr = &$arr;
    for (i, val) in arr.diag().iter().enumerate() {
      if *val > 0.0 {
        return Err(proptest::test_runner::TestCaseError::fail(format!(
          "prop_assert_array_diag_nonpositive!({}): diag[{i}] = {val} > 0",
          stringify!($arr),
        )));
      }
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
