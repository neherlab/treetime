use crate::testing::map_like::MapLike;
use regex::Regex;
use std::fmt::Debug;
use std::sync::LazyLock;

static NDARRAY_METADATA_RE: LazyLock<Regex> = LazyLock::new(|| {
  Regex::new(r"(, shape=\[[^\]]*\], strides=\[[^\]]*\], layout=\w+ \(0x\w+\))|(, const ndim=\d+)").unwrap()
});

#[macro_export]
macro_rules! pretty_assert_eq {
  ($left:expr, $right:expr) => {{
    pretty_assertions::assert_eq!(
      $crate::testing::assert::format_newlines(format!("{:#?}", $left)),
      $crate::testing::assert::format_newlines(format!("{:#?}", $right))
    );
  }};
  ($left:expr, $right:expr,) => {{
    pretty_assertions::assert_eq!(
      $crate::testing::assert::format_newlines(format!("{:#?}", $left)),
      $crate::testing::assert::format_newlines(format!("{:#?}", $right))
    );
  }};
  ($left:expr, $right:expr, $($arg:tt)+) => {{
    pretty_assertions::assert_eq!(
      $crate::testing::assert::format_newlines(format!("{:#?}", $left)),
      $crate::testing::assert::format_newlines(format!("{:#?}", $right)),
      $($arg)+
    );
  }};
}

pub fn maps_eq<K, V, M, F>(expected: &M, actual: &M, cmp: F) -> bool
where
  M: MapLike<K, V>,
  F: Fn(&V, &V) -> bool,
{
  expected.map_len() == actual.map_len()
    && expected
      .map_iter()
      .all(|(k, exp)| actual.map_get(k).is_some_and(|act| cmp(exp, act)))
}

#[macro_export]
macro_rules! pretty_assert_map_abs_diff_eq {
  ($expected:expr, $actual:expr, epsilon = $eps:expr $(,)?) => {{
    let expected = &$expected;
    let actual = &$actual;
    let eps = $eps;
    if !$crate::testing::assert::maps_eq(expected, actual, |a, b| approx::abs_diff_eq!(a, b, epsilon = eps)) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_map(expected),
        $crate::testing::assert::format_map(actual),
      );
    }
  }};
}

#[macro_export]
macro_rules! pretty_assert_map_ulps_eq {
  ($expected:expr, $actual:expr, max_ulps = $ulps:expr $(,)?) => {{
    let expected = &$expected;
    let actual = &$actual;
    let ulps = $ulps;
    if !$crate::testing::assert::maps_eq(expected, actual, |a, b| approx::ulps_eq!(a, b, max_ulps = ulps)) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_map(expected),
        $crate::testing::assert::format_map(actual),
      );
    }
  }};
}

pub fn format_map<M: Debug>(map: &M) -> String {
  let s = format!("{map:#?}");
  strip_ndarray_metadata(&s)
}

fn strip_ndarray_metadata(s: &str) -> String {
  NDARRAY_METADATA_RE.replace_all(s, "").into_owned()
}

#[macro_export]
macro_rules! pretty_assert_ulps_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::ulps_eq!(lhs, rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{lhs:#?}")),
        $crate::testing::assert::format_array(format!("{rhs:#?}")),
      );
    }
  }};
  ($lhs:expr, $rhs:expr, $($opt:ident = $val:expr,)* $msg:literal $(, $arg:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::ulps_eq!(lhs, rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{lhs:#?}")),
        $crate::testing::assert::format_array(format!("{rhs:#?}")),
        $msg $(, $arg)*
      );
    }
  }};
}

#[macro_export]
macro_rules! pretty_assert_abs_diff_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::abs_diff_eq!(lhs, rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{lhs:#?}")),
        $crate::testing::assert::format_array(format!("{rhs:#?}")),
      );
    }
  }};
  ($lhs:expr, $rhs:expr, $($opt:ident = $val:expr,)* $msg:literal $(, $arg:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if !approx::abs_diff_eq!(lhs, rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{lhs:#?}")),
        $crate::testing::assert::format_array(format!("{rhs:#?}")),
        $msg $(, $arg)*
      );
    }
  }};
}

#[macro_export]
macro_rules! pretty_assert_neg_inf {
  ($actual:expr $(,)?) => {{
    let actual: f64 = $actual;
    if !$crate::testing::assert::is_neg_inf(actual) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_scalar(format!("{actual:#?}")),
        $crate::testing::assert::format_scalar(format!("{:#?}", f64::NEG_INFINITY)),
      );
    }
  }};
  ($actual:expr, $msg:literal $(, $arg:expr)* $(,)?) => {{
    let actual: f64 = $actual;
    if !$crate::testing::assert::is_neg_inf(actual) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_scalar(format!("{actual:#?}")),
        $crate::testing::assert::format_scalar(format!("{:#?}", f64::NEG_INFINITY)),
        $msg $(, $arg)*
      );
    }
  }};
}

/// Replace newlines with NEL (U+0085) to keep multi-line Debug output on one diff line in pretty_assertions.
pub fn format_newlines(s: impl AsRef<str>) -> String {
  s.as_ref().replace('\n', "\u{0085}")
}

pub fn format_array(s: impl AsRef<str>) -> String {
  strip_ndarray_metadata(&format_newlines(s))
}

pub fn format_scalar(s: impl AsRef<str>) -> String {
  format_newlines(s)
}

pub fn is_neg_inf(actual: f64) -> bool {
  actual.is_infinite() && actual.is_sign_negative()
}

#[macro_export]
macro_rules! pretty_assert_array_eq {
  ($lhs:expr, $rhs:expr $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if lhs != rhs {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{lhs:#?}")),
        $crate::testing::assert::format_array(format!("{rhs:#?}")),
      );
    }
  }};
  ($lhs:expr, $rhs:expr, $msg:literal $(, $arg:expr)* $(,)?) => {{
    let lhs = &$lhs;
    let rhs = &$rhs;
    if lhs != rhs {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{lhs:#?}")),
        $crate::testing::assert::format_array(format!("{rhs:#?}")),
        $msg $(, $arg)*
      );
    }
  }};
}

#[macro_export]
macro_rules! assert_error {
  ($result:expr, $expected_message:expr) => {{
    let error = match $result {
      Ok(_) => panic!("expected Err, got Ok"),
      Err(e) => e,
    };
    let actual_message = $crate::error::report_to_string(&error);
    pretty_assertions::assert_eq!(actual_message, $expected_message);
  }};
}

#[cfg(test)]
mod tests {
  use std::panic::{AssertUnwindSafe, catch_unwind, set_hook, take_hook};
  use std::sync::{LazyLock, Mutex};

  use crate::testing::assert::is_neg_inf;

  static PANIC_HOOK_LOCK: LazyLock<Mutex<()>> = LazyLock::new(|| Mutex::new(()));

  #[test]
  fn test_is_neg_inf_accepts_only_negative_infinity() {
    assert!(is_neg_inf(f64::NEG_INFINITY));
    assert!(!is_neg_inf(f64::INFINITY));
    assert!(!is_neg_inf(f64::NAN));
    assert!(!is_neg_inf(-1.0));
    assert!(!is_neg_inf(-0.0));
    assert!(!is_neg_inf(0.0));
  }

  #[test]
  fn test_pretty_assert_neg_inf_accepts_negative_infinity() {
    let result = catch_unwind(|| pretty_assert_neg_inf!(f64::NEG_INFINITY));
    result.unwrap();
  }

  #[test]
  fn test_pretty_assert_neg_inf_rejects_positive_infinity() {
    let result = catch_panic_silent(|| pretty_assert_neg_inf!(f64::INFINITY));
    assert!(result.is_err());
  }

  #[test]
  fn test_pretty_assert_neg_inf_rejects_nan() {
    let result = catch_panic_silent(|| pretty_assert_neg_inf!(f64::NAN));
    assert!(result.is_err());
  }

  fn catch_panic_silent(f: impl FnOnce()) -> std::thread::Result<()> {
    let _guard = match PANIC_HOOK_LOCK.lock() {
      Ok(guard) => guard,
      Err(err) => err.into_inner(),
    };
    let original_hook = take_hook();
    set_hook(Box::new(|_| {}));
    let result = catch_unwind(AssertUnwindSafe(f));
    set_hook(original_hook);
    result
  }
}
