use crate::testing::map_like::MapLike;
use std::fmt::Debug;

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
    let eps = $eps;
    if !$crate::testing::assert::maps_eq(&$expected, &$actual, |a, b| approx::abs_diff_eq!(a, b, epsilon = eps)) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_map(&$expected),
        $crate::testing::assert::format_map(&$actual),
      );
    }
  }};
}

#[macro_export]
macro_rules! pretty_assert_map_ulps_eq {
  ($expected:expr, $actual:expr, max_ulps = $ulps:expr $(,)?) => {{
    let ulps = $ulps;
    if !$crate::testing::assert::maps_eq(&$expected, &$actual, |a, b| approx::ulps_eq!(a, b, max_ulps = ulps)) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_map(&$expected),
        $crate::testing::assert::format_map(&$actual),
      );
    }
  }};
}

pub fn format_map<M: Debug>(map: &M) -> String {
  let s = format!("{map:#?}");
  strip_ndarray_metadata(&s)
}

fn strip_ndarray_metadata(s: &str) -> String {
  let re =
    regex::Regex::new(r"(, shape=\[[^\]]*\], strides=\[[^\]]*\], layout=\w+ \(0x\w+\))|(, const ndim=\d+)").unwrap();
  re.replace_all(s, "").into_owned()
}

#[macro_export]
macro_rules! pretty_assert_ulps_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    if ! approx::ulps_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{:#?}", $lhs)),
        $crate::testing::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
  ($lhs:expr, $rhs:expr, $($opt:ident = $val:expr,)* $msg:literal $(, $arg:expr)* $(,)?) => {{
    if ! approx::ulps_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{:#?}", $lhs)),
        $crate::testing::assert::format_array(format!("{:#?}", $rhs)),
        $msg $(, $arg)*
      );
    }
  }};
}

#[macro_export]
macro_rules! pretty_assert_abs_diff_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)* $(,)?) => {{
    if ! approx::abs_diff_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{:#?}", $lhs)),
        $crate::testing::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
  ($lhs:expr, $rhs:expr, $($opt:ident = $val:expr,)* $msg:literal $(, $arg:expr)* $(,)?) => {{
    if ! approx::abs_diff_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::testing::assert::format_array(format!("{:#?}", $lhs)),
        $crate::testing::assert::format_array(format!("{:#?}", $rhs)),
        $msg $(, $arg)*
      );
    }
  }};
}

pub fn format_newlines(s: impl AsRef<str>) -> String {
  s.as_ref().replace('\n', "\u{0085}")
}

pub fn format_array(s: impl AsRef<str>) -> String {
  strip_ndarray_metadata(&format_newlines(s))
}

#[macro_export]
macro_rules! assert_error {
  ($result:expr, $expected_message:expr) => {{
    let error = $result.unwrap_err();
    let actual_message = $crate::error::report_to_string(&error);
    pretty_assertions::assert_eq!(actual_message, $expected_message);
  }};
}
