#[macro_export]
macro_rules! pretty_assert_eq {
  ($left:expr, $right:expr) => {{
    pretty_assertions::assert_eq!(
      format_newlines(format!("{:#?}", $left))
      format_newlines(format!("{:#?}", $right))
    );
  }};
  ($left:expr, $right:expr,) => {{
    pretty_assertions::assert_eq!(
      format_newlines(format!("{:#?}", $left))
      format_newlines(format!("{:#?}", $right))
    );
  }};
}

#[macro_export]
macro_rules! pretty_assert_ulps_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*) => {{
    if ! approx::ulps_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::assert::format_array(format!("{:#?}", $lhs)),
        $crate::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {{
    if ! approx::ulps_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::assert::format_array(format!("{:#?}", $lhs)),
        $crate::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
}

#[macro_export]
macro_rules! pretty_assert_abs_diff_eq {
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*) => {{
    if ! approx::abs_diff_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::assert::format_array(format!("{:#?}", $lhs)),
        $crate::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {{
    if ! approx::abs_diff_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::assert::format_array(format!("{:#?}", $lhs)),
        $crate::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
}

pub fn format_newlines(s: impl AsRef<str>) -> String {
  s.as_ref().replace('\n', "\u{0085}")
}

pub fn format_array(s: impl AsRef<str>) -> String {
  let re = regex::Regex::new(r", shape=.*").unwrap();
  re.replace(&format_newlines(s), "").into_owned()
}

#[macro_export]
macro_rules! assert_error {
  ($result:expr, $expected_message:expr) => {{
    let error = $result.unwrap_err();
    let actual_message = $crate::error::report_to_string(&error);
    pretty_assertions::assert_eq!(actual_message, $expected_message);
  }};
}
