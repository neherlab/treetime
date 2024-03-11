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
        $crate::utils::assert::format_array(format!("{:#?}", $lhs)),
        $crate::utils::assert::format_array(format!("{:#?}", $rhs)),
      );
    }
  }};
  ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {{
    if ! approx::ulps_eq!($lhs, $rhs, $($opt = $val,)*) {
      pretty_assertions::assert_eq!(
        $crate::utils::assert::format_array(format!("{:#?}", $lhs)),
        $crate::utils::assert::format_array(format!("{:#?}", $rhs)),
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
