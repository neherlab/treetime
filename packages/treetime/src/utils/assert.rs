#[macro_export]
macro_rules! pretty_assert_eq {
  ($left:expr, $right:expr) => {{
    pretty_assertions::assert_eq!(
      format!("{:#?}", $left).replace("\n", "\u{0085}"),
      format!("{:#?}", $right).replace("\n", "\u{0085}")
    );
  }};
}
