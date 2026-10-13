macro_rules! impl_from_mirror {
  ($cli:ident => $core:ident { $($variant:ident),+ $(,)? }) => {
    impl From<$cli> for $core {
      fn from(value: $cli) -> Self {
        match value {
          $($cli::$variant => Self::$variant),+
        }
      }
    }
  };
}

pub(crate) use impl_from_mirror;
