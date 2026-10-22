#[cfg(feature = "clap")]
use clap::{Arg, ValueHint};

#[cfg(feature = "clap")]
pub fn value_name_from_hint(arg: Arg) -> Arg {
  let derived = arg.get_id().as_str().to_uppercase();
  let named_after_field = arg
    .get_value_names()
    .is_some_and(|names| names.iter().map(ToString::to_string).eq([derived]));
  if !named_after_field {
    return arg;
  }
  match arg.get_value_hint() {
    ValueHint::FilePath => arg.value_name("FILE"),
    ValueHint::DirPath => arg.value_name("DIR"),
    ValueHint::AnyPath => arg.value_name("PATH"),
    _ => arg,
  }
}
