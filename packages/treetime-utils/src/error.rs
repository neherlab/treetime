use color_eyre::Report;
use std::any::Any;

pub fn panic_message(payload: &(dyn Any + Send)) -> Option<String> {
  payload
    .downcast_ref::<&str>()
    .map(|message| (*message).to_owned())
    .or_else(|| payload.downcast_ref::<String>().cloned())
}

pub struct ReportChain {
  pub message: String,
  pub causes: Vec<String>,
}

impl ReportChain {
  pub fn of(report: &Report) -> Self {
    let mut chain = report.chain().map(ToString::to_string);
    Self {
      message: chain.next().unwrap_or_default(),
      causes: chain.collect(),
    }
  }
}

pub fn report_to_string(report: &Report) -> String {
  let strings: Vec<String> = report.chain().map(ToString::to_string).collect();
  strings.join(": ")
}

pub fn to_eyre_error<T, E: Into<eyre::Error>>(val_or_err: Result<T, E>) -> Result<T, Report> {
  val_or_err.map_err(|report| make_report!(report))
}

#[macro_export(local_inner_macros)]
macro_rules! make_error {
  ($($arg:tt)*) => {
    {
      Err(eyre::eyre!(std::format!($($arg)*)))
    }
  };
}

pub use make_error;

#[macro_export(local_inner_macros)]
macro_rules! make_report {
  ($($arg:tt)*) => {
    {
      eyre::eyre!($($arg)*)
    }
  };
}

pub use make_report;

#[macro_export(local_inner_macros)]
macro_rules! make_internal_error {
  ($($arg:tt)*) => {
    {
      let msg_external = std::format!($($arg)*);
      let msg = std::format!("{msg_external}. This is an internal error. Please report it to developers.");
      Err(eyre::eyre!(msg))
    }
  };
}

pub use make_internal_error;

#[macro_export(local_inner_macros)]
macro_rules! make_internal_report {
  ($($arg:tt)*) => {
    {
      let msg_external = std::format!($($arg)*);
      let msg = std::format!("{msg_external}. This is an internal error. Please report it to developers.");
      eyre::eyre!(msg)
    }
  };
}

pub use make_internal_report;

#[cfg(test)]
mod tests {
  use crate::error::panic_message;
  use pretty_assertions::assert_eq;
  use std::panic::{catch_unwind, panic_any};

  #[test]
  fn test_error_panic_message_of_a_literal_a_formatted_and_a_foreign_payload() {
    let literal = catch_unwind(|| panic!("literal")).unwrap_err();
    let formatted = catch_unwind(|| panic!("formatted {}", 1)).unwrap_err();
    let foreign = catch_unwind(|| panic_any(7_u8)).unwrap_err();
    assert_eq!(
      (Some("literal".to_owned()), Some("formatted 1".to_owned()), None),
      (
        panic_message(literal.as_ref()),
        panic_message(formatted.as_ref()),
        panic_message(foreign.as_ref())
      )
    );
  }
}
