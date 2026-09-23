use eyre::{Report, WrapErr};
use std::env::{self, VarError};

pub fn env_var_optional(name: &str) -> Result<Option<String>, Report> {
  env_var_result_to_optional(name, env::var(name))
}

fn env_var_result_to_optional(name: &str, result: Result<String, VarError>) -> Result<Option<String>, Report> {
  match result {
    Ok(value) => Ok(Some(value)),
    Err(VarError::NotPresent) => Ok(None),
    Err(error @ VarError::NotUnicode(_)) => {
      Err(error).wrap_err_with(|| format!("When reading environment variable '{name}'"))
    },
  }
}

#[cfg(test)]
mod tests {
  use super::env_var_result_to_optional;
  use crate::assert_error;
  use pretty_assertions::assert_eq;
  use std::env::VarError;
  use std::ffi::OsString;

  #[test]
  fn test_env_var_result_to_optional_present() {
    let actual = env_var_result_to_optional("NAME", Ok("value".to_owned())).unwrap();
    assert_eq!(Some("value".to_owned()), actual);
  }

  #[test]
  fn test_env_var_result_to_optional_unset_is_none() {
    let actual = env_var_result_to_optional("NAME", Err(VarError::NotPresent)).unwrap();
    assert_eq!(None, actual);
  }

  #[test]
  fn test_env_var_result_to_optional_not_unicode_is_error() {
    let actual = env_var_result_to_optional("NAME", Err(VarError::NotUnicode(OsString::from("x"))));
    assert_error!(
      actual,
      "When reading environment variable 'NAME': environment variable was not valid unicode: \"x\""
    );
  }
}
