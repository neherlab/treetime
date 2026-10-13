#[cfg(test)]
pub(crate) mod tests {
  use app_commands::json_value::SparseConfig;
  use serde_json::Value;

  pub(crate) fn sparse(value: Value) -> SparseConfig {
    let Value::Object(settings) = value else {
      panic!("a test config is a mapping of settings");
    };
    SparseConfig(settings)
  }
}
