#[cfg(test)]
mod tests {
  use crate::job::module_log;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::io::json::to_json_value;

  #[test]
  fn test_job_log_names_the_module_and_the_thread_count_of_the_run() {
    assert_eq!(
      json!({
        "type": "log",
        "data": { "level": "info", "message": "The run uses the single-threaded WebAssembly module with 1 thread" },
      }),
      to_json_value(&module_log(1)).unwrap()
    );
  }
}
