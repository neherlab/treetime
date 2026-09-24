#[cfg(test)]
mod tests {
  use app_commands::command::AppCommand;
  use helpers::{config_for, output_files, run_cli, run_napi, run_server};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( AppCommand::Timetree)]
  #[case::clock(    AppCommand::Clock)]
  #[case::ancestral(AppCommand::Ancestral)]
  #[case::mugration(AppCommand::Mugration)]
  #[case::optimize( AppCommand::Optimize)]
  #[case::prune(    AppCommand::Prune)]
  #[trace]
  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_transport_parity_cli_server_and_napi_write_identical_outputs(#[case] command: AppCommand) {
    let config = config_for(command);
    let work = tempdir().unwrap();

    let cli_dir = work.path().join("cli");
    run_cli(command, &config, work.path(), &cli_dir);
    let server_dir = run_server(command, &config, &work.path().join("server")).await;
    let napi_dir = work.path().join("napi");
    run_napi(command, &config, &napi_dir);

    let cli = output_files(&cli_dir);
    assert!(!cli.is_empty(), "the CLI wrote no outputs");
    assert_eq!(cli, output_files(&server_dir), "server outputs differ from CLI outputs");
    assert_eq!(cli, output_files(&napi_dir), "N-API outputs differ from CLI outputs");
  }

  mod helpers {
    use crate::cli::treetime_cli::treetime_parse_cli_args;
    use crate::run::run_command;
    use app_commands::command::AppCommand;
    use app_commands::job::JobRegistry;
    use app_napi::jobs::start_job;
    use app_server::create_router;
    use app_server::state::ServerConfig;
    use axum::body::Body;
    use axum::http::Request;
    use serde_json::{Value, json};
    use std::collections::BTreeMap;
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use tower::ServiceExt;
    use treetime::progress::NoopProgress;

    pub(super) fn data_dir() -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .canonicalize()
        .unwrap()
    }

    pub(super) fn config_for(command: AppCommand) -> Value {
      let data = data_dir();
      let zika = data.join("zika/20");
      let flu = data.join("flu/h3n2/20");
      match command {
        AppCommand::Timetree => json!({
          "tree": zika.join("tree.nwk"),
          "metadata": zika.join("metadata.tsv"),
          "alignment": [zika.join("aln.fasta.xz")],
          "max_iter": 2,
          "seed": 7,
        }),
        AppCommand::Clock => json!({
          "tree": zika.join("tree.nwk"),
          "metadata": zika.join("metadata.tsv"),
        }),
        AppCommand::Ancestral => json!({
          "tree": zika.join("tree.nwk"),
          "alignment": [zika.join("aln.fasta.xz")],
        }),
        AppCommand::Mugration => json!({
          "tree": zika.join("tree.nwk"),
          "metadata": zika.join("metadata.tsv"),
          "attribute": "country",
        }),
        AppCommand::Optimize => json!({
          "tree": flu.join("tree.nwk"),
          "alignment": [flu.join("aln.fasta.xz")],
        }),
        AppCommand::Prune => json!({
          "tree": zika.join("tree.nwk"),
          "alignment": [zika.join("aln.fasta.xz")],
          "prune_short": 1e-6,
          "prune_empty": true,
          "merge_shared_mutations": true,
        }),
      }
    }

    pub(super) fn run_cli(command: AppCommand, config: &Value, work: &Path, output_dir: &Path) {
      let config_path = work.join("config.json");
      fs::write(&config_path, config.to_string()).unwrap();
      let argv = [
        "treetime".to_owned(),
        command.to_string(),
        "--config".to_owned(),
        config_path.to_string_lossy().into_owned(),
        "--output-all".to_owned(),
        output_dir.to_string_lossy().into_owned(),
      ];
      let args = treetime_parse_cli_args(argv).unwrap();
      run_command(args.command, &NoopProgress).unwrap();
    }

    pub(super) async fn run_server(command: AppCommand, config: &Value, out_dir: &Path) -> PathBuf {
      let router = create_router(
        ServerConfig {
          data_dir: data_dir(),
          out_dir: out_dir.to_path_buf(),
        },
        None,
      )
      .unwrap();
      let request = Request::post(format!("/api/{command}"))
        .header("content-type", "application/json")
        .body(Body::from(config.to_string()))
        .unwrap();
      let response = router.oneshot(request).await.unwrap();
      let body = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
      let text = String::from_utf8(body.to_vec()).unwrap();
      let terminal: Value = text
        .split("\n\n")
        .filter(|block| block.contains("event: terminal"))
        .filter_map(|block| block.lines().find_map(|line| line.strip_prefix("data: ")))
        .map(|data| serde_json::from_str(data).unwrap())
        .next()
        .unwrap();
      assert_eq!(Some("ok"), terminal["status"].as_str(), "server job failed: {terminal}");
      out_dir.join(terminal["job_id"].as_str().unwrap())
    }

    pub(super) fn run_napi(command: AppCommand, config: &Value, output_dir: &Path) {
      let mut config = config.clone();
      config["output_all"] = json!(output_dir);
      let registry = Arc::new(JobRegistry::default());
      let job = start_job(&registry, "parity", &command.to_string(), &config.to_string()).unwrap();
      let terminal = serde_json::to_value(job.run(|_event| {})).unwrap();
      assert_eq!(Some("ok"), terminal["status"].as_str(), "N-API job failed: {terminal}");
    }

    pub(super) fn output_files(dir: &Path) -> BTreeMap<String, String> {
      fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_file())
        .map(|path| {
          let name = path.file_name().unwrap().to_string_lossy().into_owned();
          let content = normalized_content(&path);
          (name, content)
        })
        .collect()
    }

    fn normalized_content(path: &Path) -> String {
      let name = path.file_name().unwrap().to_string_lossy();
      let bytes = fs::read(path).unwrap();
      if name.ends_with(".auspice.json") {
        let mut value: Value = serde_json::from_slice(&bytes).unwrap();
        remove_key(&mut value, &["meta", "updated"]);
        return value.to_string();
      }
      if name.ends_with(".augur-node-data.json") {
        let mut value: Value = serde_json::from_slice(&bytes).unwrap();
        remove_key(&mut value, &["generated_by", "version"]);
        return value.to_string();
      }
      String::from_utf8_lossy(&bytes).into_owned()
    }

    fn remove_key(value: &mut Value, key_path: &[&str]) {
      let Some((last, parents)) = key_path.split_last() else {
        return;
      };
      let parent = parents
        .iter()
        .try_fold(value, |value, key| value.get_mut(*key))
        .and_then(Value::as_object_mut);
      if let Some(parent) = parent {
        parent.remove(*last);
      }
    }
  }
}
