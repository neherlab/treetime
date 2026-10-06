#[cfg(test)]
mod tests {
  use crate::backend::{DesktopService, PortScope, stream_response};
  use crate::port::PortReply;
  use app_commands::app_paths::AppFolderEnv;
  use app_output::output_plan::OutputSelection;
  use axum::body::{Body, Bytes};
  use axum::response::Response;
  use helpers::{
    Ended, ancestral_request, archive_contents, event_ids, fetch, fetch_bytes, fetch_json, fetch_json_in, header,
    open_fetch, wait_for_terminal,
  };
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::fs;
  use std::io;
  use std::sync::atomic::AtomicBool;
  use std::sync::mpsc::RecvTimeoutError;
  use std::time::Duration;
  use tempfile::tempdir;
  use treetime_schema::version_info;
  use treetime_utils::io::json::to_json_value;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_backend_created_run_starts_at_once_and_writes_its_outputs() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(false)));
    let id = record["id"].as_str().unwrap().to_owned();
    let terminal = wait_for_terminal(&service, &id);

    let (_, record) = fetch_json(&service, "GET", &format!("/api/runs/{id}"), None);
    let (_, files) = fetch_json(&service, "GET", &format!("/api/runs/{id}/files"), None);
    let auspice = files
      .as_array()
      .unwrap()
      .iter()
      .filter(|file| file["path"] == "ancestral.auspice.json")
      .map(|file| file["kind"].clone())
      .collect::<Vec<_>>();
    assert_eq!(
      (
        json!("ok"),
        json!("ok"),
        vec![to_json_value(&OutputSelection::Auspice).unwrap()]
      ),
      (terminal["data"]["status"].clone(), record["status"].clone(), auspice)
    );
  }

  #[test]
  fn test_backend_deferred_run_waits_for_start() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(true)));
    let id = record["id"].as_str().unwrap().to_owned();
    let (_, started) = fetch_json(&service, "POST", &format!("/api/runs/{id}/start"), Some(&json!({})));
    let terminal = wait_for_terminal(&service, &id);
    assert_eq!(
      (json!("created"), json!("running"), json!("ok")),
      (
        record["status"].clone(),
        started["status"].clone(),
        terminal["data"]["status"].clone()
      )
    );
  }

  #[test]
  fn test_backend_run_cannot_start_twice() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(false)));
    let id = record["id"].as_str().unwrap().to_owned();
    let (status, error) = fetch_json(&service, "POST", &format!("/api/runs/{id}/start"), Some(&json!({})));
    wait_for_terminal(&service, &id);
    assert_eq!(
      (
        409,
        json!({ "code": "conflict", "message": format!("run `{id}` has already started"), "causes": [] })
      ),
      (status, error)
    );
  }

  #[test]
  fn test_backend_host_port_saves_a_run_file_and_the_archive_to_chosen_paths() {
    let root = tempdir().unwrap();
    let target = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(false)));
    let id = record["id"].as_str().unwrap().to_owned();
    wait_for_terminal(&service, &id);

    let file = target.path().join("tree.nwk");
    let archive = target.path().join("run.zip");
    let save = format!("/api/runs/{id}/save");
    let saved_file = fetch_json_in(
      &service,
      PortScope::Host,
      "POST",
      &save,
      Some(&json!({ "path": "ancestral.nwk", "destination": file })),
    );
    let saved_archive = fetch_json_in(
      &service,
      PortScope::Host,
      "POST",
      &save,
      Some(&json!({ "destination": archive })),
    );

    assert_eq!(
      (
        (200, Value::Null),
        (200, Value::Null),
        fetch_bytes(&service, &format!("/api/runs/{id}/file?path=ancestral.nwk")),
        archive_contents(fetch_bytes(&service, &format!("/api/runs/{id}/archive")))
      ),
      (
        saved_file,
        saved_archive,
        fs::read(&file).unwrap(),
        archive_contents(fs::read(&archive).unwrap())
      )
    );
  }

  #[test]
  fn test_backend_renderer_port_has_no_save_route() {
    let root = tempdir().unwrap();
    let target = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(true)));
    let id = record["id"].as_str().unwrap().to_owned();
    let (status, error) = fetch_json(
      &service,
      "POST",
      &format!("/api/runs/{id}/save"),
      Some(&json!({ "destination": target.path().join("run.zip") })),
    );
    assert_eq!(
      (404, json!("not_found"), 0),
      (
        status,
        error["code"].clone(),
        fs::read_dir(target.path()).unwrap().count()
      )
    );
  }

  #[test]
  fn test_backend_fetch_answers_a_json_response_through_the_router() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let exchange = fetch(&service, PortScope::Renderer, "GET", "/api/version", vec![], None);
    let body: Value = serde_json::from_slice(&exchange.body).unwrap();
    assert_eq!(
      (
        200,
        Some(o!("application/json")),
        json!(version_info().version),
        Ended::End
      ),
      (
        exchange.status,
        header(&exchange.headers, "content-type"),
        body["version"].clone(),
        exchange.ended
      )
    );
  }

  #[test]
  fn test_backend_fetch_streams_run_events_and_resumes_after_the_last_event_id() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (status, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(false)));
    assert_eq!(200, status);
    let id = record["id"].as_str().unwrap().to_owned();

    let all = fetch(
      &service,
      PortScope::Renderer,
      "GET",
      &format!("/api/runs/{id}/events"),
      vec![],
      None,
    );
    let resumed = fetch(
      &service,
      PortScope::Renderer,
      "GET",
      &format!("/api/runs/{id}/events"),
      vec![(o!("last-event-id"), o!("1"))],
      None,
    );

    let all_ids = event_ids(&all.body);
    assert!(all_ids.len() > 3, "a run sends more than three events: {all_ids:?}");
    assert_eq!(
      (
        Some(o!("text/event-stream")),
        (0..all_ids.len()).collect::<Vec<_>>(),
        all_ids[2..].to_vec(),
        (Ended::End, Ended::End)
      ),
      (
        header(&all.headers, "content-type"),
        all_ids.clone(),
        event_ids(&resumed.body),
        (all.ended, resumed.ended)
      )
    );
  }

  #[test]
  fn test_backend_fetch_abort_ends_an_open_event_stream_without_an_end_reply() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(true)));
    let id = record["id"].as_str().unwrap().to_owned();

    let (abort, replies) = open_fetch(
      &service,
      PortScope::Renderer,
      "GET",
      &format!("/api/runs/{id}/events"),
      vec![],
      None,
    );
    let head = replies.recv_timeout(Duration::from_secs(60)).unwrap();
    abort.abort();
    let after_abort = loop {
      match replies.recv_timeout(Duration::from_secs(60)) {
        Ok(PortReply::End { .. } | PortReply::Reset { .. }) => break Err("the exchange ended with a reply"),
        Ok(_) => {},
        Err(RecvTimeoutError::Disconnected) => break Ok(()),
        Err(RecvTimeoutError::Timeout) => break Err("the exchange kept its reply callback after the abort"),
      }
    };
    assert_eq!(
      (true, Ok(())),
      (matches!(head, PortReply::Head { status: 200, .. }), after_abort)
    );
  }

  #[test]
  fn test_backend_fetch_streams_app_events_of_run_changes() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (abort, replies) = open_fetch(&service, PortScope::Renderer, "GET", "/api/events", vec![], None);
    let head = replies.recv_timeout(Duration::from_secs(60)).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(true)));

    let mut body = vec![];
    while !String::from_utf8_lossy(&body).contains("\n\n") {
      match replies.recv_timeout(Duration::from_secs(60)).unwrap() {
        PortReply::Chunk { data, .. } => body.extend_from_slice(&data),
        _ => panic!("the app event stream ended"),
      }
    }
    abort.abort();
    let text = String::from_utf8(body).unwrap();
    let data = text
      .lines()
      .find_map(|line| line.strip_prefix("data: "))
      .map(|data| serde_json::from_str::<Value>(data).unwrap())
      .unwrap();
    assert_eq!(
      (true, true, json!("run-created"), record["id"].clone()),
      (
        matches!(head, PortReply::Head { status: 200, .. }),
        text.contains("event: run-created"),
        data["kind"].clone(),
        data["run"]["id"].clone()
      )
    );
  }

  #[test]
  fn test_backend_fetch_answers_a_missing_run_with_a_typed_error_body() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (status, body) = fetch_json(&service, "GET", "/api/runs/missing", None);
    assert_eq!((404, json!("not_found")), (status, body["code"].clone()));
  }

  #[test]
  fn test_backend_fetch_answers_a_malformed_request_with_an_invalid_request_response() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let exchange = fetch(
      &service,
      PortScope::Renderer,
      "NOT A METHOD",
      "/api/version",
      vec![],
      None,
    );
    let body: Value = serde_json::from_slice(&exchange.body).unwrap();
    assert_eq!(
      (
        400,
        json!({
          "code": "invalid_request",
          "message": "the request `NOT A METHOD /api/version` is malformed: invalid HTTP method",
          "causes": []
        }),
        Ended::End
      ),
      (exchange.status, body, exchange.ended)
    );
  }

  #[test]
  fn test_backend_fetch_passes_a_binary_body_unchanged() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(true)));
    let id = record["id"].as_str().unwrap().to_owned();
    let bytes = vec![0xff, 0xfe, 0x00, 0xc3, 0x28, 0x80];
    let upload = fetch(
      &service,
      PortScope::Renderer,
      "PUT",
      &format!("/api/runs/{id}/inputs/blob.bin"),
      vec![(o!("content-type"), o!("application/octet-stream"))],
      Some(bytes.clone()),
    );
    let uploaded: Value = serde_json::from_slice(&upload.body).unwrap();
    let stored = fs::read(uploaded["path"].as_str().unwrap()).unwrap();
    assert_eq!((200, bytes), (upload.status, stored));
  }

  #[test]
  fn test_backend_reject_request_answers_a_complete_invalid_request_response() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let replies = service.reject_request(3, o!("the request has no method")).unwrap();
    let summary = replies
      .into_iter()
      .map(|reply| match reply {
        PortReply::Head { seq, status, headers } => json!({
          "head": [seq, status, header(&headers.into_iter().map(|header| (header.name, header.value)).collect::<Vec<_>>(), "content-type")]
        }),
        PortReply::Chunk { seq, data } => json!({ "chunk": [seq, serde_json::from_slice::<Value>(&data).unwrap()] }),
        PortReply::End { seq } => json!({ "end": seq }),
        PortReply::Reset { seq, message } => json!({ "reset": [seq, message] }),
      })
      .collect::<Vec<_>>();
    assert_eq!(
      vec![
        json!({ "head": [3, 400, "application/json"] }),
        json!({ "chunk": [3, { "code": "invalid_request", "message": "the request has no method", "causes": [] }] }),
        json!({ "end": 3 }),
      ],
      summary
    );
  }

  #[tokio::test]
  async fn test_backend_body_failure_after_the_head_resets_the_exchange() {
    let chunks: Vec<Result<Bytes, io::Error>> =
      vec![Ok(Bytes::from_static(b"partial")), Err(io::Error::other("disk gone"))];
    let response = Response::new(Body::from_stream(tokio_stream::iter(chunks)));
    let replies = Mutex::new(vec![]);
    stream_response(5, response, &AtomicBool::new(false), &|reply| {
      replies.lock().push(reply);
      true
    })
    .await;
    let summary = replies
      .into_inner()
      .into_iter()
      .map(|reply| match reply {
        PortReply::Head { status, .. } => json!({ "head": status }),
        PortReply::Chunk { data, .. } => json!({ "chunk": String::from_utf8(data.to_vec()).unwrap() }),
        PortReply::End { .. } => json!("end"),
        PortReply::Reset { message, .. } => json!({ "reset": message }),
      })
      .collect::<Vec<_>>();
    assert_eq!(
      vec![
        json!({ "head": 200 }),
        json!({ "chunk": "partial" }),
        json!({ "reset": "When reading the response body: disk gone" }),
      ],
      summary
    );
  }

  #[test]
  fn test_backend_scope_names_one_of_the_two_routers() {
    assert_eq!(
      (PortScope::Host, PortScope::Renderer),
      (PortScope::parse("host").unwrap(), PortScope::parse("renderer").unwrap())
    );
  }

  #[test]
  fn test_backend_scope_rejects_an_unknown_name() {
    assert_error!(
      PortScope::parse("window"),
      "unknown port scope 'window': Matching variant not found. This is an internal error. Please report it to \
       developers."
    );
  }

  #[test]
  fn test_backend_opens_the_runs_folder_named_in_the_settings() {
    let root = tempdir().unwrap();
    let elsewhere = root.path().join("elsewhere");
    let first = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (status, _) = fetch_json(&first, "PUT", "/api/workspace", Some(&json!({ "path": elsewhere })));
    drop(first);
    let second = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (_, record) = fetch_json(&second, "POST", "/api/runs", Some(&ancestral_request(true)));
    let id = record["id"].as_str().unwrap().to_owned();
    let (_, workspace) = fetch_json(&second, "GET", "/api/workspace", None);
    assert_eq!(
      (
        200,
        json!({ "path": elsewhere, "default_path": root.path().join("runs") }),
        true
      ),
      (status, workspace, elsewhere.join(&id).join("run.json").is_file())
    );
  }

  #[test]
  fn test_backend_keeps_the_ui_preferences_in_the_settings_file() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let ui = json!({ "theme": "light", "sidebar_width": 300 });
    fetch_json(&service, "PUT", "/api/app-settings/ui", Some(&ui));
    let reopened = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (status, settings) = fetch_json(&reopened, "GET", "/api/app-settings", None);
    assert_eq!((200, json!({ "ui": ui })), (status, settings));
  }

  #[test]
  fn test_backend_lists_the_examples_of_the_examples_folder() {
    let root = tempdir().unwrap();
    let dataset = root.path().join("examples").join("zika").join("20");
    fs::create_dir_all(&dataset).unwrap();
    fs::write(dataset.join("tree.nwk"), "(A:1,B:1);\n").unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (status, catalog) = fetch_json(&service, "GET", "/api/datasets", None);
    let names = catalog["datasets"]
      .as_array()
      .unwrap()
      .iter()
      .map(|dataset| dataset["name"].clone())
      .collect::<Vec<_>>();
    assert_eq!((200, vec![json!("zika/20")]), (status, names));
  }

  #[test]
  fn test_backend_takes_a_runs_folder_from_the_environment() {
    let root = tempdir().unwrap();
    let scratch = tempdir().unwrap();
    let env = AppFolderEnv {
      runs: Some(scratch.path().to_path_buf()),
      ..AppFolderEnv::default()
    };
    let service = DesktopService::open(root.path(), &env).unwrap();
    let (_, workspace) = fetch_json(&service, "GET", "/api/workspace", None);
    assert_eq!(
      json!({
        "path": scratch.path(),
        "default_path": root.path().join("runs"),
        "fixed_by": "TREETIME_RUNS_DIR"
      }),
      workspace
    );
  }

  #[test]
  fn test_backend_opens_the_default_runs_folder_when_the_settings_name_an_unopenable_one() {
    let root = tempdir().unwrap();
    let blocked = root.path().join("blocked");
    fs::write(&blocked, "").unwrap();
    let settings = root.path().join("settings.yaml");
    fs::write(&settings, format!("paths:\n  runs: {}\n", blocked.display())).unwrap();
    let service = DesktopService::open(root.path(), &AppFolderEnv::default()).unwrap();
    let (status, workspace) = fetch_json(&service, "GET", "/api/workspace", None);
    let error = format!(
      "the runs folder '{}', named in '{}' cannot be opened: When creating the runs directory '{}': File exists (os \
       error 17). The default runs folder is in use.",
      blocked.display(),
      settings.display(),
      blocked.display()
    );
    assert_eq!(
      (
        200,
        json!({ "path": root.path().join("runs"), "default_path": root.path().join("runs"), "error": error }),
        true
      ),
      (status, workspace, root.path().join("runs").is_dir())
    );
  }

  #[test]
  fn test_backend_fails_when_the_runs_folder_of_the_environment_cannot_be_opened() {
    let root = tempdir().unwrap();
    let blocked = root.path().join("blocked");
    fs::write(&blocked, "").unwrap();
    let env = AppFolderEnv {
      runs: Some(blocked.clone()),
      ..AppFolderEnv::default()
    };
    assert_error!(
      DesktopService::open(root.path(), &env),
      format!(
        "When opening the runs folder '{}', set by TREETIME_RUNS_DIR: When creating the runs directory '{}': File \
         exists (os error 17)",
        blocked.display(),
        blocked.display()
      )
    );
  }

  #[test]
  fn test_backend_names_both_folders_when_the_default_runs_folder_cannot_be_opened_either() {
    let root = tempdir().unwrap();
    let blocked = root.path().join("blocked");
    fs::write(&blocked, "").unwrap();
    fs::write(root.path().join("runs"), "").unwrap();
    let settings = root.path().join("settings.yaml");
    fs::write(&settings, format!("paths:\n  runs: {}\n", blocked.display())).unwrap();
    assert_error!(
      DesktopService::open(root.path(), &AppFolderEnv::default()),
      format!(
        "the runs folder '{blocked}', named in '{settings}' cannot be opened: When creating the runs directory \
         '{blocked}': File exists (os error 17); the default runs folder '{default}' cannot be opened either: When \
         creating the runs directory '{default}': File exists (os error 17)",
        blocked = blocked.display(),
        settings = settings.display(),
        default = root.path().join("runs").display()
      )
    );
  }

  #[test]
  fn test_backend_names_the_settings_file_that_it_cannot_read() {
    let root = tempdir().unwrap();
    let file = root.path().join("settings.yaml");
    fs::write(&file, "paths: [1]\n").unwrap();
    let expected = "Unexpected: unexpected sequence, expected AppPathSettings at line 1 column 8 (path: paths)";
    assert_error!(
      DesktopService::open(root.path(), &AppFolderEnv::default()),
      format!("When reading the settings file '{}': {expected}", file.display())
    );
  }

  mod helpers {
    use crate::backend::DesktopService;
    use crate::backend::PortScope;
    use crate::port::{PortHeader, PortReply, PortRequest};
    use napi::bindgen_prelude::Uint8Array;
    use serde_json::{Value, json};
    use std::io::{Cursor, Read};
    use std::path::Path;
    use std::sync::mpsc;
    use std::time::Duration;
    use tokio::task::AbortHandle;
    use treetime_utils::o;
    use zip::ZipArchive;

    #[derive(Debug, PartialEq, Eq)]
    pub(super) enum Ended {
      End,
      Reset(String),
    }

    pub(super) struct Exchange {
      pub status: u16,
      pub headers: Vec<(String, String)>,
      pub body: Vec<u8>,
      pub ended: Ended,
    }

    pub(super) fn open_fetch(
      service: &DesktopService,
      scope: PortScope,
      method: &str,
      url: &str,
      headers: Vec<(String, String)>,
      body: Option<Vec<u8>>,
    ) -> (AbortHandle, mpsc::Receiver<PortReply>) {
      let (send, receive) = mpsc::sync_channel(1024);
      let request = PortRequest {
        seq: 7,
        method: method.to_owned(),
        url: url.to_owned(),
        headers: headers
          .into_iter()
          .map(|(name, value)| PortHeader { name, value })
          .collect(),
        body: body.map(Uint8Array::new),
      };
      let abort = service.fetch(request, scope, move |reply| send.send(reply).is_ok());
      (abort, receive)
    }

    pub(super) fn fetch(
      service: &DesktopService,
      scope: PortScope,
      method: &str,
      url: &str,
      headers: Vec<(String, String)>,
      body: Option<Vec<u8>>,
    ) -> Exchange {
      let (_abort, replies) = open_fetch(service, scope, method, url, headers, body);
      let mut exchange = Exchange {
        status: 0,
        headers: vec![],
        body: vec![],
        ended: Ended::End,
      };
      loop {
        let reply = replies.recv_timeout(Duration::from_secs(120)).unwrap();
        assert_eq!(7, reply.seq());
        match reply {
          PortReply::Head { status, headers, .. } => {
            exchange.status = status;
            exchange.headers = headers.into_iter().map(|header| (header.name, header.value)).collect();
          },
          PortReply::Chunk { data, .. } => exchange.body.extend_from_slice(&data),
          PortReply::End { .. } => return exchange,
          PortReply::Reset { message, .. } => {
            exchange.ended = Ended::Reset(message);
            return exchange;
          },
        }
      }
    }

    pub(super) fn fetch_json(service: &DesktopService, method: &str, url: &str, body: Option<&Value>) -> (u16, Value) {
      fetch_json_in(service, PortScope::Renderer, method, url, body)
    }

    pub(super) fn fetch_json_in(
      service: &DesktopService,
      scope: PortScope,
      method: &str,
      url: &str,
      body: Option<&Value>,
    ) -> (u16, Value) {
      let exchange = fetch(
        service,
        scope,
        method,
        url,
        vec![(o!("content-type"), o!("application/json"))],
        body.map(|body| body.to_string().into_bytes()),
      );
      assert_eq!(Ended::End, exchange.ended);
      let body = if exchange.body.is_empty() {
        Value::Null
      } else {
        serde_json::from_slice(&exchange.body).unwrap()
      };
      (exchange.status, body)
    }

    pub(super) fn archive_contents(bytes: Vec<u8>) -> Vec<(String, Vec<u8>)> {
      let mut archive = ZipArchive::new(Cursor::new(bytes)).unwrap();
      (0..archive.len())
        .map(|index| {
          let mut file = archive.by_index(index).unwrap();
          let mut content = vec![];
          file.read_to_end(&mut content).unwrap();
          (file.name().to_owned(), content)
        })
        .collect()
    }

    pub(super) fn header(headers: &[(String, String)], name: &str) -> Option<String> {
      headers
        .iter()
        .find(|(header, _)| header == name)
        .map(|(_, value)| value.clone())
    }

    pub(super) fn event_ids(body: &[u8]) -> Vec<usize> {
      String::from_utf8(body.to_vec())
        .unwrap()
        .lines()
        .filter_map(|line| line.strip_prefix("id: "))
        .map(|id| id.parse().unwrap())
        .collect()
    }

    pub(super) fn fetch_bytes(service: &DesktopService, url: &str) -> Vec<u8> {
      let exchange = fetch(service, PortScope::Renderer, "GET", url, vec![], None);
      assert_eq!((200, Ended::End), (exchange.status, exchange.ended));
      exchange.body
    }

    pub(super) fn wait_for_terminal(service: &DesktopService, id: &str) -> Value {
      let exchange = fetch(
        service,
        PortScope::Renderer,
        "GET",
        &format!("/api/runs/{id}/events"),
        vec![],
        None,
      );
      String::from_utf8(exchange.body)
        .unwrap()
        .lines()
        .filter_map(|line| line.strip_prefix("data: "))
        .map(|data| serde_json::from_str::<Value>(data).unwrap())
        .find(|event| event["type"] == "terminal")
        .unwrap()
    }

    pub(super) fn ancestral_request(defer_start: bool) -> Value {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({
        "command": "ancestral",
        "config": {
          "tree": zika.join("tree.nwk"),
          "alignment": [zika.join("aln.fasta.xz")],
        },
        "defer_start": defer_start,
      })
    }
  }
}
