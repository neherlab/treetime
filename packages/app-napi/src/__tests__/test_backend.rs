#[cfg(test)]
mod tests {
  use crate::backend::DesktopService;
  use crate::port::PortReply;
  use app_commands::bridge::operations::OperationRequest;
  use app_commands::runs::record::{RunRecord, RunStatus};
  use app_output::output_plan::OutputSelection;
  use helpers::{Ended, ancestral_request, call, event_ids, fetch, fetch_json, header, open_fetch, wait_for_terminal};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::fs;
  use std::sync::mpsc::RecvTimeoutError;
  use std::time::Duration;
  use tempfile::tempdir;
  use treetime_schema::version_info;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_backend_created_run_starts_at_once_and_writes_its_outputs() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let record: RunRecord = call(
      &service,
      &json!({ "operation": "create-run", "args": { "request": ancestral_request(false) } }),
    );
    let terminal = wait_for_terminal(&service, &record);

    let record: RunRecord = call(
      &service,
      &json!({ "operation": "get-run", "args": { "id": record.id } }),
    );
    let files: Vec<Value> = call(
      &service,
      &json!({ "operation": "run-files", "args": { "id": record.id } }),
    );
    let auspice = files
      .iter()
      .filter(|file| file["path"] == "ancestral.auspice.json")
      .map(|file| file["kind"].clone())
      .collect::<Vec<_>>();
    assert_eq!(
      (
        json!("ok"),
        RunStatus::Ok,
        vec![serde_json::to_value(OutputSelection::Auspice).unwrap()]
      ),
      (terminal["data"]["status"].clone(), record.status, auspice)
    );
  }

  #[test]
  fn test_backend_deferred_run_waits_for_start() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let record: RunRecord = call(
      &service,
      &json!({ "operation": "create-run", "args": { "request": ancestral_request(true) } }),
    );
    let started: RunRecord = call(
      &service,
      &json!({ "operation": "start-run", "args": { "id": record.id, "request": {} } }),
    );
    let terminal = wait_for_terminal(&service, &record);
    assert_eq!(
      (RunStatus::Created, RunStatus::Running, json!("ok")),
      (record.status, started.status, terminal["data"]["status"].clone())
    );
  }

  #[test]
  fn test_backend_run_cannot_start_twice() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let record: RunRecord = call(
      &service,
      &json!({ "operation": "create-run", "args": { "request": ancestral_request(false) } }),
    );
    let request: OperationRequest =
      serde_json::from_value(json!({ "operation": "start-run", "args": { "id": record.id, "request": {} } })).unwrap();
    assert_error!(
      request.handle(service.app()),
      format!("run `{}` has already started", record.id.as_str())
    );
    wait_for_terminal(&service, &record);
  }

  #[test]
  fn test_backend_saves_a_run_file_and_the_archive_to_chosen_paths() {
    let root = tempdir().unwrap();
    let target = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let record: RunRecord = call(
      &service,
      &json!({ "operation": "create-run", "args": { "request": ancestral_request(false) } }),
    );
    wait_for_terminal(&service, &record);

    let file = target.path().join("tree.nwk");
    let archive = target.path().join("run.zip");
    service.save_run_file(&record.id, "ancestral.nwk", &file).unwrap();
    service.save_run_archive(&record.id, &archive).unwrap();

    let source = service.app().runs().file_path(&record.id, "ancestral.nwk").unwrap();
    assert_eq!(
      (
        fs::read(&source).unwrap(),
        service.app().runs().zip(&record.id).unwrap()
      ),
      (fs::read(&file).unwrap(), fs::read(&archive).unwrap())
    );
  }

  #[test]
  fn test_backend_refuses_to_save_a_file_outside_the_run_folder() {
    let root = tempdir().unwrap();
    let target = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let record: RunRecord = call(
      &service,
      &json!({ "operation": "create-run", "args": { "request": ancestral_request(true) } }),
    );
    assert_error!(
      service.save_run_file(&record.id, "../run.json", &target.path().join("run.json")),
      "file path `../run.json` must name a file inside the run's output folder"
    );
    assert_eq!(0, fs::read_dir(target.path()).unwrap().count());
  }

  #[test]
  fn test_backend_fetch_answers_a_json_response_through_the_router() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let exchange = fetch(&service, "GET", "/api/version", vec![], None);
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
    let service = DesktopService::open(root.path()).unwrap();
    let (status, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(false)));
    assert_eq!(200, status);
    let id = record["id"].as_str().unwrap().to_owned();

    let all = fetch(&service, "GET", &format!("/api/runs/{id}/events"), vec![], None);
    let resumed = fetch(
      &service,
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
    let service = DesktopService::open(root.path()).unwrap();
    let (_, record) = fetch_json(&service, "POST", "/api/runs", Some(&ancestral_request(true)));
    let id = record["id"].as_str().unwrap().to_owned();

    let (abort, replies) = open_fetch(&service, "GET", &format!("/api/runs/{id}/events"), vec![], None);
    let head = replies.recv_timeout(Duration::from_secs(60)).unwrap();
    abort.abort();
    let after_abort = loop {
      match replies.recv_timeout(Duration::from_secs(60)) {
        Ok(PortReply::End { .. } | PortReply::Error { .. }) => break Err("the exchange ended with a reply"),
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
    let service = DesktopService::open(root.path()).unwrap();
    let (abort, replies) = open_fetch(&service, "GET", "/api/events", vec![], None);
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
    let service = DesktopService::open(root.path()).unwrap();
    let (status, body) = fetch_json(&service, "GET", "/api/runs/missing", None);
    assert_eq!((404, json!("not_found")), (status, body["code"].clone()));
  }

  #[test]
  fn test_backend_fetch_reports_a_malformed_request_as_an_invalid_request_error() {
    let root = tempdir().unwrap();
    let service = DesktopService::open(root.path()).unwrap();
    let exchange = fetch(&service, "NOT A METHOD", "/api/version", vec![], None);
    let Ended::Error(error) = exchange.ended else {
      panic!("the exchange did not end with an error reply");
    };
    assert_eq!(
      (
        0,
        o!("invalid_request"),
        o!("the request `NOT A METHOD /api/version` is malformed: invalid HTTP method")
      ),
      (exchange.status, error.code, error.message)
    );
  }

  mod helpers {
    use crate::backend::DesktopService;
    use crate::port::{PortError, PortHeader, PortReply, PortRequest};
    use app_commands::bridge::operations::OperationRequest;
    use app_commands::runs::record::RunRecord;
    use serde::de::DeserializeOwned;
    use serde_json::{Value, json};
    use std::path::Path;
    use std::sync::mpsc;
    use std::time::Duration;
    use tokio::task::AbortHandle;
    use treetime_utils::o;

    #[derive(Debug, PartialEq, Eq)]
    pub(super) enum Ended {
      End,
      Error(PortError),
    }

    pub(super) struct Exchange {
      pub status: u16,
      pub headers: Vec<(String, String)>,
      pub body: Vec<u8>,
      pub ended: Ended,
    }

    pub(super) fn open_fetch(
      service: &DesktopService,
      method: &str,
      url: &str,
      headers: Vec<(String, String)>,
      body: Option<String>,
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
        body,
      };
      let abort = service.fetch(request, move |reply| send.send(reply).is_ok());
      (abort, receive)
    }

    pub(super) fn fetch(
      service: &DesktopService,
      method: &str,
      url: &str,
      headers: Vec<(String, String)>,
      body: Option<String>,
    ) -> Exchange {
      let (_abort, replies) = open_fetch(service, method, url, headers, body);
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
          PortReply::Error { error, .. } => {
            exchange.ended = Ended::Error(error);
            return exchange;
          },
        }
      }
    }

    pub(super) fn fetch_json(service: &DesktopService, method: &str, url: &str, body: Option<&Value>) -> (u16, Value) {
      let exchange = fetch(
        service,
        method,
        url,
        vec![(o!("content-type"), o!("application/json"))],
        body.map(Value::to_string),
      );
      assert_eq!(Ended::End, exchange.ended);
      (exchange.status, serde_json::from_slice(&exchange.body).unwrap())
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

    pub(super) fn call<T: DeserializeOwned>(service: &DesktopService, request: &Value) -> T {
      let request: OperationRequest = serde_json::from_value(request.clone()).unwrap();
      serde_json::from_str(&request.handle(service.app()).unwrap()).unwrap()
    }

    pub(super) fn wait_for_terminal(service: &DesktopService, record: &RunRecord) -> Value {
      let (send, receive) = mpsc::sync_channel(1024);
      service
        .app()
        .runs()
        .subscribe(
          &record.id,
          0,
          Box::new(move |event| send.send(serde_json::to_value(event).unwrap()).is_ok()),
        )
        .unwrap();
      receive.iter().find(|event| event["type"] == "terminal").unwrap()
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
