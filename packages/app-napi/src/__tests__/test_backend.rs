#[cfg(test)]
mod tests {
  use crate::backend::DesktopService;
  use app_commands::bridge::operations::OperationRequest;
  use app_commands::runs::record::{RunRecord, RunStatus};
  use app_output::output_plan::OutputSelection;
  use helpers::{ancestral_request, call, wait_for_terminal};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

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

  mod helpers {
    use crate::backend::DesktopService;
    use app_commands::bridge::operations::OperationRequest;
    use app_commands::runs::record::RunRecord;
    use serde::de::DeserializeOwned;
    use serde_json::{Value, json};
    use std::path::Path;
    use std::sync::mpsc;

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
