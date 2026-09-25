use app_commands::runs::events::{RunEvent, Subscriber};
use log::error;
use parking_lot::Mutex;
use std::sync::Arc;

pub(crate) struct EventForwarder<S> {
  sink: Arc<Mutex<Option<S>>>,
}

impl<S: Fn(String) -> bool + Send + 'static> EventForwarder<S> {
  pub(crate) fn new(sink: S) -> Self {
    Self {
      sink: Arc::new(Mutex::new(Some(sink))),
    }
  }

  pub(crate) fn subscriber(&self) -> Subscriber {
    let sink = Arc::clone(&self.sink);
    Box::new(move |event: &RunEvent| {
      let mut sink = sink.lock();
      let Some(send) = sink.as_ref() else {
        return false;
      };
      let delivered = match serde_json::to_string(event) {
        Ok(json) => send(json),
        Err(err) => {
          error!("When serializing a run event: {err}");
          false
        },
      };
      if !delivered {
        sink.take();
      }
      delivered
    })
  }

  pub(crate) fn close(&self) {
    self.sink.lock().take();
  }
}
