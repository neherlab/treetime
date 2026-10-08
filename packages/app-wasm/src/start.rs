use js_sys::Function;
use log::{Level, error};
use std::cell::RefCell;
use std::panic::{self, PanicHookInfo};
use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

thread_local! {
  static PANIC_SINK: RefCell<Option<Function>> = const { RefCell::new(None) };
}

#[wasm_bindgen(start)]
pub fn start() -> Result<(), JsError> {
  panic::set_hook(Box::new(on_panic));
  console_log::init_with_level(Level::Warn)
    .map_err(|err| JsError::new(&format!("When installing the console logger of the back end: {err}")))
}

#[wasm_bindgen(js_name = setPanicSink)]
pub fn set_panic_sink(#[wasm_bindgen(unchecked_param_type = "(text: string) => void")] sink: &Function) {
  PANIC_SINK.with_borrow_mut(|current| *current = Some(sink.clone()));
}

fn on_panic(info: &PanicHookInfo<'_>) {
  console_error_panic_hook::hook(info);
  let text = JsValue::from(info.to_string());
  PANIC_SINK.with_borrow(|sink| {
    if let Some(Err(err)) = sink.as_ref().map(|sink| sink.call1(&JsValue::UNDEFINED, &text)) {
      error!("The panic sink of the back end failed: {err:?}");
    }
  });
}
