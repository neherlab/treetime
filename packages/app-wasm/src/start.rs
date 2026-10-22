use log::Level;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

#[wasm_bindgen(start)]
pub fn start() -> Result<(), JsError> {
  console_log::init_with_level(Level::Warn)
    .map_err(|err| JsError::new(&format!("When installing the console logger of the module: {err}")))
}
