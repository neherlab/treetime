use app_backend::runs::files::FileChunks;
use treetime_utils::error::report_to_string;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

#[wasm_bindgen]
pub struct ChunkReader {
  chunks: FileChunks,
}

#[wasm_bindgen]
impl ChunkReader {
  pub fn chunk(&mut self) -> Result<Option<Vec<u8>>, JsError> {
    self
      .chunks
      .next()
      .transpose()
      .map_err(|report| JsError::new(&report_to_string(&report)))
  }
}

impl From<FileChunks> for ChunkReader {
  fn from(chunks: FileChunks) -> Self {
    Self { chunks }
  }
}
