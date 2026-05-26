use std::path::PathBuf;

#[derive(Clone, Debug)]
pub struct ServerConfig {
  pub data_dir: PathBuf,
  pub out_dir: PathBuf,
}
