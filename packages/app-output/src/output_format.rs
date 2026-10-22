use std::path::Path;
use treetime_utils::io::compression::remove_compression_ext;
use treetime_utils::io::file::is_path_stdout;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChartFormat {
  Svg,
  Png,
  Pdf,
}

impl ChartFormat {
  pub fn from_path(path: &Path) -> Option<Self> {
    if is_path_stdout(path) {
      return Some(Self::Svg);
    }
    match lowercase_extension(path)?.as_str() {
      "svg" => Some(Self::Svg),
      "png" => Some(Self::Png),
      "pdf" => Some(Self::Pdf),
      _ => None,
    }
  }
}

fn lowercase_extension(path: &Path) -> Option<String> {
  remove_compression_ext(path)
    .extension()
    .and_then(|ext| ext.to_str())
    .map(str::to_ascii_lowercase)
}
