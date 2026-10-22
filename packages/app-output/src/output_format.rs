use crate::output_plan::OutputSelection;
use derive_more::{Display, Error};
use std::path::{Path, PathBuf};
use treetime_io::csv::TableFormat;
use treetime_utils::io::compression::remove_compression_ext;
use treetime_utils::io::file::is_path_stdout;

pub fn table_format(path: &Path) -> Option<TableFormat> {
  if is_path_stdout(path) {
    return Some(TableFormat::Tsv);
  }
  TableFormat::from_path(path)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FileFormat {
  Table,
  Chart,
  Fixed,
}

impl FileFormat {
  pub fn check(self, selection: OutputSelection, path: &Path) -> Result<(), UnsupportedExtension> {
    let known = match self {
      Self::Table => table_format(path).is_some(),
      Self::Chart => ChartFormat::from_path(path).is_some(),
      Self::Fixed => true,
    };
    if known {
      return Ok(());
    }
    Err(UnsupportedExtension {
      selection,
      path: path.to_path_buf(),
      extensions: self.extensions(),
    })
  }

  const fn extensions(self) -> &'static str {
    match self {
      Self::Table => ".tsv (tab-separated) or .csv (comma-separated)",
      Self::Chart => ".svg, .png or .pdf",
      Self::Fixed => "",
    }
  }
}

#[derive(Clone, Debug, Eq, PartialEq, Display, Error)]
#[display(
  "The file name '{}' of {} must end in {extensions}, optionally followed by a compression suffix such as .gz, or be '-' for standard output",
  path.display(),
  selection.flag_name()
)]
pub struct UnsupportedExtension {
  pub selection: OutputSelection,
  pub path: PathBuf,
  pub extensions: &'static str,
}

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
