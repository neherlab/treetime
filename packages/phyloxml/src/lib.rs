pub mod types;

pub use crate::types::*;
use quick_xml::DeError;
use quick_xml::de::from_reader;

pub fn phyloxml_read(reader: impl std::io::Read) -> Result<Phyloxml, DeError> {
  let reader = std::io::BufReader::new(reader);
  from_reader(reader)
}

pub fn phyloxml_write(writer: impl std::io::Write, phyloxml: &Phyloxml) -> std::io::Result<()> {
  details::to_writer_pretty(writer, "phyloxml", phyloxml)
}

mod details {
  use quick_xml::se::to_string_with_root;
  use serde::Serialize;
  use std::io::Cursor;

  pub fn to_writer_pretty<W, T>(writer: W, root_tag: &str, data: &T) -> std::io::Result<()>
  where
    W: std::io::Write,
    T: Serialize,
  {
    let s = to_string_with_root(root_tag, data).map_err(std::io::Error::other)?;

    let reader = xml::ParserConfig::new()
      .trim_whitespace(true)
      .ignore_comments(false)
      .create_reader(Cursor::new(s));

    let mut writer = xml::EmitterConfig::new()
      .perform_indent(true)
      .normalize_empty_elements(false)
      .autopad_comments(false)
      .create_writer(writer);

    for event in reader {
      if let Some(event) = event.map_err(to_io)?.as_writer_event() {
        writer.write(event).map_err(to_io)?;
      }
    }
    Ok(())
  }

  fn to_io<E>(e: E) -> std::io::Error
  where
    E: Into<Box<dyn std::error::Error + Send + Sync>>,
  {
    std::io::Error::other(e)
  }
}
