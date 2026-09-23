use crate::env::env_var_optional;
use crate::error::report_to_string;
use crate::io::fs::extension;
use color_eyre::{Help, SectionExt};
use eyre::{Report, WrapErr};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use log::{debug, error};
use num::Integer;
use num_traits::NumCast;
use std::error::Error;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

#[cfg(not(target_arch = "wasm32"))]
use bzip2::read::MultiBzDecoder;
#[cfg(not(target_arch = "wasm32"))]
use bzip2::write::BzEncoder;

#[cfg(not(target_arch = "wasm32"))]
use xz2::read::XzDecoder;
#[cfg(not(target_arch = "wasm32"))]
use xz2::write::XzEncoder;

pub fn remove_compression_ext(filepath: impl AsRef<Path>) -> PathBuf {
  let compressed_exts = ["bz2", "xz", "zst", "gz"];
  let path = filepath.as_ref();

  path
    .extension()
    .and_then(|ext| ext.to_str())
    .filter(|ext| compressed_exts.iter().any(|&e| e.eq_ignore_ascii_case(ext)))
    .map_or_else(|| path.to_path_buf(), |_| path.with_extension(""))
}

pub struct Decompressor<'r> {
  decompressor: Box<dyn Read + 'r>,
  compression_type: CompressionType,
  filepath: Option<String>,
}

impl<'r> Decompressor<'r> {
  pub fn new<R: 'r + Read>(reader: R, compression_type: &CompressionType) -> Result<Self, Report> {
    let decompressor: Box<dyn Read> = match compression_type {
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Bzip2 => Box::new(MultiBzDecoder::new(reader)),
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Xz => Box::new(XzDecoder::new_multi_decoder(reader)),
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Zstandard => {
        Box::new(zstd::Decoder::new(reader).wrap_err("When creating the Zstandard decoder")?)
      },
      CompressionType::Gzip => Box::new(MultiGzDecoder::new(reader)),
      CompressionType::None => Box::new(reader),
    };

    Ok(Self {
      decompressor,
      compression_type: compression_type.clone(),
      filepath: None,
    })
  }

  pub fn from_str_and_path(content: &'r str, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let reader = content.as_bytes();
    let (compression_type, ext) = guess_compression_from_filepath(filepath);
    Self::new(reader, &compression_type)
  }

  pub fn from_path<R: 'r + Read>(reader: R, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let (compression_type, ext) = guess_compression_from_filepath(filepath);
    Self::new(reader, &compression_type)
  }
}

impl Read for Decompressor<'_> {
  fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
    self
      .decompressor
      .read(buf)
      .wrap_err_with(|| "While decompressing file")
      .with_section(|| {
        self
          .filepath
          .clone()
          .unwrap_or_else(|| "None".to_owned())
          .header("Filename")
      })
      .with_section(|| self.compression_type.clone().header("Decompressor"))
      .map_err(|report| io::Error::other(report_to_string(&report)))
  }
}

pub struct Compressor<'w> {
  compressor: Box<dyn Write + Send + 'w>,
  compression_type: CompressionType,
  filepath: Option<String>,
}

impl<'w> Compressor<'w> {
  pub fn new<W: 'w + Write + Send>(writer: W, compression_type: &CompressionType) -> Result<Self, Report> {
    let compressor: Box<dyn Write + Send + 'w> = match compression_type {
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Bzip2 => Box::new(BzEncoder::new(writer, bzip2::Compression::new(get_comp_level("BZ2")?))),
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Xz => Box::new(XzEncoder::new(writer, get_comp_level("XZ")?)),
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Zstandard => Box::new(
        zstd::Encoder::new(writer, get_comp_level("ZST")?)
          .wrap_err("When creating the Zstandard encoder")?
          .auto_finish(),
      ),
      CompressionType::Gzip => Box::new(GzEncoder::new(writer, flate2::Compression::new(get_comp_level("GZ")?))),
      CompressionType::None => Box::new(writer),
    };

    Ok(Self {
      compressor,
      compression_type: compression_type.clone(),
      filepath: None,
    })
  }

  pub fn from_path<W: 'w + Write + Send>(writer: W, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let (compression_type, ext) = guess_compression_from_filepath(filepath);
    Self::new(writer, &compression_type)
  }
}

impl Write for Compressor<'_> {
  fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
    self
      .compressor
      .write(buf)
      .wrap_err_with(|| "While compressing file")
      .with_section(|| {
        self
          .filepath
          .clone()
          .unwrap_or_else(|| "None".to_owned())
          .header("Filename")
      })
      .with_section(|| self.compression_type.clone().header("Compressor"))
      .map_err(|report| io::Error::other(report_to_string(&report)))
  }

  fn flush(&mut self) -> io::Result<()> {
    self
      .compressor
      .flush()
      .wrap_err_with(|| "While flushing compressed file")
      .with_section(|| {
        self
          .filepath
          .clone()
          .unwrap_or_else(|| "None".to_owned())
          .header("Filename")
      })
      .with_section(|| self.compression_type.clone().header("Compressor"))
      .map_err(|report| io::Error::other(report_to_string(&report)))
  }
}

impl Drop for Compressor<'_> {
  fn drop(&mut self) {
    if let Err(e) = self.flush() {
      error!(
        "Failed to flush compressor on drop: {e}{}",
        self
          .filepath
          .as_ref()
          .map_or(String::new(), |p| format!(" (file: {p})"))
      );
    }
  }
}

pub fn guess_compression_from_filepath(filepath: impl AsRef<Path>) -> (CompressionType, String) {
  let filepath = filepath.as_ref();

  match extension(filepath).map(|ext| ext.to_lowercase()) {
    None => (CompressionType::None, "".to_owned()),
    Some(ext) => {
      let compression_type: CompressionType = match ext.as_str() {
        #[cfg(not(target_arch = "wasm32"))]
        "bz2" => CompressionType::Bzip2,
        #[cfg(not(target_arch = "wasm32"))]
        "xz" => CompressionType::Xz,
        #[cfg(not(target_arch = "wasm32"))]
        "zst" => CompressionType::Zstandard,
        "gz" => CompressionType::Gzip,
        _ => CompressionType::None,
      };

      debug!(
        "When processing '{}': detected file extension '{ext}'. \
        Will be using compression algorithm: '{compression_type}'",
        filepath.display()
      );

      (compression_type, ext)
    },
  }
}

#[derive(strum_macros::Display, Clone)]
pub enum CompressionType {
  #[cfg(not(target_arch = "wasm32"))]
  Bzip2,
  #[cfg(not(target_arch = "wasm32"))]
  Xz,
  #[cfg(not(target_arch = "wasm32"))]
  Zstandard,

  Gzip,
  None,
}

#[allow(
  clippy::unwrap_used,
  reason = "default compression level 2 is representable in every NumCast integer target"
)]
fn get_comp_level<I>(ext: &str) -> Result<I, Report>
where
  I: FromStr + Integer + NumCast,
  I::Err: Error + Send + Sync + 'static,
{
  let var_name = format!("{}_COMPRESSION", ext.to_uppercase());
  match env_var_optional(&var_name)? {
    Some(value) => value
      .parse::<I>()
      .wrap_err_with(|| format!("When parsing compression level '{value}' from environment variable '{var_name}'")),
    None => Ok(NumCast::from(2).unwrap()),
  }
}
