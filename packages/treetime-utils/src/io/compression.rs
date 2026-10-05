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

// TODO: keep an eye on alternative crates that don't rely on C, and replace when they are stable enough.
// TODO: keep an eye on efforts of bringing (pieces of) libc to `wasm32-unknown-unknown`, and enable `bzip2`, `xz2`
// and `zstd` crates for wasm builds when it's compatible enough: https://github.com/rustwasm/team/issues/291
// Crates `bzip2`, `xz2` and `zstd` depend on corresponding C libraries and require libc in order to build.
// libc is not present for `wasm32-unknown-unknown` target, so we disable these crates.
#[cfg(not(target_arch = "wasm32"))]
use bzip2::read::MultiBzDecoder;
#[cfg(not(target_arch = "wasm32"))]
use bzip2::write::BzEncoder;

#[cfg(not(target_arch = "wasm32"))]
use xz2::read::XzDecoder;
#[cfg(not(target_arch = "wasm32"))]
use xz2::write::XzEncoder;

pub const COMPRESSION_EXTENSIONS: [&str; 4] = ["bz2", "xz", "zst", "gz"];

pub fn remove_compression_ext(filepath: impl AsRef<Path>) -> PathBuf {
  let path = filepath.as_ref();

  path
    .extension()
    .and_then(|ext| ext.to_str())
    .filter(|ext| COMPRESSION_EXTENSIONS.iter().any(|&e| e.eq_ignore_ascii_case(ext)))
    .map_or_else(|| path.to_path_buf(), |_| path.with_extension(""))
}

pub struct Decompressor<'r> {
  decompressor: Box<dyn Read + 'r>,
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

    Ok(Self { decompressor })
  }

  pub fn from_path<R: 'r + Read>(reader: R, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let (compression_type, _) = guess_compression_from_filepath(filepath);
    Self::new(reader, &compression_type)
  }
}

impl Read for Decompressor<'_> {
  fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
    self
      .decompressor
      .read(buf)
      .map_err(|error| io::Error::new(error.kind(), format!("While decompressing file: {error}")))
  }
}

trait Encoder: Write + Send {
  fn finish_encoding(self: Box<Self>) -> io::Result<()>;
}

#[cfg(not(target_arch = "wasm32"))]
impl<W: Write + Send> Encoder for BzEncoder<W> {
  fn finish_encoding(self: Box<Self>) -> io::Result<()> {
    self.finish()?.flush()
  }
}

#[cfg(not(target_arch = "wasm32"))]
impl<W: Write + Send> Encoder for XzEncoder<W> {
  fn finish_encoding(self: Box<Self>) -> io::Result<()> {
    self.finish()?.flush()
  }
}

#[cfg(not(target_arch = "wasm32"))]
impl<W: Write + Send> Encoder for zstd::Encoder<'_, W> {
  fn finish_encoding(self: Box<Self>) -> io::Result<()> {
    self.finish()?.flush()
  }
}

impl<W: Write + Send> Encoder for GzEncoder<W> {
  fn finish_encoding(self: Box<Self>) -> io::Result<()> {
    self.finish()?.flush()
  }
}

struct Uncompressed<W>(W);

impl<W: Write> Write for Uncompressed<W> {
  fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
    self.0.write(buf)
  }

  fn flush(&mut self) -> io::Result<()> {
    self.0.flush()
  }
}

impl<W: Write + Send> Encoder for Uncompressed<W> {
  fn finish_encoding(mut self: Box<Self>) -> io::Result<()> {
    self.0.flush()
  }
}

pub struct Compressor<'w> {
  encoder: Option<Box<dyn Encoder + 'w>>,
  compression_type: CompressionType,
  filepath: Option<String>,
}

impl<'w> Compressor<'w> {
  pub fn new<W: 'w + Write + Send>(writer: W, compression_type: &CompressionType) -> Result<Self, Report> {
    let encoder: Box<dyn Encoder + 'w> = match compression_type {
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Bzip2 => Box::new(BzEncoder::new(writer, bzip2::Compression::new(get_comp_level("BZ2")?))),
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Xz => Box::new(XzEncoder::new(writer, get_comp_level("XZ")?)),
      #[cfg(not(target_arch = "wasm32"))]
      CompressionType::Zstandard => {
        Box::new(zstd::Encoder::new(writer, get_comp_level("ZST")?).wrap_err("When creating the Zstandard encoder")?)
      },
      CompressionType::Gzip => Box::new(GzEncoder::new(writer, flate2::Compression::new(get_comp_level("GZ")?))),
      CompressionType::None => Box::new(Uncompressed(writer)),
    };

    Ok(Self {
      encoder: Some(encoder),
      compression_type: compression_type.clone(),
      filepath: None,
    })
  }

  pub fn from_path<W: 'w + Write + Send>(writer: W, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let (compression_type, _) = guess_compression_from_filepath(filepath);
    let mut compressor = Self::new(writer, &compression_type)?;
    compressor.filepath = Some(filepath.display().to_string());
    Ok(compressor)
  }

  pub fn finish(mut self) -> Result<(), Report> {
    match self.encoder.take() {
      Some(encoder) => self.with_context(encoder.finish_encoding(), "While finishing compressed file"),
      None => Ok(()),
    }
  }

  fn encoder(&mut self) -> io::Result<&mut Box<dyn Encoder + 'w>> {
    self
      .encoder
      .as_mut()
      .ok_or_else(|| io::Error::other("Compressed file was already finished"))
  }

  fn with_context<T>(&self, result: io::Result<T>, message: &'static str) -> Result<T, Report> {
    result
      .wrap_err(message)
      .with_section(|| {
        self
          .filepath
          .clone()
          .unwrap_or_else(|| "None".to_owned())
          .header("Filename")
      })
      .with_section(|| self.compression_type.clone().header("Compressor"))
  }
}

impl Write for Compressor<'_> {
  fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
    let result = self.encoder().and_then(|encoder| encoder.write(buf));
    self
      .with_context(result, "While compressing file")
      .map_err(|report| io::Error::other(report_to_string(&report)))
  }

  fn flush(&mut self) -> io::Result<()> {
    let result = self.encoder().and_then(|encoder| encoder.flush());
    self
      .with_context(result, "While flushing compressed file")
      .map_err(|report| io::Error::other(report_to_string(&report)))
  }
}

impl Drop for Compressor<'_> {
  fn drop(&mut self) {
    if let Some(encoder) = self.encoder.take()
      && let Err(e) = self.with_context(encoder.finish_encoding(), "While finishing compressed file")
    {
      error!("Failed to finish compressed file on drop: {}", report_to_string(&e));
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
