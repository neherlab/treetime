use crate::io::compression::{CompressionType, Compressor, Decompressor};
use crate::io::fs::ensure_dir;
use eyre::{Report, WrapErr};
use log::info;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write, stdin, stdout};
use std::path::{Path, PathBuf};

pub const DEFAULT_FILE_BUF_SIZE: usize = 256 * 1024;

pub fn open_file_or_stdin(filepath: impl AsRef<Path>) -> Result<BufReader<Decompressor<'static>>, Report> {
  let filepath = filepath.as_ref();
  if is_path_stdin(filepath) {
    return open_stdin();
  }
  let file = File::open(filepath).wrap_err_with(|| format!("When opening file '{}'", filepath.display()))?;
  let buf_file = BufReader::with_capacity(DEFAULT_FILE_BUF_SIZE, file);
  let decompressor = Decompressor::from_path(buf_file, filepath)?;
  Ok(BufReader::with_capacity(DEFAULT_FILE_BUF_SIZE, decompressor))
}

fn open_stdin() -> Result<BufReader<Decompressor<'static>>, Report> {
  info!("Reading from standard input");

  #[cfg(not(target_arch = "wasm32"))]
  non_wasm::warn_if_tty();

  let decompressor = Decompressor::new(stdin(), &CompressionType::None)?;
  Ok(BufReader::with_capacity(DEFAULT_FILE_BUF_SIZE, decompressor))
}

pub fn read_file_with<T>(
  filepath: impl AsRef<Path>,
  read: impl FnOnce(BufReader<Decompressor<'static>>) -> Result<T, Report>,
) -> Result<T, Report> {
  let filepath = filepath.as_ref();
  let reader = open_file_or_stdin(filepath)?;
  read(reader).wrap_err_with(|| format!("When reading file '{}'", filepath.display()))
}

pub struct FileWriter {
  inner: BufWriter<Compressor<'static>>,
  filepath: PathBuf,
}

impl FileWriter {
  pub fn filepath(&self) -> &Path {
    &self.filepath
  }

  pub fn finish(self) -> Result<(), Report> {
    let filepath = self.filepath;
    self
      .inner
      .into_inner()
      .map_err(|error| Report::new(error.into_error()))
      .wrap_err("While flushing the output buffer")
      .and_then(Compressor::finish)
      .wrap_err_with(|| format!("When writing file '{}'", filepath.display()))
  }
}

impl Write for FileWriter {
  fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
    self.inner.write(buf)
  }

  fn flush(&mut self) -> std::io::Result<()> {
    self.inner.flush()
  }
}

pub fn create_file_or_stdout(filepath: impl AsRef<Path>) -> Result<FileWriter, Report> {
  let filepath = filepath.as_ref();

  let file: Box<dyn Write + Sync + Send> = if is_path_stdout(filepath) {
    info!("File path is '{}'. Writing to standard output.", filepath.display());
    Box::new(BufWriter::with_capacity(DEFAULT_FILE_BUF_SIZE, stdout()))
  } else {
    ensure_dir(filepath)?;
    Box::new(File::create(filepath).wrap_err_with(|| format!("When creating file '{}'", filepath.display()))?)
  };

  let buf_file = BufWriter::with_capacity(DEFAULT_FILE_BUF_SIZE, file);
  let compressor = Compressor::from_path(buf_file, filepath)?;
  let buf_compressor = BufWriter::with_capacity(DEFAULT_FILE_BUF_SIZE, compressor);
  Ok(FileWriter {
    inner: buf_compressor,
    filepath: filepath.to_owned(),
  })
}

pub fn write_file_with(
  filepath: impl AsRef<Path>,
  write: impl FnOnce(&mut FileWriter) -> Result<(), Report>,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut writer = create_file_or_stdout(filepath)?;
  write(&mut writer).wrap_err_with(|| format!("When writing file '{}'", filepath.display()))?;
  writer.finish()
}

pub fn is_path_stdin(filepath: impl AsRef<Path>) -> bool {
  let filepath = filepath.as_ref();
  filepath == "-" || filepath == "/dev/stdin"
}

pub fn is_path_stdout(filepath: impl AsRef<Path>) -> bool {
  let filepath = filepath.as_ref();
  filepath == "-" || filepath == "/dev/stdout"
}

#[cfg(not(target_arch = "wasm32"))]
mod non_wasm {
  use log::warn;
  use std::io::{IsTerminal, stdin};

  const TTY_WARNING: &str = r#"Reading from standard input which is a TTY (e.g. an interactive terminal). This is likely not what you meant. Instead:

 - if you want to read from the output of another program, pipe it in and pass '-' as the file path:

    cat /path/to/file | treetime <command> --alignment - <your other flags>

 - if you want to read from a file, pass its path instead of '-'
"#;

  pub(super) fn warn_if_tty() {
    if stdin().is_terminal() {
      warn!("{TTY_WARNING}");
    }
  }
}
