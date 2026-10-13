use crate::examples_download::DownloadProgress;
use eyre::{Report, WrapErr};
use reqwest::blocking::Client;
use rustls::crypto::CryptoProvider;
use rustls::crypto::ring::default_provider;
use std::fs::File;
use std::io::{Read, Write};
use std::time::Duration;
use treetime_utils::{make_error, make_report};

const CHUNK_SIZE: usize = 1 << 16;

const CONNECT_TIMEOUT: Duration = Duration::from_secs(30);

pub(crate) fn fetch(url: &str, file: &mut File, progress: &mut impl FnMut(DownloadProgress)) -> Result<(), Report> {
  install_crypto_provider()?;
  let client = Client::builder()
    .user_agent(concat!("treetime/", env!("CARGO_PKG_VERSION")))
    .connect_timeout(CONNECT_TIMEOUT)
    .timeout(None)
    .build()
    .wrap_err("When creating the HTTP client")?;
  let mut response = client.get(url).send()?.error_for_status()?;
  let total = response.content_length();
  let mut received = 0;
  let mut buffer = vec![0; CHUNK_SIZE];
  progress(DownloadProgress { received, total });
  loop {
    let read = response
      .read(&mut buffer)
      .wrap_err_with(|| format!("When reading the response of '{url}'"))?;
    if read == 0 {
      break;
    }
    file
      .write_all(&buffer[..read])
      .wrap_err("When writing the downloaded data")?;
    received += u64::try_from(read)?;
    progress(DownloadProgress { received, total });
  }
  file.flush().wrap_err("When writing the downloaded data")?;
  if let Some(total) = total.filter(|total| *total != received) {
    return make_error!("the download ended after {received} of {total} bytes");
  }
  Ok(())
}

fn install_crypto_provider() -> Result<(), Report> {
  if CryptoProvider::get_default().is_some() || default_provider().install_default().is_ok() {
    return Ok(());
  }
  CryptoProvider::get_default()
    .map(|_| ())
    .ok_or_else(|| make_report!("the TLS crypto provider could not be installed"))
}
