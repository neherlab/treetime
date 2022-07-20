use crate::*;
use curl::easy::Easy;
use std::fs;

impl Config {
    /// Download archive from AWS S3, and expand into `${out_dir}/*.so`
    pub fn download<P: AsRef<Path>>(&self, out_dir: P) -> Result<()> {
        let out_dir = out_dir.as_ref();
        if out_dir.exists() {
            fs::create_dir_all(&out_dir)?;
        }
        let data = read_from_url(&format!("{}/{}.tar.zst", s3_addr(), self.name()))?;
        let zstd = zstd::stream::read::Decoder::new(data.as_slice())?;
        let mut arc = tar::Archive::new(zstd);
        arc.unpack(&out_dir)?;
        Ok(())
    }
}

/// Helper for download file from URL
///
/// - This function expands obtained data into memory space
///
fn read_from_url(url: &str) -> Result<Vec<u8>> {
    let mut data = Vec::new();
    let mut handle = Easy::new();
    handle.fail_on_error(true)?;
    handle.url(url)?;
    {
        let mut transfer = handle.transfer();
        transfer
            .write_function(|new_data| {
                data.extend_from_slice(new_data);
                Ok(new_data.len())
            })
            .unwrap();
        transfer.perform().unwrap();
    }
    Ok(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! impl_test_download {
        ($name:expr) => {
            paste::item! {
                #[test]
                fn [<download_$name>]() -> Result<()> {
                    let name = $name;
                    let cfg = Config::from_str(name)?;
                    cfg.download(format!("test_download/{}", name))?;
                    Ok(())
                }
            }
        };
    }

    #[cfg(target_os = "windows")]
    mod macos {
        use super::*;
        impl_test_download!("mkl-dynamic-lp64-seq");
        impl_test_download!("mkl-dynamic-ilp64-seq");
        impl_test_download!("mkl-static-lp64-seq");
        impl_test_download!("mkl-static-ilp64-seq");
    }

    #[cfg(target_os = "macos")]
    mod macos {
        use super::*;
        impl_test_download!("mkl-dynamic-lp64-seq");
        impl_test_download!("mkl-dynamic-lp64-iomp");
        impl_test_download!("mkl-dynamic-ilp64-seq");
        impl_test_download!("mkl-dynamic-ilp64-iomp");
    }

    #[cfg(target_os = "linux")]
    mod linux {
        use super::*;
        impl_test_download!("mkl-dynamic-lp64-seq");
        impl_test_download!("mkl-dynamic-lp64-iomp");
        impl_test_download!("mkl-dynamic-ilp64-seq");
        impl_test_download!("mkl-dynamic-ilp64-iomp");
        impl_test_download!("mkl-static-lp64-seq");
        impl_test_download!("mkl-static-lp64-iomp");
        impl_test_download!("mkl-static-ilp64-seq");
        impl_test_download!("mkl-static-ilp64-iomp");
    }
}
