use crate::*;
use std::{fs, io};

impl Entry {
    pub fn package(&self, out_dir: &Path) -> Result<PathBuf> {
        fs::create_dir_all(out_dir)?;
        let out = out_dir.join(format!("{}.tar.zst", self.name()));
        if out.exists() {
            bail!("Output archive already exits: {}", out.display());
        }
        let f = fs::File::create(&out)?;
        let buf = io::BufWriter::new(f);
        let zstd = zstd::stream::write::Encoder::new(buf, 6)?;
        let mut ar = tar::Builder::new(zstd);
        ar.mode(tar::HeaderMode::Deterministic);
        for (path, name) in self.found_files() {
            let lib = path.join(&name);
            ar.append_path_with_name(lib, name)?;
        }
        let zstd = ar.into_inner()?;
        zstd.finish()?;
        Ok(out)
    }
}
