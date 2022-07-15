use crate::io::fs::read_file_to_string;
use bio::io::newick;
use bio_types::phylogeny::Tree;
use eyre::{Report, WrapErr};
use std::path::Path;

pub fn read_nwk(nwk_file_path: impl AsRef<Path>) -> Result<Tree, Report> {
  let nwk_file_path = nwk_file_path.as_ref();
  let nwk_str = read_file_to_string(&nwk_file_path)?;
  let nwk_tree =
    newick::read(nwk_str.as_bytes()).wrap_err_with(|| format!("When parsing Newick file {nwk_file_path:#?}"))?;
  Ok(nwk_tree)
}
