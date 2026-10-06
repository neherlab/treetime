use crate::nwk::{NwkParse, graph_from_newick, log_read_warnings, tree_read_options};
use eyre::{Report, WrapErr};
use std::io::Read;
use std::path::Path;
use treetime_utils::io::file::read_file_with;
use treetime_utils::make_error;
use util_newick::{is_nexus, newick_from_str, nexus_from_str};

pub const TREE_EXTENSIONS: [&str; 6] = ["nwk", "newick", "tree", "tre", "nex", "nexus"];

pub fn tree_read_file(filepath: impl AsRef<Path>) -> Result<NwkParse, Report> {
  read_file_with(filepath, tree_read)
}

pub fn tree_read(mut reader: impl Read) -> Result<NwkParse, Report> {
  let mut input = String::new();
  reader.read_to_string(&mut input).wrap_err("When reading the tree")?;
  let options = tree_read_options();
  if !is_nexus(&input) {
    let tree = newick_from_str(&input, &options).wrap_err("When reading Newick")?;
    log_read_warnings(&tree.warnings);
    return graph_from_newick(&tree.graph).wrap_err("When reading Newick");
  }
  let mut file = nexus_from_str(&input, &options).wrap_err("When reading Nexus")?;
  log_read_warnings(&file.warnings);
  match (file.trees.pop(), file.trees.len()) {
    (Some(tree), 0) => {
      log_read_warnings(&tree.tree.warnings);
      graph_from_newick(&tree.tree.graph).wrap_err_with(|| format!("When reading Nexus tree '{}'", tree.name))
    },
    (None, _) => make_error!("The Nexus file contains no tree"),
    (Some(_), others) => make_error!(
      "The Nexus file contains {} trees, but TreeTime reads exactly one tree",
      others + 1
    ),
  }
}
