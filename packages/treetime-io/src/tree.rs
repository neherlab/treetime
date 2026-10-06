use crate::nwk::{NwkParse, graph_from_newick};
use eyre::{Report, WrapErr};
use std::io::Read;
use std::path::Path;
use treetime_utils::io::file::read_file_with;
use treetime_utils::make_error;
use util_newick::{NewickReadOptions, is_nexus, newick_from_string, nexus_from_string};

pub const TREE_EXTENSIONS: [&str; 6] = ["nwk", "newick", "tree", "tre", "nex", "nexus"];

pub fn tree_read_file(filepath: impl AsRef<Path>) -> Result<NwkParse, Report> {
  read_file_with(filepath, tree_read)
}

pub fn tree_read(mut reader: impl Read) -> Result<NwkParse, Report> {
  let mut input = String::new();
  reader.read_to_string(&mut input).wrap_err("When reading the tree")?;
  let options = NewickReadOptions::default();
  if !is_nexus(&input) {
    let graph = newick_from_string(&input, &options).wrap_err("When reading Newick")?;
    return graph_from_newick(&graph).wrap_err("When reading Newick");
  }
  let mut trees = nexus_from_string(&input, &options).wrap_err("When reading Nexus")?;
  match (trees.pop(), trees.len()) {
    (Some(tree), 0) => {
      graph_from_newick(&tree.graph).wrap_err_with(|| format!("When reading Nexus tree '{}'", tree.name))
    },
    (None, _) => make_error!("The Nexus file contains no tree"),
    (Some(_), others) => make_error!(
      "The Nexus file contains {} trees, but TreeTime reads exactly one tree",
      others + 1
    ),
  }
}
