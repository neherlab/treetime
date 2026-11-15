use crate::nwk::{NwkParse, graph_from_newick, log_read_warnings, tree_read_options};
use eyre::{Report, WrapErr};
use std::io::Read;
use std::path::Path;
use treetime_graph::merge_by_name::{MergedGraph, NamedTree, merge_trees_by_name};
use treetime_utils::io::file::read_file_with;
use treetime_utils::make_error;
use util_newick::{NewickDialect, is_nexus, newick_from_str, nexus_from_str};

pub const TREE_EXTENSIONS: [&str; 6] = ["nwk", "newick", "tree", "tre", "nex", "nexus"];

pub fn tree_read_file(filepath: impl AsRef<Path>, dialect: NewickDialect) -> Result<NwkParse, Report> {
  read_file_with(filepath, |reader| tree_read(reader, dialect))
}

pub fn tree_read(mut reader: impl Read, dialect: NewickDialect) -> Result<NwkParse, Report> {
  let mut input = String::new();
  reader.read_to_string(&mut input).wrap_err("When reading the tree")?;
  let options = tree_read_options(dialect);
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

pub fn segment_trees_read_files(paths: &[impl AsRef<Path>], dialect: NewickDialect) -> Result<MergedGraph, Report> {
  let parsed = paths
    .iter()
    .map(|path| {
      tree_read_file(path, dialect).wrap_err_with(|| format!("When reading tree '{}'", path.as_ref().display()))
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let labels = paths
    .iter()
    .map(|path| path.as_ref().display().to_string())
    .collect::<Vec<_>>();
  let trees = parsed
    .iter()
    .zip(&labels)
    .map(|(tree, label)| NamedTree {
      label,
      graph: &tree.graph,
      names: &tree.names,
      unnamed_nodes: &tree.unnamed_nodes,
      branch_lengths: &tree.branch_lengths,
    })
    .collect::<Vec<_>>();
  merge_trees_by_name(&trees).wrap_err("When merging the segment trees by node name")
}
