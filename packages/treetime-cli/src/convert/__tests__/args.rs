#[cfg(test)]
mod tests {
  use crate::convert::args::{TreeFormat, guess_tree_format_from_filename};
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  // Compound extensions
  #[case::auspice_json(    "tree.auspice.json",    Some(TreeFormat::Auspice))]
  #[case::graph_json(      "tree.graph.json",      Some(TreeFormat::PhyloGraph))]
  #[case::mat_json(        "tree.mat.json",         Some(TreeFormat::MatJson))]
  #[case::mat_pb(          "tree.mat.pb",           Some(TreeFormat::MatPb))]
  #[case::phylo_xml(       "tree.phylo.xml",        Some(TreeFormat::Phyloxml))]
  #[case::phyloxml_json(   "tree.phyloxml.json",    Some(TreeFormat::PhyloxmlJson))]
  // Simple extensions
  #[case::nwk(             "tree.nwk",              Some(TreeFormat::Newick))]
  #[case::newick(          "tree.newick",            Some(TreeFormat::Newick))]
  #[case::nex(             "tree.nex",               Some(TreeFormat::Nexus))]
  #[case::nexus(           "tree.nexus",             Some(TreeFormat::Nexus))]
  // Case insensitivity
  #[case::upper_nwk(       "tree.NWK",              Some(TreeFormat::Newick))]
  #[case::mixed_auspice(   "tree.Auspice.JSON",     Some(TreeFormat::Auspice))]
  // With directory path
  #[case::path_compound(   "/data/trees/tree.auspice.json", Some(TreeFormat::Auspice))]
  #[case::path_simple(     "/data/trees/tree.nwk",          Some(TreeFormat::Newick))]
  // Unknown or no extension
  #[case::unknown_ext(     "tree.xyz",              None)]
  #[case::no_extension(    "tree",                  None)]
  #[case::bare_json(       "tree.json",             None)]
  #[case::bare_xml(        "tree.xml",              None)]
  #[case::bare_pb(         "tree.pb",               None)]
  #[trace]
  fn test_args_guess_tree_format(#[case] filename: &str, #[case] expected: Option<TreeFormat>) {
    let actual = guess_tree_format_from_filename(filename);
    assert_eq!(expected, actual, "filename: {filename}");
  }
}
