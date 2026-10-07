use crate::error::{NewickWriteError, WriteContext};
use crate::nexus::grammar::{Rule, matches};
use crate::nexus::types::{NexusTreeRef, NexusWriteOptions};
use crate::write::comments::encode_comment;
use crate::write::newick::write_tree;
use crate::write::options::Spaces;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;
use std::io;

pub fn nexus_to_writer(
  writer: &mut impl io::Write,
  trees: &[NexusTreeRef<'_>],
  options: &NexusWriteOptions,
) -> Result<(), NewickWriteError> {
  let text = nexus_to_string(trees, options)?;
  writer.write_all(text.as_bytes()).context("When writing NEXUS")?;
  Ok(())
}

pub fn nexus_to_string(trees: &[NexusTreeRef<'_>], options: &NexusWriteOptions) -> Result<String, NewickWriteError> {
  let mut text = String::new();
  write_nexus(&mut text, trees, options).context("When writing NEXUS")?;
  Ok(text)
}

fn write_nexus(
  out: &mut String,
  trees: &[NexusTreeRef<'_>],
  options: &NexusWriteOptions,
) -> Result<(), NewickWriteError> {
  let taxa = taxon_labels(trees);
  out.push_str("#NEXUS\n");
  if !taxa.is_empty() {
    out.push_str("\nBegin Taxa;\n");
    writeln!(out, "  Dimensions NTax={};", taxa.len())?;
    out.push_str("  TaxLabels\n");
    for label in &taxa {
      out.push_str("    ");
      push_word(out, label, options);
      out.push('\n');
    }
    out.push_str("  ;\nEnd;\n");
  }
  out.push_str("\nBegin Trees;\n");
  let translate: Option<BTreeMap<String, String>> = options.translate.then(|| {
    taxa
      .iter()
      .enumerate()
      .map(|(idx, label)| ((*label).to_owned(), (idx + 1).to_string()))
      .collect()
  });
  if options.translate && !taxa.is_empty() {
    out.push_str("  Translate\n");
    for (idx, label) in taxa.iter().enumerate() {
      write!(out, "    {} ", idx + 1)?;
      push_word(out, label, options);
      out.push_str(if idx + 1 < taxa.len() { ",\n" } else { "\n" });
    }
    out.push_str("  ;\n");
  }
  for tree in trees {
    write_tree_command(out, tree, options, translate.as_ref())
      .with_context(|| format!("When writing the tree {:?}", tree.name))?;
  }
  out.push_str("End;\n");
  Ok(())
}

fn write_tree_command(
  out: &mut String,
  tree: &NexusTreeRef<'_>,
  options: &NexusWriteOptions,
  translate: Option<&BTreeMap<String, String>>,
) -> Result<(), NewickWriteError> {
  out.push_str("  Tree ");
  push_word(out, tree.name, options);
  for comment in tree.comments {
    if let Some(text) = encode_comment(comment, options.newick.dialect)? {
      out.push(' ');
      out.push_str(&text);
    }
  }
  out.push_str(" = ");
  write_tree(out, tree.graph, &options.newick, translate)?;
  out.push('\n');
  Ok(())
}

fn taxon_labels<'t>(trees: &[NexusTreeRef<'t>]) -> Vec<&'t str> {
  let mut seen = BTreeSet::new();
  let mut labels = Vec::new();
  for tree in trees {
    for node in tree.graph.preorder() {
      if !tree.graph.is_leaf(node) {
        continue;
      }
      if let Some(name) = tree.graph.node(node).name()
        && seen.insert(name)
      {
        labels.push(name);
      }
    }
  }
  labels
}

fn push_word(out: &mut String, word: &str, options: &NexusWriteOptions) {
  let written = match options.newick.spaces {
    Spaces::Quote => word.to_owned(),
    Spaces::Underscore => word.replace(' ', "_"),
  };
  if matches(Rule::nexus_word_exact, &written) {
    out.push_str(&written);
  } else {
    out.push('\'');
    out.push_str(&word.replace('\'', "''"));
    out.push('\'');
  }
}
