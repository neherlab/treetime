//! Nexus container format reader and writer.

use crate::parse::newick_from_string;
use crate::types::{NewickGraph, NewickWriteOptions, NexusTree};
use crate::write::newick_to_writer;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::{Read, Write};

/// Parse a Nexus string, returning all trees from TREES blocks.
pub fn nexus_from_string(input: &str) -> Result<Vec<NexusTree>, Report> {
  let input = input.trim();
  let header_end = input.find(|c: char| c.is_ascii_whitespace()).unwrap_or(input.len());
  if !input[..header_end].eq_ignore_ascii_case("#nexus") {
    return Err(eyre::eyre!("Input does not start with #NEXUS header"));
  }

  let blocks = parse_nexus_blocks(input);
  let translate_table = extract_translate_table(&blocks);
  extract_trees(&blocks, &translate_table)
}

/// Parse Nexus from an `impl Read` source.
pub fn nexus_from_reader(mut reader: impl Read) -> Result<Vec<NexusTree>, Report> {
  let mut input = String::new();
  reader.read_to_string(&mut input)?;
  nexus_from_string(&input)
}

/// Write trees as a Nexus string with Taxa and Trees blocks.
pub fn nexus_to_string(trees: &[NexusTree], options: &NewickWriteOptions) -> Result<String, Report> {
  let mut buf = Vec::new();
  nexus_to_writer(&mut buf, trees, options)?;
  Ok(String::from_utf8(buf).expect("Nexus output should be valid UTF-8"))
}

/// Write trees in Nexus format to an `impl Write` destination.
pub fn nexus_to_writer(
  writer: &mut impl Write,
  trees: &[NexusTree],
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  writeln!(writer, "#NEXUS")?;
  writeln!(writer)?;

  // Collect all leaf names across all trees
  let mut taxa: Vec<String> = Vec::new();
  for tree in trees {
    collect_leaf_names(&tree.graph, &mut taxa);
  }
  taxa.sort();
  taxa.dedup();

  // Taxa block
  writeln!(writer, "Begin Taxa;")?;
  writeln!(writer, "  Dimensions ntax={};", taxa.len())?;
  writeln!(writer, "  TaxLabels")?;
  for name in &taxa {
    if needs_nexus_quoting(name) {
      writeln!(writer, "    '{}'", name.replace('\'', "''"))?;
    } else {
      writeln!(writer, "    {name}")?;
    }
  }
  writeln!(writer, "  ;")?;
  writeln!(writer, "End;")?;
  writeln!(writer)?;

  // Trees block
  writeln!(writer, "Begin Trees;")?;
  for tree in trees {
    let name = &tree.name;
    if needs_nexus_quoting(name) {
      write!(writer, "  Tree '{}' = ", name.replace('\'', "''"))?;
    } else {
      write!(writer, "  Tree {name} = ")?;
    }
    newick_to_writer(writer, &tree.graph, options)?;
    writeln!(writer)?;
  }
  writeln!(writer, "End;")?;

  Ok(())
}

fn collect_leaf_names(graph: &NewickGraph, taxa: &mut Vec<String>) {
  for node in &graph.nodes {
    if node.children.is_empty() {
      if let Some(name) = &node.name {
        taxa.push(name.clone());
      }
    }
  }
}

fn needs_nexus_quoting(s: &str) -> bool {
  s.contains(|c: char| matches!(c, '(' | ')' | '[' | ']' | ',' | ';' | ':' | '\'') || c.is_whitespace())
}

struct NexusBlock {
  name: String,
  content: String,
}

#[allow(clippy::string_slice)]
fn parse_nexus_blocks(input: &str) -> Vec<NexusBlock> {
  let mut blocks = Vec::new();
  let lower = input.to_ascii_lowercase();
  let mut pos = 0;

  while pos < input.len() {
    // Find "begin <name>;"
    let Some(begin_pos) = lower[pos..].find("begin ") else {
      break;
    };
    let begin_start = pos + begin_pos + 6; // after "begin "

    let Some(semi_pos) = input[begin_start..].find(';') else {
      break;
    };
    let block_name = input[begin_start..begin_start + semi_pos].trim().to_ascii_lowercase();
    let content_start = begin_start + semi_pos + 1;

    // Find "end;"
    let Some(end_pos) = lower[content_start..].find("end;") else {
      break;
    };
    let content = input[content_start..content_start + end_pos].to_owned();

    blocks.push(NexusBlock {
      name: block_name,
      content,
    });

    pos = content_start + end_pos + 4;
  }

  blocks
}

#[allow(clippy::string_slice)]
fn extract_translate_table(blocks: &[NexusBlock]) -> BTreeMap<String, String> {
  let mut table = BTreeMap::new();

  for block in blocks {
    if block.name != "trees" {
      continue;
    }

    let lower = block.content.to_ascii_lowercase();
    let Some(translate_pos) = lower.find("translate") else {
      continue;
    };

    let after_translate = &block.content[translate_pos + 9..];
    let Some(semi_pos) = after_translate.find(';') else {
      continue;
    };
    let translate_body = &after_translate[..semi_pos];

    for entry in split_translate_entries(translate_body) {
      let entry = entry.trim();
      if entry.is_empty() {
        continue;
      }
      let mut parts = entry.splitn(2, char::is_whitespace);
      let Some(key) = parts.next() else { continue };
      let Some(value) = parts.next() else { continue };
      let key = key.trim().to_owned();
      let value = value.trim().to_owned();
      let value = strip_nexus_quotes(&value);
      table.insert(key, value);
    }
  }

  table
}

fn split_translate_entries(body: &str) -> Vec<String> {
  let mut entries = Vec::new();
  let mut current = String::new();
  let mut in_single_quote = false;

  for ch in body.chars() {
    match ch {
      '\'' => {
        in_single_quote = !in_single_quote;
        current.push(ch);
      },
      ',' if !in_single_quote => {
        entries.push(std::mem::take(&mut current));
      },
      _ => {
        current.push(ch);
      },
    }
  }
  if !current.is_empty() {
    entries.push(current);
  }
  entries
}

fn strip_nexus_quotes(s: &str) -> String {
  if let Some(inner) = s.strip_prefix('\'').and_then(|i| i.strip_suffix('\'')) {
    return inner.replace("''", "'");
  }
  if let Some(inner) = s.strip_prefix('"').and_then(|i| i.strip_suffix('"')) {
    return inner.replace("\"\"", "\"");
  }
  s.to_owned()
}

#[allow(clippy::string_slice)]
fn extract_trees(blocks: &[NexusBlock], translate_table: &BTreeMap<String, String>) -> Result<Vec<NexusTree>, Report> {
  let mut trees = Vec::new();

  for block in blocks {
    if block.name != "trees" {
      continue;
    }

    // Find all "tree <name> = <newick>;" entries
    let lower = block.content.to_ascii_lowercase();
    let mut search_pos = 0;

    while search_pos < block.content.len() {
      // Skip past translate block if present
      let remaining_lower = &lower[search_pos..];
      let Some(tree_pos) = find_tree_command(remaining_lower) else {
        break;
      };

      let abs_pos = search_pos + tree_pos + 4; // after "tree"
      let after_tree = &block.content[abs_pos..];

      // Find "=" separating name from newick
      let Some(eq_pos) = after_tree.find('=') else {
        break;
      };
      let tree_name = after_tree[..eq_pos].trim();
      let tree_name = strip_nexus_quotes(tree_name);

      // Find ";" terminating the newick string
      let after_eq = &after_tree[eq_pos + 1..];
      let Some(semi_pos) = find_newick_semicolon(after_eq) else {
        break;
      };
      let nwk_str = after_eq[..=semi_pos].trim();

      let graph = parse_newick_with_translate(nwk_str, translate_table)
        .wrap_err_with(|| format!("In Nexus tree '{tree_name}'"))?;

      trees.push(NexusTree { name: tree_name, graph });

      search_pos = abs_pos + eq_pos + 1 + semi_pos + 1;
    }
  }

  Ok(trees)
}

#[allow(clippy::string_slice)]
fn find_tree_command(lower: &str) -> Option<usize> {
  let mut pos = 0;
  while pos < lower.len() {
    let found = lower[pos..].find("tree")?;
    let abs = pos + found;

    // Must be preceded by whitespace or start of string
    let before_ok = abs == 0 || lower.as_bytes()[abs - 1].is_ascii_whitespace() || lower.as_bytes()[abs - 1] == b';';
    // Must be followed by whitespace
    let after_ok = abs + 4 < lower.len() && lower.as_bytes()[abs + 4].is_ascii_whitespace();

    // Must not be inside "translate" keyword
    if before_ok && after_ok {
      let prefix = &lower[..abs].trim_end();
      if !prefix.ends_with("translat") {
        return Some(found + pos);
      }
    }

    pos = abs + 4;
  }
  None
}

fn find_newick_semicolon(s: &str) -> Option<usize> {
  let mut in_single_quote = false;
  let mut bracket_depth = 0_u32;
  for (i, c) in s.char_indices() {
    match c {
      '\'' if bracket_depth == 0 => in_single_quote = !in_single_quote,
      '[' if !in_single_quote => bracket_depth += 1,
      ']' if !in_single_quote && bracket_depth > 0 => bracket_depth -= 1,
      ';' if !in_single_quote && bracket_depth == 0 => return Some(i),
      _ => {},
    }
  }
  None
}

fn parse_newick_with_translate(
  nwk_str: &str,
  translate_table: &BTreeMap<String, String>,
) -> Result<NewickGraph, Report> {
  let mut graph = newick_from_string(nwk_str)?;

  if translate_table.is_empty() {
    return Ok(graph);
  }

  // Replace integer tokens with full taxon names
  for node in &mut graph.nodes {
    if let Some(name) = &node.name {
      if let Some(translated) = translate_table.get(name) {
        node.name = Some(translated.clone());
      }
    }
  }

  Ok(graph)
}
