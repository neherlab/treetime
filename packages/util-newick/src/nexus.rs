use crate::parse::newick_from_string;
use crate::types::{NewickGraph, NewickLabel, NewickReadOptions, NewickWriteOptions, NexusTree};
use crate::write::write_newick;
use eyre::{Report, WrapErr, eyre};
use pest::Parser;
use pest::error::Error as PestError;
use pest::iterators::Pair;
use pest_derive::Parser;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Read};

pub fn nexus_from_reader(mut reader: impl Read, options: &NewickReadOptions) -> Result<Vec<NexusTree>, Report> {
  let mut input = String::new();
  reader.read_to_string(&mut input).wrap_err("When reading Nexus input")?;
  nexus_from_string(&input, options)
}

pub fn nexus_from_string(input: &str, options: &NewickReadOptions) -> Result<Vec<NexusTree>, Report> {
  if !is_nexus(input) {
    return Err(eyre!("Input does not start with #NEXUS header"));
  }
  let blocks = NexusParser::parse(Rule::nexus, input)
    .map_err(|error| Report::new(rename_rules(error)))
    .wrap_err("Failed to parse Nexus input")?;
  let mut trees = Vec::new();
  for block in blocks.filter(|pair| pair.as_rule() == Rule::block) {
    read_block(block, options, &mut trees)?;
  }
  Ok(trees)
}

pub fn is_nexus(input: &str) -> bool {
  input
    .trim_start_matches('\u{feff}')
    .trim_start()
    .get(..6)
    .is_some_and(|header| header.eq_ignore_ascii_case("#nexus"))
}

pub fn nexus_to_writer(
  writer: &mut impl io::Write,
  trees: &[NexusTree],
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  write_nexus(writer, trees, options).wrap_err("When writing Nexus")
}

pub fn nexus_to_string(trees: &[NexusTree], options: &NewickWriteOptions) -> Result<String, Report> {
  let mut buffer = Vec::new();
  nexus_to_writer(&mut buffer, trees, options)?;
  Ok(String::from_utf8(buffer)?)
}

fn write_nexus_word(writer: &mut impl io::Write, word: &str) -> Result<(), Report> {
  if needs_nexus_quoting(word) {
    write!(writer, "'{}'", word.replace('\'', "''"))?;
  } else {
    writer.write_all(word.as_bytes())?;
  }
  Ok(())
}

#[derive(Parser)]
#[grammar_inline = r##"
nexus             = _{ SOI ~ "\u{FEFF}"? ~ header ~ block* ~ EOI }
header            = @{ "#" ~ ^"nexus" ~ boundary }
block             =  { begin_keyword ~ word ~ ";" ~ command* ~ block_end }
block_end         =  { end_keyword ~ ";" }
begin_keyword     = @{ ^"begin" ~ boundary }
end_keyword       = @{ (^"endblock" | ^"end") ~ boundary }
tree_keyword      = @{ (^"utree" | ^"tree") ~ boundary }
translate_keyword = @{ ^"translate" ~ boundary }
command           = _{ tree_command | translate_command | other_command }
tree_command      = ${ tree_head ~ newick_text ~ ";" }
tree_head         = !{ tree_keyword ~ "*"? ~ word ~ "=" }
translate_command =  { translate_keyword ~ (translate_pair ~ ("," ~ translate_pair)*)? ~ ","? ~ ";" }
translate_pair    =  { word ~ word }
other_command     =  { !block_end ~ word ~ command_text ~ ";" }
command_text      = @{ (quoted_word | double_quoted | comment | !";" ~ ANY)* }
newick_text       = @{ (quoted_word | newick_comment | !";" ~ ANY)* }
word              = ${ quoted_word | double_quoted | plain_word }
quoted_word       = @{ "'" ~ ("''" | !"'" ~ ANY)* ~ "'" }
double_quoted     = @{ "\"" ~ ("\"\"" | !"\"" ~ ANY)* ~ "\"" }
plain_word        = @{ (!(";" | "=" | "," | "[" | "]" | "'" | "\"" | WHITESPACE) ~ ANY)+ }
boundary          = _{ !(ASCII_ALPHANUMERIC | "_") }
comment           = @{ "[" ~ (PUSH("[") | "]" ~ DROP | !"]" ~ ANY)* ~ "]" }
newick_comment    = @{ "[" ~ ("\"" ~ (!"\"" ~ ANY)* ~ "\"" | PUSH("[") | "]" ~ DROP | !"]" ~ ANY)* ~ "]" }
COMMENT           = _{ comment }
WHITESPACE        = _{ " " | "\t" | NEWLINE }
"##]
struct NexusParser;

fn rename_rules(error: PestError<Rule>) -> PestError<Rule> {
  error.renamed_rules(|rule| {
    match rule {
      Rule::header => "'#NEXUS' header",
      Rule::block | Rule::begin_keyword => "'Begin' block",
      Rule::block_end | Rule::end_keyword => "'End;'",
      Rule::tree_keyword => "'Tree' command",
      Rule::translate_keyword => "'Translate' command",
      Rule::tree_command | Rule::tree_head => "'Tree' command",
      Rule::translate_command | Rule::translate_pair => "'Translate' entry",
      Rule::other_command | Rule::command | Rule::command_text => "command",
      Rule::newick_text => "Newick tree",
      Rule::word | Rule::quoted_word | Rule::double_quoted | Rule::plain_word => "word",
      Rule::comment | Rule::newick_comment | Rule::COMMENT => "comment",
      Rule::EOI => "end of input",
      Rule::nexus | Rule::boundary | Rule::WHITESPACE => "input",
    }
    .to_owned()
  })
}

fn read_block(block: Pair<'_, Rule>, options: &NewickReadOptions, trees: &mut Vec<NexusTree>) -> Result<(), Report> {
  let mut parts = block.into_inner().filter(|part| part.as_rule() != Rule::begin_keyword);
  let is_trees_block = parts
    .next()
    .is_some_and(|name| word_value(name).eq_ignore_ascii_case("trees"));
  if !is_trees_block {
    return Ok(());
  }
  let mut translate = BTreeMap::new();
  for command in parts {
    match command.as_rule() {
      Rule::translate_command => {
        for pair in command
          .into_inner()
          .filter(|part| part.as_rule() == Rule::translate_pair)
        {
          let mut words = pair.into_inner();
          if let (Some(key), Some(value)) = (words.next(), words.next()) {
            translate.insert(word_value(key), word_value(value));
          }
        }
      },
      Rule::tree_command => trees.push(read_tree(command, &translate, options)?),
      _ => {},
    }
  }
  Ok(())
}

fn read_tree(
  command: Pair<'_, Rule>,
  translate: &BTreeMap<String, String>,
  options: &NewickReadOptions,
) -> Result<NexusTree, Report> {
  let (line, _) = command.line_col();
  let mut parts = command.into_inner();
  let name = parts
    .next()
    .and_then(|head| head.into_inner().find(|part| part.as_rule() == Rule::word))
    .map(word_value)
    .unwrap_or_default();
  let newick = parts.next().map(|text| text.as_str()).unwrap_or_default();
  let mut graph =
    newick_from_string(newick, options).wrap_err_with(|| format!("In Nexus tree '{name}' at line {line}"))?;
  translate_names(&mut graph, translate);
  Ok(NexusTree { name, graph })
}

fn translate_names(graph: &mut NewickGraph, translate: &BTreeMap<String, String>) {
  for node in &mut graph.nodes {
    if let Some(NewickLabel::Name(name)) = &mut node.label {
      if let Some(translated) = translate.get(name.as_str()) {
        name.clone_from(translated);
      }
    }
  }
}

fn word_value(word: Pair<'_, Rule>) -> String {
  let text = word.as_str();
  match word.into_inner().next().map(|inner| inner.as_rule()) {
    Some(Rule::quoted_word) => unquote(text, '\''),
    Some(Rule::double_quoted) => unquote(text, '"'),
    _ => text.to_owned(),
  }
}

fn unquote(text: &str, quote: char) -> String {
  text
    .strip_prefix(quote)
    .and_then(|inner| inner.strip_suffix(quote))
    .unwrap_or(text)
    .replace(&format!("{quote}{quote}"), &quote.to_string())
}

fn write_nexus(writer: &mut impl io::Write, trees: &[NexusTree], options: &NewickWriteOptions) -> Result<(), Report> {
  let taxa: BTreeSet<&str> = trees
    .iter()
    .flat_map(|tree| tree.graph.nodes.iter())
    .filter(|node| node.children.is_empty())
    .filter_map(|node| node.name())
    .collect();

  writeln!(writer, "#NEXUS")?;
  writeln!(writer)?;
  writeln!(writer, "Begin Taxa;")?;
  writeln!(writer, "  Dimensions ntax={};", taxa.len())?;
  writeln!(writer, "  TaxLabels")?;
  for name in taxa {
    write!(writer, "    ")?;
    write_nexus_word(writer, name)?;
    writeln!(writer)?;
  }
  writeln!(writer, "  ;")?;
  writeln!(writer, "End;")?;
  writeln!(writer)?;
  writeln!(writer, "Begin Trees;")?;
  for tree in trees {
    write!(writer, "  Tree ")?;
    write_nexus_word(writer, &tree.name)?;
    write!(writer, " = ")?;
    write_newick(writer, &tree.graph, options).wrap_err_with(|| format!("In Nexus tree '{}'", tree.name))?;
    writeln!(writer)?;
  }
  writeln!(writer, "End;")?;
  Ok(())
}

fn needs_nexus_quoting(word: &str) -> bool {
  word.is_empty()
    || word.contains(|c: char| {
      c.is_whitespace()
        || matches!(
          c,
          '('
            | ')'
            | '['
            | ']'
            | '{'
            | '}'
            | '/'
            | '\\'
            | ','
            | ';'
            | ':'
            | '='
            | '*'
            | '\''
            | '"'
            | '`'
            | '+'
            | '-'
            | '<'
            | '>'
        )
    })
}
