use pest::Parser;
use pest::iterators::Pairs;
use pest_derive::Parser;

#[derive(Parser)]
#[grammar = "common.pest"]
#[grammar = "nexus.pest"]
pub(crate) struct NexusParser;

pub(crate) fn parse(rule: Rule, input: &str) -> Result<Pairs<'_, Rule>, Box<pest::error::Error<Rule>>> {
  NexusParser::parse(rule, input).map_err(Box::new)
}

pub(crate) fn matches(rule: Rule, input: &str) -> bool {
  NexusParser::parse(rule, input).is_ok()
}

pub(crate) const fn rule_name(rule: Rule) -> &'static str {
  match rule {
    Rule::header | Rule::nexus_header => "'#NEXUS' header",
    Rule::begin_keyword | Rule::begin_command => "'Begin' command",
    Rule::end_keyword | Rule::end_command => "'End' command",
    Rule::tree_keyword
    | Rule::tree_command_ws
    | Rule::tree_head_strict
    | Rule::tree_head_loose
    | Rule::tree_tail
    | Rule::tree_strict
    | Rule::tree_tolerant => "'Tree' command",
    Rule::newick_text => "Newick tree",
    Rule::translate_keyword
    | Rule::translate_strict
    | Rule::translate_tolerant
    | Rule::translate_pair
    | Rule::translate_pair_loose => "'Translate' entry",
    Rule::taxlabels_keyword | Rule::taxlabels_command => "'TaxLabels' command",
    Rule::dimensions_keyword | Rule::dimensions_command => "'Dimensions' command",
    Rule::properties_keyword | Rule::properties_command | Rule::property => "'Properties' command",
    Rule::other_command | Rule::command_text | Rule::known_keyword => "command",
    Rule::word
    | Rule::plain_word
    | Rule::loose_text
    | Rule::quoted_label
    | Rule::nexus_word_safe
    | Rule::nexus_word_exact => "word",
    Rule::double_quoted => "quoted text",
    Rule::plain_comment
    | Rule::scan_comment
    | Rule::scan_comment_strict
    | Rule::scan_annotation
    | Rule::scan_annotation_strict
    | Rule::annotation_start
    | Rule::scan_string
    | Rule::scan_quote_start
    | Rule::COMMENT => "comment",
    Rule::equals => "'='",
    Rule::separator => "','",
    Rule::terminator => "';'",
    Rule::EOI => "end of input",
    Rule::bom => "byte order mark",
    Rule::WHITESPACE
    | Rule::boundary
    | Rule::scan_word
    | Rule::command_extent
    | Rule::command_extent_final
    | Rule::trivia_only
    | Rule::command_strict
    | Rule::command_tolerant => "input",
  }
}
