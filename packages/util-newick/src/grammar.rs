use crate::dialect::{NewickAnnotations, NewickDialect, NewickStructure};
use crate::read::options::ReadMode;
use pest::Parser;
use pest::iterators::Pairs;
use pest_derive::Parser;

#[derive(Parser)]
#[grammar = "common.pest"]
#[grammar = "newick.pest"]
pub(crate) struct NewickParser;

pub(crate) fn parse(rule: Rule, input: &str) -> Result<Pairs<'_, Rule>, Box<pest::error::Error<Rule>>> {
  NewickParser::parse(rule, input).map_err(Box::new)
}

pub(crate) fn matches(rule: Rule, input: &str) -> bool {
  NewickParser::parse(rule, input).is_ok()
}

pub(crate) const fn start_rule(dialect: NewickDialect, mode: ReadMode) -> Rule {
  match (dialect.structure, dialect.annotations, mode) {
    (NewickStructure::Classic, NewickAnnotations::Plain, ReadMode::Strict) => Rule::classic_plain_strict,
    (NewickStructure::Classic, NewickAnnotations::Plain, ReadMode::Tolerant) => Rule::classic_plain_tolerant,
    (NewickStructure::Classic, NewickAnnotations::Beast, ReadMode::Strict) => Rule::classic_beast_strict,
    (NewickStructure::Classic, NewickAnnotations::Beast, ReadMode::Tolerant) => Rule::classic_beast_tolerant,
    (NewickStructure::Classic, NewickAnnotations::Nhx, ReadMode::Strict) => Rule::classic_nhx_strict,
    (NewickStructure::Classic, NewickAnnotations::Nhx, ReadMode::Tolerant) => Rule::classic_nhx_tolerant,
    (NewickStructure::Classic, NewickAnnotations::MrBayes, ReadMode::Strict) => Rule::classic_mrbayes_strict,
    (NewickStructure::Classic, NewickAnnotations::MrBayes, ReadMode::Tolerant) => Rule::classic_mrbayes_tolerant,
    (NewickStructure::ENewick, NewickAnnotations::Plain, ReadMode::Strict) => Rule::enewick_plain_strict,
    (NewickStructure::ENewick, NewickAnnotations::Plain, ReadMode::Tolerant) => Rule::enewick_plain_tolerant,
    (NewickStructure::ENewick, NewickAnnotations::Beast, ReadMode::Strict) => Rule::enewick_beast_strict,
    (NewickStructure::ENewick, NewickAnnotations::Beast, ReadMode::Tolerant) => Rule::enewick_beast_tolerant,
    (NewickStructure::ENewick, NewickAnnotations::Nhx, ReadMode::Strict) => Rule::enewick_nhx_strict,
    (NewickStructure::ENewick, NewickAnnotations::Nhx, ReadMode::Tolerant) => Rule::enewick_nhx_tolerant,
    (NewickStructure::ENewick, NewickAnnotations::MrBayes, ReadMode::Strict) => Rule::enewick_mrbayes_strict,
    (NewickStructure::ENewick, NewickAnnotations::MrBayes, ReadMode::Tolerant) => Rule::enewick_mrbayes_tolerant,
    (NewickStructure::Rich, NewickAnnotations::Plain, ReadMode::Strict) => Rule::rich_plain_strict,
    (NewickStructure::Rich, NewickAnnotations::Plain, ReadMode::Tolerant) => Rule::rich_plain_tolerant,
    (NewickStructure::Rich, NewickAnnotations::Beast, ReadMode::Strict) => Rule::rich_beast_strict,
    (NewickStructure::Rich, NewickAnnotations::Beast, ReadMode::Tolerant) => Rule::rich_beast_tolerant,
    (NewickStructure::Rich, NewickAnnotations::Nhx, ReadMode::Strict) => Rule::rich_nhx_strict,
    (NewickStructure::Rich, NewickAnnotations::Nhx, ReadMode::Tolerant) => Rule::rich_nhx_tolerant,
    (NewickStructure::Rich, NewickAnnotations::MrBayes, ReadMode::Strict) => Rule::rich_mrbayes_strict,
    (NewickStructure::Rich, NewickAnnotations::MrBayes, ReadMode::Tolerant) => Rule::rich_mrbayes_tolerant,
  }
}

pub(crate) const fn comment_rule(annotations: NewickAnnotations) -> Rule {
  match annotations {
    NewickAnnotations::Plain => Rule::plain_comment_only,
    NewickAnnotations::Beast => Rule::beast_comment_only,
    NewickAnnotations::Nhx => Rule::nhx_comment_only,
    NewickAnnotations::MrBayes => Rule::mrbayes_comment_only,
  }
}

pub(crate) const fn rule_name(rule: Rule) -> &'static str {
  match rule {
    Rule::open => "'('",
    Rule::close => "')'",
    Rule::comma => "','",
    Rule::colon => "':'",
    Rule::end => "';'",
    Rule::number | Rule::number_exact | Rule::integer_exact => "number",
    Rule::label
    | Rule::network_label
    | Rule::network_label_exact
    | Rule::unquoted_label
    | Rule::unquoted_label_plain
    | Rule::support_label
    | Rule::support_label_exact
    | Rule::safe_name
    | Rule::safe_name_exact
    | Rule::label_char
    | Rule::label_end => "label",
    Rule::quoted_label => "quoted label",
    Rule::hybrid_tag | Rule::hybrid_marker | Rule::hybrid_kind | Rule::hybrid_index => "hybrid tag",
    Rule::plain_comment
    | Rule::plain_comment_exact
    | Rule::scan_comment
    | Rule::scan_comment_strict
    | Rule::scan_annotation
    | Rule::scan_annotation_strict
    | Rule::annotation_start
    | Rule::scan_quote_start
    | Rule::scan_string
    | Rule::malformed_annotation
    | Rule::scan_plain_comment
    | Rule::plain_any_comment
    | Rule::plain_comment_only
    | Rule::beast_comment_only
    | Rule::nhx_comment_only
    | Rule::mrbayes_comment_only => "comment",
    Rule::beast_comment | Rule::beast_any_comment => "BEAST annotation",
    Rule::beast_pair | Rule::beast_key | Rule::beast_bare_key | Rule::beast_bare_key_exact => "annotation key",
    Rule::beast_value
    | Rule::beast_scalar
    | Rule::beast_value_end
    | Rule::double_quoted
    | Rule::single_quoted
    | Rule::color
    | Rule::boolean
    | Rule::bare_string => "annotation value",
    Rule::array_open | Rule::array_close | Rule::beast_array | Rule::beast_array_first | Rule::beast_array_next => {
      "annotation array"
    },
    Rule::nhx_comment
    | Rule::nhx_any_comment
    | Rule::nhx_tag
    | Rule::nhx_key
    | Rule::nhx_key_exact
    | Rule::nhx_value
    | Rule::nhx_part
    | Rule::nhx_part_exact
    | Rule::nhx_color_exact
    | Rule::color_channel => "NHX annotation",
    Rule::mrbayes_comment
    | Rule::mrbayes_any_comment
    | Rule::mrbayes_kind
    | Rule::mrbayes_token
    | Rule::mrbayes_token_exact => "MrBayes comment",
    Rule::rooting | Rule::rooting_value | Rule::rooting_exact => "rooting comment",
    Rule::weight | Rule::weight_exact => "tree weight comment",
    Rule::EOI => "end of input",
    Rule::bom => "byte order mark",
    Rule::WHITESPACE
    | Rule::delimiter
    | Rule::beast_ws
    | Rule::plain_open
    | Rule::plain_field
    | Rule::beast_open
    | Rule::beast_field
    | Rule::nhx_open
    | Rule::nhx_field
    | Rule::mrbayes_open
    | Rule::mrbayes_field
    | Rule::classic_plain_tail
    | Rule::classic_plain_items
    | Rule::classic_plain_strict
    | Rule::classic_plain_tolerant
    | Rule::classic_beast_tail
    | Rule::classic_beast_items
    | Rule::classic_beast_preamble
    | Rule::classic_beast_strict
    | Rule::classic_beast_tolerant
    | Rule::classic_nhx_tail
    | Rule::classic_nhx_items
    | Rule::classic_nhx_strict
    | Rule::classic_nhx_tolerant
    | Rule::classic_mrbayes_tail
    | Rule::classic_mrbayes_items
    | Rule::classic_mrbayes_preamble
    | Rule::classic_mrbayes_strict
    | Rule::classic_mrbayes_tolerant
    | Rule::enewick_plain_tail
    | Rule::enewick_plain_items
    | Rule::enewick_plain_strict
    | Rule::enewick_plain_tolerant
    | Rule::enewick_beast_tail
    | Rule::enewick_beast_items
    | Rule::enewick_beast_preamble
    | Rule::enewick_beast_strict
    | Rule::enewick_beast_tolerant
    | Rule::enewick_nhx_tail
    | Rule::enewick_nhx_items
    | Rule::enewick_nhx_strict
    | Rule::enewick_nhx_tolerant
    | Rule::enewick_mrbayes_tail
    | Rule::enewick_mrbayes_items
    | Rule::enewick_mrbayes_preamble
    | Rule::enewick_mrbayes_strict
    | Rule::enewick_mrbayes_tolerant
    | Rule::rich_plain_tail
    | Rule::rich_plain_items
    | Rule::rich_plain_preamble
    | Rule::rich_plain_strict
    | Rule::rich_plain_tolerant
    | Rule::rich_beast_tail
    | Rule::rich_beast_items
    | Rule::rich_beast_preamble
    | Rule::rich_beast_strict
    | Rule::rich_beast_tolerant
    | Rule::rich_nhx_tail
    | Rule::rich_nhx_items
    | Rule::rich_nhx_preamble
    | Rule::rich_nhx_strict
    | Rule::rich_nhx_tolerant
    | Rule::rich_mrbayes_tail
    | Rule::rich_mrbayes_items
    | Rule::rich_mrbayes_preamble
    | Rule::rich_mrbayes_strict
    | Rule::rich_mrbayes_tolerant
    | Rule::scan_word
    | Rule::scan_text
    | Rule::scan_text_final
    | Rule::tree_scan
    | Rule::tree_scan_final
    | Rule::tree_scan_plain
    | Rule::tree_scan_plain_final
    | Rule::trivia_only => "tree",
  }
}
