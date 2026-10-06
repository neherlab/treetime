use crate::dialect::NewickDialect;
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
  match (dialect, mode) {
    (NewickDialect::Classic, ReadMode::Strict) => Rule::classic_strict,
    (NewickDialect::Classic, ReadMode::Tolerant) => Rule::classic_tolerant,
    (NewickDialect::Beast, ReadMode::Strict) => Rule::beast_strict,
    (NewickDialect::Beast, ReadMode::Tolerant) => Rule::beast_tolerant,
    (NewickDialect::MrBayes, ReadMode::Strict) => Rule::mrbayes_strict,
    (NewickDialect::MrBayes, ReadMode::Tolerant) => Rule::mrbayes_tolerant,
    (NewickDialect::Nhx, ReadMode::Strict) => Rule::nhx_strict,
    (NewickDialect::Nhx, ReadMode::Tolerant) => Rule::nhx_tolerant,
    (NewickDialect::ENewick, ReadMode::Strict) => Rule::enewick_strict,
    (NewickDialect::ENewick, ReadMode::Tolerant) => Rule::enewick_tolerant,
    (NewickDialect::Rich, ReadMode::Strict) => Rule::rich_strict,
    (NewickDialect::Rich, ReadMode::Tolerant) => Rule::rich_tolerant,
  }
}

pub(crate) const fn comment_rule(dialect: NewickDialect) -> Rule {
  match dialect {
    NewickDialect::Classic => Rule::classic_comment_only,
    NewickDialect::Beast => Rule::beast_comment_only,
    NewickDialect::MrBayes => Rule::mrbayes_comment_only,
    NewickDialect::Nhx => Rule::nhx_comment_only,
    NewickDialect::ENewick => Rule::enewick_comment_only,
    NewickDialect::Rich => Rule::rich_comment_only,
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
    | Rule::classic_comment
    | Rule::enewick_comment
    | Rule::rich_comment
    | Rule::classic_comment_only
    | Rule::beast_comment_only
    | Rule::mrbayes_comment_only
    | Rule::nhx_comment_only
    | Rule::enewick_comment_only
    | Rule::rich_comment_only => "comment",
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
    Rule::rooting | Rule::rooting_value => "rooting comment",
    Rule::weight => "tree weight comment",
    Rule::EOI => "end of input",
    Rule::bom => "byte order mark",
    Rule::WHITESPACE
    | Rule::delimiter
    | Rule::beast_ws
    | Rule::classic_open
    | Rule::classic_field
    | Rule::classic_tail
    | Rule::classic_items
    | Rule::beast_open
    | Rule::beast_field
    | Rule::beast_tail
    | Rule::beast_items
    | Rule::beast_preamble
    | Rule::mrbayes_open
    | Rule::mrbayes_field
    | Rule::mrbayes_tail
    | Rule::mrbayes_items
    | Rule::mrbayes_preamble
    | Rule::nhx_open
    | Rule::nhx_field
    | Rule::nhx_tail
    | Rule::nhx_items
    | Rule::enewick_open
    | Rule::enewick_field
    | Rule::enewick_tail
    | Rule::enewick_items
    | Rule::rich_open
    | Rule::rich_field
    | Rule::rich_tail
    | Rule::rich_items
    | Rule::rich_preamble
    | Rule::classic_strict
    | Rule::classic_tolerant
    | Rule::beast_strict
    | Rule::beast_tolerant
    | Rule::mrbayes_strict
    | Rule::mrbayes_tolerant
    | Rule::nhx_strict
    | Rule::nhx_tolerant
    | Rule::enewick_strict
    | Rule::enewick_tolerant
    | Rule::rich_strict
    | Rule::rich_tolerant
    | Rule::scan_word
    | Rule::tree_extent
    | Rule::tree_extent_final
    | Rule::trivia_only => "tree",
  }
}
