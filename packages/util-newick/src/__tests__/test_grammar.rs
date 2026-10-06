#[cfg(test)]
mod tests {
  use crate::grammar::{Rule, matches};
  use helpers::{prefix, tokens};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::unnamed_hybrid(       "#H1",        "network_label(hybrid_tag(hybrid_marker:'#' hybrid_kind:'H' hybrid_index:'1'))")]
  #[case::hash_inside_name(     "A#x#H1",     "network_label(unquoted_label:'A#x' hybrid_tag(hybrid_marker:'#' hybrid_kind:'H' hybrid_index:'1'))")]
  #[case::quoted_name_tag(      "'x y'#H1",   "network_label(quoted_label:''x y'' hybrid_tag(hybrid_marker:'#' hybrid_kind:'H' hybrid_index:'1'))")]
  #[case::acceptor(             "x##LGT2",    "network_label(unquoted_label:'x' hybrid_tag(hybrid_marker:'##' hybrid_kind:'LGT' hybrid_index:'2'))")]
  #[case::tag_without_kind(     "A#1",        "network_label(unquoted_label:'A' hybrid_tag(hybrid_marker:'#' hybrid_index:'1'))")]
  #[case::digits_then_letters(  "A#H1B",      "network_label(unquoted_label:'A#H1B')")]
  #[case::support(              "80.5/95",    "network_label(support_label(number:'80.5' number:'95'))")]
  #[case::support_then_text(    "80.5/95x",   "network_label(unquoted_label:'80.5/95x')")]
  #[case::quoted_number(        "'95'",       "network_label(quoted_label:''95'')")]
  #[case::apostrophe_inside(    "it's",       "network_label(unquoted_label:'it's')")]
  #[trace]
  fn test_grammar_network_label(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(Rule::network_label, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::hash_is_a_name_char(  "A#H1",       "label(unquoted_label_plain:'A#H1')")]
  #[case::gisaid_name(          "EPI_ISL#402124", "label(unquoted_label_plain:'EPI_ISL#402124')")]
  #[case::support(              "80.5/95",    "label(support_label(number:'80.5' number:'95'))")]
  #[case::quoted_number(        "'95'",       "label(quoted_label:''95'')")]
  #[case::escaped_quote(        "'it''s'",    "label(quoted_label:''it''s'')")]
  #[case::exponent(             "1e-5",       "label(support_label(number:'1e-5'))")]
  #[case::infinity_is_a_name(   "inf",        "label(unquoted_label_plain:'inf')")]
  #[trace]
  fn test_grammar_label(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(Rule::label, input));
  }

  #[test]
  fn test_grammar_label_ends_before_whitespace() {
    assert_eq!(Some("A#H1".to_owned()), prefix(Rule::network_label, "A#H1 ,"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nested(              "[a[b]c]",                    "plain_comment:'[a[b]c]'")]
  #[case::quote_is_plain_text( "[5\" tall]",                 "plain_comment:'[5\" tall]'")]
  #[case::empty(               "[]",                         "plain_comment:'[]'")]
  #[trace]
  fn test_grammar_plain_comment(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(Rule::plain_comment, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::number(            "[&a=1]",                "beast_comment(beast_pair(beast_bare_key:'a' number:'1'))")]
  #[case::bare_key(          "[&flag]",               "beast_comment(beast_pair(beast_bare_key:'flag'))")]
  #[case::double_quoted(     "[&a=\"x]y\"]",          "beast_comment(beast_pair(beast_bare_key:'a' double_quoted:'\"x]y\"'))")]
  #[case::single_quoted(     "[&a='x,y']",            "beast_comment(beast_pair(beast_bare_key:'a' single_quoted:''x,y''))")]
  #[case::quoted_key(        "[&\"a b\"=1]",          "beast_comment(beast_pair(double_quoted:'\"a b\"' number:'1'))")]
  #[case::color(             "[&c=#FF0000]",          "beast_comment(beast_pair(beast_bare_key:'c' color:'#FF0000'))")]
  #[case::boolean(           "[&b=True]",             "beast_comment(beast_pair(beast_bare_key:'b' boolean:'True'))")]
  #[case::bare_string(       "[&s=Cote d'Ivoire]",    "beast_comment(beast_pair(beast_bare_key:'s' bare_string:'Cote d'Ivoire'))")]
  #[case::number_then_text(  "[&s=1.5abc]",           "beast_comment(beast_pair(beast_bare_key:'s' bare_string:'1.5abc'))")]
  #[case::spaces(            "[& a = 1 , b = x ]",    "beast_comment(beast_pair(beast_bare_key:'a' number:'1') beast_pair(beast_bare_key:'b' bare_string:'x'))")]
  #[case::empty_array(       "[&a={}]",               "beast_comment(beast_pair(beast_bare_key:'a' array_open:'{' array_close:'}'))")]
  #[case::nested_array(      "[&a={1,{2,{}}}]",       "beast_comment(beast_pair(beast_bare_key:'a' array_open:'{' number:'1' array_open:'{' number:'2' array_open:'{' array_close:'}' array_close:'}' array_close:'}'))")]
  #[case::empty_comment(     "[&]",                   "beast_comment")]
  #[trace]
  fn test_grammar_beast_comment(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(Rule::beast_comment, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::missing_value(     "[&a=]")]
  #[case::trailing_comma(    "[&a={1,}]")]
  #[case::leading_comma(     "[&a={,1}]")]
  #[case::unclosed_array(    "[&a={1]")]
  #[case::unclosed_string(   "[&a=\"x]")]
  #[case::nhx(               "[&&NHX:S=x]")]
  #[case::mrbayes(           "[&B name 1]")]
  #[trace]
  fn test_grammar_beast_comment_rejects(#[case] input: &str) {
    assert_eq!((None, true), (tokens(Rule::beast_comment, input), matches(Rule::malformed_annotation, input)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tags(        "[&&NHX:S=human:B=90]",      "nhx_comment(nhx_tag(nhx_key:'S' nhx_value(nhx_part:'human')) nhx_tag(nhx_key:'B' nhx_value(nhx_part:'90')))")]
  #[case::lowercase(   "[&&nhx:D]",                 "nhx_comment(nhx_tag(nhx_key:'D'))")]
  #[case::parts(       "[&&NHX:Ev=1>2>dup]",        "nhx_comment(nhx_tag(nhx_key:'Ev' nhx_value(nhx_part:'1' nhx_part:'2' nhx_part:'dup')))")]
  #[case::comma_value( "[&&NHX:m=A1T,C2G]",         "nhx_comment(nhx_tag(nhx_key:'m' nhx_value(nhx_part:'A1T,C2G')))")]
  #[case::empty(       "[&&NHX]",                   "nhx_comment")]
  #[trace]
  fn test_grammar_nhx_comment(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(Rule::nhx_comment, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::branch(      "[&B TK02Brlens 0.1]",       "mrbayes_comment(mrbayes_kind:'B' mrbayes_token:'TK02Brlens' mrbayes_token:'0.1')")]
  #[case::event(       "[&E ibr 2: 0.1 0.2]",       "mrbayes_comment(mrbayes_kind:'E' mrbayes_token:'ibr' mrbayes_token:'2:' mrbayes_token:'0.1' mrbayes_token:'0.2')")]
  #[case::node(        "[&N height 3.5]",           "mrbayes_comment(mrbayes_kind:'N' mrbayes_token:'height' mrbayes_token:'3.5')")]
  #[trace]
  fn test_grammar_mrbayes_comment(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(Rule::mrbayes_comment, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rooted(       Rule::rooting, "[&R]",      "rooting(rooting_value:'R')")]
  #[case::unrooted(     Rule::rooting, "[&u]",      "rooting(rooting_value:'u')")]
  #[case::weight(       Rule::weight,  "[&W 0.25]", "weight(number:'0.25')")]
  #[case::weight_lower( Rule::weight,  "[&w  1]",   "weight(number:'1')")]
  #[trace]
  fn test_grammar_tree_comments(#[case] rule: Rule, #[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(rule, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::quote_in_plain_comments( Rule::classic_strict, "(A[5\" tall],B[say \"hi\"]);",  "classic_strict(open:'(' label(unquoted_label_plain:'A') plain_comment:'[5\" tall]' comma:',' label(unquoted_label_plain:'B') plain_comment:'[say \"hi\"]' close:')' end:';' EOI:'')")]
  #[case::rich_fields(             Rule::rich_strict,    "[&U](A:1::0.4)#H1;",             "rich_strict(rooting(rooting_value:'U') open:'(' network_label(unquoted_label:'A') colon:':' number:'1' colon:':' colon:':' number:'0.4' close:')' network_label(hybrid_tag(hybrid_marker:'#' hybrid_kind:'H' hybrid_index:'1')) end:';' EOI:'')")]
  #[case::beast_branch(            Rule::beast_strict,   "(A:[&r=1]2[c]);",                "beast_strict(open:'(' label(unquoted_label_plain:'A') colon:':' beast_comment(beast_pair(beast_bare_key:'r' number:'1')) number:'2' plain_comment:'[c]' close:')' end:';' EOI:'')")]
  #[case::comment_before_open(     Rule::classic_strict, "[c](A);",                        "classic_strict(plain_comment:'[c]' open:'(' label(unquoted_label_plain:'A') close:')' end:';' EOI:'')")]
  #[case::enewick_annotation(      Rule::enewick_strict, "(A[&a=1]);",                     "enewick_strict(open:'(' network_label(unquoted_label:'A') malformed_annotation:'[&a=1]' close:')' end:';' EOI:'')")]
  #[trace]
  fn test_grammar_tree_tokens(#[case] rule: Rule, #[case] input: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), tokens(rule, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::classic_two_fields(   Rule::classic_strict,  "(A:1:2);")]
  #[case::beast_second_label(   Rule::beast_strict,    "(A B);")]
  #[case::strict_bom(           Rule::classic_strict,  "\u{feff}(A);")]
  #[case::strict_no_semicolon(  Rule::classic_strict,  "(A)")]
  #[case::rich_four_fields(     Rule::rich_strict,     "(A:1:2:3:4);")]
  #[trace]
  fn test_grammar_tree_rejects(#[case] rule: Rule, #[case] input: &str) {
    assert_eq!(None, tokens(rule, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tolerant_bom(           Rule::classic_tolerant, "\u{feff}(A)")]
  #[case::tolerant_no_semicolon(  Rule::beast_tolerant,   "(A[&a=1])")]
  #[trace]
  fn test_grammar_tolerant_accepts(#[case] rule: Rule, #[case] input: &str) {
    assert_ne!(None, tokens(rule, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::safe_name(          Rule::safe_name_exact,      "A_b-1",     true)]
  #[case::safe_leading_space( Rule::safe_name_exact,      " A",        false)]
  #[case::safe_trailing_space(Rule::safe_name_exact,      "A ",        false)]
  #[case::safe_quote(         Rule::safe_name_exact,      "it's",      false)]
  #[case::safe_hybrid_tail(   Rule::safe_name_exact,      "A#1",       false)]
  #[case::safe_hash(          Rule::safe_name_exact,      "A#x",       true)]
  #[case::safe_empty(         Rule::safe_name_exact,      "",          false)]
  #[case::number(             Rule::number_exact,         "-1.5e3",    true)]
  #[case::number_space(       Rule::number_exact,         " 1",        false)]
  #[case::number_inf(         Rule::number_exact,         "inf",       false)]
  #[case::support(            Rule::support_label_exact,  "80/95",     true)]
  #[case::bare_key(           Rule::beast_bare_key_exact, "height_95%_HPD", true)]
  #[case::bare_key_space(     Rule::beast_bare_key_exact, "a b",       false)]
  #[case::bare_key_ampersand( Rule::beast_bare_key_exact, "&NHX",      false)]
  #[case::nhx_key(            Rule::nhx_key_exact,        "GN",        true)]
  #[case::nhx_key_colon(      Rule::nhx_key_exact,        "a:b",       false)]
  #[case::nhx_part_greater(   Rule::nhx_part_exact,       "a>b",       false)]
  #[case::nhx_part_comma(     Rule::nhx_part_exact,       "a,b c",     true)]
  #[case::integer(            Rule::integer_exact,        "+12",       true)]
  #[case::integer_decimal(    Rule::integer_exact,        "1.0",       false)]
  #[case::nhx_color(          Rule::nhx_color_exact,      "255.0.12",  true)]
  #[case::nhx_color_short(    Rule::nhx_color_exact,      "255.0",     false)]
  #[case::mrbayes_token(      Rule::mrbayes_token_exact,  "a]",        false)]
  #[case::plain_balanced(     Rule::plain_comment_exact,  "[a[b]]",    true)]
  #[case::plain_unbalanced(   Rule::plain_comment_exact,  "[a]b]",     false)]
  #[trace]
  fn test_grammar_exact_rules(#[case] rule: Rule, #[case] input: &str, #[case] expected: bool) {
    assert_eq!(expected, matches(rule, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::quoted_semicolon(   Rule::tree_scan,       "(A,'B;C');(D);",           Some("(A,'B;C')"))]
  #[case::comment_semicolon(  Rule::tree_scan,       "(A[x;y]);",                Some("(A[x;y])"))]
  #[case::string_bracket(     Rule::tree_scan,       "(A[&a=\"x];y\"]);(B);",   Some("(A[&a=\"x];y\"])"))]
  #[case::nhx_quote(          Rule::tree_scan,       "(A[&&NHX:S=a\"b];(B\");", Some("(A[&&NHX:S=a\"b]"))]
  #[case::apostrophe_in_name( Rule::tree_scan,       "(it's,B);(C,D'x);",        Some("(it's,B)"))]
  #[case::open_comment(       Rule::tree_scan,       "(A[x;",                    Some("(A"))]
  #[case::open_quote(         Rule::tree_scan,       "(A,'x;",                   Some("(A,"))]
  #[case::open_string(        Rule::tree_scan,       "(A[&a=\"x];",             Some("(A"))]
  #[case::final_open_comment( Rule::tree_scan_final, "(A[x;",                    Some("(A[x"))]
  #[case::final_no_semicolon( Rule::tree_scan_final, "(A,B)",                    Some("(A,B)"))]
  #[trace]
  fn test_grammar_tree_scan(#[case] rule: Rule, #[case] input: &str, #[case] expected: Option<&str>) {
    assert_eq!(expected.map(str::to_owned), prefix(rule, input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::comments(  " [c]\n[d [e]] ",  true)]
  #[case::empty(     "",                true)]
  #[case::bom(       "\u{feff}\n",      true)]
  #[case::tree(      "[c] (A)",         false)]
  #[trace]
  fn test_grammar_trivia_only(#[case] input: &str, #[case] expected: bool) {
    assert_eq!(expected, matches(Rule::trivia_only, input));
  }

  mod helpers {
    use crate::grammar::{Rule, parse};
    use pest::iterators::{Pair, Pairs};

    pub(super) fn tokens(rule: Rule, input: &str) -> Option<String> {
      let pairs = parse(rule, input).ok()?;
      let end = pairs.clone().next_back().map(|pair| pair.as_span().end());
      (end == Some(input.len())).then(|| render(pairs))
    }

    pub(super) fn prefix(rule: Rule, input: &str) -> Option<String> {
      let mut pairs = parse(rule, input).ok()?;
      pairs.next().map(|pair| pair.as_str().to_owned())
    }

    fn render(pairs: Pairs<'_, Rule>) -> String {
      pairs.map(render_pair).collect::<Vec<_>>().join(" ")
    }

    fn render_pair(pair: Pair<'_, Rule>) -> String {
      let rule = format!("{:?}", pair.as_rule());
      let text = pair.as_str().to_owned();
      let inner = pair.into_inner();
      if inner.peek().is_none() {
        let is_composite = matches!(rule.as_str(), "beast_comment" | "nhx_comment");
        return if is_composite { rule } else { format!("{rule}:'{text}'") };
      }
      format!("{rule}({})", render(inner))
    }
  }
}
