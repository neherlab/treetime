#[cfg(test)]
mod tests {
  use crate::dialect::NewickDialect;
  use crate::read::error::NewickErrorKind;
  use crate::read::options::{NewickReadOptions, ReadMode};
  use crate::read::stream::newick_from_str;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const BEAST: &str = "((A[&rate=1.5]:1,B:1)[&posterior=0.9]:1,C:2);";
  const NHX: &str = "((A[&&NHX:S=human]:1,B:1):1,C:2);";
  const MRBAYES: &str = "[&U]((A:1[&B TK02Brlens 0.1],B:1):1,C:2);";
  const ENEWICK: &str = "((A,(B)x#H1),(x#H1,C));";
  const RICH: &str = "[&U]((A:1:90,(B)#H1:::0.3),(#H1:::0.7,C));";
  const CLASSIC: &str = "((A:1,B:1)95:1,C:2);";
  const GISAID: &str = "(EPI_ISL#402124:1,B:1);";

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_all(          BEAST,    NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::Beast))]
  #[case::beast_own(          BEAST,    vec![NewickDialect::Beast],                       Ok(NewickDialect::Beast))]
  #[case::beast_as_classic(   BEAST,    vec![NewickDialect::Classic],                     Ok(NewickDialect::Classic))]
  #[case::beast_as_nhx(       BEAST,    vec![NewickDialect::Nhx],                         Err(NewickErrorKind::Annotation))]
  #[case::nhx_all(            NHX,      NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::Nhx))]
  #[case::nhx_own(            NHX,      vec![NewickDialect::Nhx],                         Ok(NewickDialect::Nhx))]
  #[case::nhx_as_beast(       NHX,      vec![NewickDialect::Beast],                       Err(NewickErrorKind::Annotation))]
  #[case::mrbayes_all(        MRBAYES,  NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::MrBayes))]
  #[case::mrbayes_as_beast(   MRBAYES,  vec![NewickDialect::Beast],                       Err(NewickErrorKind::Annotation))]
  #[case::enewick_all(        ENEWICK,  NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::Rich))]
  #[case::enewick_own(        ENEWICK,  vec![NewickDialect::ENewick],                     Ok(NewickDialect::ENewick))]
  #[case::enewick_as_classic( ENEWICK,  vec![NewickDialect::Classic],                     Ok(NewickDialect::Classic))]
  #[case::rich_all(           RICH,     NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::Rich))]
  #[case::rich_as_enewick(    RICH,     vec![NewickDialect::ENewick],                     Err(NewickErrorKind::Syntax))]
  #[case::classic_all(        CLASSIC,  NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::Rich))]
  #[case::classic_own(        CLASSIC,  vec![NewickDialect::Classic],                     Ok(NewickDialect::Classic))]
  #[case::order_given(        CLASSIC,  vec![NewickDialect::Classic, NewickDialect::Rich], Ok(NewickDialect::Classic))]
  #[case::gisaid_all(         GISAID,   NewickDialect::ALL.to_vec(),                      Ok(NewickDialect::Rich))]
  #[case::gisaid_treetime(    GISAID,   vec![NewickDialect::Beast, NewickDialect::Nhx, NewickDialect::MrBayes, NewickDialect::Classic], Ok(NewickDialect::Beast))]
  #[trace]
  fn test_read_dialects_selection(#[case] input: &str, #[case] dialects: Vec<NewickDialect>, #[case] expected: Result<NewickDialect, NewickErrorKind>) {
    let options = NewickReadOptions {
      dialects,
      ..NewickReadOptions::default()
    };

    let actual = newick_from_str(input, &options).map(|tree| tree.dialect).map_err(|error| error.kind);

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_read_dialects_mapping_error_moves_to_next_dialect() {
    let options = NewickReadOptions::all_dialects();

    let tree = newick_from_str("((C)#H1,(#H1)#H1);", &options).unwrap();

    assert_eq!(NewickDialect::Nhx, tree.dialect);
  }

  #[test]
  fn test_read_dialects_strict_forms_before_tolerant_forms() {
    let options = NewickReadOptions {
      dialects: vec![NewickDialect::Beast, NewickDialect::MrBayes],
      mode: ReadMode::Tolerant,
      ..NewickReadOptions::default()
    };

    let tree = newick_from_str(MRBAYES, &options).unwrap();

    assert_eq!((NewickDialect::MrBayes, 0), (tree.dialect, tree.warnings.len()));
  }

  #[test]
  fn test_read_dialects_tolerant_form_after_all_strict_forms_fail() {
    let options = NewickReadOptions {
      dialects: vec![NewickDialect::Nhx, NewickDialect::Beast],
      mode: ReadMode::Tolerant,
      ..NewickReadOptions::default()
    };

    let tree = newick_from_str("(A[&a=],B[&&NHX:S=x]);", &options).unwrap();

    assert_eq!((NewickDialect::Nhx, 1), (tree.dialect, tree.warnings.len()));
  }

  #[test]
  fn test_read_dialects_trees_choose_dialects_independently() {
    let input = format!("{BEAST}\n{NHX}\n{ENEWICK}\n");
    let options = NewickReadOptions::all_dialects();

    let dialects: Vec<NewickDialect> = crate::read::stream::newick_trees(input.as_bytes(), options)
      .map(|tree| tree.unwrap().dialect)
      .collect();

    assert_eq!(
      vec![NewickDialect::Beast, NewickDialect::Nhx, NewickDialect::Rich],
      dialects
    );
  }
}
