#[cfg(test)]
mod tests {
  use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
  use crate::commands::shared::topology_order_args::{
    TopologyOrderArg, TopologyOrderArgs, TopologyOrderTargetSourceArg,
  };
  use crate::commands::shared::tree_input::{TreeDialectArg, TreeDialectArgs, read_input_tree};
  use crate::runs::warnings::WarningCollector;
  use clap::Parser;
  use helpers::{clock_dialect, write_tree};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use schemars::{JsonSchema, SchemaGenerator};
  use serde_json::json;
  use std::collections::BTreeMap;
  use tempfile::tempdir;
  use treetime::progress::NoopProgress;
  use treetime_io::nwk::{NewickDialect, TREE_DIALECT_DEFAULT, nwk_read};
  use treetime_utils::io::json::{from_json_value, to_json_value};

  const COALRE: &str = "((A[&segments={0,1}]:1,(B[&segments={0,1}]:1)#H0[&segments={0}]:0.5)[&segments={0,1}]:1,(#H0[&segments={1}]:0.7,C[&segments={0,1}]:1)[&segments={0,1}]:1);";

  #[rustfmt::skip]
  #[rstest]
  #[case::default(        &[],                                  NewickDialect::BEAST)]
  #[case::enewick_beast(  &["--tree-dialect=enewick,beast"],    NewickDialect::ENEWICK_BEAST)]
  #[case::classic_nhx(    &["--tree-dialect", "classic,nhx"],   NewickDialect::NHX)]
  #[case::rich_mrbayes(   &["--tree-dialect=rich,mrbayes"],     NewickDialect { structure: NewickDialect::RICH.structure, annotations: NewickDialect::MRBAYES.annotations })]
  #[trace]
  fn test_tree_input_dialect_flag(#[case] flags: &[&str], #[case] expected: NewickDialect) {
    assert_eq!(expected, clock_dialect(flags).unwrap());
  }

  #[test]
  fn test_tree_input_dialect_flag_rejects_a_half_pair() {
    let actual = clock_dialect(&["--tree-dialect=classic"]).unwrap_err();

    assert_eq!(
      "error: invalid value 'classic' for '--tree-dialect <TREE_DIALECT>'\n  [possible values: classic,plain, classic,beast, classic,nhx, classic,mrbayes, enewick,plain, enewick,beast, enewick,nhx, enewick,mrbayes, rich,plain, rich,beast, rich,nhx, rich,mrbayes]\n\n  tip: a similar value exists: 'classic,nhx'\n\nFor more information, try '--help'.\n",
      actual
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::given(    json!({ "tree_dialect": "enewick,beast" }), Ok(NewickDialect::ENEWICK_BEAST))]
  #[case::absent(   json!({}),                                  Ok(TREE_DIALECT_DEFAULT))]
  #[case::invalid(  json!({ "tree_dialect": "beast" }),         Err("When converting a JSON value: Unexpected: invalid value: \"beast\" is not a Newick dialect: expected <structure>,<annotations> with structure one of classic, enewick, rich and annotations one of plain, beast, nhx, mrbayes".to_owned()))]
  #[trace]
  fn test_tree_input_dialect_config(#[case] config: serde_json::Value, #[case] expected: Result<NewickDialect, String>) {
    let actual = from_json_value::<TreeDialectArgs>(&config).map(|args| args.dialect()).map_err(|error| format!("{error:#}"));

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_tree_input_dialect_config_writes_the_pair_text() {
    let actual = to_json_value(&TreeDialectArgs::default()).unwrap();

    assert_eq!(json!({ "tree_dialect": "classic,beast" }), actual);
  }

  #[test]
  fn test_tree_input_dialect_schema_lists_every_pair() {
    let schema = TreeDialectArg::json_schema(&mut SchemaGenerator::default());

    let expected: Vec<String> = NewickDialect::pairs().map(|dialect| dialect.to_string()).collect();
    assert_eq!(json!({ "type": "string", "enum": expected }), schema.as_value().clone());
  }

  #[test]
  fn test_tree_input_network_is_error() {
    let dir = tempdir().unwrap();
    let path = write_tree(dir.path(), COALRE);
    let log = WarningCollector::new(&NoopProgress);

    let actual = format!(
      "{:#}",
      read_input_tree(&path, NewickDialect::ENEWICK_BEAST, &log).unwrap_err()
    );

    assert_eq!(
      format!(
        "When reading file '{}': When reading Newick: The tree contains the network node '#H0' with 2 parents; TreeTime reads trees only",
        path.display()
      ),
      actual
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::enewick_splits_tags(  NewickDialect::ENEWICK, btreemap! { "A" => 1, "B" => 0 })]
  #[case::beast_keeps_tags(     NewickDialect::BEAST,   btreemap! {})]
  #[trace]
  fn test_tree_input_reference_topology_uses_the_dialect(
    #[case] dialect: NewickDialect,
    #[case] expected: BTreeMap<&str, usize>,
  ) {
    let dir = tempdir().unwrap();
    let reference = write_tree(dir.path(), "(B#1:1,A#2:1);");
    let tree = nwk_read(b"(A:1,B:1)root;".as_slice()).unwrap();
    let names = tree.names();
    let args = TopologyOrderArgs {
      topology_order: Some(TopologyOrderArg::TargetOrder),
      topology_order_target_source: Some(TopologyOrderTargetSourceArg::ReferenceTopology),
      topology_order_target_file: Some(reference),
      ..TopologyOrderArgs::default()
    };

    let spec = args.resolve_topology_order(&tree.graph, &names, None, dialect).unwrap();

    let actual: BTreeMap<&str, usize> = spec
      .target_order
      .iter()
      .map(|(key, &position)| (names[key].as_deref().unwrap(), position))
      .collect();
    assert_eq!(expected, actual);
  }

  mod helpers {
    use super::{Parser, TreetimeClockArgs, TreetimeClockArgsRaw};
    use std::fs;
    use std::path::{Path, PathBuf};
    use treetime_io::nwk::NewickDialect;

    pub(super) fn clock_dialect(flags: &[&str]) -> Result<NewickDialect, String> {
      let args = ["clock", "--tree=tree.nwk", "--metadata=metadata.tsv"]
        .iter()
        .chain(flags);
      let raw = TreetimeClockArgsRaw::try_parse_from(args).map_err(|error| error.to_string())?;
      let args = TreetimeClockArgs::try_from(raw).map_err(|error| error.to_string())?;
      Ok(args.tree_dialect.dialect())
    }

    pub(super) fn write_tree(dir: &Path, text: &str) -> PathBuf {
      let path = dir.join("tree.nwk");
      fs::write(&path, text).unwrap();
      path
    }
  }
}
