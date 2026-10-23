#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{MapRun, run};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::progress::{RunWarning, RunWarningKind};
  use treetime_io::auspice_types::AuspiceGeoResolution;
  use treetime_io::lat_longs_tsv::Coordinates;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_lat_longs_give_the_auspice_tree_a_map_of_the_located_trait_values() -> Result<(), Report> {
    let actual = run("country\tUSA\t38\t-97\ncountry\tgermany\t51\t10\n")?;

    let expected = MapRun {
      lat_longs: actual.lat_longs.clone(),
      geo_resolutions: vec![AuspiceGeoResolution {
        key: o!("country"),
        title: None,
        demes: btreemap! {
          o!("germany") => Coordinates { latitude: 51.0, longitude: 10.0 },
          o!("usa") => Coordinates { latitude: 38.0, longitude: -97.0 },
        },
      }],
      panels: vec_of_owned!["tree", "map"],
      warnings: vec![RunWarning {
        kind: RunWarningKind::LocationsWithoutCoordinates,
        message: format!(
          "The coordinates file '{}' has no coordinates for the country values peru. The map of the Auspice tree leaves them out.",
          actual.lat_longs
        ),
        names: vec_of_owned!["peru"],
      }],
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rstest]
  #[case::other_resolution("region\tusa\t38\t-97\n")]
  #[case::other_values("country\tfrance\t46\t2\n")]
  #[trace]
  fn test_lat_longs_without_trait_values_give_no_map(#[case] lat_longs: &str) -> Result<(), Report> {
    let actual = run(lat_longs)?;

    let expected = MapRun {
      lat_longs: actual.lat_longs.clone(),
      geo_resolutions: vec![],
      panels: vec_of_owned!["tree"],
      warnings: vec![RunWarning {
        kind: RunWarningKind::LocationsWithoutCoordinates,
        message: format!(
          "The coordinates file '{}' has no coordinates for the country values germany, peru, usa. The Auspice tree has no map.",
          actual.lat_longs
        ),
        names: vec_of_owned!["germany", "peru", "usa"],
      }],
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_support::tests::RecordingLog;
    use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
    use crate::commands::mugration::run::run_mugration;
    use crate::commands::shared::metadata::MetadataArgs;
    use crate::commands::shared::output_args::OutputCoreArgs;
    use crate::commands::shared::tree_input::TreeArgs;
    use eyre::{Report, WrapErr};
    use std::fs;
    use tempfile::tempdir;
    use treetime::cancel::NoopCancel;
    use treetime::progress::{NoopProgress, RunWarning};
    use treetime_io::auspice_types::{AuspiceGeoResolution, AuspiceTree};
    use treetime_utils::io::json::json_read_file;
    use treetime_utils::o;

    const TREE: &str = "((A:0.1,B:0.1)AB:0.1,C:0.1)root;";
    const METADATA: &str = "strain\tcountry\nA\tusa\nB\tgermany\nC\tperu\n";

    #[derive(Debug, PartialEq)]
    pub(super) struct MapRun {
      pub(super) lat_longs: String,
      pub(super) geo_resolutions: Vec<AuspiceGeoResolution>,
      pub(super) panels: Vec<String>,
      pub(super) warnings: Vec<RunWarning>,
    }

    pub(super) fn run(lat_longs: &str) -> Result<MapRun, Report> {
      let dir = tempdir().wrap_err("When creating a temporary directory")?;
      fs::write(dir.path().join("tree.nwk"), TREE).wrap_err("When writing the tree fixture")?;
      fs::write(dir.path().join("metadata.tsv"), METADATA).wrap_err("When writing the metadata fixture")?;
      let lat_longs_path = dir.path().join("lat_longs.tsv");
      fs::write(&lat_longs_path, lat_longs).wrap_err("When writing the coordinates fixture")?;
      let auspice = dir.path().join("tree.json");
      let args = TreetimeMugrationArgs::try_from(TreetimeMugrationArgsRaw {
        tree: TreeArgs {
          tree: Some(dir.path().join("tree.nwk")),
        },
        metadata: MetadataArgs {
          metadata: Some(dir.path().join("metadata.tsv")),
        },
        attribute: Some(o!("country")),
        lat_longs: Some(lat_longs_path.clone()),
        output: OutputCoreArgs {
          output_tree_auspice: Some(auspice.clone()),
          ..OutputCoreArgs::default()
        },
        ..TreetimeMugrationArgsRaw::default()
      })?;
      let log = RecordingLog::default();

      run_mugration(&args, &NoopCancel, &NoopProgress, &log)?;

      let tree: AuspiceTree = json_read_file(&auspice)?;
      Ok(MapRun {
        lat_longs: lat_longs_path.display().to_string(),
        geo_resolutions: tree.data.meta.geo_resolutions,
        panels: tree.data.meta.panels,
        warnings: log.into_warnings(),
      })
    }
  }
}
