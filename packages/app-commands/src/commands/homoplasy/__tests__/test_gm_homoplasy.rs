#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{V0Report, load_case, run_v1};
  use pretty_assertions::assert_eq;

  #[test]
  fn test_gm_homoplasy_zika_86_matches_v0() -> Result<(), Report> {
    let (inputs, expected) = load_case("zika_86_jc")?;
    let actual = V0Report::from_v1(&run_v1(&inputs)?);
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
    use crate::commands::homoplasy::result::{HomoplasyResult, MutationTable};
    use crate::commands::homoplasy::run::run_homoplasy;
    use clap::Parser;
    use eyre::Report;
    use itertools::Itertools;
    use serde::Deserialize;
    use std::cmp::Reverse;
    use std::collections::BTreeMap;
    use std::path::{Path, PathBuf};
    use std::str::FromStr;
    use tempfile::tempdir;
    use treetime::cancel::NoopCancel;
    use treetime::progress::NoopProgress;
    use treetime::seq::mutation::Sub;
    use treetime_utils::fmt::float::float_to_exponential;
    use treetime_utils::io::json::json_read_file;

    const V0_ROOT_BRANCH_LENGTH: f64 = 0.001;

    #[derive(Debug, Deserialize)]
    pub(super) struct GmInputs {
      tree: String,
      aln: String,
      n: usize,
    }

    #[derive(Debug, PartialEq, Eq, Deserialize)]
    pub(super) struct V0Report {
      total_branch_length: String,
      mutations: usize,
      multiplicities: BTreeMap<String, usize>,
      terminal_branch_length: String,
      terminal_mutations: usize,
      terminal_multiplicities: BTreeMap<String, usize>,
      genome_length: usize,
      site_hits: Vec<V0SiteHits>,
      log_likelihood_difference: String,
      top_mutations: Vec<V0Row>,
      top_terminal_mutations: Vec<V0Row>,
      taxa: Vec<V0Row>,
    }

    impl V0Report {
      pub(super) fn from_v1(result: &HomoplasyResult) -> Self {
        let substitutions = &result.substitutions;
        Self {
          total_branch_length: float_to_exponential(substitutions.total_branch_length + V0_ROOT_BRANCH_LENGTH, 3),
          mutations: substitutions.all.mutations,
          multiplicities: multiplicities(&substitutions.all),
          terminal_branch_length: float_to_exponential(substitutions.terminal_branch_length, 3),
          terminal_mutations: substitutions.terminal.mutations,
          terminal_multiplicities: multiplicities(&substitutions.terminal),
          genome_length: substitutions.genome_length,
          site_hits: substitutions
            .site_hits
            .iter()
            .filter(|row| row.sites > 0)
            .map(|row| V0SiteHits {
              hits: row.hits,
              sites: row.sites,
              expected: format!("{:.2}", row.expected),
            })
            .collect(),
          log_likelihood_difference: float_to_exponential(substitutions.log_likelihood_difference, 3),
          top_mutations: recurrent(&substitutions.all),
          top_terminal_mutations: recurrent(&substitutions.terminal),
          taxa: result
            .taxa
            .iter()
            .filter(|taxon| !taxon.homoplasic_mutations.is_empty())
            .map(|taxon| V0Row {
              name: taxon.name.clone(),
              count: taxon.homoplasic_mutations.len(),
            })
            .collect(),
        }
      }

      fn with_v1_tie_order(mut self) -> Result<Self, Report> {
        self.top_mutations = sorted_by_position(self.top_mutations)?;
        self.top_terminal_mutations = sorted_by_position(self.top_terminal_mutations)?;
        self.taxa = self
          .taxa
          .into_iter()
          .sorted_by(|a, b| b.count.cmp(&a.count).then_with(|| a.name.cmp(&b.name)))
          .collect();
        Ok(self)
      }
    }

    #[derive(Debug, PartialEq, Eq, Deserialize)]
    struct V0SiteHits {
      hits: usize,
      sites: usize,
      expected: String,
    }

    #[derive(Debug, PartialEq, Eq, Deserialize)]
    struct V0Row {
      name: String,
      count: usize,
    }

    pub(super) fn load_case(name: &str) -> Result<(GmInputs, V0Report), Report> {
      let mut inputs: BTreeMap<String, GmInputs> = json_read_file(fixtures().join("gm_homoplasy_inputs.json"))?;
      let mut outputs: BTreeMap<String, V0Report> = json_read_file(fixtures().join("gm_homoplasy_outputs.json"))?;
      Ok((inputs.remove(name).unwrap(), outputs.remove(name).unwrap().with_v1_tie_order()?))
    }

    pub(super) fn run_v1(inputs: &GmInputs) -> Result<HomoplasyResult, Report> {
      let out = tempdir()?;
      let root = workspace_root();
      let raw = TreetimeHomoplasyArgsRaw::try_parse_from([
        "homoplasy",
        &format!("--tree={}", root.join(&inputs.tree).display()),
        &format!("--alignment={}", root.join(&inputs.aln).display()),
        "--model=jc69",
        "--dense=true",
        "--zero-based",
        "--detailed",
        &format!("-n={}", inputs.n),
        &format!("--output-homoplasy-stats={}", out.path().join("stats.json").display()),
      ])?;
      let args = TreetimeHomoplasyArgs::try_from(raw)?;
      run_homoplasy(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }

    fn multiplicities(table: &MutationTable) -> BTreeMap<String, usize> {
      table
        .multiplicities
        .iter()
        .map(|row| (row.branches.to_string(), row.mutations))
        .collect()
    }

    fn recurrent(table: &MutationTable) -> Vec<V0Row> {
      table
        .ranked
        .iter()
        .take_while(|mutation| mutation.multiplicity > 1)
        .map(|mutation| V0Row {
          name: mutation.mutation.clone(),
          count: mutation.multiplicity,
        })
        .collect()
    }

    fn sorted_by_position(rows: Vec<V0Row>) -> Result<Vec<V0Row>, Report> {
      let keyed: Vec<(Reverse<usize>, usize, V0Row)> = rows
        .into_iter()
        .map(|row| Ok((Reverse(row.count), Sub::from_str(&row.name)?.pos(), row)))
        .collect::<Result<_, Report>>()?;
      Ok(
        keyed
          .into_iter()
          .sorted_by(|a, b| (a.0, a.1, &a.2.name).cmp(&(b.0, b.1, &b.2.name)))
          .map(|(_, _, row)| row)
          .collect(),
      )
    }

    fn workspace_root() -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    fn fixtures() -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("src/commands/homoplasy/__tests__/__fixtures__")
    }
  }
}
