#[cfg(test)]
mod tests {
  use crate::homoplasy::pipeline::{IndelKey, SiteBranches};
  use crate::homoplasy::recurrence::Recurrence;
  use crate::seq::indel::InDelKind;
  use crate::test_utils::deletion;
  use eyre::Report;
  use helpers::{Scenario, sub};
  use maplit::btreemap;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use treetime_primitives::Seq;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[test]
  fn test_homoplasy_ranks_substitutions_by_branches_then_position() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_substitutions(&[
      ("root", "AB", &["A10G", "C20T"]),
      ("CD", "C", &["A10G"]),
      ("CD", "D", &["C20T", "G30A"]),
      ("AB", "A", &["T40C"]),
    ])?;
    let output = scenario.run(0)?;

    let all = &output.substitutions.all;
    assert_eq!(6, all.count);
    assert_eq!(btreemap! {1 => 2, 2 => 2}, all.histogram);
    assert_eq!(
      vec![
        Recurrence {
          mutation: sub("A10G"),
          branches: scenario.nodes(&["AB", "C"]),
        },
        Recurrence {
          mutation: sub("C20T"),
          branches: scenario.nodes(&["AB", "D"]),
        },
        Recurrence {
          mutation: sub("G30A"),
          branches: scenario.nodes(&["D"]),
        },
        Recurrence {
          mutation: sub("T40C"),
          branches: scenario.nodes(&["A"]),
        },
      ],
      all.ranked
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_terminal_table_counts_leaf_branches_only() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_substitutions(&[
      ("root", "AB", &["A10G", "C20T"]),
      ("CD", "C", &["A10G"]),
      ("CD", "D", &["C20T", "G30A"]),
      ("AB", "A", &["T40C"]),
    ])?;
    let output = scenario.run(0)?;

    let terminal = &output.substitutions.terminal;
    assert_eq!(4, terminal.count);
    assert_eq!(btreemap! {1 => 4}, terminal.histogram);
    let ranked: Vec<_> = terminal
      .ranked
      .iter()
      .map(|recurrence| recurrence.mutation.clone())
      .collect();
    assert_eq!(vec![sub("A10G"), sub("C20T"), sub("G30A"), sub("T40C")], ranked);
    Ok(())
  }

  #[test]
  fn test_homoplasy_site_histogram_matches_poisson_formula() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_substitutions(&[
      ("root", "AB", &["A10G", "C20T"]),
      ("CD", "C", &["A10G"]),
      ("CD", "D", &["C20T", "G30A"]),
      ("AB", "A", &["T40C"]),
    ])?;
    let sites = scenario.run(0)?.substitutions.sites;

    assert_eq!(50, sites.genome_length);
    assert_eq!(
      vec![46, 2, 2],
      sites.rows.iter().map(|row| row.sites).collect::<Vec<_>>()
    );
    pretty_assert_abs_diff_eq!(
      array![44.346_021_835_857_876, 5.321_522_620_302_945, 0.319_291_357_218_176_7],
      sites.rows.iter().map(|row| row.expected).collect::<Array1<f64>>(),
      epsilon = 1e-10
    );
    pretty_assert_abs_diff_eq!(
      -1.140_831_793_246_832_1,
      sites.log_likelihood_difference,
      epsilon = 1e-10
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_constant_sites_extend_the_genome() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_substitutions(&[
      ("root", "AB", &["A10G", "C20T"]),
      ("CD", "C", &["A10G"]),
      ("CD", "D", &["C20T", "G30A"]),
      ("AB", "A", &["T40C"]),
    ])?;
    let sites = scenario.run(10)?.substitutions.sites;

    assert_eq!(60, sites.genome_length);
    assert_eq!(
      vec![56, 2, 2],
      sites.rows.iter().map(|row| row.sites).collect::<Vec<_>>()
    );
    pretty_assert_abs_diff_eq!(
      array![54.290_245_082_157_57, 5.429_024_508_215_757_5, 0.271_451_225_410_787_9],
      sites.rows.iter().map(|row| row.expected).collect::<Array1<f64>>(),
      epsilon = 1e-10
    );
    pretty_assert_abs_diff_eq!(
      -1.181_185_129_015_297_3,
      sites.log_likelihood_difference,
      epsilon = 1e-10
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_expected_counts_for_rate_one_half() -> Result<(), Report> {
    let scenario = Scenario::with_length(10)?.with_substitutions(&[
      ("root", "AB", &["A1G"]),
      ("root", "CD", &["A2G"]),
      ("AB", "A", &["A3G"]),
      ("AB", "B", &["A4G"]),
      ("CD", "C", &["A5G"]),
    ])?;
    let sites = scenario.run(0)?.substitutions.sites;

    assert_eq!(vec![5, 5], sites.rows.iter().map(|row| row.sites).collect::<Vec<_>>());
    pretty_assert_abs_diff_eq!(
      array![10.0 * (-0.5_f64).exp(), 5.0 * (-0.5_f64).exp()],
      sites.rows.iter().map(|row| row.expected).collect::<Array1<f64>>(),
      epsilon = 1e-10
    );
    pretty_assert_abs_diff_eq!(
      0.809_018_614_972_192_3,
      sites.log_likelihood_difference,
      epsilon = 1e-10
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_without_substitutions_has_all_sites_unhit() -> Result<(), Report> {
    let output = Scenario::new()?.run(5)?;

    assert_eq!(0, output.substitutions.all.count);
    assert_eq!(1, output.substitutions.sites.rows.len());
    assert_eq!(55, output.substitutions.sites.rows[0].sites);
    pretty_assert_abs_diff_eq!(55.0, output.substitutions.sites.rows[0].expected, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(
      0.0,
      output.substitutions.sites.log_likelihood_difference,
      epsilon = 1e-10
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_rejects_an_empty_genome() -> Result<(), Report> {
    let result = Scenario::with_length(0)?.run(0);

    assert_error!(
      result,
      "The homoplasy statistics need at least one site, but the alignment is empty and no constant sites were given"
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_tree_lengths() -> Result<(), Report> {
    let output = Scenario::new()?.run(0)?;

    pretty_assert_abs_diff_eq!(2.1, output.substitutions.total_branch_length, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(
      0.846_870_397_798_857_4,
      output.substitutions.terminal_branch_length,
      epsilon = 1e-10
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_leaf_lists_terminal_substitutions_at_sites_hit_elsewhere() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_substitutions(&[
      ("root", "AB", &["A10G", "C20T"]),
      ("CD", "C", &["A10G"]),
      ("CD", "D", &["C20T", "G30A"]),
      ("AB", "A", &["T40C"]),
    ])?;
    let output = scenario.run(0)?;

    let expected = btreemap! {
      scenario.node("C") => vec![sub("A10G")],
      scenario.node("D") => vec![sub("C20T")],
    };
    assert_eq!(expected, output.substitutions.leaves);
    Ok(())
  }

  #[test]
  fn test_homoplasy_leaf_counts_a_reversion() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_substitutions(&[("root", "AB", &["A10G"]), ("AB", "A", &["G10A"])])?;
    let output = scenario.run(0)?;

    assert_eq!(
      btreemap! {scenario.node("A") => vec![sub("G10A")]},
      output.substitutions.leaves
    );
    assert_eq!(btreemap! {1 => 2}, output.substitutions.all.histogram);
    Ok(())
  }

  #[test]
  fn test_homoplasy_puts_ambiguous_changes_in_their_own_lists() -> Result<(), Report> {
    let scenario = Scenario::new()?
      .with_raw_substitutions(&[("AB", "B", &["G5R", "A6N"]), ("CD", "C", &["G5R"])])?
      .with_bridged_substitutions(&[("AB", "B", &["G5R"]), ("CD", "C", &["G5R"])])?;
    let output = scenario.run(0)?;

    assert_eq!(0, output.substitutions.all.count);
    assert_eq!(3, output.ambiguous.all.count);
    assert_eq!(
      vec![
        Recurrence {
          mutation: sub("G5R"),
          branches: scenario.nodes(&["B", "C"]),
        },
        Recurrence {
          mutation: sub("A6N"),
          branches: scenario.nodes(&["B"]),
        },
      ],
      output.ambiguous.all.ranked
    );
    assert_eq!(
      vec![
        SiteBranches {
          position: 4,
          branches: 2,
        },
        SiteBranches {
          position: 5,
          branches: 1,
        },
      ],
      output.ambiguous.sites
    );
    assert_eq!(
      btreemap! {scenario.node("B") => 2, scenario.node("C") => 1},
      output.ambiguous.leaves
    );
    Ok(())
  }

  #[test]
  fn test_homoplasy_counts_a_change_through_an_unknown_ancestor_once() -> Result<(), Report> {
    let scenario = Scenario::new()?
      .with_raw_substitutions(&[("root", "CD", &["A7N"]), ("CD", "C", &["N7G"])])?
      .with_bridged_substitutions(&[("CD", "C", &["A7G"])])?;
    let output = scenario.run(0)?;

    assert_eq!(
      vec![Recurrence {
        mutation: sub("A7G"),
        branches: scenario.nodes(&["C"]),
      }],
      output.substitutions.all.ranked
    );
    let ambiguous: Vec<_> = output
      .ambiguous
      .all
      .ranked
      .iter()
      .map(|recurrence| recurrence.mutation.clone())
      .collect();
    assert_eq!(vec![sub("A7N"), sub("N7G")], ambiguous);
    Ok(())
  }

  #[test]
  fn test_homoplasy_counts_a_recurrent_deletion() -> Result<(), Report> {
    let scenario = Scenario::new()?.with_raw_indels(&[
      ("root", "AB", deletion((2, 4), Seq::try_from_slice(b"AC")?)),
      ("CD", "D", deletion((2, 4), Seq::try_from_slice(b"AC")?)),
      ("CD", "C", deletion((6, 7), Seq::try_from_slice(b"G")?)),
    ])?;
    let output = scenario.run(0)?;

    let recurrent = IndelKey {
      kind: InDelKind::Deletion,
      range: (2, 4),
      sequence: Seq::try_from_slice(b"AC")?,
    };
    let single = IndelKey {
      kind: InDelKind::Deletion,
      range: (6, 7),
      sequence: Seq::try_from_slice(b"G")?,
    };
    assert_eq!(
      vec![
        Recurrence {
          mutation: recurrent,
          branches: scenario.nodes(&["AB", "D"]),
        },
        Recurrence {
          mutation: single,
          branches: scenario.nodes(&["C"]),
        },
      ],
      output.indels.all.ranked
    );
    assert_eq!(btreemap! {1 => 1, 2 => 1}, output.indels.all.histogram);
    assert_eq!(2, output.indels.terminal.count);
    assert_eq!(btreemap! {scenario.node("D") => 1}, output.indels.leaves);
    assert_eq!(0, output.substitutions.all.count);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::homoplasy::pipeline::{HomoplasyInput, HomoplasyOutput, HomoplasyParams, run};
    use crate::seq::indel::InDel;
    use crate::seq::mutation::{Mutation, MutationTrack, Sub};
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::{OptionExt, Report};
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use std::str::FromStr;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::{NwkParse, nwk_read};

    const NEWICK: &str = "((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)root;";

    pub(super) fn sub(text: &str) -> Sub {
      Sub::from_str(text).unwrap()
    }

    pub(super) struct Scenario {
      parse: NwkParse,
      names: BTreeMap<GraphNodeKey, Option<String>>,
      raw: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
      bridged: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
      sequence_length: usize,
    }

    impl Scenario {
      pub(super) fn new() -> Result<Self, Report> {
        Self::with_length(50)
      }

      pub(super) fn with_length(sequence_length: usize) -> Result<Self, Report> {
        let parse = nwk_read(NEWICK.as_bytes())?;
        let names = parse.names();
        let empty: BTreeMap<GraphEdgeKey, Vec<Mutation>> =
          parse.graph.get_edges().map(|edge| (edge.key(), vec![])).collect();
        Ok(Self {
          parse,
          names,
          raw: empty.clone(),
          bridged: empty,
          sequence_length,
        })
      }

      pub(super) fn with_substitutions(self, edges: &[(&str, &str, &[&str])]) -> Result<Self, Report> {
        self.with_raw_substitutions(edges)?.with_bridged_substitutions(edges)
      }

      pub(super) fn with_raw_substitutions(mut self, edges: &[(&str, &str, &[&str])]) -> Result<Self, Report> {
        for &(parent, child, subs) in edges {
          let key = self.edge(parent, child)?;
          self.raw.insert(key, substitutions(subs));
        }
        Ok(self)
      }

      pub(super) fn with_bridged_substitutions(mut self, edges: &[(&str, &str, &[&str])]) -> Result<Self, Report> {
        for &(parent, child, subs) in edges {
          let key = self.edge(parent, child)?;
          self.bridged.insert(key, substitutions(subs));
        }
        Ok(self)
      }

      pub(super) fn with_raw_indels(mut self, edges: &[(&str, &str, InDel)]) -> Result<Self, Report> {
        for (parent, child, indel) in edges {
          let key = self.edge(parent, child)?;
          let mutation = Mutation::indel(MutationTrack::Nucleotide, indel)?;
          self.raw.entry(key).or_default().push(mutation.clone());
          self.bridged.entry(key).or_default().push(mutation);
        }
        Ok(self)
      }

      pub(super) fn run(&self, constant_sites: usize) -> Result<HomoplasyOutput, Report> {
        let alphabet = Alphabet::new(AlphabetName::Nuc)?;
        let branch_lengths: BTreeMap<GraphEdgeKey, f64> = self
          .parse
          .branch_lengths
          .iter()
          .map(|(&key, length)| (key, length.unwrap_or_default()))
          .collect();
        let input = HomoplasyInput {
          graph: &self.parse.graph,
          bridged_mutations: &self.bridged,
          raw_mutations: &self.raw,
          branch_lengths: &branch_lengths,
          alphabet: &alphabet,
          sequence_length: self.sequence_length,
        };
        run(&HomoplasyParams { constant_sites }, &input).map_err(|err| err.into_report())
      }

      pub(super) fn node(&self, name: &str) -> GraphNodeKey {
        find_node_key_by_name(&self.parse.graph, &self.names, name).unwrap()
      }

      pub(super) fn nodes(&self, names: &[&str]) -> Vec<GraphNodeKey> {
        names.iter().map(|name| self.node(name)).sorted().collect()
      }

      fn edge(&self, parent: &str, child: &str) -> Result<GraphEdgeKey, Report> {
        find_edge_key(&self.parse.graph, &self.names, parent, child).ok_or_eyre("edge not found")
      }
    }

    fn substitutions(subs: &[&str]) -> Vec<Mutation> {
      subs
        .iter()
        .map(|text| Mutation::substitution(MutationTrack::Nucleotide, sub(text)))
        .collect()
    }
  }
}
