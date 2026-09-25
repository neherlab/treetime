#[cfg(test)]
mod tests {
  use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
  use crate::commands::prune::run::run_prune;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output_args::{OutputCoreArgs, PruneOutputSelection};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::path::{Path, PathBuf};
  use tempfile::TempDir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_io::usher_mat::UsherTree;
  use treetime_utils::io::json::json_read_file;

  use helpers::*;

  #[test]
  fn test_prune_auspice_and_mat_carry_the_same_per_edge_mutations() -> Result<(), Report> {
    let outdir = TempDir::new()?;
    run_flu_20_prune(outdir.path(), true)?;
    let auspice = auspice_mutations(&read_auspice(outdir.path())?);
    let mat = mat_mutations(&json_read_file::<UsherTree, _>(outdir.path().join("prune.mat.json"))?)?;
    assert_eq!(auspice, mat);
    Ok(())
  }

  #[test]
  fn test_prune_mutation_count_equals_fitch_parsimony_score() -> Result<(), Report> {
    let outdir = TempDir::new()?;
    run_flu_20_prune(outdir.path(), true)?;
    let auspice = auspice_mutations(&read_auspice(outdir.path())?);
    let total: usize = auspice.values().map(BTreeSet::len).sum();
    assert_eq!(FLU_H3N2_20_FITCH_SCORE, total);
    Ok(())
  }

  #[test]
  fn test_prune_auspice_without_alignment_has_no_mutations() -> Result<(), Report> {
    let outdir = TempDir::new()?;
    run_flu_20_prune(outdir.path(), false)?;
    let tree = read_auspice(outdir.path())?;
    let mutated = auspice_mutations(&tree)
      .into_values()
      .filter(|mutations| !mutations.is_empty())
      .count();
    assert_eq!((0, None), (mutated, tree.data.root_sequence));
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) const FLU_H3N2_20_FITCH_SCORE: usize = 179;

    const MAT_NUCLEOTIDES: [char; 4] = ['A', 'C', 'G', 'T'];

    pub(super) fn run_flu_20_prune(outdir: &Path, with_alignment: bool) -> Result<(), Report> {
      let root = project_root();
      let alignment = if with_alignment {
        vec![root.join("data/flu/h3n2/20/aln.fasta.xz")]
      } else {
        vec![]
      };
      let args = TreetimePruneArgs::try_from(TreetimePruneArgsRaw {
        alignment: AlignmentArgs { alignment },
        tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
        prune_empty: with_alignment,
        prune_short: (!with_alignment).then_some(1e-3),
        output: OutputCoreArgs {
          output_all: Some(outdir.to_path_buf()),
          ..OutputCoreArgs::default()
        },
        output_selection: vec![
          PruneOutputSelection::Auspice,
          PruneOutputSelection::MatJson,
          PruneOutputSelection::MatPb,
        ],
        ..TreetimePruneArgsRaw::default()
      })?;
      run_prune(&args, &NoopCancel, &NoopProgress)?;
      Ok(())
    }

    pub(super) fn read_auspice(outdir: &Path) -> Result<AuspiceTree, Report> {
      json_read_file(outdir.join("prune.auspice.json"))
    }

    pub(super) fn auspice_mutations(tree: &AuspiceTree) -> BTreeMap<String, BTreeSet<String>> {
      let mut mutations = BTreeMap::new();
      collect_auspice_mutations(&tree.tree, &mut mutations);
      mutations
    }

    pub(super) fn mat_mutations(tree: &UsherTree) -> Result<BTreeMap<String, BTreeSet<String>>, Report> {
      eyre::ensure!(
        tree.condensed_nodes.len() == tree.node_mutations.len(),
        "MAT has {} nodes but {} mutation lists",
        tree.condensed_nodes.len(),
        tree.node_mutations.len()
      );
      tree
        .condensed_nodes
        .iter()
        .zip(&tree.node_mutations)
        .map(|(node, list)| {
          let mutations = list
            .mutation
            .iter()
            .map(|mutation| {
              let [mut_nuc] = mutation.mut_nuc.as_slice() else {
                eyre::bail!(
                  "MAT mutation at {} has {} target states",
                  mutation.position,
                  mutation.mut_nuc.len()
                );
              };
              Ok(format!(
                "{}{}{}",
                mat_nucleotide(mutation.par_nuc)?,
                mutation.position,
                mat_nucleotide(*mut_nuc)?
              ))
            })
            .collect::<Result<BTreeSet<_>, Report>>()?;
          Ok((node.node_name.clone(), mutations))
        })
        .collect()
    }

    fn collect_auspice_mutations(node: &AuspiceTreeNode, mutations: &mut BTreeMap<String, BTreeSet<String>>) {
      let nuc = node
        .branch_attrs
        .mutations
        .get("nuc")
        .map(|list| list.iter().cloned().collect())
        .unwrap_or_default();
      mutations.insert(node.name.clone(), nuc);
      for child in &node.children {
        collect_auspice_mutations(child, mutations);
      }
    }

    fn mat_nucleotide(index: i32) -> Result<char, Report> {
      usize::try_from(index)
        .ok()
        .and_then(|index| MAT_NUCLEOTIDES.get(index).copied())
        .ok_or_else(|| eyre::eyre!("MAT nucleotide index {index} is out of range"))
    }

    fn project_root() -> PathBuf {
      PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .map(PathBuf::from)
        .expect("project has workspace root")
    }
  }
}
