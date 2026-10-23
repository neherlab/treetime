#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{project_root, write_column_blocks};
  use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
  use crate::commands::prune::run::run_prune;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output_args::{OutputCoreArgs, PruneOutputSelection};
  use eyre::Report;
  use helpers::*;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::path::{Path, PathBuf};
  use tempfile::TempDir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_io::usher_mat::UsherTree;
  use treetime_utils::io::json::json_read_file;

  #[test]
  fn test_prune_auspice_and_mat_carry_the_same_per_edge_mutations() -> Result<(), Report> {
    let outdir = TempDir::new()?;
    run_flu_20_prune(outdir.path(), true)?;
    let auspice = auspice_mutations(&read_auspice(outdir.path())?);
    let mat = mat_mutations(&json_read_file::<UsherTree>(outdir.path().join("prune.mat.json"))?)?;
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

  #[test]
  fn test_prune_two_alignments_write_the_mutations_of_each_alignment_to_its_own_files() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let alignment = project_root().join("data/zika/20/aln.fasta.xz");
    let blocks = write_column_blocks(&alignment, ZIKA_CUT, dir.path())?;
    let whole_out = dir.path().join("whole");
    let split_out = dir.path().join("split");
    run_zika_prune(vec![alignment], &whole_out)?;
    run_zika_prune(blocks.to_vec(), &split_out)?;

    let whole = auspice_mutations(&json_read_file(whole_out.join("prune.auspice.json"))?);
    let [first, second] = [1, 2].map(|alignment| {
      let auspice = auspice_mutations(&json_read_file(
        split_out.join(format!("prune.{alignment}.auspice.json")),
      )?);
      let mat = mat_mutations(&json_read_file(split_out.join(format!("prune.{alignment}.mat.json")))?)?;
      eyre::ensure!(
        auspice == mat,
        "Auspice and MAT mutations of alignment {alignment} differ"
      );
      Ok::<_, Report>(auspice)
    });
    let (first, second) = (first?, second?);
    let joined = first
      .into_iter()
      .map(|(name, mutations)| {
        let shifted = second[&name].iter().map(|mutation| shift_position(mutation, ZIKA_CUT));
        let mutations = mutations
          .into_iter()
          .map(Ok)
          .chain(shifted)
          .collect::<Result<BTreeSet<_>, Report>>()?;
        Ok((name, mutations))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?;

    assert_eq!(whole, joined);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::output_args::PruneOutputSelectionArgs;
    use crate::commands::shared::tree_input::TreeArgs;

    pub(super) const ZIKA_CUT: usize = 5000;

    pub(super) fn run_zika_prune(alignment: Vec<PathBuf>, outdir: &Path) -> Result<(), Report> {
      let root = project_root();
      let args = TreetimePruneArgs::try_from(TreetimePruneArgsRaw {
        alignment: AlignmentArgs { alignment },
        tree: TreeArgs {
          tree: Some(root.join("data/zika/20/tree.nwk")),
        },
        collapse_empty: true,
        output: OutputCoreArgs {
          output_all: Some(outdir.to_path_buf()),
          ..OutputCoreArgs::default()
        },
        output_selection: PruneOutputSelectionArgs {
          output_selection: vec![PruneOutputSelection::Auspice, PruneOutputSelection::MatJson],
        },
        ..TreetimePruneArgsRaw::default()
      })?;
      run_prune(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }

    pub(super) fn shift_position(mutation: &str, offset: usize) -> Result<String, Report> {
      let (from, rest) = mutation.split_at(1);
      let (position, to) = rest.split_at(rest.len() - 1);
      Ok(format!("{from}{}{to}", position.parse::<usize>()? + offset))
    }

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
        tree: TreeArgs {
          tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
        },
        collapse_empty: with_alignment,
        collapse_short: (!with_alignment).then_some(1e-3),
        output: OutputCoreArgs {
          output_all: Some(outdir.to_path_buf()),
          ..OutputCoreArgs::default()
        },
        output_selection: PruneOutputSelectionArgs {
          output_selection: vec![
            PruneOutputSelection::Auspice,
            PruneOutputSelection::MatJson,
            PruneOutputSelection::MatPb,
          ],
        },
        ..TreetimePruneArgsRaw::default()
      })?;
      run_prune(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;
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
  }
}
