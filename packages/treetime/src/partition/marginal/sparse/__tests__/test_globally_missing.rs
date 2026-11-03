#[cfg(test)]
mod tests {
  use crate::test_utils::{NUC_ALPHABET, make_nonuniform_gtr, run_sparse_marginal};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::{NwkParse, nwk_read};
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::o;

  #[test]
  fn test_internal_profiles_keep_entries_only_at_clade_level_no_data_sites() -> Result<(), Report> {
    let NwkParse {
      graph,
      branch_lengths,
      names,
      ..
    } = nwk_read(b"((A:0.1,B:0.2)X:0.1,(C:0.1,D:0.2)Y:0.3)R;".as_slice())?;
    let alignment = fasta_read(
      b">A\nAN-NCN\n>B\nGN-NCN\n>C\nAN-TC-\n>D\nAN-TG-\n".as_slice(),
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect_vec();
    let (_, reconstruction) = run_sparse_marginal(&graph, &branch_lengths, &names, &alignment, make_nonuniform_gtr()?)?;

    let no_data_entries: BTreeMap<String, Vec<usize>> = graph
      .get_nodes()
      .filter(|node| !node.is_leaf())
      .map(|node| {
        let positions = reconstruction.node_states[&node.key()]
          .profile
          .variable
          .keys()
          .copied()
          .filter(|pos| [1, 2, 3, 5].contains(pos))
          .collect_vec();
        (names[node.key()].to_owned(), positions)
      })
      .collect();

    assert_eq!(
      btreemap! {
        o!("R") => vec![3],
        o!("X") => vec![3],
        o!("Y") => vec![3],
      },
      no_data_entries,
      "columns 1, 2 and 5 hold no data in any tip; column 3, without data below X only, stays an explicit site of every internal node"
    );
    Ok(())
  }
}
