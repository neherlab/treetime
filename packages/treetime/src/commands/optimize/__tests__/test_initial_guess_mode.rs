#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::InitialGuessMode;
  use crate::commands::optimize::run::should_run_initial_guess;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  const TREE_WITH_LENGTHS: &str = "((A:0.1,B:0.2)AB:0.05,C:0.3)root:0.01;";

  /// Bio crate newick parser returns `f32::NAN` for missing branch lengths,
  /// which becomes `Some(NaN)` after cast and wrapping.
  const TREE_WITHOUT_LENGTHS: &str = "((A,B)AB,C)root;";

  #[test]
  fn test_initial_guess_mode_default_is_auto() {
    assert_eq!(InitialGuessMode::Auto, InitialGuessMode::default());
  }

  #[test]
  fn test_initial_guess_mode_always_runs() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    assert!(should_run_initial_guess(InitialGuessMode::Always, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_skips() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    assert!(!should_run_initial_guess(InitialGuessMode::Never, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_skips_when_all_lengths_present() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    assert!(!should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_runs_when_one_length_none() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(None);
    assert!(should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_runs_when_all_lengths_none() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(None);
    }
    assert!(should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_detects_nan_from_newick() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITHOUT_LENGTHS)?;
    // Bio crate stores missing weights as Some(NaN), not None.
    // Auto mode must detect NaN as a missing branch length.
    assert!(should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_detects_explicit_nan() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(f64::NAN));
    assert!(should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }
}
