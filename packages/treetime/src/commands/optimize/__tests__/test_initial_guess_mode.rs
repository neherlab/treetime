#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::InitialGuessMode;
  use crate::commands::optimize::run::should_run_initial_guess;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NWK: &str = "((A:0.1,B:0.2)AB:0.05,C:0.3)root:0.01;";

  #[test]
  fn test_initial_guess_mode_always_returns_true() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NWK)?;
    assert!(should_run_initial_guess(InitialGuessMode::Always, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_returns_false() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NWK)?;
    assert!(!should_run_initial_guess(InitialGuessMode::Never, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_skips_when_all_lengths_present() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NWK)?;
    assert!(!should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_runs_when_one_length_missing() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NWK)?;
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(None);
    assert!(should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_runs_when_all_lengths_missing() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NWK)?;
    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(None);
    }
    assert!(should_run_initial_guess(InitialGuessMode::Auto, &graph));
    Ok(())
  }
}
