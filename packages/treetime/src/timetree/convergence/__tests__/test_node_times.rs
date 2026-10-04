#[cfg(test)]
mod tests {
  use crate::test_utils::{empty_time_inference, find_node_key_by_name};
  use crate::timetree::convergence::node_times::capture_node_times;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read;

  #[test]
  fn test_capture_node_times_keeps_only_finite_committed_times() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:1.0,B:1.0)I:1.0,C:1.0)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let mut inference = empty_time_inference(&graph);
    for (name, time) in [
      ("A", Some(2010.0)),
      ("B", Some(f64::NAN)),
      ("C", Some(f64::INFINITY)),
      ("I", None),
      ("root", Some(2000.0)),
    ] {
      inference
        .posterior
        .get_mut(&key(name))
        .expect("every fixture node has a posterior")
        .time = time;
    }

    let actual = capture_node_times(&graph, &inference);

    let expected = btreemap! {
      key("A") => 2010.0,
      key("root") => 2000.0,
    };
    assert_eq!(expected, actual);
    Ok(())
  }
}
