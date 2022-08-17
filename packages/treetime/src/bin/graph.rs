use ctor::ctor;
use eyre::Report;
use itertools::Itertools;
use std::borrow::Borrow;
use std::error::Error;
use std::fmt::{Debug, Display};
use std::time::Duration;
use treetime::graph::breadth_first::GraphTraversalContinuation;
use treetime::graph::examples::get_example_graph_with_multiple_roots;
use treetime::io::file::create_file;
use treetime::utils::global_init::global_init;

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let mut graph = get_example_graph_with_multiple_roots()?;

  #[cfg(debug_assertions)]
  {
    use parking_lot::deadlock;
    use std::thread;
    use std::time::Duration;

    // Create a background thread which checks for deadlocks every 10s
    thread::spawn(move || loop {
      thread::sleep(Duration::from_secs(5));
      let deadlocks = deadlock::check_deadlock();
      println!("{} deadlocks detected", deadlocks.len());

      if deadlocks.is_empty() {
        continue;
      }

      for (i, threads) in deadlocks.iter().enumerate() {
        println!("Deadlock #{}", i);
        for t in threads {
          println!("Thread Id {:#?}", t.thread_id());
          println!("{:#?}", t.backtrace());
        }
      }
    });
  }

  println!("Traverse forward:");
  println!(
    "{:^6} | {:^16} | {:^7} | {:^7}",
    "Node", "Parents", "Is leaf", "Is root"
  );
  graph.par_iter_breadth_first_forward(|node| {
    // Mutable access to node payload
    node.payload.name = format!("{}*", &node.payload.name);

    // Read-write access to list pairs of `(edge, parent)`, where `edge is the edge leading to the `parent`.
    // Note: order of writes is not specified. Write access is exclusive and blocks other threads.
    // Note: there is no access to children, because they are not guaranteed to be processed at this point yet.
    let parent_names = node
      .parents
      .iter()
      .map(|(parent, edge)| {
        let parent = parent.read();
        parent.name.clone()
      })
      .join(", ");

    println!(
      "{:<6} | {:<16} | {:<7} | {:<7}",
      node.payload.name, parent_names, node.is_leaf, node.is_root
    );

    GraphTraversalContinuation::Continue
  });

  graph.print_graph(create_file("tmp/graph2.dot")?)?;

  println!();

  println!("Traverse backward:");
  println!(
    "{:^6} | {:^16} | {:^7} | {:^7}",
    "Node", "Children", "Is leaf", "Is root"
  );
  graph.par_iter_breadth_first_backward(|node| {
    // Mutable access to node payload
    node.payload.name = format!("{}*", &node.payload.name);

    // Access to list pairs of `(edge, child)`, where `edge` is the edge leading to that `child`
    // Note: order of writes is not specified. Write access is exclusive and blocks other threads.
    // Note: there is no access to parents, because they are not guaranteed to be processed at this point yet.
    let child_names = node
      .children
      .iter()
      .map(|(child, edge)| {
        let child = child.read();
        child.name.clone()
      })
      .join(", ");

    println!(
      "{:<6} | {:<16} | {:<7} | {:<7}",
      node.payload.name, child_names, node.is_leaf, node.is_root
    );

    GraphTraversalContinuation::Continue
  });

  Ok(())
}
