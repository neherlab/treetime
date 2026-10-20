// aux-build:graph.rs
// compile-flags: --test --crate-name core_algorithms

extern crate graph;

use graph::schedule::GraphSchedule;

pub fn root_of_signature(tree: &graph::schedule::TreeSchedule<'_>) -> Option<usize> {
    tree.root()
}

#[test]
fn root_of_test_function() {
    let roots = [0];
    let tree: graph::schedule::TreeSchedule<'_> = GraphSchedule::new(&roots).into_tree().unwrap();
    assert_eq!(Some(0), root_of_signature(&tree));
}

#[cfg(test)]
mod tests {
    use graph::schedule::{GraphSchedule, TreeSchedule};

    fn root_of_helper(tree: &TreeSchedule<'_>) -> Option<usize> {
        tree.root()
    }

    #[test]
    fn root_of_test_module() {
        let roots = [0];
        let tree: TreeSchedule<'_> = GraphSchedule::new(&roots).into_tree().unwrap();
        assert_eq!(Some(0), root_of_helper(&tree));
    }
}
