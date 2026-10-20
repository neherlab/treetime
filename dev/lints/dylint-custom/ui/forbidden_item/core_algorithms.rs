// aux-build:graph.rs
#![crate_type = "lib"]

extern crate graph;

use graph::schedule::{GraphSchedule, TreeSchedule};

pub fn root_of_signature(tree: &TreeSchedule<'_>) -> Option<usize> {
    tree.root()
}

pub fn root_of_let_type(roots: &[usize]) -> Option<usize> {
    let tree: TreeSchedule<'_> = GraphSchedule::new(roots).into_tree()?;
    tree.root()
}

pub fn root_of_method_reference(roots: &[usize]) -> Option<usize> {
    let trees: Vec<_> = Some(GraphSchedule::new(roots))
        .into_iter()
        .map(GraphSchedule::into_tree)
        .collect::<Option<_>>()?;
    trees.first().and_then(Holder::root_of_tree)
}

pub fn root_count_of_graph(roots: &[usize]) -> usize {
    GraphSchedule::new(roots).roots().len()
}

pub struct Holder<'g> {
    pub tree: TreeSchedule<'g>,
    pub trees: Vec<TreeSchedule<'g>>,
}

impl Holder<'_> {
    pub fn root_of_tree(tree: &TreeSchedule<'_>) -> Option<usize> {
        tree.root()
    }
}

pub trait Rooted {
    fn rooted(&self) -> bool;
}

impl Rooted for TreeSchedule<'_> {
    fn rooted(&self) -> bool {
        self.root().is_some()
    }
}
