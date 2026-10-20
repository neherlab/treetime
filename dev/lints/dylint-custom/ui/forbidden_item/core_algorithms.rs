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

#[derive(Clone, Debug)]
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

macro_rules! tree_of {
    ($roots:expr) => {
        GraphSchedule::new($roots).into_tree()
    };
}

macro_rules! tree_type {
    () => {
        TreeSchedule<'_>
    };
}

pub fn root_of_local_macro(roots: &[usize]) -> Option<usize> {
    let tree: tree_type!() = tree_of!(roots)?;
    tree.root()
}

pub fn root_of_external_macro(roots: &[usize]) -> Option<usize> {
    graph::tree_of_roots!(roots)?.root()
}
