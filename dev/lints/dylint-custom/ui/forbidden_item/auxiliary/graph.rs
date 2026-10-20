pub mod schedule {
    pub struct GraphSchedule<'g> {
        roots: &'g [usize],
    }

    impl<'g> GraphSchedule<'g> {
        pub fn new(roots: &'g [usize]) -> Self {
            Self { roots }
        }

        pub fn roots(&self) -> &'g [usize] {
            self.roots
        }

        pub fn into_tree(self) -> Option<TreeSchedule<'g>> {
            (self.roots.len() == 1).then_some(TreeSchedule { schedule: self })
        }
    }

    pub struct TreeSchedule<'g> {
        schedule: GraphSchedule<'g>,
    }

    impl TreeSchedule<'_> {
        pub fn root(&self) -> Option<usize> {
            self.schedule.roots.first().copied()
        }
    }
}
