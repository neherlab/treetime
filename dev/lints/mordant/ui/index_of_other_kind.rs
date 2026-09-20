// A place indexed by names of one kind must not also be indexed by a name of another kind.

struct Lockfile {
    dependencies: Vec<&'static str>,
    // One slot per dependency, holding the package it resolved to.
    resolutions: Vec<u32>,
    packages: Vec<&'static str>,
}

fn crossed_table_is_flagged(l: &Lockfile, dep_id: u32, pkg_id: u32) -> bool {
    let dep = l.dependencies[dep_id as usize];
    let resolved = l.resolutions[dep_id as usize];
    let name = l.packages[pkg_id as usize];
    // Flagged: `resolutions` is the per-dependency table and `pkg_id` is what indexes `packages`.
    let stale = l.resolutions[pkg_id as usize];
    dep == name && resolved == stale
}

fn crossed_get_is_flagged(
    sources: &[&str],
    parts: &[u32],
    source_index: usize,
    part_index: usize,
) -> bool {
    let path = sources[source_index];
    let part = parts.get(part_index);
    let sibling = parts[part_index.saturating_sub(1)];
    // Flagged: `parts` is indexed by `part_index`; `source_index` is what indexes `sources`.
    let crossed = parts.get(source_index);
    path.is_empty() && part == crossed && sibling == 0
}

struct Graph {
    files: Vec<&'static str>,
    parts: Vec<Vec<u32>>,
    flags: Vec<bool>,
}

impl Graph {
    fn crossed_field_is_flagged(&mut self, source_index: u32, part_index: u32) {
        let _ = self.files[source_index as usize];
        let _ = self.parts[source_index as usize][part_index as usize];
        // Flagged: `self.files` is the per-source table; `part_index` indexes `self.parts[..]`.
        let _ = self.files.get_mut(part_index as usize);
    }

    // Fine: a two-level table is a different place at each level.
    fn nested_table_is_fine(&self, source_index: u32, part_index: u32) -> u32 {
        let _ = self.files[source_index as usize];
        self.parts[source_index as usize][part_index as usize]
            + self.parts[source_index as usize].as_slice()[part_index as usize]
            + u32::from(self.flags[part_index as usize])
    }

    // Fine: parallel columns share one index kind.
    fn parallel_columns_are_fine(&self, source_index: u32) -> bool {
        self.files[source_index as usize].is_empty()
            && self.parts[source_index as usize].is_empty()
            && self.flags[source_index as usize]
    }
}

struct Dep {
    entry_id: u32,
    dep_id: u32,
}

// Flagged, naming the index the home table does use: `dep_idx` and `d.dep_id` are one kind.
fn renamed_crossing_is_flagged(
    states: &[u32],
    dependencies: &[u32],
    d: &Dep,
    root_id: usize,
) -> u32 {
    let unvisited = states[root_id];
    let name = dependencies[d.dep_id as usize];
    let dep_idx = d.entry_id as usize;
    unvisited + name + states[dep_idx]
}

// Flagged though the crossing comes first: `parts` is named after `part_index`.
fn crossing_first_is_flagged(
    sources: &[u32],
    parts: &[u32],
    source_index: usize,
    part_index: usize,
) -> u32 {
    let crossed = parts[source_index];
    crossed + parts[part_index] + sources[source_index]
}

// Flagged once: `dep_id` on the table `pkg_id` indexes; `package_id` spells `pkg_id` out and stays quiet.
fn third_kind_is_flagged(
    table: &[u32],
    deps: &[u32],
    packages: &[u32],
    pkg_id: usize,
    package_id: usize,
    dep_id: usize,
) -> u32 {
    table[pkg_id]
        + table[pkg_id]
        + table[dep_id]
        + table[package_id]
        + deps[dep_id]
        + packages[package_id]
}

#[derive(Clone, Copy)]
struct FileId(usize);
#[derive(Clone, Copy)]
struct FnId(usize);
struct Db {
    files: Vec<u32>,
    fns: Vec<u32>,
}
struct Files(Vec<u32>);

impl std::ops::Index<FileId> for Db {
    type Output = u32;
    fn index(&self, i: FileId) -> &u32 {
        &self.files[i.0]
    }
}
impl std::ops::Index<FnId> for Db {
    type Output = u32;
    fn index(&self, i: FnId) -> &u32 {
        &self.fns[i.0]
    }
}
impl std::ops::Index<FileId> for Files {
    type Output = u32;
    fn index(&self, i: FileId) -> &u32 {
        &self.0[i.0]
    }
}

// Fine: typed indices already tell the kinds apart.
fn typed_indices_are_fine(db: &Db, files: &Files, file_id: FileId, fn_id: FnId) -> u32 {
    db[fn_id] + db[fn_id] + files[file_id] + db[file_id]
}

// Fine: a role name for the same kind has no table of its own name.
fn role_name_is_fine(nodes: &[u32], depths: &[u32], node_idx: usize, parent_idx: usize) -> u32 {
    nodes[node_idx] + depths[node_idx] + nodes[parent_idx] + depths[parent_idx]
}

// Fine: the two prefixes abbreviate one word.
fn abbreviation_is_fine(
    l: &Lockfile,
    dep_id: u32,
    dependency_id: u32,
    pkg_id: u32,
    package_id: u32,
) {
    let _ = l.dependencies[dep_id as usize];
    let _ = l.resolutions[dep_id as usize];
    let _ = l.resolutions[dependency_id as usize];
    let _ = l.packages[pkg_id as usize];
    let _ = l.packages[package_id as usize];
}

// Fine: two locals of one name are two places.
fn shadowed_local_is_fine(l: &Lockfile, per_package: &[u32], dep_id: u32, pkg_id: u32) -> u32 {
    let resolutions = &l.resolutions;
    let a = resolutions[dep_id as usize] + l.dependencies[dep_id as usize].len() as u32;
    let resolutions = per_package;
    a + resolutions[pkg_id as usize] + l.packages[pkg_id as usize].len() as u32
}

// Fine: names without a kind suffix claim nothing.
fn unsuffixed_names_are_fine(v: &[u8], w: &[u8], i: usize, at: usize, n: usize) -> u8 {
    v[i] + v[at] + w[at] + v[n] + w[0]
}

fn main() {
    let l = Lockfile {
        dependencies: vec!["a"],
        resolutions: vec![0],
        packages: vec!["a"],
    };
    let _ = crossed_table_is_flagged(&l, 0, 0);
    let _ = crossed_get_is_flagged(&["a"], &[0], 0, 0);
    let mut g = Graph {
        files: vec!["a"],
        parts: vec![vec![0]],
        flags: vec![false],
    };
    g.crossed_field_is_flagged(0, 0);
    let _ = g.nested_table_is_fine(0, 0);
    let _ = g.parallel_columns_are_fine(0);
    let d = Dep {
        entry_id: 0,
        dep_id: 0,
    };
    let _ = renamed_crossing_is_flagged(&[0], &[0], &d, 0);
    let _ = crossing_first_is_flagged(&[0], &[0], 0, 0);
    let _ = third_kind_is_flagged(&[0], &[0], &[0], 0, 0, 0);
    let db = Db {
        files: vec![0],
        fns: vec![0],
    };
    let _ = typed_indices_are_fine(&db, &Files(vec![0]), FileId(0), FnId(0));
    let _ = role_name_is_fine(&[0], &[0], 0, 0);
    abbreviation_is_fine(&l, 0, 0, 0, 0);
    let _ = shadowed_local_is_fine(&l, &[0], 0, 0);
    let _ = unsuffixed_names_are_fine(&[0], &[0], 0, 0, 0);
}
