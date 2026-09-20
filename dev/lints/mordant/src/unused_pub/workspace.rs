//! Which crate prints the findings for which library. `cargo check
//! --workspace` compiles every member as its own rustc process. A library
//! is compiled before the crates that use it, so its unused items are only
//! known once those crates are done. The report for a library therefore
//! comes from a root: a member no other member depends on, which cargo
//! compiles after everything below it. `cargo metadata` names the members
//! and the edges between them.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use serde::Deserialize;

#[derive(Deserialize)]
struct Metadata {
    packages: Vec<Package>,
    workspace_root: PathBuf,
}

#[derive(Deserialize)]
struct Package {
    name: String,
    targets: Vec<Target>,
    dependencies: Vec<Dep>,
}

#[derive(Deserialize)]
struct Target {
    name: String,
    kind: Vec<String>,
}

#[derive(Deserialize)]
struct Dep {
    name: String,
    kind: Option<String>,
}

const LIB_KINDS: &[&str] = &["lib", "rlib", "dylib", "cdylib", "staticlib", "proc-macro"];

impl Package {
    fn lib_crate(&self) -> Option<String> {
        self.targets
            .iter()
            .find(|t| t.kind.iter().any(|k| LIB_KINDS.contains(&k.as_str())))
            .map(|t| crate_name(&t.name))
    }

    fn bin_crates(&self) -> Vec<String> {
        self.targets
            .iter()
            .filter(|t| t.kind.iter().any(|k| k == "bin"))
            .map(|t| crate_name(&t.name))
            .collect()
    }
}

fn crate_name(target: &str) -> String {
    target.replace('-', "_")
}

/// What cargo is compiling this crate as.
#[derive(Clone, Copy, PartialEq)]
pub enum Kind {
    Lib,
    Bin,
}

pub struct Workspace {
    pub root: PathBuf,
    pub kind: Kind,
    /// `<target>/mordant/unused_pub/`, where every member leaves its files.
    pub dir: PathBuf,
    /// This crate's package.
    pub package: String,
    /// Package name of every member, to recognise files of removed ones.
    pub members: BTreeSet<String>,
    /// Library crate names whose findings this crate prints. Empty unless
    /// this crate is the reporting target of a root package.
    pub reports: BTreeSet<String>,
    /// Library crate name -> its package, for every member with a library.
    pub libraries: BTreeMap<String, String>,
}

/// `None` when rustc is not running under cargo in a workspace that has
/// this crate as a member target; the lint then judges the crate alone.
pub fn locate(local_crate: &str, kind: Kind) -> Option<Workspace> {
    let package = std::env::var("CARGO_PKG_NAME").ok()?;
    let manifest_dir = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR")?);
    let meta = metadata(&manifest_dir)?;
    let me = meta.packages.iter().find(|p| p.name == package)?;
    let member = match kind {
        Kind::Lib => me.lib_crate().as_deref() == Some(local_crate),
        Kind::Bin => me.bin_crates().iter().any(|b| b == local_crate),
    };
    if !member {
        // A build script, or a crate compiled outside its package.
        return None;
    }
    let dir = crate::baseline::target_dir(&meta.workspace_root)
        .join("mordant")
        .join("unused_pub");
    let graph = Graph::new(&meta.packages);
    let reports = if graph.reporting_crate(me) == Some(local_crate.to_string()) {
        graph.owned_by(&package)
    } else {
        BTreeSet::new()
    };
    Some(Workspace {
        root: meta.workspace_root,
        kind,
        dir,
        package,
        members: meta.packages.iter().map(|p| p.name.clone()).collect(),
        reports,
        libraries: meta
            .packages
            .iter()
            .filter_map(|p| Some((p.lib_crate()?, p.name.clone())))
            .collect(),
    })
}

fn metadata(manifest_dir: &Path) -> Option<Metadata> {
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
    let out = Command::new(cargo)
        .args([
            "metadata",
            "--format-version",
            "1",
            "--no-deps",
            "--offline",
        ])
        .current_dir(manifest_dir)
        .env_remove("RUSTC_WRAPPER")
        .env_remove("RUSTC_WORKSPACE_WRAPPER")
        .stdin(Stdio::null())
        .stderr(Stdio::null())
        .output()
        .ok()?;
    out.status
        .success()
        .then(|| serde_json::from_slice(&out.stdout).ok())?
}

/// Member packages and the normal and build dependency edges between them.
struct Graph<'a> {
    packages: BTreeMap<&'a str, &'a Package>,
    /// Package -> the member packages it depends on directly.
    deps: BTreeMap<&'a str, BTreeSet<&'a str>>,
}

impl<'a> Graph<'a> {
    fn new(packages: &'a [Package]) -> Self {
        let by_name: BTreeMap<&str, &Package> =
            packages.iter().map(|p| (p.name.as_str(), p)).collect();
        let deps = packages
            .iter()
            .map(|p| {
                let edges = p
                    .dependencies
                    .iter()
                    .filter(|d| d.kind.as_deref() != Some("dev"))
                    .filter_map(|d| by_name.get_key_value(d.name.as_str()).map(|(k, _)| *k))
                    .collect();
                (p.name.as_str(), edges)
            })
            .collect();
        Graph {
            packages: by_name,
            deps,
        }
    }

    /// Every member reachable from `root` through dependency edges,
    /// `root` itself included.
    fn below(&self, root: &str) -> BTreeSet<&'a str> {
        let mut seen = BTreeSet::new();
        let mut stack: Vec<&str> = self
            .packages
            .get_key_value(root)
            .map(|(k, _)| *k)
            .into_iter()
            .collect();
        while let Some(p) = stack.pop() {
            if seen.insert(p) {
                stack.extend(self.deps.get(p).into_iter().flatten().copied());
            }
        }
        seen
    }

    /// Members that no member depends on.
    fn roots(&self) -> Vec<&'a str> {
        let depended_on: BTreeSet<&str> = self.deps.values().flatten().copied().collect();
        self.packages
            .keys()
            .filter(|p| !depended_on.contains(*p))
            .copied()
            .collect()
    }

    /// The one target of a root package that prints findings: its last
    /// binary by name if it has any, since a binary is compiled after the
    /// package's own library and may use it, else the library. `None` for
    /// a package that is not a root.
    fn reporting_crate(&self, package: &Package) -> Option<String> {
        if !self.roots().contains(&package.name.as_str()) {
            return None;
        }
        package
            .bin_crates()
            .into_iter()
            .max()
            .or_else(|| package.lib_crate())
    }

    /// Library crate names of the members whose report belongs to `root`.
    /// A library below several roots belongs to the root with the most
    /// members below it, then the first by name, so exactly one prints it.
    fn owned_by(&self, root: &str) -> BTreeSet<String> {
        let roots = self.roots();
        let below: BTreeMap<&str, BTreeSet<&str>> =
            roots.iter().map(|r| (*r, self.below(r))).collect();
        let owner = |lib: &str| {
            roots
                .iter()
                .filter(|r| below[*r].contains(lib))
                .max_by(|a, b| below[*a].len().cmp(&below[*b].len()).then_with(|| b.cmp(a)))
                .copied()
        };
        self.packages
            .values()
            .filter(|p| owner(&p.name) == Some(root))
            .filter_map(|p| p.lib_crate())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pkg(name: &str, kinds: &[(&str, &str)], deps: &[&str]) -> Package {
        Package {
            name: name.into(),
            targets: kinds
                .iter()
                .map(|(n, k)| Target {
                    name: (*n).into(),
                    kind: vec![(*k).into()],
                })
                .collect(),
            dependencies: deps
                .iter()
                .map(|d| Dep {
                    name: (*d).into(),
                    kind: None,
                })
                .collect(),
        }
    }

    #[test]
    fn the_larger_root_owns_a_shared_library() {
        let packages = vec![
            pkg("core", &[("core", "lib")], &[]),
            pkg("parse", &[("parse", "lib")], &["core"]),
            pkg("app", &[("app", "staticlib")], &["parse", "core"]),
            pkg("shim", &[("shim", "bin")], &["core"]),
        ];
        let g = Graph::new(&packages);
        assert_eq!(g.roots(), vec!["app", "shim"]);
        assert_eq!(
            g.owned_by("app"),
            BTreeSet::from(["app".to_string(), "core".into(), "parse".into()])
        );
        assert!(g.owned_by("shim").is_empty());
        assert_eq!(g.reporting_crate(&packages[2]), Some("app".into()));
        assert_eq!(g.reporting_crate(&packages[3]), Some("shim".into()));
        assert_eq!(g.reporting_crate(&packages[0]), None);
    }

    #[test]
    fn a_root_with_binaries_reports_from_the_last_binary() {
        let p = pkg(
            "tool",
            &[("tool", "lib"), ("b-two", "bin"), ("a-one", "bin")],
            &[],
        );
        let packages = vec![p];
        let g = Graph::new(&packages);
        assert_eq!(g.reporting_crate(&packages[0]), Some("b_two".into()));
    }

    #[test]
    fn dev_dependencies_are_not_edges() {
        let mut user = pkg("user", &[("user", "lib")], &[]);
        user.dependencies.push(Dep {
            name: "helper".into(),
            kind: Some("dev".into()),
        });
        let packages = vec![pkg("helper", &[("helper", "lib")], &[]), user];
        let g = Graph::new(&packages);
        assert_eq!(g.roots(), vec!["helper", "user"]);
    }
}
