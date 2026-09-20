// A value declared under one integer alias must not flow into a place declared under another.

pub type PackageId = u32;
pub type DependencyId = u32;
// An alias of an alias names the same kind, not a new one.
pub type PkgId = PackageId;
pub type NameHash = u64;

pub const INVALID_PACKAGE: PackageId = PackageId::MAX;
// Flagged: a dependency sentinel derived from the package sentinel.
pub const ROOT_DEPENDENCY: DependencyId = INVALID_PACKAGE - 1;
// Fine: same alias.
pub const LAST_PACKAGE: PackageId = INVALID_PACKAGE - 1;

pub struct Resolution {
    pub package: PackageId,
    pub dependency: DependencyId,
    pub name_hash: NameHash,
    pub count: u32,
}

pub fn link(package: PackageId, dependency: DependencyId) -> u32 {
    package ^ dependency
}

pub fn first_dependency(r: &Resolution) -> DependencyId {
    r.dependency
}

pub fn swapped_call_is_flagged(pkg: PackageId, dep: DependencyId) -> u32 {
    // Flagged: both arguments cross.
    link(dep, pkg)
}

pub fn struct_literal_is_flagged(pkg: PackageId, dep: DependencyId) -> Resolution {
    Resolution {
        // Flagged: the field is a package id, the value a dependency id.
        package: dep,
        // Fine: same alias.
        dependency: dep,
        // Fine: a `u32` alias into a `u64` alias cannot compile without a cast.
        name_hash: 0,
        // Fine: the field has no alias.
        count: pkg,
    }
}

pub fn assignment_is_flagged(r: &mut Resolution, other: &Resolution) {
    // Flagged: a field declared `DependencyId` written into one declared `PackageId`.
    r.package = other.dependency;
    // Fine: same field kind.
    r.dependency = other.dependency;
}

pub fn let_annotation_is_flagged(r: &Resolution) -> PackageId {
    // Flagged: the call returns `DependencyId`.
    let pkg: PackageId = first_dependency(r);
    pkg
}

pub fn return_is_flagged(dep: DependencyId) -> PackageId {
    // Flagged: the tail is declared `DependencyId`, the signature `PackageId`.
    let _seen = dep;
    dep
}

pub fn branch_tail_is_flagged(r: &Resolution, pick: u8) -> PackageId {
    // Flagged: each branch is a value the function may return; two cross.
    if pick == 0 {
        r.dependency
    } else {
        match pick {
            1 => r.package,
            _ => first_dependency(r),
        }
    }
}

pub fn explicit_return_is_flagged(r: &Resolution, early: bool) -> PackageId {
    if early {
        // Flagged.
        return r.dependency;
    }
    r.package
}

pub struct Pinned(pub PackageId, pub u32);

pub fn tuple_constructor_is_flagged(pkg: PackageId, dep: DependencyId) -> (Pinned, Pinned) {
    // Flagged: the first positional field is declared `PackageId`. Fine: the
    // second has no alias, and `pkg` matches.
    (Pinned(dep, dep), Pinned(pkg, dep))
}

pub fn comparison_is_flagged(r: &Resolution, dep: DependencyId) -> bool {
    // Flagged: a dependency id tested against the package sentinel.
    if dep == INVALID_PACKAGE {
        return false;
    }
    // Fine: same kind; a literal has no kind.
    dep != r.dependency && r.package < LAST_PACKAGE && dep > 3
}

// One representation per platform, like `c_int`: not an identity.
#[cfg(not(windows))]
type Backing = u32;
#[cfg(windows)]
type Backing = u64;

pub struct Handle(Backing);

impl Handle {
    pub fn cfg_selected_alias_is_fine(self) -> PackageId {
        self.0
    }
}

pub fn same_alias_is_fine(a: PackageId, r: &Resolution, dep: DependencyId) -> u32 {
    let b: PackageId = a;
    link(b, dep) + link(r.package, r.dependency)
}

pub fn alias_of_alias_is_fine(p: PkgId, dep: DependencyId) -> PackageId {
    let _ = link(p, dep);
    p
}

pub fn literal_and_arithmetic_are_fine(pkg: PackageId, dep: DependencyId) -> u32 {
    let _ = link(0, 1);
    // Arithmetic between two ids computes a number, not an id of either kind.
    link(pkg + dep, dep - pkg)
}

pub fn cast_is_fine(dep: DependencyId) -> PackageId {
    let _ = link(dep as PackageId, dep);
    dep as u32
}

pub fn plain_u32_is_fine(n: u32, r: &Resolution) -> PackageId {
    let _ = link(n, r.count);
    let m: u32 = r.dependency;
    m
}

fn main() {}
