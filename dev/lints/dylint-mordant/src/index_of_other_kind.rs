use std::collections::{HashMap, HashSet};
use std::ops::ControlFlow;

use clippy_utils::visitors::for_each_expr_without_closures;
use rustc_hir::def::Res;
use rustc_hir::{Body, Expr, ExprKind, QPath, UnOp};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::{Span, Symbol};

use crate::baseline::emit_with_note;
use crate::hir_shapes::value_name;

rustc_session::declare_lint! {
    /// Flags an index whose name claims one kind (`source_index`, `pkg_id`:
    /// the non-empty prefix before `_index`, `_idx`, `_id` or `_i`) landing
    /// on a place that everywhere else in the function is indexed by names
    /// of another kind, when the function also shows the crossing name
    /// indexing a table named after it: `sources[source_index]`,
    /// `parts[part_index]` twice, then `parts[source_index]`. This looks
    /// like the wrong table, and since both indices are plain integers it
    /// compiles. `[]`, `.get`, `.get_mut` and `get_unchecked*` with an
    /// integer operand all count as indexing (a newtyped index is already
    /// told apart by the compiler); a place is the binding or item at the
    /// root plus the
    /// fields, zero-argument accessors and earlier indices on the way
    /// (`self.graph.parts`, `lockfile.packages.names()`, `parts[..]`), so a
    /// two-level table is a different place at each level and two locals
    /// both called `resolutions` are two places. The place's own kind is the
    /// one it is named after, else the one indexing it more often, else the
    /// one with no table of its own in the function; when nothing tells them
    /// apart the crossing is still there and is reported at the kind used
    /// second.
    ///
    /// The claim is read off names only, never types: silent when either
    /// name carries no kind suffix (`i`, `n`, `at`), when the two prefixes
    /// share or abbreviate a word (`dep_id`/`dependency_id`,
    /// `pkg_id`/`package_id`, `other_chunk_index`/`chunk_index`), when the
    /// place is itself named after the crossing kind, and when the crossing
    /// name has no table of its own name in the function that the other name
    /// leaves alone, which is how a role name for the same kind reads
    /// (`nodes[parent_idx]` beside `nodes[node_idx]`, `symbols[existing_id]`
    /// beside `symbols[symbol_id]`).
    pub INDEX_OF_OTHER_KIND,
    Warn,
    "a place indexed by names of two different index kinds within one function"
}

rustc_session::declare_lint_pass!(IndexOfOtherKind => [INDEX_OF_OTHER_KIND]);

/// The index kind a name claims: the non-empty prefix before `_index`,
/// `_idx`, `_id` or `_i`.
fn claimed_kind(name: &str) -> Option<&str> {
    ["_index", "_idx", "_id", "_i"]
        .iter()
        .find_map(|suffix| name.strip_suffix(suffix))
        .filter(|prefix| !prefix.is_empty())
}

/// `short` spells `long`: equal, a prefix of it, or its letters in order
/// from the same initial (`pkg`/`package`, `func`/`function`), three letters
/// at least so `id`-sized fragments do not unify everything.
fn abbreviates(short: &str, long: &str) -> bool {
    if short == long {
        return true;
    }
    if short.len() < 3 || short.len() > long.len() || short.as_bytes()[0] != long.as_bytes()[0] {
        return false;
    }
    let mut rest = long.chars();
    short.chars().all(|c| rest.any(|l| l == c))
}

/// `entries` -> `entry`, `hashes` -> `hash`, `parts` -> `part`; anything
/// else unchanged.
fn singular(word: &str) -> &str {
    if let Some(stem) = word.strip_suffix("ies") {
        // `y` is not in the stem, but `abbreviates` only needs the prefix.
        return stem;
    }
    for tail in ["shes", "ches", "xes", "sses"] {
        if word.ends_with(tail) {
            return &word[..word.len() - 2];
        }
    }
    word.strip_suffix('s')
        .filter(|stem| !stem.ends_with('s'))
        .unwrap_or(word)
}

/// Two `_`-separated names name one kind when any word of one is, or
/// abbreviates, a word of the other, plurals aside: `unresolved_dep` and
/// `dependencies`, `other_chunk` and `chunk`, `pkg` and `package_names`.
fn same_kind(a: &str, b: &str) -> bool {
    fn words(s: &str) -> impl Iterator<Item = &str> {
        s.split('_').filter(|w| !w.is_empty()).map(singular)
    }
    words(a).any(|x| words(b).any(|y| abbreviates(x, y) || abbreviates(y, x)))
}

/// Zero-argument methods that expose the receiver's own elements rather
/// than select a different table, so `v.as_slice()[i]` indexes `v`.
const VIEWS: &[&str] = &[
    "as_slice",
    "as_mut_slice",
    "slice",
    "slice_mut",
    "as_ref",
    "as_mut",
    "borrow",
    "borrow_mut",
    "deref",
    "deref_mut",
    "iter",
    "iter_mut",
    "unwrap",
];

/// Where an index lands, as an identity (`key`, roots told apart by
/// resolution), as it reads in the source (`shown`), and the names along
/// it (`root`, each field and accessor) that may say what kind indexes it.
struct Place {
    key: String,
    shown: String,
    names: Vec<Symbol>,
}

impl Place {
    /// Some name along the place is the table of `kind`: `sources` or
    /// `self.graph.input_files` for `source`, `lockfile.packages.names()`
    /// for `pkg`.
    fn named_after(&self, kind: &str) -> bool {
        self.names.iter().any(|n| same_kind(n.as_str(), kind))
    }
}

/// The place `e` denotes: a local, parameter or item at the root, then
/// fields, non-view zero-argument accessors and earlier indices. `&`, `*`
/// and HIR temporaries are transparent. Anything else at the root (a call
/// with arguments, a literal) is no place.
fn place_of(cx: &LateContext<'_>, mut e: &Expr<'_>) -> Option<Place> {
    let mut path = Vec::new();
    let mut names = Vec::new();
    loop {
        match e.kind {
            ExprKind::AddrOf(_, _, inner)
            | ExprKind::Unary(UnOp::Deref, inner)
            | ExprKind::DropTemps(inner) => e = inner,
            ExprKind::Field(inner, ident) => {
                path.push(format!(".{}", ident.name));
                names.push(ident.name);
                e = inner;
            }
            ExprKind::Index(inner, _, _) => {
                path.push("[..]".to_string());
                e = inner;
            }
            ExprKind::MethodCall(seg, recv, [], _) => {
                if !VIEWS.contains(&seg.ident.name.as_str()) {
                    path.push(format!(".{}()", seg.ident.name));
                    names.push(seg.ident.name);
                }
                e = recv;
            }
            ExprKind::Path(ref qpath) => {
                let root = match cx.qpath_res(qpath, e.hir_id) {
                    Res::Local(id) => format!("{id:?}"),
                    Res::Def(_, did) => format!("{did:?}"),
                    _ => return None,
                };
                let name = match qpath {
                    QPath::Resolved(_, p) => p.segments.last()?.ident.name,
                    QPath::TypeRelative(_, seg) => seg.ident.name,
                };
                names.push(name);
                path.reverse();
                let tail = path.concat();
                return Some(Place {
                    key: format!("{root}{tail}"),
                    shown: format!("{name}{tail}"),
                    names,
                });
            }
            _ => return None,
        }
    }
}

/// `base[idx]`, `base.get(idx)`, `base.get_mut(idx)`,
/// `base.get_unchecked(idx)`, `base.get_unchecked_mut(idx)`.
fn index_parts<'h>(e: &'h Expr<'h>) -> Option<(&'h Expr<'h>, &'h Expr<'h>)> {
    match e.kind {
        ExprKind::Index(base, idx, _) => Some((base, idx)),
        ExprKind::MethodCall(seg, recv, [idx], _)
            if matches!(
                seg.ident.name.as_str(),
                "get" | "get_mut" | "get_unchecked" | "get_unchecked_mut"
            ) =>
        {
            Some((recv, idx))
        }
        _ => None,
    }
}

struct Site {
    place: usize,
    kind: Symbol,
    name: Symbol,
    span: Span,
}

/// One index kind's sites on one place, names that share or abbreviate a
/// word counted as one kind (`pkg_id` with `package_id`).
struct KindHere<'a> {
    kinds: Vec<Symbol>,
    first: &'a Site,
    sites: Vec<&'a Site>,
}

impl KindHere<'_> {
    fn has(&self, kind: Symbol) -> bool {
        self.kinds
            .iter()
            .any(|k| same_kind(k.as_str(), kind.as_str()))
    }

    /// Some name of this kind indexes `place` somewhere in the body.
    fn reaches(&self, reach: &HashMap<Symbol, HashSet<usize>>, place: usize) -> bool {
        self.kinds.iter().any(|k| reach[k].contains(&place))
    }
}

/// One body's index sites and the places they land on.
struct Sites {
    places: Vec<Place>,
    sites: Vec<Site>,
}

impl Sites {
    fn collect<'tcx>(cx: &LateContext<'tcx>, body: &Body<'tcx>) -> Self {
        let mut places: Vec<Place> = Vec::new();
        let mut place_ids: HashMap<String, usize> = HashMap::new();
        let mut sites = Vec::new();
        for_each_expr_without_closures(body.value, |e: &'tcx Expr<'tcx>| {
            if !e.span.from_expansion()
                && let Some((base, idx)) = index_parts(e)
                && cx.typeck_results().expr_ty(idx).peel_refs().is_integral()
                && let Some(name) = value_name(idx)
                && let Some(kind) = claimed_kind(name.name.as_str())
                && let Some(place) = place_of(cx, base)
            {
                let next = places.len();
                let id = *place_ids.entry(place.key.clone()).or_insert(next);
                if id == next {
                    places.push(place);
                }
                sites.push(Site {
                    place: id,
                    kind: Symbol::intern(kind),
                    name: name.name,
                    span: e.span,
                });
            }
            ControlFlow::<()>::Continue(())
        });
        Self { places, sites }
    }

    /// For each place indexed by two kinds, every site of the visiting kind,
    /// with the place's own kind's first site and count and the visiting
    /// kind's first site on a table of its own. A kind is visiting when the
    /// place is not named after it, the body indexes a table named after it
    /// that the other kind never touches, and the other kind is not the
    /// lesser claim: named on the place, else used there more often, else
    /// (one use each) without such a table of its own, else merely used
    /// first, which decides only where the one report of a symmetric
    /// crossing goes, not whether there is one.
    fn crossings(&self) -> Vec<(&Site, &Site, usize, &Site)> {
        let mut reach: HashMap<Symbol, HashSet<usize>> = HashMap::new();
        let mut by_place: HashMap<usize, Vec<KindHere<'_>>> = HashMap::new();
        for site in &self.sites {
            reach.entry(site.kind).or_default().insert(site.place);
            let kinds = by_place.entry(site.place).or_default();
            match kinds.iter_mut().find(|k| k.has(site.kind)) {
                Some(k) => {
                    if !k.kinds.contains(&site.kind) {
                        k.kinds.push(site.kind);
                    }
                    if site.span.lo() < k.first.span.lo() {
                        k.first = site;
                    }
                    k.sites.push(site);
                }
                None => kinds.push(KindHere {
                    kinds: vec![site.kind],
                    first: site,
                    sites: vec![site],
                }),
            }
        }
        // `k`'s first site on a table named after it that `other` never
        // indexes.
        let home_of = |k: &KindHere<'_>, other: &KindHere<'_>| {
            self.sites
                .iter()
                .filter(|s| {
                    k.kinds.contains(&s.kind)
                        && !other.reaches(&reach, s.place)
                        && self.places[s.place].named_after(s.kind.as_str())
                })
                .min_by_key(|s| (s.place, s.span.lo()))
        };
        let mut out = Vec::new();
        for (&place, kinds) in &by_place {
            let named = |k: &KindHere<'_>| {
                k.kinds
                    .iter()
                    .any(|s| self.places[place].named_after(s.as_str()))
            };
            let rank = |k: &KindHere<'_>| (named(k), k.sites.len());
            for (i, k) in kinds.iter().enumerate() {
                if named(k) {
                    continue;
                }
                let mut majors: Vec<&KindHere<'_>> = kinds
                    .iter()
                    .enumerate()
                    .filter(|&(j, m)| j != i && rank(m) >= rank(k))
                    .map(|(_, m)| m)
                    .collect();
                majors.sort_by_key(|m| (std::cmp::Reverse(rank(m)), m.first.span.lo()));
                let Some((major, home)) = majors.into_iter().find_map(|m| {
                    let h = home_of(k, m)?;
                    (rank(m) > rank(k)
                        || home_of(m, k).is_none()
                        || m.first.span.lo() < k.first.span.lo())
                    .then_some((m, h))
                }) else {
                    continue;
                };
                for &site in &k.sites {
                    out.push((site, major.first, major.sites.len(), home));
                }
            }
        }
        out.sort_by_key(|(site, ..)| site.span.lo());
        out
    }
}

impl<'tcx> LateLintPass<'tcx> for IndexOfOtherKind {
    fn check_body(&mut self, cx: &LateContext<'tcx>, body: &Body<'tcx>) {
        let sites = Sites::collect(cx, body);
        for (site, major_first, major_n, home) in sites.crossings() {
            // The table of the visiting kind may be indexed under another
            // name of that kind (`dep_id` there, `dep_idx` here).
            let by = if home.name == site.name {
                format!("`{}`", site.name)
            } else {
                format!("`{}` (the same kind as `{}`)", home.name, site.name)
            };
            let place = &sites.places[site.place].shown;
            let home_place = &sites.places[home.place].shown;
            let name = site.name;
            emit_with_note(
                cx,
                INDEX_OF_OTHER_KIND,
                site.span,
                format!(
                    "`{place}` is indexed by `{name}` here. Elsewhere in this function it is \
                     indexed by `{major}` ({major_n} place{s}), and {by} indexes `{home_place}`. \
                     This looks like the wrong table",
                    major = major_first.name,
                    s = if major_n == 1 { "" } else { "s" },
                ),
                major_first.span,
                format!("`{place}` indexed the usual way"),
                format!(
                    "give `{place}` and `{home_place}` their own index newtypes so \
                     `{place}[{name}]` cannot compile. If this line is intended, rename the \
                     index to say so"
                ),
            );
        }
    }
}
