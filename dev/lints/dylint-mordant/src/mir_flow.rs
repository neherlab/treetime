//! Helpers for reading a fn's MIR, shared by the lints that analyse fn
//! bodies. Nothing here decides what to report; each caller does that.
//!
//! * `mir_for`: the MIR body before optimization (`?` is still `Try::branch`).
//! * `assert_panics`: panics from `Assert` terminators, which have no HIR call.
//! * `place_info`: a place as an [`Atom`] (a local plus leading field
//!   indices) and how precisely the atom matches the place ([`Exactness`]).
//! * `switch_operand_atoms`: the atoms a `SwitchInt` terminator reads.
//! * [`FlowGraph`]: successors, predecessors, post-dominators, reachability,
//!   and which branch blocks decide whether a block runs, over every block
//!   except cleanup and the empty `unreachable` one, plus an added EXIT node.

use std::collections::VecDeque;

use rustc_data_structures::graph::dominators::{Dominators, dominators};
use rustc_data_structures::graph::{DirectedGraph, Predecessors, StartNode, Successors};
use rustc_hir::LangItem;
use rustc_hir::def_id::LocalDefId;
use rustc_index::bit_set::DenseBitSet;
use rustc_index::{IndexSlice, IndexVec};
use rustc_middle::mir::visit::{PlaceContext, Visitor};
use rustc_middle::mir::{
    AssertKind, BasicBlock, BasicBlockData, Body, Local, Location, Operand, Place, ProjectionElem,
    TerminatorKind, VarDebugInfoContents,
};
use rustc_middle::ty::TyCtxt;
use rustc_span::{Span, Symbol};

// ── MIR access ───────────────────────────────────────────────────────────────

pub(crate) fn mir_for<'tcx>(tcx: TyCtxt<'tcx>, def: LocalDefId) -> Option<MirRef<'tcx>> {
    if !tcx.def_kind(def).is_fn_like() || !tcx.is_mir_available(def.to_def_id()) {
        return None;
    }
    // The same at every opt level. `optimized_mir` consumes it
    // (`Steal::steal`), normally only at codegen; fall back if another driver
    // already ran it.
    let steal = tcx.mir_drops_elaborated_and_const_checked(def);
    if steal.is_stolen() {
        Some(MirRef::Opt(tcx.optimized_mir(def.to_def_id())))
    } else {
        Some(MirRef::Steal(steal.borrow()))
    }
}

pub(crate) enum MirRef<'tcx> {
    Steal(rustc_data_structures::sync::MappedReadGuard<'tcx, Body<'tcx>>),
    Opt(&'tcx Body<'tcx>),
}

impl<'tcx> std::ops::Deref for MirRef<'tcx> {
    type Target = Body<'tcx>;
    fn deref(&self) -> &Body<'tcx> {
        match self {
            MirRef::Steal(g) => g,
            MirRef::Opt(b) => b,
        }
    }
}

// ── panics with no call ──────────────────────────────────────────────────────

/// The panic function rustc calls when this assertion fails, as its lang
/// item.
///
/// [`AssertKind::panic_function`] is the same mapping codegen uses, but it
/// raises an internal compiler error on `BoundsCheck` and
/// `MisalignedPointerDereference`, whose panics take runtime arguments.
/// Those two are matched here first; `BoundsCheck` is the most common kind.
fn assert_panic_lang_item(kind: &AssertKind<Operand<'_>>) -> LangItem {
    match kind {
        AssertKind::BoundsCheck { .. } => LangItem::PanicBoundsCheck,
        AssertKind::MisalignedPointerDereference { .. } => {
            LangItem::PanicMisalignedPointerDereference
        }
        other => other.panic_function(),
    }
}

/// Every panic `body` reaches through an `Assert` terminator rather than a
/// call, as (lang item called, span that caused it). rustc adds these while
/// building MIR -- bounds checks, arithmetic overflow, division or remainder
/// by zero -- and they have no call in HIR, so a HIR visitor does not find
/// them.
/// **Overflow asserts exist only with `-C overflow-checks` on** (the debug
/// default, and what `#[rustc_inherit_overflow_checks]` propagates); bounds
/// checks and division by zero are always emitted. A missing overflow assert
/// means the profile did not ask for one, not that the code cannot overflow.
pub(crate) fn assert_panics<'a>(body: &'a Body<'_>) -> impl Iterator<Item = (LangItem, Span)> + 'a {
    body.basic_blocks.iter().filter_map(|data| {
        let term = data.terminator.as_ref()?;
        let TerminatorKind::Assert { msg, .. } = &term.kind else {
            return None;
        };
        Some((assert_panic_lang_item(msg), term.source_info.span))
    })
}

// ── places as atoms ──────────────────────────────────────────────────────────

/// A local plus its leading field indices: `_3.1.0`. Indexing adds
/// [`ANY_ELEM`], a deref is skipped, and any other projection ends the path.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub(crate) struct Atom {
    pub(crate) local: Local,
    pub(crate) path: Vec<u32>,
}

/// Path entry for "some element of": all indices count as one, but a check
/// on the container (a slice's length) is not a check on its elements.
pub(crate) const ANY_ELEM: u32 = u32::MAX;

impl Atom {
    pub(crate) fn whole(local: Local) -> Self {
        Atom {
            local,
            path: Vec::new(),
        }
    }
    pub(crate) fn extended(&self, tail: &[u32]) -> Self {
        // `node = node.next` in a loop would grow the path without limit;
        // past a few levels the distinction stops mattering, so cap it.
        const MAX_PATH: usize = 6;
        let mut path = self.path.clone();
        let room = MAX_PATH.saturating_sub(path.len());
        path.extend_from_slice(&tail[..tail.len().min(room)]);
        Atom {
            local: self.local,
            path,
        }
    }
    pub(crate) fn overlaps(&self, other: &Atom) -> bool {
        self.local == other.local && self.path.iter().zip(&other.path).all(|(a, b)| a == b)
    }
    /// `self` (an atom a branch reads) is `stored` or part of it. A branch on
    /// the whole (`lexer.next()?`) says nothing about a stored part (`.log`).
    pub(crate) fn inspects(&self, stored: &Atom) -> bool {
        self.local == stored.local && self.path.starts_with(&stored.path)
    }
}

/// How precisely `PlaceInfo::atom` matches the place that was read. Once a
/// `Downcast` is seen the state stays `VariantPayload`, so there is no
/// "payload but inexact" state; callers check the payload case first.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Exactness {
    /// Pure field access: the atom is exact.
    Exact,
    /// A `Downcast` appeared: this reads a variant payload of the atom.
    VariantPayload,
    /// Some other projection: the atom names more than was read.
    Inexact,
}

pub(crate) struct PlaceInfo {
    pub(crate) atom: Atom,
    pub(crate) exactness: Exactness,
    pub(crate) index_locals: Vec<Local>,
}

pub(crate) fn place_info(place: Place<'_>) -> PlaceInfo {
    let mut path = Vec::new();
    let mut exactness = Exactness::Exact;
    let mut index_locals = Vec::new();
    for elem in place.projection.iter() {
        match elem {
            ProjectionElem::Field(f, _) if exactness == Exactness::Exact => path.push(f.as_u32()),
            ProjectionElem::Field(..) => {}
            // `(*r).f`: treat the reference as the value, so a deref
            // neither ends the field path nor makes it inexact.
            ProjectionElem::Deref => {}
            ProjectionElem::Downcast(..) => exactness = Exactness::VariantPayload,
            ProjectionElem::Index(v) => {
                index_locals.push(v);
                if exactness == Exactness::Exact {
                    path.push(ANY_ELEM);
                }
            }
            ProjectionElem::ConstantIndex { .. } | ProjectionElem::Subslice { .. }
                if exactness == Exactness::Exact =>
            {
                path.push(ANY_ELEM);
            }
            // Any other projection makes the atom inexact, unless a
            // `Downcast` was already seen.
            _ if exactness == Exactness::VariantPayload => {}
            _ => exactness = Exactness::Inexact,
        }
    }
    PlaceInfo {
        atom: Atom {
            local: place.local,
            path,
        },
        exactness,
        index_locals,
    }
}

pub(crate) fn switch_operand_atoms(body: &Body<'_>, bb: BasicBlock) -> Vec<Atom> {
    match &body.basic_blocks[bb].terminator().kind {
        TerminatorKind::SwitchInt { discr, .. } => discr
            .place()
            .map(|p| {
                let info = place_info(p);
                let mut v: Vec<Atom> = info.index_locals.into_iter().map(Atom::whole).collect();
                v.push(info.atom);
                v
            })
            .unwrap_or_default(),
        _ => Vec::new(),
    }
}

// ── the flow graph ───────────────────────────────────────────────────────────
//
// `FlowGraph` is the body's basic blocks, minus unwind cleanup and the empty
// `unreachable` block, with the normal edges between them and one extra node,
// EXIT, that every returning block leads to. Blocks that cannot return (they
// panic, abort or never finish) are kept: a caller that moves a run of
// blocks out of a fn moves its panicking blocks too and must check them like
// any other. Three rules differ from the raw MIR edges:
// * A returning block's one successor is EXIT, so "every returning path from
//   `b` passes through `x`" is the question "does `x` post-dominate `b`".
// * A diverging block (a call that never returns, `abort`, `Unreachable`
//   after statements) has no successors, so a block no path returns from has
//   no post-dominator. Linking these to EXIT instead would make EXIT the
//   nearest post-dominator of every block that can reach a `panic!`.
// * Unwind cleanup blocks and the shared empty `unreachable` block (the
//   `otherwise` target of every enum `match` after `SimplifyCfg`) are not
//   nodes. Cleanup only leads to cleanup, so leaving it out changes no
//   dominance; the shared block would give every `match` a common successor.
// Only post-dominators are computed here, by running rustc's dominator
// routine from EXIT over the reversed edges. Dominators from the entry block
// are rustc's own `body.basic_blocks.dominators()`: the blocks and edges
// this graph leaves out never lead back to a node, so both agree on every
// node.

type Bits = DenseBitSet<BasicBlock>;

/// The terminator returns to the caller: `return`, or `become` (tail call).
fn leaves_fn(data: &BasicBlockData<'_>) -> bool {
    matches!(
        data.terminator.as_ref().map(|t| &t.kind),
        Some(TerminatorKind::Return | TerminatorKind::TailCall { .. })
    )
}

/// Normal (non-unwind) successors, sorted, deduplicated; none for cleanup.
fn raw_successors(body: &Body<'_>, data: &BasicBlockData<'_>) -> Vec<BasicBlock> {
    let Some(term) = data.terminator.as_ref().filter(|_| !data.is_cleanup) else {
        return Vec::new();
    };
    // Unwind targets are always cleanup blocks, so this keeps only normal edges.
    let mut v: Vec<_> = term
        .successors()
        .filter(|b| !body.basic_blocks[*b].is_cleanup)
        .collect();
    v.sort();
    v.dedup();
    v
}

/// See the note above. EXIT is numbered one past the last block.
pub(crate) struct FlowGraph {
    /// Indexed by block and by EXIT; empty for a block that is not a node.
    succs: IndexVec<BasicBlock, Vec<BasicBlock>>,
    /// `succs` inverted; `preds[exit]` is the returning blocks.
    preds: IndexVec<BasicBlock, Vec<BasicBlock>>,
    /// The blocks that are nodes, and EXIT.
    nodes: Bits,
    exit: BasicBlock,
    /// Dominators over the reversed edges from EXIT: post-dominators.
    pdom: Dominators<BasicBlock>,
}

/// The input rustc's dominator routine takes: an edge list, its inverse, a
/// start node. Swapping the lists and starting at EXIT gives post-dominators.
struct Edges<'g> {
    out: &'g IndexSlice<BasicBlock, Vec<BasicBlock>>,
    into: &'g IndexSlice<BasicBlock, Vec<BasicBlock>>,
    root: BasicBlock,
}

impl DirectedGraph for Edges<'_> {
    type Node = BasicBlock;
    fn num_nodes(&self) -> usize {
        self.out.len()
    }
}

impl StartNode for Edges<'_> {
    fn start_node(&self) -> BasicBlock {
        self.root
    }
}

impl Successors for Edges<'_> {
    fn successors(&self, node: BasicBlock) -> impl Iterator<Item = BasicBlock> {
        self.out[node].iter().copied()
    }
}

impl Predecessors for Edges<'_> {
    fn predecessors(&self, node: BasicBlock) -> impl Iterator<Item = BasicBlock> {
        self.into[node].iter().copied()
    }
}

impl FlowGraph {
    pub(crate) fn new(body: &Body<'_>) -> Self {
        let exit = BasicBlock::from_usize(body.basic_blocks.len());
        let size = exit.as_usize() + 1;
        let mut nodes = Bits::new_empty(size);
        for (b, data) in body.basic_blocks.iter_enumerated() {
            let shared_unreachable = data.terminator.is_some() && data.is_empty_unreachable();
            if !data.is_cleanup && !shared_unreachable {
                nodes.insert(b);
            }
        }
        nodes.insert(exit);
        let mut succs: IndexVec<BasicBlock, Vec<BasicBlock>> =
            IndexVec::from_elem_n(Vec::new(), size);
        let mut preds: IndexVec<BasicBlock, Vec<BasicBlock>> =
            IndexVec::from_elem_n(Vec::new(), size);
        for (b, data) in body.basic_blocks.iter_enumerated() {
            if !nodes.contains(b) {
                continue;
            }
            let out = &mut succs[b];
            if leaves_fn(data) {
                out.push(exit);
            } else {
                out.extend(
                    raw_successors(body, data)
                        .into_iter()
                        .filter(|s| nodes.contains(*s)),
                );
            }
            for &s in out.iter() {
                preds[s].push(b);
            }
        }
        let pdom = dominators(&Edges {
            out: &preds,
            into: &succs,
            root: exit,
        });
        FlowGraph {
            succs,
            preds,
            nodes,
            exit,
            pdom,
        }
    }

    /// The added node every `return` leads to; one past the last block index.
    pub(crate) fn exit(&self) -> BasicBlock {
        self.exit
    }

    /// `b` is EXIT or a non-cleanup block other than the empty `unreachable`.
    pub(crate) fn contains(&self, b: BasicBlock) -> bool {
        self.nodes.contains(b)
    }

    /// Normal successors of `b`; `[EXIT]` if `b` returns; empty if it diverges.
    pub(crate) fn succs(&self, b: BasicBlock) -> &[BasicBlock] {
        &self.succs[b]
    }

    /// The nodes `b` is a successor of; for EXIT, the returning blocks.
    pub(crate) fn preds(&self, b: BasicBlock) -> &[BasicBlock] {
        &self.preds[b]
    }

    /// Some path of normal edges leads from `b` to a `return`.
    pub(crate) fn can_return(&self, b: BasicBlock) -> bool {
        self.nodes.contains(b) && self.pdom.is_reachable(b)
    }

    /// Nearest node, other than `b`, on every returning path from `b`
    /// (immediate post-dominator). `None` for EXIT, non-returning, non-node.
    pub(crate) fn ipdom(&self, b: BasicBlock) -> Option<BasicBlock> {
        if self.nodes.contains(b) {
            self.pdom.immediate_dominator(b)
        } else {
            None
        }
    }

    /// `ipdom` applied repeatedly from `b`, ending at EXIT; empty if none.
    pub(crate) fn ipdom_chain(&self, b: BasicBlock) -> impl Iterator<Item = BasicBlock> + '_ {
        std::iter::successors(self.ipdom(b), |&x| self.ipdom(x))
    }

    /// The blocks that run or not depending on which successor `b` takes
    /// (those directly control-dependent on `b`), in block order, never EXIT:
    /// `ipdom` followed from each successor up to but excluding `ipdom(b)`,
    /// stopping where a block has none. `b` itself is included when it
    /// post-dominates a successor (a loop test). Empty if < 2 successors.
    pub(crate) fn decides(&self, b: BasicBlock) -> Vec<BasicBlock> {
        let succs = self.succs(b);
        if succs.len() < 2 {
            return Vec::new();
        }
        let join = self.ipdom(b);
        let mut decided = Bits::new_empty(self.succs.len());
        for &s in succs {
            let mut node = Some(s);
            while let Some(n) = node
                && node != join
                && n != self.exit
                && decided.insert(n)
            {
                node = self.ipdom(n);
            }
        }
        decided.iter().collect()
    }

    /// The branch blocks that directly decide whether `t` runs, in block
    /// order: every returning path from one of the block's successors passes
    /// through `t`, but not every returning path from the block itself,
    /// unless the block is `t` (a loop test). In other words `t`
    /// post-dominates a successor but not strictly the block. Blocks that
    /// cannot return are ignored, so a panic or abort counts as "does not
    /// happen"; otherwise every block after `assert!(x)` would depend on `x`.
    /// A branch that only decides whether such a branch is reached (an early
    /// `return Ok` before it) is not included; `control_deps` adds those.
    pub(crate) fn direct_control_deps(&self, t: BasicBlock) -> Vec<BasicBlock> {
        let t_post_dominates =
            |b: BasicBlock| self.pdom.is_reachable(b) && self.pdom.dominates(t, b);
        self.succs
            .iter_enumerated()
            .filter(|(a, succs)| {
                succs.len() >= 2
                    && (*a == t || !t_post_dominates(*a))
                    && succs.iter().any(|&s| t_post_dominates(s))
            })
            .map(|(a, _)| a)
            .collect()
    }

    /// Every branch block that decides whether `target` runs, directly or
    /// through another such branch, nearest first.
    pub(crate) fn control_deps(&self, target: BasicBlock) -> Vec<BasicBlock> {
        let mut seen = Bits::new_empty(self.succs.len());
        let mut queue = VecDeque::from([target]);
        let mut deps = Vec::new();
        while let Some(t) = queue.pop_front() {
            for a in self.direct_control_deps(t) {
                if seen.insert(a) {
                    deps.push(a);
                    queue.push_back(a);
                }
            }
        }
        deps
    }
}

// ── Uses, reach, names ───────────────────────────────────────────────────────

/// Records whether any visited place uses a local in `of`. Storage markers
/// do not count; callers that ignore `return`'s read of `_0` skip it.
pub(crate) struct UsesAny<'a> {
    of: &'a DenseBitSet<Local>,
    found: bool,
}

impl<'tcx> Visitor<'tcx> for UsesAny<'_> {
    fn visit_local(&mut self, local: Local, context: PlaceContext, _: Location) {
        if context.is_use() && self.of.contains(local) {
            self.found = true;
        }
    }
}

/// Whether what `visit` walks uses a local in `of`.
pub(crate) fn reads_any(of: &DenseBitSet<Local>, visit: impl FnOnce(&mut UsesAny<'_>)) -> bool {
    let mut uses = UsesAny { of, found: false };
    visit(&mut uses);
    uses.found
}

/// Every block with a path to `to` along any edge, `to` included.
pub(crate) fn reaching(body: &Body<'_>, to: BasicBlock) -> DenseBitSet<BasicBlock> {
    let preds = body.basic_blocks.predecessors();
    let mut seen = DenseBitSet::new_empty(body.basic_blocks.len());
    seen.insert(to);
    let mut pending = vec![to];
    while let Some(block) = pending.pop() {
        for &pred in &preds[block] {
            if seen.insert(pred) {
                pending.push(pred);
            }
        }
    }
    seen
}

/// The name the user wrote for `local`, if any. A debug-info entry does not
/// count if it names only a field of `local`, if `local` holds only part of
/// the named variable, or if it comes from a macro or desugaring (`iter`,
/// `residual`).
pub(crate) fn local_name(body: &Body<'_>, local: Local) -> Option<Symbol> {
    body.var_debug_info
        .iter()
        .find_map(|info| match info.value {
            VarDebugInfoContents::Place(place)
                if place.as_local() == Some(local)
                    && info.composite.is_none()
                    && !info.source_info.span.from_expansion() =>
            {
                Some(info.name)
            }
            VarDebugInfoContents::Place(_) | VarDebugInfoContents::Const(_) => None,
        })
}
