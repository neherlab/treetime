//! Searches the body for the largest group of statements that could move into
//! a separate non-generic function.

use std::collections::VecDeque;

use rustc_index::bit_set::DenseBitSet;
use rustc_index::{IndexSlice, IndexVec};
use rustc_middle::mir::visit::Visitor;
use rustc_middle::mir::{self, BasicBlock, Local, Location, Rvalue, StatementKind, TerminatorKind};
use rustc_middle::ty::TyCtxt;

use super::classify::{BlockFacts, locals_using_param, per_copy_consts};
use super::region_io::{LocalFacts, inner_signature};
use crate::mir_flow::{FlowGraph, local_name, reaching, reads_any};

/// A set of usable blocks with one entry, one exit, enough hand-written
/// statements and an `inner_signature`: it could become a separate fn.
#[derive(Clone, Debug)]
pub(super) struct SharedPart {
    pub(super) blocks: DenseBitSet<BasicBlock>,
    /// Counted items, the printed size. `hand_written` ranks parts.
    pub(super) size: usize,
    hand_written: usize,
    pub(super) params: Vec<Local>,
    pub(super) returns: Vec<Local>,
    pub(super) other_parts: usize,
}

const MAX_OTHER_PARTS: usize = 3;

// Each check is cheap. A limit of 32 dropped a real finding (ui test `tail_after_refused_head`).
const MAX_IO_CHECKS_PER_ROUND: usize = 256;

/// In the flow graph, nothing dependent, no local whose type has a
/// parameter, and no tail call (which needs the caller's own signature).
fn usable_blocks<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    facts: &IndexSlice<BasicBlock, BlockFacts>,
    flow: &FlowGraph,
) -> DenseBitSet<BasicBlock> {
    let generic = locals_using_param(tcx, body);
    let mut usable = DenseBitSet::new_empty(body.basic_blocks.len());
    for (block, data) in body.basic_blocks.iter_enumerated() {
        if !flow.contains(block) || facts[block].dependent != 0 {
            continue;
        }
        if let Some(terminator) = &data.terminator
            && matches!(terminator.kind, TerminatorKind::TailCall { .. })
        {
            continue;
        }
        let generic_use = reads_any(&generic, |uses| {
            for statement in &data.statements {
                uses.visit_statement(statement, Location::START);
            }
            if let Some(terminator) = &data.terminator
                && !matches!(terminator.kind, TerminatorKind::Return)
            {
                uses.visit_terminator(terminator, Location::START);
            }
        });
        if !generic_use {
            usable.insert(block);
        }
    }
    usable
}

/// Per block, an upper bound on the hand-written size of a part entered there: the total over
/// the blocks reached only through it (its dominator subtree), each path cut at an unusable block.
fn size_bounds(
    facts: &IndexSlice<BasicBlock, BlockFacts>,
    body: &mir::Body<'_>,
    usable: &DenseBitSet<BasicBlock>,
    n: usize,
) -> IndexVec<BasicBlock, usize> {
    let mut children: IndexVec<BasicBlock, Vec<BasicBlock>> = IndexVec::from_elem_n(Vec::new(), n);
    let dominators = body.basic_blocks.dominators();
    for block in body.basic_blocks.indices() {
        if let Some(parent) = dominators.immediate_dominator(block) {
            children[parent].push(block);
        }
    }
    let mut bound: IndexVec<BasicBlock, usize> = IndexVec::from_elem_n(0, n);
    let mut stack: Vec<(BasicBlock, bool)> = vec![(mir::START_BLOCK, false)];
    while let Some((block, summed_children)) = stack.pop() {
        if summed_children {
            if usable.contains(block) {
                bound[block] = facts[block].hand_written as usize
                    + children[block].iter().map(|&c| bound[c]).sum::<usize>();
            }
        } else {
            stack.push((block, true));
            stack.extend(children[block].iter().map(|&c| (c, false)));
        }
    }
    bound
}

/// The first `members` blocks collected from an entry. `inner_signature` not yet run.
struct Candidate {
    exit: BasicBlock,
    members: usize,
    size: usize,
    hand_written: usize,
}

/// Breadth-first search from one entry. `entries`: edges from outside into
/// a member other than the entry. `exits`: edges from a member to outside.
struct EntryCandidates<'a> {
    flow: &'a FlowGraph,
    facts: &'a IndexSlice<BasicBlock, BlockFacts>,
    usable: &'a DenseBitSet<BasicBlock>,
    in_part: DenseBitSet<BasicBlock>,
    members: Vec<BasicBlock>,
    size: usize,
    hand_written: usize,
    entries: usize,
    exits: usize,
    queue: VecDeque<BasicBlock>,
}

impl<'a> EntryCandidates<'a> {
    fn new(
        flow: &'a FlowGraph,
        facts: &'a IndexSlice<BasicBlock, BlockFacts>,
        usable: &'a DenseBitSet<BasicBlock>,
        n: usize,
    ) -> Self {
        EntryCandidates {
            flow,
            facts,
            usable,
            in_part: DenseBitSet::new_empty(n + 1),
            members: Vec::new(),
            size: 0,
            hand_written: 0,
            entries: 0,
            exits: 0,
            queue: VecDeque::new(),
        }
    }

    fn reset(&mut self) {
        for &block in &self.members {
            self.in_part.remove(block);
        }
        self.members.clear();
        self.size = 0;
        self.hand_written = 0;
        self.entries = 0;
        self.exits = 0;
        self.queue.clear();
    }

    fn entry(&self) -> BasicBlock {
        self.members[0]
    }

    fn add(&mut self, block: BasicBlock) {
        debug_assert!(self.usable.contains(block) && !self.in_part.contains(block));
        let is_entry = self.members.is_empty();
        self.in_part.insert(block);
        self.members.push(block);
        self.size += self.facts[block].counted as usize;
        self.hand_written += self.facts[block].hand_written as usize;
        for &pred in self.flow.preds(block) {
            if pred == block {
                continue;
            }
            if self.in_part.contains(pred) {
                self.exits -= 1;
            } else if !is_entry {
                self.entries += 1;
            }
        }
        let entry = self.entry();
        for &succ in self.flow.succs(block) {
            if succ == block {
                continue;
            }
            if self.in_part.contains(succ) {
                if succ != entry {
                    self.entries -= 1;
                }
            } else {
                self.exits += 1;
            }
        }
        self.queue.push_back(block);
    }

    /// Adds usable blocks reachable from the queue without entering `stop_at`
    /// or EXIT. `false` on meeting a block no candidate from here may contain.
    fn add_reachable_except(
        &mut self,
        stop_at: BasicBlock,
        excluded: &DenseBitSet<BasicBlock>,
    ) -> bool {
        let exit = self.flow.exit();
        while let Some(block) = self.queue.pop_front() {
            for &succ in self.flow.succs(block) {
                if succ == stop_at || succ == exit || self.in_part.contains(succ) {
                    continue;
                }
                if !self.usable.contains(succ) || excluded.contains(succ) {
                    return false;
                }
                self.add(succ);
            }
        }
        true
    }

    fn single_entry(&self) -> bool {
        let entry = self.entry();
        self.entries == 0
            && (entry == mir::START_BLOCK
                || self
                    .flow
                    .preds(entry)
                    .iter()
                    .any(|&p| !self.in_part.contains(p)))
    }

    /// Every edge out goes to `exit`. Only correct while `exit` is not a member.
    fn single_exit(&self, exit: BasicBlock) -> bool {
        let exit_preds_inside = self
            .flow
            .preds(exit)
            .iter()
            .filter(|&&p| self.in_part.contains(p))
            .count();
        self.exits == exit_preds_inside
    }

    fn candidate(&self, exit: BasicBlock, min_statements: usize) -> Option<Candidate> {
        (self.hand_written >= min_statements
            && !self.in_part.contains(exit)
            && self.single_entry()
            && self.single_exit(exit))
        .then_some(Candidate {
            exit,
            members: self.members.len(),
            size: self.size,
            hand_written: self.hand_written,
        })
    }

    /// Every candidate from `entry`, smallest first. A single exit is on every path from the
    /// entry to EXIT (a post-dominator), so the exits tried are the nearest one, then its own.
    /// An exit with a path back to the entry is refused: the new call would run in a loop.
    fn collect(
        &mut self,
        body: &mir::Body<'_>,
        entry: BasicBlock,
        excluded: &DenseBitSet<BasicBlock>,
        min_statements: usize,
    ) -> Vec<Candidate> {
        self.reset();
        self.add(entry);
        let mut passed = Vec::new();
        let mut previous: Option<BasicBlock> = None;
        let mut reaches_entry: Option<DenseBitSet<BasicBlock>> = None;
        let chain: Vec<BasicBlock> = self.flow.ipdom_chain(entry).collect();
        for exit in chain {
            if let Some(previous) = previous
                && !self.in_part.contains(previous)
            {
                if !self.usable.contains(previous) || excluded.contains(previous) {
                    break;
                }
                self.add(previous);
            }
            if !self.add_reachable_except(exit, excluded) {
                break;
            }
            let returns_pair = self
                .flow
                .preds(exit)
                .iter()
                .any(|&pred| self.in_part.contains(pred) && ends_in_overflow_check(body, pred));
            let in_loop = exit != self.flow.exit()
                && reaches_entry
                    .get_or_insert_with(|| reaching(body, entry))
                    .contains(exit);
            if !returns_pair && !in_loop {
                passed.extend(self.candidate(exit, min_statements));
            }
            previous = Some(exit);
        }
        passed
    }

    fn blocks(&self, members: usize) -> DenseBitSet<BasicBlock> {
        let mut blocks = DenseBitSet::new_empty(self.facts.len());
        for &block in &self.members[..members] {
            blocks.insert(block);
        }
        blocks
    }
}

/// Whether `block` ends in the overflow `assert` of a `count += 1` whose
/// `(sum, overflowed)` pair is an unnamed local that the next block reads.
fn ends_in_overflow_check(body: &mir::Body<'_>, block: BasicBlock) -> bool {
    let data = &body.basic_blocks[block];
    let TerminatorKind::Assert { cond, msg, .. } = &data.terminator().kind else {
        return false;
    };
    let Some(pair) = cond.place().map(|tested| tested.local) else {
        return false;
    };
    matches!(**msg, mir::AssertKind::Overflow(..))
        && local_name(body, pair).is_none()
        && data.statements.iter().any(|statement| {
            matches!(
                &statement.kind,
                StatementKind::Assign(assign)
                    if assign.0.as_local() == Some(pair)
                        && matches!(assign.1, Rvalue::BinaryOp(op, _) if op.is_overflowing())
            )
        })
}

struct Search<'a, 'mir, 'tcx> {
    tcx: TyCtxt<'tcx>,
    body: &'mir mir::Body<'tcx>,
    flow: &'a FlowGraph,
    builder: EntryCandidates<'a>,
    size_bound: IndexVec<BasicBlock, usize>,
    per_copy_consts: DenseBitSet<Local>,
    blocks_under_const_branch: DenseBitSet<BasicBlock>,
    locals: Option<LocalFacts>,
    min_statements: usize,
}

impl Search<'_, '_, '_> {
    /// One search over usable blocks not in `excluded`: the best valid part
    /// by hand-written size, then earlier entry. May miss the largest.
    fn round(&mut self, excluded: &DenseBitSet<BasicBlock>) -> Option<SharedPart> {
        let mut best: Option<SharedPart> = None;
        let mut io_checks_left = MAX_IO_CHECKS_PER_ROUND;
        for &entry in self.body.basic_blocks.reverse_postorder() {
            if io_checks_left == 0 {
                break;
            }
            if !self.builder.usable.contains(entry)
                || excluded.contains(entry)
                || self.blocks_under_const_branch.contains(entry)
                || !self.flow.can_return(entry)
            {
                continue;
            }
            let bound = self.size_bound[entry];
            if bound < self.min_statements || best.as_ref().is_some_and(|b| bound <= b.hand_written)
            {
                continue;
            }
            let passed = self
                .builder
                .collect(self.body, entry, excluded, self.min_statements);
            let best_hand_written = best.as_ref().map_or(0, |b| b.hand_written);
            let from = passed.partition_point(|c| c.hand_written <= best_hand_written);
            if let Some(found) = self.largest_passing(entry, &passed[from..], &mut io_checks_left) {
                best = Some(found);
            }
        }
        best
    }

    /// Runs `inner_signature` on `candidates` (ascending, all above the best
    /// so far): the largest first, then a binary search down if refused.
    fn largest_passing(
        &mut self,
        entry: BasicBlock,
        candidates: &[Candidate],
        io_checks_left: &mut usize,
    ) -> Option<SharedPart> {
        let mut found = None;
        let (mut low, mut high) = (0, candidates.len());
        let mut index = high.checked_sub(1)?;
        while low < high && *io_checks_left > 0 {
            *io_checks_left -= 1;
            let candidate = &candidates[index];
            let blocks = self.builder.blocks(candidate.members);
            let exit = (candidate.exit != self.flow.exit()).then_some(candidate.exit);
            let locals = self
                .locals
                .get_or_insert_with(|| LocalFacts::new(self.body, &self.per_copy_consts));
            match inner_signature(self.tcx, self.body, locals, &blocks, entry, exit) {
                Some(io) => {
                    found = Some(SharedPart {
                        blocks,
                        size: candidate.size,
                        hand_written: candidate.hand_written,
                        params: io.params,
                        returns: io.returns,
                        other_parts: 0,
                    });
                    low = index + 1;
                }
                None => high = index,
            }
            index = low + (high - low) / 2;
        }
        found
    }
}

/// The largest valid part of `body`, then a count of further disjoint ones.
pub(super) fn best_shared_part<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    facts: &IndexSlice<BasicBlock, BlockFacts>,
    min_statements: usize,
) -> Option<SharedPart> {
    let n = body.basic_blocks.len();
    let total: usize = facts
        .iter()
        .filter(|f| f.dependent == 0)
        .map(|f| f.hand_written as usize)
        .sum();
    if total < min_statements {
        return None;
    }
    let flow = FlowGraph::new(body);
    let usable = usable_blocks(tcx, body, facts, &flow);
    let usable_total: usize = usable.iter().map(|b| facts[b].hand_written as usize).sum();
    if usable_total < min_statements {
        return None;
    }
    let (consts, branch_blocks) = per_copy_consts(tcx, body, &flow);
    let mut search = Search {
        tcx,
        body,
        flow: &flow,
        builder: EntryCandidates::new(&flow, facts, &usable, n),
        size_bound: size_bounds(facts, body, &usable, n),
        per_copy_consts: consts,
        blocks_under_const_branch: branch_blocks,
        locals: None,
        min_statements,
    };
    let mut excluded = DenseBitSet::new_empty(n);
    let mut part = search.round(&excluded)?;
    excluded.union(&part.blocks);
    while part.other_parts < MAX_OTHER_PARTS
        && let Some(other) = search.round(&excluded)
    {
        part.other_parts += 1;
        excluded.union(&other.blocks);
    }
    Some(part)
}
