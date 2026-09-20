// compile-flags: --edition=2021
//
// A fact about a field of `self` (a length, a flag, a pointer), then a
// statement that runs code this function cannot see and which can reach the
// field, then the fact used against the field as if nothing had happened.
// The ui config names `Vm::run_callback`, `dispatch*`, and the trait-impl
// methods `Worker::run_job` and `Runner::schedule` as this project's own
// re-entry points; those count without any argument analysis.

use std::cell::RefCell;
use std::ptr::NonNull;

trait Sink {
    fn notify(&self);
    fn drain(&self, items: &mut Vec<u32>);
}

struct Stdout;

impl Sink for Stdout {
    fn notify(&self) {}
    fn drain(&self, items: &mut Vec<u32>) {
        items.clear();
    }
}

struct Vm;

impl Vm {
    fn run_callback(&self) {}
    fn tick(&self) {}
}

trait Runner {
    fn run_job(&self);
    fn schedule(&self);
    fn poll(&self);
}

struct Worker;

impl Runner for Worker {
    fn run_job(&self) {}
    fn schedule(&self) {}
    fn poll(&self) {}
}

/// A collection of this crate: its methods may mutate through `&self`, so a
/// fact about it can go stale in a `&self` method.
struct Log {
    lines: RefCell<Vec<u32>>,
}

impl Log {
    fn len(&self) -> usize {
        self.lines.borrow().len()
    }
    fn at(&self, i: usize) -> u32 {
        self.lines.borrow()[i]
    }
    fn remove(&self, i: usize) -> u32 {
        self.lines.borrow_mut().remove(i)
    }
}

/// The same, with the cell behind a pointer: `Freeze` on its own, but a
/// shared reference still gets to the cell.
struct Journal {
    lines: Box<RefCell<Vec<u32>>>,
}

impl Journal {
    fn len(&self) -> usize {
        self.lines.borrow().len()
    }
    fn remove(&self, i: usize) -> u32 {
        self.lines.borrow_mut().remove(i)
    }
}

/// A collection of this crate that nothing can change through `&`.
struct Ring(Vec<u32>);

impl Ring {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::ops::Index<usize> for Ring {
    type Output = u32;
    fn index(&self, i: usize) -> &u32 {
        &self.0[i]
    }
}

async fn yield_now() {}

struct Host {
    items: Vec<u32>,
    buf: Vec<u8>,
    pending: Option<u32>,
    log: Log,
    journal: Journal,
    ring: Ring,
    raw: NonNull<u8>,
    borrowed: &'static [u8],
    on_event: fn(&mut Host),
    sink: Box<dyn Sink>,
    vm: Vm,
    worker: Worker,
    tasks: Vec<Box<dyn Fn(&Vm)>>,
}

fn no_op(_: &mut Host) {}

fn peek(_: &Host) {}

fn peek_items(_: &Vec<u32>) {}

impl Host {
    fn dispatch_event(&self) {}

    fn plain(&self) {}

    fn count_across_fn_pointer_field_given_self_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        (self.on_event)(self);
        self.items[n - 1]
    }

    fn position_across_dyn_given_the_field_is_flagged(&mut self) -> u32 {
        let last = self.items.len() - 1;
        self.sink.drain(&mut self.items);
        self.items.remove(last)
    }

    fn gated_count_across_closure_given_self_is_flagged<F: FnMut(&mut Host)>(
        &mut self,
        mut f: F,
    ) -> u32 {
        let n = self.items.len();
        f(self);
        if n > 0 { self.items[0] } else { 0 }
    }

    fn loop_bound_across_closure_given_self_is_flagged<F: FnMut(&mut Host)>(
        &mut self,
        mut f: F,
    ) -> u32 {
        let n = self.items.len();
        f(self);
        let mut total = 0;
        for i in 0..n {
            total += self.items[i];
        }
        total
    }

    fn flag_across_fn_pointer_param_given_self_is_flagged(&mut self, f: fn(&mut Host)) -> u32 {
        let had = self.pending.is_some();
        f(self);
        if had { self.pending.take().unwrap() } else { 0 }
    }

    fn pointer_across_dyn_given_self_is_flagged(&mut self, f: &dyn Fn(&mut Host)) -> u8 {
        let p = self.buf.as_ptr();
        f(self);
        unsafe { *p }
    }

    fn crate_collection_across_dyn_in_shared_method_is_flagged(&self) -> u32 {
        let n = self.log.len();
        self.sink.notify();
        self.log.remove(n - 1)
    }

    async fn crate_collection_across_await_in_shared_method_is_flagged(&self) -> u32 {
        let n = self.log.len();
        yield_now().await;
        self.log.remove(n - 1)
    }

    fn count_across_configured_callee_is_flagged(&mut self) -> u32 {
        let n = self.items.len() as u32;
        self.vm.run_callback();
        self.items[n as usize - 1]
    }

    fn count_across_configured_prefix_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        self.dispatch_event();
        self.items[n - 1]
    }

    fn closure_not_given_self_is_fine<F: FnMut(&[u32])>(&mut self, mut f: F) -> u32 {
        let n = self.items.len();
        f(&self.items);
        self.items[n - 1]
    }

    fn dyn_not_given_the_field_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        self.sink.notify();
        self.items[n - 1]
    }

    fn std_collection_in_shared_method_is_fine(&self) -> u32 {
        let n = self.items.len();
        self.sink.notify();
        self.items[n - 1]
    }

    async fn await_in_exclusive_method_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        yield_now().await;
        self.items[n - 1]
    }

    fn crate_collection_read_after_reentry_is_fine(&self) -> u32 {
        let n = self.log.len();
        self.sink.notify();
        self.log.at(0) + n as u32
    }

    fn pointer_copied_out_of_a_wrapper_is_fine(&mut self, f: fn(&mut Host)) -> u8 {
        let p = self.raw.as_ptr();
        f(self);
        unsafe { *p }
    }

    fn fact_about_a_reference_field_is_fine(&mut self, f: fn(&mut Host)) -> u8 {
        let n = self.borrowed.len();
        f(self);
        self.borrowed[n - 1]
    }

    fn count_across_plain_method_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        self.plain();
        self.items[n - 1]
    }

    fn count_across_unconfigured_method_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        self.vm.tick();
        self.items[n - 1]
    }

    fn requeried_after_reentry_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        (self.on_event)(self);
        if n <= self.items.len() { self.items[n - 1] } else { 0 }
    }

    fn used_only_before_reentry_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        let v = self.items[n - 1];
        (self.on_event)(self);
        v
    }

    fn fact_about_a_local_is_fine(&mut self) -> u32 {
        let local = vec![1, 2, 3];
        let n = local.len();
        (self.on_event)(self);
        local[n - 1]
    }

    fn non_panicking_reuse_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        (self.on_event)(self);
        self.items.truncate(n);
        self.items.get(n).copied().unwrap_or(0)
    }

    fn field_replaced_after_reentry_is_fine(&mut self) -> usize {
        let n = self.items.len();
        (self.on_event)(self);
        self.items = Vec::new();
        self.items.get(n).map_or(n, |v| *v as usize)
    }

    fn binding_reassigned_after_reentry_is_fine(&mut self) -> u32 {
        let mut n = self.items.len();
        (self.on_event)(self);
        n = 0;
        self.items[n]
    }

    fn other_field_is_fine(&mut self) -> u8 {
        let n = self.items.len();
        (self.on_event)(self);
        self.buf[n]
    }

    fn local_closure_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        let bump = |_: &mut Host| 1;
        let extra = bump(self);
        self.items[n - 1] + extra
    }

    fn pending_is_read_properly(&self) -> u32 {
        if let Some(p) = self.pending { p } else { 0 }
    }

    fn index_in_a_condition_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        (self.on_event)(self);
        if self.items[n - 1] > 0 { 1 } else { 0 }
    }

    fn removal_in_a_scrutinee_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        (self.on_event)(self);
        match self.items.remove(n - 1) {
            0 => 0,
            v => v,
        }
    }

    fn reentry_before_a_break_inside_the_statement_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        for i in 0..2 {
            if i == 1 {
                (self.on_event)(self);
                break;
            }
        }
        self.items[n - 1]
    }

    fn configured_impl_method_named_by_its_type_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        self.worker.run_job();
        self.items[n - 1]
    }

    fn configured_impl_method_named_by_its_trait_is_flagged(&mut self) -> u32 {
        let n = self.items.len();
        self.worker.schedule();
        self.items[n - 1]
    }

    fn crate_type_with_the_cell_behind_a_box_is_flagged(&self) -> u32 {
        let n = self.journal.len();
        self.sink.notify();
        self.journal.remove(n - 1)
    }

    fn frozen_crate_type_given_self_is_flagged(&mut self) -> u32 {
        let n = self.ring.len();
        (self.on_event)(self);
        self.ring[n - 1]
    }

    fn crate_collection_given_shared_self_is_flagged(&mut self, f: fn(&Host)) -> u32 {
        let n = self.log.len();
        f(self);
        self.log.remove(n - 1)
    }

    fn reentry_in_a_returning_branch_is_fine(&mut self, early: bool) -> u32 {
        let n = self.items.len();
        if early {
            (self.on_event)(self);
            return 0;
        }
        self.items[n - 1]
    }

    fn returned_reentry_is_fine(&mut self, early: bool, f: fn(&mut Host) -> u32) -> u32 {
        let n = self.items.len();
        if early {
            return f(self);
        }
        self.items[n - 1]
    }

    fn reentry_before_a_break_out_of_the_block_is_fine(&mut self, early: bool) -> u32 {
        let mut total = 0;
        for _ in 0..2 {
            let n = self.items.len();
            if early {
                (self.on_event)(self);
                break;
            }
            total += self.items[n - 1];
        }
        total
    }

    fn configured_callee_inside_a_stored_closure_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        self.tasks.push(Box::new(|vm: &Vm| vm.run_callback()));
        self.items[n - 1]
    }

    fn frozen_crate_type_in_shared_method_is_fine(&self) -> u32 {
        let n = self.ring.len();
        self.sink.notify();
        self.ring[n - 1]
    }

    fn fn_pointer_taking_shared_self_is_fine(&mut self, f: fn(&Host)) -> u32 {
        let n = self.items.len();
        f(self);
        self.items[n - 1]
    }

    fn field_handed_out_shared_by_coercion_is_fine(&mut self, f: fn(&Vec<u32>)) -> u32 {
        let n = self.items.len();
        f(&mut self.items);
        self.items[n - 1]
    }

    fn unconfigured_impl_method_is_fine(&mut self) -> u32 {
        let n = self.items.len();
        self.worker.poll();
        self.items[n - 1]
    }
}

fn main() {
    static BYTES: [u8; 2] = [1, 2];
    let mut byte = 7u8;
    let mut h = Host {
        items: vec![1, 2, 3],
        buf: vec![4],
        pending: Some(5),
        log: Log {
            lines: RefCell::new(vec![6]),
        },
        journal: Journal {
            lines: Box::new(RefCell::new(vec![7])),
        },
        ring: Ring(vec![8]),
        raw: NonNull::from(&mut byte),
        borrowed: &BYTES,
        on_event: no_op,
        sink: Box::new(Stdout),
        vm: Vm,
        worker: Worker,
        tasks: Vec::new(),
    };
    let _ = h.count_across_fn_pointer_field_given_self_is_flagged()
        + h.gated_count_across_closure_given_self_is_flagged(|_| ())
        + h.loop_bound_across_closure_given_self_is_flagged(|_| ())
        + h.flag_across_fn_pointer_param_given_self_is_flagged(no_op)
        + u32::from(h.pointer_across_dyn_given_self_is_flagged(&|_| ()))
        + h.crate_collection_across_dyn_in_shared_method_is_flagged()
        + h.count_across_configured_callee_is_flagged()
        + h.count_across_configured_prefix_is_flagged()
        + h.closure_not_given_self_is_fine(|_| ())
        + h.dyn_not_given_the_field_is_fine()
        + h.std_collection_in_shared_method_is_fine()
        + h.crate_collection_read_after_reentry_is_fine()
        + u32::from(h.pointer_copied_out_of_a_wrapper_is_fine(no_op))
        + u32::from(h.fact_about_a_reference_field_is_fine(no_op))
        + h.count_across_plain_method_is_fine()
        + h.count_across_unconfigured_method_is_fine()
        + h.requeried_after_reentry_is_fine()
        + h.used_only_before_reentry_is_fine()
        + h.fact_about_a_local_is_fine()
        + h.local_closure_is_fine()
        + h.pending_is_read_properly()
        + h.position_across_dyn_given_the_field_is_flagged();
    let _ = h.crate_collection_across_await_in_shared_method_is_flagged();
    let _ = h.await_in_exclusive_method_is_fine();
    let _ = h.non_panicking_reuse_is_fine();
    let _ = h.binding_reassigned_after_reentry_is_fine();
    let _ = u32::from(h.other_field_is_fine());
    let _ = h.field_replaced_after_reentry_is_fine();
    let _ = h.index_in_a_condition_is_flagged()
        + h.removal_in_a_scrutinee_is_flagged()
        + h.reentry_before_a_break_inside_the_statement_is_flagged()
        + h.configured_impl_method_named_by_its_type_is_flagged()
        + h.configured_impl_method_named_by_its_trait_is_flagged()
        + h.crate_type_with_the_cell_behind_a_box_is_flagged()
        + h.frozen_crate_type_given_self_is_flagged()
        + h.crate_collection_given_shared_self_is_flagged(peek)
        + h.reentry_in_a_returning_branch_is_fine(false)
        + h.returned_reentry_is_fine(false, |_| 0)
        + h.reentry_before_a_break_out_of_the_block_is_fine(false)
        + h.configured_callee_inside_a_stored_closure_is_fine()
        + h.frozen_crate_type_in_shared_method_is_fine()
        + h.fn_pointer_taking_shared_self_is_fine(peek)
        + h.field_handed_out_shared_by_coercion_is_fine(peek_items)
        + h.unconfigured_impl_method_is_fine();
    for task in &h.tasks {
        task(&h.vm);
    }
}
