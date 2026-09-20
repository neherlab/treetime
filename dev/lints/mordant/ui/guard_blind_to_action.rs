// A permission-flavored guard must read everything its gated action touches.

struct Sched {
    queue: Vec<u32>,
    conns: Vec<u32>,
}

impl Sched {
    fn can_donate(&self) -> bool {
        self.queue.is_empty()
    }

    fn can_flush(&self) -> bool {
        self.queue.is_empty() && self.conns.is_empty()
    }

    // Flagged: `detach` touches `conns`, which `can_donate` never reads.
    fn donate(&mut self) {
        if !self.can_donate() {
            return;
        }
        self.detach();
    }

    // Fine: `can_flush` reads both fields `detach` touches.
    fn flush(&mut self) {
        if self.can_flush() {
            self.detach();
        }
    }

    fn detach(&mut self) {
        self.queue.clear();
        self.conns = Vec::new();
    }

    // Flagged: the gated call's result is let-bound, the shape the real
    // donation bug used.
    fn donate_binding(&mut self) -> usize {
        if !self.can_donate() {
            return 0;
        }
        let n = self.drain();
        n
    }

    fn drain(&mut self) -> usize {
        self.conns = Vec::new();
        self.conns.len()
    }

    // Fine: the gated call does not mutate.
    fn report(&mut self) {
        if !self.can_donate() {
            return;
        }
        let _ = self.count();
    }

    fn count(&self) -> usize {
        self.queue.len() + self.conns.len()
    }

    // Fine: `detach` sits under a condition of its own; `can_donate` opened the
    // method, it did not approve this call.
    fn donate_when_idle(&mut self) {
        if !self.can_donate() {
            return;
        }
        if self.conns.is_empty() {
            self.detach();
        }
    }

    // Fine: same, in the then-branch form.
    fn donate_if_idle(&mut self) {
        if self.can_donate() {
            if self.conns.is_empty() {
                self.detach();
            }
        }
    }

    fn can_touch_conns(&self) -> bool {
        self.conns.is_empty()
    }

    // Fine: two stacked guards approve `detach`, and between them they read
    // both fields it touches.
    fn donate_checked(&mut self) {
        if !self.can_donate() {
            return;
        }
        if !self.can_touch_conns() {
            return;
        }
        self.detach();
    }

    // Fine: same, one guard nested in the other's then-branch.
    fn donate_checked_nested(&mut self) {
        if self.can_donate() {
            if self.can_touch_conns() {
                self.detach();
            }
        }
    }

    fn can_log(&self) -> bool {
        self.queue.len() < 8
    }

    // Flagged once, naming both guards: together they still never read
    // `conns`.
    fn donate_logged(&mut self) {
        if !self.can_donate() {
            return;
        }
        if self.can_log() {
            self.detach();
        }
    }

    // Fine: `can_shed` hands `self` to a free function, so what it reads is
    // unknown; a guard whose coverage cannot be computed accuses nothing.
    fn can_shed(&self) -> bool {
        policy_allows(self)
    }

    fn shed(&mut self) {
        if !self.can_shed() {
            return;
        }
        self.detach();
    }

    // Fine: same one call away — `can_shed_now` is followed into `can_shed`,
    // and the escape there makes this coverage unknown too.
    fn can_shed_now(&self) -> bool {
        self.queue.is_empty() && self.can_shed()
    }

    fn shed_now(&mut self) {
        if self.can_shed_now() {
            self.detach();
        }
    }
}

fn policy_allows(s: &Sched) -> bool {
    s.queue.is_empty()
}

fn main() {
    let mut s = Sched {
        queue: vec![1],
        conns: vec![2],
    };
    let _ = s.can_donate();
    let _ = s.can_flush();
    s.donate();
    let _ = s.donate_binding();
    s.flush();
    s.report();
    s.donate_when_idle();
    s.donate_if_idle();
    s.donate_checked();
    s.donate_checked_nested();
    s.donate_logged();
    let _ = s.can_touch_conns() && s.can_log();
    s.shed();
    s.shed_now();
    s.detach();
}
