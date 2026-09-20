// Bool fields only ever assigned together encode one state.

// Flagged: running/done are co-assigned in two methods and nowhere else.
struct Task {
    running: bool,
    done: bool,
    retries: u32,
}

impl Task {
    fn start(&mut self) {
        self.running = true;
        self.done = false;
    }

    fn finish(&mut self) {
        self.running = false;
        self.done = true;
    }
}

// Fine: `a` has a lone write, so the fields are independent.
struct Independent {
    a: bool,
    b: bool,
}

impl Independent {
    fn set_a(&mut self) {
        self.a = true;
    }

    fn set_both(&mut self) {
        self.a = true;
        self.b = true;
    }
}

// Fine: the lone write to `open` goes through a Box, and is still a write.
struct Boxed {
    open: bool,
    ready: bool,
}

impl Boxed {
    fn up(&mut self) {
        self.open = true;
        self.ready = true;
    }

    fn down(&mut self) {
        self.open = false;
        self.ready = false;
    }
}

fn close_boxed(b: &mut Box<Boxed>) {
    b.open = false;
}

// Fine: the lone write to `dirty` is a compound assignment.
struct Ored {
    dirty: bool,
    seen: bool,
}

impl Ored {
    fn reset(&mut self) {
        self.dirty = false;
        self.seen = false;
    }

    fn mark(&mut self) {
        self.dirty = true;
        self.seen = true;
    }

    fn touch(&mut self, changed: bool) {
        self.dirty |= changed;
    }
}

fn main() {
    let mut t = Task {
        running: false,
        done: false,
        retries: 0,
    };
    t.start();
    t.finish();
    t.retries += 1;

    let mut i = Independent { a: false, b: false };
    i.set_a();
    i.set_both();

    let mut b = Box::new(Boxed {
        open: false,
        ready: false,
    });
    b.up();
    b.down();
    close_boxed(&mut b);

    let mut o = Ored {
        dirty: false,
        seen: false,
    };
    o.reset();
    o.mark();
    o.touch(true);
}
