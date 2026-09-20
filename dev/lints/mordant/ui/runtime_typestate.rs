// A bool field guarded at the entry of 2+ methods is a runtime ordering
// invariant; one guard alone is not a pattern.

struct Conn {
    ready: bool,
    sent: u32,
}

impl Conn {
    fn send(&mut self) {
        if !self.ready {
            return;
        }
        self.sent += 1;
    }

    fn flush(&mut self) {
        if !self.ready {
            return;
        }
        self.sent = 0;
    }

    fn reset(&mut self) {
        self.sent = 0;
    }

    // The transition that makes `ready` a state and not a setting.
    fn connect(&mut self) {
        self.ready = true;
    }
}

struct Once {
    armed: bool,
    fired: u32,
}

impl Once {
    fn fire(&mut self) {
        if self.armed {
            return;
        }
        self.fired += 1;
    }
}

struct Printer {
    minify: bool,
    out: String,
}

impl Printer {
    // Fine: `minify` is set once, in the literal, and never written again;
    // bailing on it is a mode, not an order of operations.
    fn newline(&mut self) {
        if self.minify {
            return;
        }
        self.out.push('\n');
    }

    fn space(&mut self) {
        if self.minify {
            return;
        }
        self.out.push(' ');
    }
}

fn main() {
    let mut c = Conn {
        ready: true,
        sent: 0,
    };
    c.connect();
    c.send();
    c.flush();
    c.reset();

    let mut o = Once {
        armed: false,
        fired: 0,
    };
    o.fire();

    let mut p = Printer {
        minify: true,
        out: String::new(),
    };
    p.newline();
    p.space();
}
