// Several `bool` parameters filled with bare `true`/`false` at a call site:
// the flag names live only in the signature, so a swapped call compiles.

// Flagged: two of the three calls pass both flags as bare literals.
fn render(text: &str, wrap: bool, color: bool) -> usize {
    text.len() + usize::from(wrap) + usize::from(color)
}

struct File {
    len: usize,
}

impl File {
    // Flagged: a method; the receiver is not a flag, the three bools are, and
    // both the method-call and the path-call form leave two of them unnamed.
    fn open(&mut self, create: bool, truncate: bool, append: bool) -> usize {
        self.len + usize::from(create) + usize::from(truncate) + usize::from(append)
    }

    // Fine: every call names what it passes, by a binding or by a comment.
    fn lock(&mut self, shared: bool, blocking: bool) -> usize {
        self.len + usize::from(shared) + usize::from(blocking)
    }

    // Fine: one bool parameter has nothing to be swapped with.
    fn sync(&mut self, data_only: bool) -> usize {
        self.len + usize::from(data_only)
    }
}

// Fine: exported, so its signature is not the crate's alone to change.
pub fn spawn(cmd: &str, inherit_stdout: bool, inherit_stderr: bool) -> usize {
    cmd.len() + usize::from(inherit_stdout) + usize::from(inherit_stderr)
}

trait Sink {
    fn write(&mut self, flush: bool, sync: bool) -> usize;
}

// Fine: a trait method's signature is dictated by the trait.
impl Sink for File {
    fn write(&mut self, flush: bool, sync: bool) -> usize {
        self.len + usize::from(flush) + usize::from(sync)
    }
}

// Fine: one bool among other parameters, called with a literal everywhere.
fn add(list: &mut Vec<u8>, enabled: bool, item: u8) {
    if enabled {
        list.push(item);
    }
}

// Fine: two bool parameters, but at most one is ever a bare literal per call.
fn connect(host: &str, tls: bool, keepalive: bool) -> usize {
    host.len() + usize::from(tls) + usize::from(keepalive)
}

fn main() {
    let wrap = true;
    let color = false;
    let _ = render("a", true, false);
    let _ = render("b", false, false);
    let _ = render("c", wrap, color);

    let mut f = File { len: 0 };
    let append = false;
    let _ = f.open(true, false, append);
    let _ = File::open(&mut f, false, false, false);

    let shared = true;
    let blocking = false;
    let _ = f.lock(shared, blocking);
    let _ = f.lock(/* shared */ true, /* blocking */ false);
    let _ = f.lock(
        true, // shared
        false,
    );

    let _ = f.sync(true);
    let _ = f.sync(false);

    let _ = spawn("ls", true, false);
    let _ = f.write(true, true);
    let _ = Sink::write(&mut f, false, false);

    let mut list = Vec::new();
    add(&mut list, true, 1);
    add(&mut list, true, 2);
    add(&mut list, false, 3);

    let keepalive = true;
    let _ = connect("h", true, keepalive);
    let _ = connect("h", cfg!(unix), cfg!(windows));
}
