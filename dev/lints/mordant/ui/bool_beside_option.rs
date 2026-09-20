// A bool field written only beside an Option field, `true` with `Some`, is that Option's `is_some()`.

// Flagged: `connected` is `peer.is_some()`: false/None at construction,
// true/Some in `open`, false/None in `close`, and nowhere else.
struct Conn {
    peer: Option<String>,
    connected: bool,
    bytes: u64,
}

impl Conn {
    fn new() -> Self {
        Conn {
            peer: None,
            connected: false,
            bytes: 0,
        }
    }

    fn open(&mut self, peer: String) {
        self.peer = Some(peer);
        self.connected = true;
    }

    fn close(&mut self) {
        self.connected = false;
        self.peer = None;
        self.bytes = 0;
    }

    fn peer(&self) -> Option<&str> {
        self.peer.as_deref()
    }

    fn is_up(&self) -> bool {
        self.connected
    }
}

// Flagged: the reverse polarity, `stale` is `fresh.is_none()`; reading the
// Option through `as_mut` cannot change that.
struct Cache {
    fresh: Option<u32>,
    stale: bool,
}

impl Cache {
    fn filled(v: u32) -> Self {
        Cache {
            fresh: Some(v),
            stale: false,
        }
    }

    fn invalidate(&mut self) {
        self.fresh = None;
        self.stale = true;
    }

    fn bump(&mut self) {
        if let Some(v) = self.fresh.as_mut() {
            *v += 1;
        }
    }
}

// Flagged: a derived `Default` builds `None`/`false`, one more site that
// agrees with `is_some()`; a derived `Clone` copies the pair as it is; a loop
// that only breaks out of itself does not part the two writes around it.
#[derive(Clone, Default)]
struct Link {
    up: Option<u32>,
    linked: bool,
}

impl Link {
    fn connect(&mut self, v: u32) {
        self.up = Some(v);
        for _ in 0..v {
            if self.up.is_some() {
                break;
            }
        }
        self.linked = true;
    }

    fn drop_link(&mut self) {
        self.up = None;
        self.linked = false;
    }
}

// Fine: the derived `Default` builds `fresh: None, stale: false`, which is
// not `fresh.is_none()`.
#[derive(Default)]
struct Memo {
    fresh: Option<u32>,
    stale: bool,
}

impl Memo {
    fn filled(v: u32) -> Self {
        Memo {
            fresh: Some(v),
            stale: false,
        }
    }

    fn invalidate(&mut self) {
        self.fresh = None;
        self.stale = true;
    }
}

// Fine: `done` is also set on its own, so it is not `error.is_some()`.
struct Scan {
    error: Option<String>,
    done: bool,
}

impl Scan {
    fn new() -> Self {
        Scan {
            error: None,
            done: false,
        }
    }

    fn fail(&mut self, e: String) {
        self.error = Some(e);
        self.done = true;
    }

    fn finish(&mut self) {
        self.done = true;
    }
}

// Fine: the Option is drained with `take()`, a write the flag does not follow.
struct Exit {
    status: Option<i32>,
    exited: bool,
}

impl Exit {
    fn new() -> Self {
        Exit {
            status: None,
            exited: false,
        }
    }

    fn on_exit(&mut self, code: i32) {
        self.exited = true;
        self.status = Some(code);
    }

    fn reap(&mut self) -> Option<i32> {
        self.status.take()
    }
}

// Fine: the Option is drained alone through a `&mut` binding split off `self`.
struct Split {
    slot: Option<i32>,
    held: bool,
}

impl Split {
    fn new() -> Self {
        Split {
            slot: None,
            held: false,
        }
    }

    fn hold(&mut self, v: i32) {
        self.slot = Some(v);
        self.held = true;
    }

    fn drain(&mut self) -> Option<i32> {
        let Split { slot, .. } = self;
        slot.take()
    }
}

// Fine: the same, the binding spelled `ref mut`.
struct Drain {
    item: Option<i32>,
    full: bool,
}

impl Drain {
    fn new() -> Self {
        Drain {
            item: None,
            full: false,
        }
    }

    fn fill(&mut self, v: i32) {
        self.item = Some(v);
        self.full = true;
    }

    fn empty(&mut self) {
        let Drain { ref mut item, .. } = *self;
        *item = None;
    }
}

// Fine: positional fields are built through the constructor fn, whose calls
// are not followed.
struct Pair(Option<u32>, bool);

impl Pair {
    fn new() -> Self {
        Pair(Some(1), false)
    }

    fn set(&mut self, v: u32) {
        self.0 = Some(v);
        self.1 = true;
    }

    fn clear(&mut self) {
        self.0 = None;
        self.1 = false;
    }
}

// Fine: the match arms without braces are exclusive, not beside each other.
enum Ev {
    Up(u32),
    Flag,
    Down,
}

struct Peer {
    addr: Option<u32>,
    known: bool,
}

impl Peer {
    fn new() -> Self {
        Peer {
            addr: None,
            known: false,
        }
    }

    fn on(&mut self, e: Ev) {
        match e {
            Ev::Up(a) => self.addr = Some(a),
            Ev::Flag => self.known = true,
            Ev::Down => {
                self.addr = None;
                self.known = false;
            }
        }
    }
}

// Fine: one block writes both polarities of each, `false` beside `None` and
// `true` beside `Some`, while construction is `Some` beside `false`.
struct Swap {
    cur: Option<u32>,
    empty: bool,
}

impl Swap {
    fn new(v: u32) -> Self {
        Swap {
            cur: Some(v),
            empty: false,
        }
    }

    fn cycle(&mut self, v: u32) {
        self.cur = None;
        self.empty = false;
        let _ = v.leading_zeros();
        self.cur = Some(v);
        self.empty = true;
    }
}

// Fine: the two writes in one block go to different instances.
struct Slot {
    task: Option<u32>,
    taken: bool,
}

impl Slot {
    fn new() -> Self {
        Slot {
            task: None,
            taken: false,
        }
    }

    fn park(&mut self) {
        self.task = None;
        self.taken = false;
    }
}

fn hand_off(from: &mut Slot, to: &mut Slot, v: u32) {
    to.task = Some(v);
    from.taken = true;
}

// Fine: the `?` between the two writes can leave with only the flag set.
struct Fetch {
    body: Option<u32>,
    started: bool,
}

impl Fetch {
    fn new() -> Self {
        Fetch {
            body: None,
            started: false,
        }
    }

    fn run(&mut self, v: u32) -> Option<()> {
        self.started = true;
        let got = v.checked_sub(1)?;
        self.body = Some(got);
        Some(())
    }

    fn reset(&mut self) {
        self.body = None;
        self.started = false;
    }
}

// Fine: the flag is written from a computed value once.
struct Probe {
    addr: Option<u32>,
    reachable: bool,
}

impl Probe {
    fn new() -> Self {
        Probe {
            addr: None,
            reachable: false,
        }
    }

    fn resolve(&mut self, addr: u32, ok: bool) {
        self.addr = Some(addr);
        self.reachable = ok;
    }
}

// Fine: the flag is only ever `false`; a constant is not a copy of anything.
struct Idle {
    job: Option<u32>,
    busy: bool,
}

impl Idle {
    fn new() -> Self {
        Idle {
            job: None,
            busy: false,
        }
    }

    fn park(&mut self) {
        self.job = None;
        self.busy = false;
    }
}

// Fine: the two writes are in different blocks, so one can run without the other.
struct Lazy {
    value: Option<u32>,
    loaded: bool,
}

impl Lazy {
    fn new() -> Self {
        Lazy {
            value: None,
            loaded: false,
        }
    }

    fn load(&mut self, v: Option<u32>) {
        self.loaded = true;
        if let Some(v) = v {
            self.value = Some(v);
        }
    }
}

// Fine: exported fields can be written by other crates.
pub struct Public {
    pub handle: Option<u32>,
    pub attached: bool,
}

impl Public {
    pub fn new() -> Self {
        Public {
            handle: None,
            attached: false,
        }
    }

    pub fn attach(&mut self, h: u32) {
        self.handle = Some(h);
        self.attached = true;
    }
}

fn main() {
    let mut c = Conn::new();
    c.open("p".to_owned());
    let _ = (c.peer(), c.is_up());
    c.close();
    c.bytes += 1;

    let mut k = Cache::filled(1);
    k.bump();
    k.invalidate();
    let _ = k.stale;

    let mut n = Link::default().clone();
    n.connect(1);
    n.drop_link();
    let _ = (n.linked, n.up);

    let mut m = Memo::default();
    let _ = Memo::filled(1);
    m.invalidate();
    let _ = (m.stale, m.fresh);

    let mut s = Scan::new();
    s.fail("e".to_owned());
    s.finish();
    let _ = (s.done, s.error.is_some());

    let mut e = Exit::new();
    e.on_exit(0);
    let _ = (e.exited, e.reap());

    let mut t = Split::new();
    t.hold(1);
    let _ = (t.held, t.drain());

    let mut d = Drain::new();
    d.fill(1);
    d.empty();
    let _ = (d.full, d.item);

    let mut two = Pair::new();
    two.set(2);
    two.clear();
    let _ = (two.1, two.0);

    let mut r = Peer::new();
    r.on(Ev::Up(1));
    r.on(Ev::Flag);
    r.on(Ev::Down);
    let _ = (r.known, r.addr);

    let mut w = Swap::new(1);
    w.cycle(2);
    let _ = (w.empty, w.cur);

    let (mut a, mut b) = (Slot::new(), Slot::new());
    hand_off(&mut a, &mut b, 1);
    a.park();
    let _ = (a.taken, b.task);

    let mut f = Fetch::new();
    let _ = f.run(1);
    f.reset();
    let _ = (f.started, f.body);

    let mut p = Probe::new();
    p.resolve(1, true);
    let _ = (p.reachable, p.addr);

    let mut i = Idle::new();
    i.park();
    let _ = (i.busy, i.job);

    let mut l = Lazy::new();
    l.load(Some(1));
    let _ = (l.loaded, l.value);

    let mut u = Public::new();
    u.attach(1);
    let _ = (u.attached, u.handle);
}
