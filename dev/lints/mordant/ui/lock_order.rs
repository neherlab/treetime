// Two locks taken in both orders across a crate is a deadlock waiting for
// the interleaving; one consistent order is fine.

use std::sync::Mutex;

struct S {
    a: Mutex<u32>,
    b: Mutex<u32>,
}

impl S {
    fn ab(&self) -> u32 {
        let ga = self.a.lock().unwrap();
        let gb = self.b.lock().unwrap();
        *ga + *gb
    }

    // Flagged together with `ab`: the reverse order.
    fn ba(&self) -> u32 {
        let gb = self.b.lock().unwrap();
        let ga = self.a.lock().unwrap();
        *ga + *gb
    }

    // Fine: the first guard is dropped before the second lock.
    fn sequential(&self) -> u32 {
        let ga = self.a.lock().unwrap();
        let v = *ga;
        drop(ga);
        let gb = self.b.lock().unwrap();
        v + *gb
    }
}

struct T {
    c: Mutex<u32>,
    d: Mutex<u32>,
}

impl T {
    // Fine: both sites agree on the order.
    fn first(&self) -> u32 {
        let gc = self.c.lock().unwrap();
        let gd = self.d.lock().unwrap();
        *gc + *gd
    }

    fn second(&self) -> u32 {
        let gc = self.c.lock().unwrap();
        let gd = self.d.lock().unwrap();
        *gc * *gd
    }
}

struct Pool {
    x: Mutex<u32>,
    y: Mutex<u32>,
}

impl Pool {
    fn xy(&self) -> u32 {
        let gx = self.x.lock().unwrap();
        let gy = self.y.lock().unwrap();
        *gx + *gy
    }
}

struct Cache {
    x: Mutex<u32>,
    y: Mutex<u32>,
}

impl Cache {
    // Fine: `Cache::y` and `Cache::x` are not `Pool::x` and `Pool::y`; a
    // field name shared by two types is not one lock.
    fn yx(&self) -> u32 {
        let gy = self.y.lock().unwrap();
        let gx = self.x.lock().unwrap();
        *gx + *gy
    }
}

struct Deferred {
    e: Mutex<u32>,
    f: Mutex<u32>,
}

impl Deferred {
    fn ef(&self) -> u32 {
        let ge = self.e.lock().unwrap();
        let gf = self.f.lock().unwrap();
        *ge + *gf
    }

    // Fine: the reverse order is only inside a closure built while `f` is
    // held; here it runs after the drop, and in general when its holder says.
    fn fe_later(&self) -> u32 {
        let gf = self.f.lock().unwrap();
        let job = || *self.e.lock().unwrap();
        let v = *gf;
        drop(gf);
        v + job()
    }
}

struct Spin {
    g: Mutex<u32>,
    h: Mutex<u32>,
}

impl Spin {
    fn gh(&self) -> u32 {
        let gg = self.g.lock().unwrap();
        let gh = self.h.lock().unwrap();
        *gg + *gh
    }

    // Flagged together with `gh`: the reverse order sits directly in a
    // `loop {}` body, which is a block of its own.
    fn hg_in_loop(&self) -> u32 {
        loop {
            let gh = self.h.lock().unwrap();
            let gg = self.g.lock().unwrap();
            if *gg > 0 {
                return *gg + *gh;
            }
        }
    }
}

fn main() {
    let s = S {
        a: Mutex::new(1),
        b: Mutex::new(2),
    };
    let _ = s.ab() + s.ba() + s.sequential();
    let t = T {
        c: Mutex::new(3),
        d: Mutex::new(4),
    };
    let _ = t.first() + t.second();
    let p = Pool {
        x: Mutex::new(5),
        y: Mutex::new(6),
    };
    let c = Cache {
        x: Mutex::new(7),
        y: Mutex::new(8),
    };
    let _ = p.xy() + c.yx();
    let d = Deferred {
        e: Mutex::new(9),
        f: Mutex::new(10),
    };
    let _ = d.ef() + d.fe_later();
    let sp = Spin {
        g: Mutex::new(11),
        h: Mutex::new(12),
    };
    let _ = sp.gh() + sp.hg_in_loop();
}
