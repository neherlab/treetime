fn main() {
    let _ = Positional::new(1, 2, 3, 4, 5, 6);
    let _ = Built::builder().finish(1);
}

pub struct Positional {
    a: u8,
    b: u8,
    c: u8,
    d: u8,
    e: u8,
    f: u8,
}

impl Positional {
    pub fn new(a: u8, b: u8, c: u8, d: u8, e: u8, f: u8) -> Self {
        Self { a, b, c, d, e, f }
    }
}

pub struct Built {
    a: u8,
    b: u8,
    c: u8,
    d: u8,
    e: u8,
    f: u8,
}

impl Built {
    pub fn builder() -> BuiltBuilder {
        BuiltBuilder
    }

    fn __orig_new(a: u8) -> Self {
        Self { a, b: a, c: a, d: a, e: a, f: a }
    }
}

pub struct BuiltBuilder;

impl BuiltBuilder {
    pub fn finish(self, a: u8) -> Built {
        Built::__orig_new(a)
    }
}
