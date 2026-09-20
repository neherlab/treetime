// A transmute or pointer cast producing a type that something already converts
// the same source into is flagged outside that type's module and impls; the
// type's own code, identity transmutes, and types nothing converts into are not.

mod level {
    #[derive(Clone, Copy, PartialEq, Debug)]
    #[repr(u8)]
    pub enum Level {
        Low = 0,
        Mid = 1,
        High = 2,
    }

    impl std::convert::TryFrom<u8> for Level {
        type Error = u8;
        fn try_from(n: u8) -> Result<Self, u8> {
            match n {
                0 => Ok(Level::Low),
                1 => Ok(Level::Mid),
                2 => Ok(Level::High),
                _ => Err(n),
            }
        }
    }

    // Fine: the type's own module decides the layout; this sits beside the check.
    pub fn trusted(n: u8) -> Level {
        unsafe { core::mem::transmute::<u8, Level>(n) }
    }
}

mod code {
    #[derive(Clone, Copy)]
    #[repr(u16)]
    pub enum Code {
        Ok = 0,
        Retry = 1,
    }

    impl Code {
        pub fn from_raw(n: u32) -> Option<Code> {
            match n {
                0 => Some(Code::Ok),
                1 => Some(Code::Retry),
                _ => None,
            }
        }
    }
}

mod fd {
    #[repr(transparent)]
    pub struct Fd(pub(crate) i32);

    impl Fd {
        // An infallible conversion still decides how an i32 becomes an Fd.
        pub fn new(raw: i32) -> Fd {
            Fd(if raw < 0 { -1 } else { raw })
        }

        // Fine as a conversion source: `unsafe fn` promises no check.
        pub unsafe fn adopt(raw: i64) -> Fd {
            Fd(raw as i32)
        }
    }
}

mod meters {
    #[derive(Clone, Copy)]
    #[repr(transparent)]
    pub struct Meters(pub(crate) u32);

    impl From<u32> for Meters {
        fn from(n: u32) -> Meters {
            Meters(n.min(40_000_000))
        }
    }
}

mod raw {
    // Nothing converts into this: a transmute is the only door there is.
    #[derive(Clone, Copy)]
    #[repr(u8)]
    pub enum Raw {
        A = 0,
        B = 1,
    }

    // Only an unsafe constructor: it promises nothing a transmute skips.
    #[repr(transparent)]
    pub struct Slot(pub(crate) u8);

    impl Slot {
        pub unsafe fn from_raw(n: u8) -> Slot {
            Slot(n)
        }
    }

    pub struct View<'a> {
        pub bytes: &'a [u8],
    }

    impl<'a> View<'a> {
        pub fn parse(bytes: &'a [u8]) -> Option<View<'a>> {
            (!bytes.is_empty()).then_some(View { bytes })
        }
    }
}

mod px {
    use std::marker::PhantomData;

    pub struct Screen;
    pub struct World;

    #[repr(transparent)]
    pub struct Px<S>(pub(crate) f32, PhantomData<S>);

    // Only screen pixels are made from a bare float; nothing converts into
    // `Px<World>`.
    impl From<f32> for Px<Screen> {
        fn from(v: f32) -> Px<Screen> {
            Px(v.max(0.0), PhantomData)
        }
    }
}

mod hdr {
    #[repr(transparent)]
    pub struct Hdr(pub(crate) [u8; 4]);

    // A borrow-taking conversion: the usual shape for buffer newtypes.
    impl<'a> From<&'a [u8; 4]> for Hdr {
        fn from(b: &'a [u8; 4]) -> Hdr {
            Hdr([b[0] & 0x7f, b[1], b[2], b[3]])
        }
    }
}

mod name {
    // An unsized view type: its checked constructor returns a reference.
    #[repr(transparent)]
    pub struct NameRef(pub(crate) [u8]);

    impl NameRef {
        pub fn new(b: &[u8]) -> Option<&NameRef> {
            if b.contains(&0) {
                return None;
            }
            Some(unsafe { &*(b as *const [u8] as *const NameRef) })
        }
    }

    // A wrapper over a borrow whose constructor takes an elided lifetime.
    #[derive(Clone, Copy)]
    #[repr(transparent)]
    pub struct Str<'a>(pub(crate) &'a [u8]);

    impl Str<'_> {
        pub fn new(b: &[u8]) -> Option<Str<'_>> {
            std::str::from_utf8(b).ok().map(|_| Str(b))
        }
    }
}

mod gate {
    // A struct whose `Option<Self>` constructor checks the field it stores:
    // `unchecked_construction` names that check, so this lint leaves a
    // `mem::transmute` into it alone.
    pub struct Gate {
        pub(crate) v: u8,
    }

    impl Gate {
        pub fn new(v: u8) -> Option<Gate> {
            if v < 2 { Some(Gate { v }) } else { None }
        }
    }
}

use std::convert::TryFrom;

use code::Code;
use fd::Fd;
use gate::Gate;
use hdr::Hdr;
use level::Level;
use meters::Meters;
use name::{NameRef, Str};
use px::{Px, Screen, World};
use raw::{Raw, Slot, View};

// Flagged: `Level::try_from` exists for exactly this pair.
fn level_from_wire(n: u8) -> Level {
    unsafe { core::mem::transmute::<u8, Level>(n) }
}

// Flagged: the conversion takes u32; a value transmute from any integer skips it.
fn code_from_wire(n: u32) -> Code {
    unsafe { core::mem::transmute::<u16, Code>(n as u16) }
}

// Flagged: transmute_copy reads through the reference and reinterprets.
fn code_copied(n: &u16) -> Code {
    unsafe { core::mem::transmute_copy::<u16, Code>(n) }
}

// Flagged: an infallible constructor is still the conversion this skips.
fn fd_from_env(raw: i32) -> Fd {
    unsafe { core::mem::transmute::<i32, Fd>(raw) }
}

// Flagged: a pointer cast between the exact pair `From<u32> for Meters` covers.
fn meters_in_place(n: &u32) -> Meters {
    unsafe { *(n as *const u32 as *const Meters) }
}

// Flagged: `.cast()` is the same reinterpretation.
fn meters_read(p: *const u32) -> Meters {
    unsafe { p.cast::<Meters>().read() }
}

// Flagged: transmuting references compares the pointees.
fn meters_ref(n: &u32) -> &Meters {
    unsafe { core::mem::transmute::<&u32, &Meters>(n) }
}

// Flagged: the transmute is `&u32 -> &Meters` whatever the argument coerced from.
fn meters_boxed(v: &Box<u32>) -> &Meters {
    unsafe { core::mem::transmute::<&u32, &Meters>(v) }
}

// Flagged: `.cast()` through an auto-deref'd `&*const u32` still casts `*const u32`.
fn meters_each(ptrs: &[*const u32]) -> u32 {
    let mut sum = 0;
    for p in ptrs.iter() {
        sum += unsafe { p.cast::<Meters>().read() }.0;
    }
    sum
}

// Flagged: `From<f32> for Px<Screen>` covers exactly this instantiation.
fn px_screen(v: f32) -> Px<Screen> {
    unsafe { core::mem::transmute::<f32, Px<Screen>>(v) }
}

// Fine: nothing converts an f32 into `Px<World>`; the `Px<Screen>` impl is not it.
fn px_world(v: f32) -> Px<World> {
    unsafe { core::mem::transmute::<f32, Px<World>>(v) }
}

// Flagged: `From<&[u8; 4]> for Hdr` is the conversion, borrow or not.
fn hdr_value(b: &[u8; 4]) -> Hdr {
    unsafe { core::mem::transmute::<[u8; 4], Hdr>(*b) }
}

// Flagged: the same conversion against a pointer view of the same bytes.
fn hdr_view(b: &[u8; 4]) -> &Hdr {
    unsafe { &*(b as *const [u8; 4] as *const Hdr) }
}

// Flagged: `NameRef::new` returns the checked reference this cast makes unchecked.
fn name_view(b: &[u8]) -> &NameRef {
    unsafe { &*(b as *const [u8] as *const NameRef) }
}

// Flagged: `Str::new(&[u8])` takes an elided lifetime and is still the conversion.
fn str_from_bytes(b: &[u8]) -> Str<'_> {
    unsafe { core::mem::transmute::<&[u8], Str<'_>>(b) }
}

// Fine here: `unchecked_construction` reports this one and names the check.
fn gate_from_wire(n: u8) -> Gate {
    unsafe { core::mem::transmute::<u8, Gate>(n) }
}

// Fine: a byte buffer viewed as the type; pointer casts need the exact source.
fn meters_from_bytes(bytes: &[u8]) -> Meters {
    unsafe { bytes.as_ptr().cast::<Meters>().read_unaligned() }
}

// Fine: nothing converts a u8 into Raw, so there is no check to skip.
fn raw_from_wire(n: u8) -> Raw {
    unsafe { core::mem::transmute::<u8, Raw>(n) }
}

// Fine: Slot's only constructor is unsafe.
fn slot_from_wire(n: u8) -> Slot {
    unsafe { core::mem::transmute::<u8, Slot>(n) }
}

// Fine: only the lifetime changes; the value already went through `parse`.
fn extend<'a>(v: View<'a>) -> View<'static> {
    unsafe { core::mem::transmute::<View<'a>, View<'static>>(v) }
}

// Fine: a trait impl of the type is its own code wherever it is written.
impl Default for Level {
    fn default() -> Self {
        unsafe { core::mem::transmute::<u8, Level>(0) }
    }
}

// Fine: an integer-to-pointer cast reinterprets no pointee.
fn from_addr(addr: usize) -> *const Meters {
    addr as *const Meters
}

// Fine: casting to the same pointee, mutability aside.
fn constness(p: *mut Meters) -> *const Meters {
    p as *const Meters
}

// Fine: a type with interior mutability is viewed in place; `From` would make
// a new cell, not a view of this one.
fn as_atomic(n: &usize) -> usize {
    use std::sync::atomic::{AtomicUsize, Ordering};
    unsafe { (*(n as *const usize).cast::<AtomicUsize>()).load(Ordering::Relaxed) }
}

fn main() {
    let _ = level_from_wire(1);
    let _ = code_from_wire(1);
    let _ = code_copied(&1);
    let _ = fd_from_env(2);
    let _ = meters_in_place(&3);
    let _ = meters_read(&3);
    let _ = meters_ref(&3);
    let _ = meters_boxed(&Box::new(3)).0;
    let _ = meters_each(&[&3u32 as *const u32]);
    let _ = px_screen(1.0).0;
    let _ = px_world(1.0).0;
    let _ = hdr_value(&[0; 4]).0;
    let _ = hdr_view(&[0; 4]).0;
    let _ = name_view(b"x").0.len();
    let _ = str_from_bytes(b"x").0;
    let _ = gate_from_wire(0).v;
    let _ = meters_from_bytes(&[0, 0, 0, 0]);
    let _ = raw_from_wire(0);
    let _ = slot_from_wire(0);
    let _ = extend(View::parse(b"x").unwrap()).bytes;
    let _ = Level::default();
    let _ = from_addr(8);
    let _ = constness(core::ptr::null_mut());
    let _ = as_atomic(&0);
    let _ = level::trusted(0);
    let _ = Code::from_raw(0);
    let _ = Fd::new(0).0;
    let _ = unsafe { Fd::adopt(0) }.0;
    let _ = Meters::from(1).0;
    let _ = Level::try_from(9);
    let _ = unsafe { Slot::from_raw(0) }.0;
    let _ = Px::<Screen>::from(1.0).0;
    let _ = Hdr::from(&[0; 4]).0;
    let _ = NameRef::new(b"x").map(|n| n.0.len());
    let _ = Str::new(b"x").map(|s| s.0);
    let _ = Gate::new(0).map(|g| g.v);
}
