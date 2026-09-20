// `pub` items that nothing uses. rustc's `dead_code` passes over them
// because `pub` alone makes them reachable from outside the crate.

pub fn never_called_is_flagged() {}

pub fn called_is_fine() {}

pub fn taken_as_value_is_fine() {}

pub fn only_calls_itself_is_flagged(n: u32) {
    if n > 0 {
        only_calls_itself_is_flagged(n - 1);
    }
}

pub struct NeverNamed;

// The impl names its type, which is not a use; the type is reported and
// its method is left out as part of it.
impl NeverNamed {
    pub fn method_of_unused_type_is_folded_in(&self) {}
}

#[derive(Clone, Debug, PartialEq)]
pub struct OnlyDerived {
    pub n: u32,
}

pub struct Built {
    pub n: u32,
}

impl Built {
    pub fn called_method_is_fine(&self) -> u32 {
        self.n
    }

    pub fn never_called_method_is_flagged(&self) -> u32 {
        self.n + 1
    }

    pub const NEVER_READ_IS_FLAGGED: u32 = 3;

    pub const READ_IS_FINE: u32 = 4;
}

pub trait NeverImplemented {
    fn required(&self);
    fn provided(&self) {}
}

pub trait Implemented {
    fn called_through_trait_is_fine(&self);
    fn never_called_trait_method_is_flagged(&self) {}
}

impl Implemented for Built {
    fn called_through_trait_is_fine(&self) {}
}

pub enum NeverMatched {
    A,
}

pub enum Matched {
    A,
    B,
}

pub const UNREAD: u32 = 1;

pub static UNREAD_STATIC: u32 = 2;

pub type NeverWritten = u32;

pub type Written = u64;

// An impl written against an alias applies to the type behind it, so the
// alias is in use even if nothing else spells it.
pub type ImplTarget = Built;

impl ImplTarget {
    pub fn through_alias_is_fine(&self) -> u32 {
        self.n
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn reached_by_symbol_is_fine() {}

#[allow(unknown_lints, unused_pub)]
pub fn allowed_is_fine() {}

pub(crate) fn crate_private_is_left_to_rustc() {}

mod closed {
    // Not reachable from outside: rustc's `dead_code` owns it.
    pub fn behind_private_module_is_left_to_rustc() {}
}

pub mod open {
    pub fn nested_never_called_is_flagged() {}

    pub fn nested_called_is_fine() {}
}

// Generated code enters a crate through `include!`; an unused item there
// is the generator's to drop, not this file's.
include!("unused_pub.included");

macro_rules! make {
    () => {
        pub fn produced_by_macro_is_fine() {}
    };
}
make!();

fn generic<T: Implemented>(t: &T) {
    t.called_through_trait_is_fine();
}

fn main() {
    called_is_fine();
    let f = taken_as_value_is_fine;
    f();
    let b = Built { n: Built::READ_IS_FINE };
    let _ = b.called_method_is_fine() + b.through_alias_is_fine();
    generic(&b);
    let _ = match Matched::A {
        Matched::A => 0,
        Matched::B => 1,
    };
    let _: Written = 0;
    open::nested_called_is_fine();
    crate_private_is_left_to_rustc();
    closed::behind_private_module_is_left_to_rustc();
}
