//! N-API command orchestrations.
//!
//! Each module owns one operation's request shape and its complete read-run-project-write workflow:
//! it converts the request into core values, calls `treetime::<op>::pipeline::run` directly, and
//! projects the core result through the shared `app-output` encoders into the default `--output-all`
//! file set. The N-API client depends on the core and the shared concern crates, never on another
//! client.

pub mod ancestral;
pub mod clock;
pub mod mugration;
pub mod optimize;
pub mod prune;
pub mod support;
pub mod timetree;
