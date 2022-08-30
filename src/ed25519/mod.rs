//! This library implements curve25519 for the purpose of benchmark.
//! This is not a pairing-friendly curve, so its instantiation will be straightforward.

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
