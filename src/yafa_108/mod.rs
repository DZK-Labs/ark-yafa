//! This library implements the Yafa-108 curve, which is a CP6 curve that embeds curve25519,
//! with a bit security level of 108.
//!
//! The name denotes that it was generated using the Cocks--Pinch method for the
//! embedding degree 6. The main feature of this curve is that the scalar field
//! equals the base field of the curve25519 curve.

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
