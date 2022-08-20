#![cfg_attr(not(feature = "std"), no_std)]
#![deny(
    warnings,
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms
)]
#![forbid(unsafe_code)]

//! This library implements the Yafa curve, which is a CP6 curve that embeds curve25519.
//!
//! The name denotes that it was generated using the Cocks--Pinch method for the
//! embedding degree 6. The main feature of this curve is that the scalar field
//! equals the base field of the curve25519 curve.
//!
//! Curve information:
//! * Base field: q = 1393127037143584247199893356983111562901745574937902324694961530768859597608968253224455327034513790598090225418769269547013027125294147720154141059771190988081
//! * Scalar field: r = 57896044618658097711785492504343953926634992332820282019728792003956564819949
//!
//! G1 curve equation: y^2 = x^3 + ax + b, where
//! * a = 11,
//! * b = 258778877286921393099430952356951813749468875819541194566590982411954936341028897879371735235881815122771684547010326496988837784319176832293682272090515709037
//!
//! G2 curve equation: y^2 = x^3 + Ax + B
//! * A = Fq3(0, 0, 11)
//! * B = Fq3(60313575868966829693953761960246825440666484139148490842577744993785104533381370224178433525672385154308079179575052372851161376922649714922222873453290823245, 0, 0)

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
