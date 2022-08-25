use ark_algebra_test_templates::{
    curves::{curve_tests, sw_tests},
    generate_bilinearity_test, generate_g1_test, generate_g2_test,
    msm::*,
};
use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::{Field, One, PrimeField};
use ark_std::{rand::Rng, test_rng};
use core::ops::MulAssign;

use crate::yafa_108::*;

generate_g1_test!(yafa; curve_tests; sw_tests;);
generate_g2_test!(yafa; curve_tests; sw_tests;);
generate_bilinearity_test!(Yafa, Fq6);