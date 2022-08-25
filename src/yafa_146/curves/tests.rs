#![allow(unused)]

use ark_algebra_test_templates::{
    curves::{curve_tests, sw_tests},
    generate_bilinearity_test, generate_g1_test, generate_g2_test,
    msm::*,
};
use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::{Field, One, PrimeField};
use ark_std::{rand::Rng, test_rng};
use core::ops::MulAssign;

use crate::yafa_146::*;

generate_g1_test!(yafa; curve_tests; sw_tests;);
generate_g2_test!(yafa; curve_tests; sw_tests;);
//generate_bilinearity_test!(Yafa, Fq12);

#[test]
fn test_bilinearity() {
    let mut rng = test_rng();
    let a: G1Projective = rng.gen();
    let b: G2Projective = rng.gen();
    let s: Fr = rng.gen();

    let mut sa = a;
    sa.mul_assign(s);
    let mut sb = b;
    sb.mul_assign(s);

    let ans1 = Yafa::pairing(sa, b);
    let ans2 = Yafa::pairing(a, sb);
    let ans3 = Yafa::pairing(a, b).pow(s.into_bigint());

    assert_eq!(ans1, ans2);
    assert_eq!(ans2, ans3);

    assert_ne!(ans1, Fq12::one());
    assert_ne!(ans2, Fq12::one());
    assert_ne!(ans3, Fq12::one());

    assert_eq!(ans1.pow(Fr::characteristic()), Fq12::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq12::one());
    assert_eq!(ans3.pow(Fr::characteristic()), Fq12::one());
}
