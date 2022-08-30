use ark_algebra_test_templates::curves::*;
use ark_ec::AffineCurve;

use super::*;

#[test]
fn test_generator() {
    let generator = EdwardsAffine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_montgomery_conversion() {
    // We want to emphasize that this Montgomery curve is not Curve25519.
    montgomery_conversion_test::<EdwardsParameters>();
}
