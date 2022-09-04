use crate::yafa_108::*;
use ark_algebra_test_templates::*;

test_group!(g1; G1Projective; sw);
test_group!(g2; G2Projective; sw);
test_group!(pairing_output; ark_ec::pairing::PairingOutput<Yafa>; msm);
test_pairing!(pairing; crate::yafa_108::Yafa);
