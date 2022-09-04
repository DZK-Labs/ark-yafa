use ark_ec::{pairing::Pairing, AffineRepr};
use ark_ff::{biginteger::BigInt, fields::Field, CyclotomicMultSubgroup, One};
use ark_std::{cfg_chunks_mut, vec::Vec};
use itertools::Itertools;

use crate::yafa_146::{Fq, Fq12, Fq2, Fr};

pub mod g1;
pub use self::g1::{G1Affine, G1Prepared, G1Projective};

pub mod g2;
pub use self::g2::{G2Affine, G2Prepared, G2Projective};

use ark_ec::pairing::{MillerLoopOutput, PairingOutput};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(test)]
mod tests;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Yafa;

impl Yafa {
    // Evaluate the line function at point p.
    fn ell(f: &mut Fq12, coeffs: &g2::EllCoeff<Fq2>, p: &G1Prepared) {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;
        let c2 = coeffs.2;
        let (px, py) = p.0.xy().unwrap();

        // This is a divisive twist
        c0.mul_assign_by_fp(py);
        c1.mul_assign_by_fp(px);
        f.mul_by_034(&c0, &c1, &c2);
    }
}

impl Pairing for Yafa {
    type ScalarField = Fr;
    type G1 = G1Projective;
    type G1Affine = G1Affine;
    type G1Prepared = G1Prepared;
    type G2 = G2Projective;
    type G2Affine = G2Affine;
    type G2Prepared = G2Prepared;
    type TargetField = Fq12;

    fn multi_miller_loop(
        a: impl IntoIterator<Item = impl Into<Self::G1Prepared>>,
        b: impl IntoIterator<Item = impl Into<Self::G2Prepared>>,
    ) -> MillerLoopOutput<Self> {
        let mut pairs = a
            .into_iter()
            .zip_eq(b)
            .filter_map(|(p, q)| {
                let (p, q) = (p.into(), q.into());
                match !p.is_zero() && !q.is_zero() {
                    true => Some((p, q.ell_coeffs.into_iter())),
                    false => None,
                }
            })
            .collect::<Vec<_>>();

        let f = cfg_chunks_mut!(pairs, 4)
            .map(|pairs| {
                let mut f = Self::TargetField::one();
                for i in TATE_LOOP_COUNT.iter().skip(1) {
                    f.square_in_place();
                    for (p, coeffs) in pairs.iter_mut() {
                        Self::ell(&mut f, &coeffs.next().unwrap(), &p);
                    }
                    match i {
                        1 | -1 => {
                            for (p, coeffs) in pairs.iter_mut() {
                                Self::ell(&mut f, &coeffs.next().unwrap(), &p);
                            }
                        }
                        0 => continue,
                        _ => unreachable!(),
                    }
                }
                f
            })
            .product::<Self::TargetField>();
        MillerLoopOutput(f)
    }

    fn final_exponentiation(f: MillerLoopOutput<Self>) -> Option<PairingOutput<Self>> {
        let f = f.0;

        // f1 = r.conjugate() = f^(p^6)
        let mut f1 = f;
        f1.conjugate_in_place();

        f.inverse().map(|mut f2| {
            // Easy part
            // f2 = f^(-1);
            // r = f^(p^6 - 1)
            let mut r = f1 * &f2;

            // f2 = f^(p^6 - 1)
            f2 = r;
            // r = f^((p^6 - 1)(p^2))
            r.frobenius_map(2);

            // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
            // r = f^((p^6 - 1)(p^2 + 1))
            r *= &f2;

            // Hard part
            let mut elt_q = r.clone();
            let mut elt_q2 = r.clone();
            let mut elt_q3 = r.clone();
            elt_q.frobenius_map(1);
            elt_q2.frobenius_map(2);
            elt_q3.frobenius_map(3);

            let w3_part = elt_q3.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_W3);
            let w2_part = elt_q2.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_W2);
            let w1_part = elt_q.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_W1);
            let w0_part = r.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_W0);

            PairingOutput(w3_part * &w2_part * &w1_part * &w0_part)
        })
    }
}

/// TWIST = (0, 1, 0)
pub const TWIST: Fq2 = Fq2::new(Fq::ZERO, Fq::ONE);

/// TATE_LOOP_COUNT =
/// 57896044618658097711785492504343953926634992332820282019728792003956564819949
pub const TATE_LOOP_COUNT: &'static [i8] = &[
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0,
    1,
];

/// FINAL_EXPONENT_LAST_CHUNK_W0 =
/// 27816679719392039590916985961594891329592754889347470862124370788722235197910030055580499179007884620636974433034081916025003956486380300643703133308392747535275498881809581
pub const FINAL_EXPONENT_LAST_CHUNK_W0: BigInt<9> = BigInt::new([
    0xbc41361f419dfcad,
    0x7e898c59e9f6e67a,
    0xa53054d2e155420a,
    0xf01d5c325a7e8db9,
    0x0062cc2d9828db4b,
    0xa9ab464e0f3efb97,
    0x4498a2a9fcde4960,
    0xf05a7a282a9da654,
    0x1ccaaec847ce3756,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W1 =
/// 19533734280962826151300710596734343137033362535294523292134692601750277351215937028705950010165895260665221949740604041087518331761534784981220098087966726908412430392012572
pub const FINAL_EXPONENT_LAST_CHUNK_W1: BigInt<9> = BigInt::new([
    0x4a8171212b68af1c,
    0x6e9fd1b1cd44b35d,
    0x10677198e7742d97,
    0x13b522f9221f979a,
    0xa9d14c2dbb4cf3db,
    0xaed043007079fa1e,
    0x7b9bdbddd59fe860,
    0x0ee6a96afcb80cfc,
    0x1437ebf93e40e187,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W2 =
/// 18787556099692390437676225874885260726423443207869755270416397806622396448189003069574413752571064036218150696957580977806429906239508024997419791567936987076281415850010664
pub const FINAL_EXPONENT_LAST_CHUNK_W2: BigInt<9> = BigInt::new([
    0x3e418200ad402828,
    0xeaac6e307c80e9ac,
    0x993fbe8fefb574c8,
    0x7e876f3ddda8ba7a,
    0xf8ff3fbf0bb49620,
    0x63883dfe09f064b0,
    0x8fd38d1af07c11e5,
    0xee91cf069f66a1eb,
    0x1372344f3c10e50b,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W3 =
/// 600746510796850913325441891520214792552399427644760978176030334322268158184417708567430611936132
pub const FINAL_EXPONENT_LAST_CHUNK_W3: BigInt<9> = BigInt::new([
    0xba0cbd2b5bed2384,
    0xdab3b5db47bc3720,
    0x6303cbb620293bf6,
    0xd93c6dabba8ba44f,
    0x4800013b93e07009,
    0x0,
    0x0,
    0x0,
    0x0,
]);
