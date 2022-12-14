use ark_ec::pairing::{MillerLoopOutput, PairingOutput};
use ark_ec::{models::short_weierstrass::SWCurveConfig, pairing::Pairing};
use ark_ff::{
    biginteger::BigInt,
    fields::{BitIteratorBE, Field},
    CyclotomicMultSubgroup, One,
};
use itertools::Itertools;

use crate::yafa_108::{Fq, Fq3, Fq6, Fr};

pub mod g1;
pub use self::g1::{G1Affine, G1Prepared, G1Projective};

pub mod g2;
pub use self::g2::{G2Affine, G2Prepared, G2Projective};

#[cfg(test)]
mod tests;

pub type GT = Fq6;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Yafa;

impl Pairing for Yafa {
    type ScalarField = Fr;
    type G1 = G1Projective;
    type G1Affine = G1Affine;
    type G1Prepared = G1Prepared;
    type G2 = G2Projective;
    type G2Affine = G2Affine;
    type G2Prepared = G2Prepared;
    type TargetField = Fq6;

    fn multi_miller_loop(
        a: impl IntoIterator<Item = impl Into<Self::G1Prepared>>,
        b: impl IntoIterator<Item = impl Into<Self::G2Prepared>>,
    ) -> MillerLoopOutput<Self> {
        let mut result = Self::TargetField::one();
        a.into_iter().zip_eq(b).for_each(|(p, q)| {
            let (p, q) = (p.into(), q.into());
            result *= &Yafa::ate_miller_loop(&p.0, &q.0);
        });

        MillerLoopOutput(result)
    }

    fn final_exponentiation(r: MillerLoopOutput<Self>) -> Option<PairingOutput<Self>> {
        Some(PairingOutput(Yafa::final_exponentiation(&r.0)))
    }
}

impl Yafa {
    pub fn ate_pairing(p: &G1Affine, q: &G2Affine) -> GT {
        Yafa::final_exponentiation(&Yafa::ate_miller_loop(p, q))
    }

    fn ate_miller_loop(p: &G1Affine, q: &G2Affine) -> GT {
        let px = p.x;
        let py = p.y;
        let qx = q.x;
        let qy = q.y;
        let mut py_twist_squared = TWIST.square();
        py_twist_squared.mul_assign_by_fp(&py);

        let mut old_rx;
        let mut old_ry;
        let mut rx = qx;
        let mut ry = qy;
        let mut f = Fq6::one();

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        // cp6_782_param_p (skipping leading zeros) in MSB to LSB order
        for bit in BitIteratorBE::without_leading_zeros(ATE_LOOP_COUNT).skip(1) {
            old_rx = rx;
            old_ry = ry;

            let old_rx_square = old_rx.square();
            let old_rx_square_3 = old_rx_square.double() + &old_rx_square;
            let old_rx_square_3_a = old_rx_square_3 + &g2::Parameters::COEFF_A;
            let old_ry_double_inverse = old_ry.double().inverse().unwrap();

            let gamma = old_rx_square_3_a * &old_ry_double_inverse;
            let gamma_twist = gamma * &TWIST;
            let gamma_old_rx = gamma * &old_rx;
            let mut gamma_twist_px = gamma_twist;
            gamma_twist_px.mul_assign_by_fp(&px);

            let x = py_twist_squared;
            let y = gamma_old_rx - &old_ry - &gamma_twist_px;
            let ell_rr_at_p: Fq6 = Fq6::new(x, y);

            rx = gamma.square() - &old_rx.double();
            ry = gamma * &(old_rx - &rx) - &old_ry;
            f = f.square() * &ell_rr_at_p;

            if bit {
                old_rx = rx;
                old_ry = ry;

                let gamma = (old_ry - &qy) * &((old_rx - &qx).inverse().unwrap());
                let gamma_twist = gamma * &TWIST;
                let gamma_qx = gamma * &qx;
                let mut gamma_twist_px = gamma_twist;
                gamma_twist_px.mul_assign_by_fp(&px);

                let x = py_twist_squared;
                let y = gamma_qx - &qy - &gamma_twist_px;
                let ell_rq_at_p: Fq6 = Fq6::new(x, y);

                rx = gamma.square() - &old_rx - &qx;
                ry = gamma * &(old_rx - &rx) - &old_ry;
                f = f * &ell_rq_at_p;
            }
        }
        f
    }

    fn final_exponentiation(value: &Fq6) -> GT {
        let value_inv = value.inverse().unwrap();
        let value_to_first_chunk = Yafa::final_exponentiation_first(value, &value_inv);
        let value_inv_to_first_chunk = Yafa::final_exponentiation_first(&value_inv, value);
        Yafa::final_exponentiation_last(&value_to_first_chunk, &value_inv_to_first_chunk)
    }

    fn final_exponentiation_first(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        let mut elt_q3 = elt.clone();
        elt_q3.frobenius_map(3);
        // elt_q3_over_elt = elt^(q^3-1)
        let elt_q3_over_elt = elt_q3 * elt_inv;
        // alpha = elt^((q^3-1) * q)
        let mut alpha = elt_q3_over_elt.clone();
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha * &elt_q3_over_elt
    }

    fn final_exponentiation_last(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        let mut elt_q = elt.clone();
        elt_q.frobenius_map(1);

        let w1_part = elt_q.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_W1);
        let w0_part = if FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            elt_inv.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)
        } else {
            elt.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)
        };

        w1_part * &w0_part
    }
}

/// TWIST = (0, 1, 0)
pub const TWIST: Fq3 = Fq3::new(Fq::ZERO, Fq::ONE, Fq::ZERO);

/// ATE_IS_LOOP_COUNT_NEG = false
pub const ATE_IS_LOOP_COUNT_NEG: bool = false;

/// ATE_LOOP_COUNT =
/// 372992428555672686789935698172304228806337918473348441569462718938730304377195952806182
pub const ATE_LOOP_COUNT: [u64; 5] = [
    0x55792745f9e71926,
    0xe5e55257bec5baaf,
    0x2678d0333aa42ad6,
    0xc7e345c9559a9a4a,
    0xc000000a,
];

/// FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG = true
pub const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;

/// FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0 =
/// 19533729308556160109794825562832707092134991522349430861538344636908041560785866966036323918363729219562215110665533950746637708295543215802259030565596621823690292271103113
pub const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInt<9> = BigInt::new([
    0x7d919d468ad5bc89,
    0x851c91204bcf98dc,
    0x1ccd844b0de55ff8,
    0x36bd1be87c106c82,
    0xcfe6ae45c3d21c71,
    0x96d3c4ddb45de345,
    0xe20009f1de3c5d1f,
    0xc534aab128f9dd7a,
    0x1437eba2e55524f8,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W1 =
/// 600746357873917083993205832519317177466334726496174638989123538389804293458724181786507272967556
pub const FINAL_EXPONENT_LAST_CHUNK_W1: BigInt<9> = BigInt::new([
    0xed29225c5506d584,
    0x4360c489b2a69d35,
    0xd497d8292a77e57a,
    0x2e5856438790676e,
    0x4800000815ea74e0,
    0x0,
    0x0,
    0x0,
    0x0,
]);
