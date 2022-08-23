use ark_ec::{models::short_weierstrass::SWCurveConfig, PairingEngine};
use ark_ff::{
    biginteger::BigInt,
    fields::{BitIteratorBE, Field},
    CyclotomicMultSubgroup, One,
};

use crate::yafa_108::{Fq, Fq3, Fq6, Fr};

pub mod g1;
pub use self::g1::{G1Affine, G1Projective};

pub mod g2;
pub use self::g2::{G2Affine, G2Projective};

#[cfg(test)]
mod tests;

pub type GT = Fq6;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Yafa;

impl PairingEngine for Yafa {
    type Fr = Fr;
    type G1Projective = G1Projective;
    type G1Affine = G1Affine;
    type G1Prepared = G1Affine;
    type G2Projective = G2Projective;
    type G2Affine = G2Affine;
    type G2Prepared = G2Affine;
    type Fq = Fq;
    type Fqe = Fq3;
    type Fqk = Fq6;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        let mut result = Self::Fqk::one();
        for &(ref p, ref q) in i {
            result *= &Yafa::ate_miller_loop(p, q);
        }
        result
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        Some(Yafa::final_exponentiation(r))
    }
}

impl Yafa {
    pub fn ate_pairing(p: &G1Affine, q: &G2Affine) -> GT {
        Yafa::final_exponentiation(&Yafa::ate_miller_loop(p, q))
    }

    fn ate_miller_loop(p: &G1Affine, q: &G2Affine) -> Fq6 {
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
/// 372992445403421670819442132301882547570428511124131210420164786679808777528556315411341
pub const ATE_LOOP_COUNT: [u64; 5] = [
    0x55792745f9e7038d,
    0xe5e55257bec5baaf,
    0x2678d0333aa42ad6,
    0x47e345c9559a9a4a,
    0xc000009c,
];

/// FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG = true
pub const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;

/// FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0 =
/// 19533731073207388012378338389467147321952799507004027602120097593540492719109989264168193066083831322222399106936053605103753388004158213616940547405614189990736827562331161
pub const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInt<9> = BigInt::new([
    0x1dc822b236540c19,
    0xb32840686d9b3d30,
    0x5fe7da3eac81af7e,
    0x7d72f866fc867962,
    0x0a1ac894fb4e064d,
    0x5b5259e8ae3c723e,
    0x8e7db93d38b0bd87,
    0x0da5a4ce11e1c642,
    0x1437ebc18a1ed9a7,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W1 =
/// 600746412144547199455308148290259497480968394648582220356252480504778992479853140369801482324344
pub const FINAL_EXPONENT_LAST_CHUNK_W1: BigInt<9> = BigInt::new([
    0xe2bb1306781fc978,
    0xdc7abc2bf0661157,
    0x17cdb96cb71d1e72,
    0x7f162be473a7d636,
    0x4800007536090574,
    0x0,
    0x0,
    0x0,
    0x0,
]);
