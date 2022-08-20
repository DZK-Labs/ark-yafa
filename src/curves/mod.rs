use ark_ec::{models::short_weierstrass::SWCurveConfig, PairingEngine};
use ark_ff::{
    biginteger::BigInt,
    fields::{BitIteratorBE, Field},
    CyclotomicMultSubgroup, One,
};

use crate::{Fq, Fq3, Fq6, Fr};

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
/// 958852482080108523503261136209276035742162691825013326170216721598135503624581
pub const ATE_LOOP_COUNT: [u64; 5] = [
    0x5579276279e71985,
    0xe5e55257bec5baaf,
    0x2678d0333aa42ad6,
    0x47e345c9559a9a4a,
    0x8,
];

/// FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG = true
pub const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;

/// FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0 =
/// 782412617055576810900917744291620294143560296654224301172315017705450234607288874792247667942436918253025244581142084110845414191684600371997508793548252111273
pub const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInt<9> = BigInt::new([
    0xc6f366e4d5e7e1a9,
    0x3ef739d45f765667,
    0xd381ac296a431b80,
    0x8f9fbc2488d2c335,
    0xae2929e3168b8148,
    0x7eb2b012cf097664,
    0xa5eec9b3ceeb34d5,
    0xffb82a96512f44dd,
    0xe3f2,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W1 =
/// 24062559822862625313599379840463328384057194034656457084274246263984232490014561516
pub const FINAL_EXPONENT_LAST_CHUNK_W1: BigInt<9> = BigInt::new([
    0xbd733d2fbf3354ec,
    0xbc296ea016247281,
    0xccc28d5cce4871f4,
    0x528f2e9cc94ee931,
    0x32bc0,
    0x0,
    0x0,
    0x0,
    0x0,
]);
