use ark_ff::{fields::*, MontFp};

use crate::yafa_146::*;

pub type Fq2 = Fp2<Fq2Config>;

pub struct Fq2Config;

impl Fp2Config for Fq2Config {
    type Fp = Fq;

    /// NONRESIDUE = 59
    const NONRESIDUE: Fq = MontFp!("59");

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[
        // NONRESIDUE**(((q^0) - 1) / 2)
        Fq::ONE,
        // NONRESIDUE**(((q^1) - 1) / 2)
        MontFp!("-1"),
    ];
}
