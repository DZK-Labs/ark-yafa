use ark_ec::{
    models::{short_weierstrass::SWCurveConfig, CurveConfig},
    short_weierstrass::{Affine, Projective},
    AffineRepr, CurveGroup,
};
use ark_ff::{Field, MontFp};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;

use crate::yafa_146::{Fq, Fr};

pub type G1Affine = Affine<Parameters>;
pub type G1Projective = Projective<Parameters>;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct G1Prepared(pub G1Affine);

impl From<G1Affine> for G1Prepared {
    fn from(other: G1Affine) -> Self {
        G1Prepared(other)
    }
}

impl From<G1Projective> for G1Prepared {
    fn from(q: G1Projective) -> Self {
        q.into_affine().into()
    }
}

impl<'a> From<&'a G1Affine> for G1Prepared {
    fn from(other: &'a G1Affine) -> Self {
        G1Prepared(*other)
    }
}

impl<'a> From<&'a G1Projective> for G1Prepared {
    fn from(q: &'a G1Projective) -> Self {
        q.into_affine().into()
    }
}

impl G1Prepared {
    pub fn is_zero(&self) -> bool {
        self.0.is_identity()
    }
}

impl Default for G1Prepared {
    fn default() -> Self {
        G1Prepared(G1Affine::generator())
    }
}

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl CurveConfig for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;

    /// COFACTOR =
    /// 600746510796850913325441891520214792552399427644760978176030334322268158184417708567424169484347
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0xba0cbd29dbed203b,
        0xdab3b5db47bc3720,
        0x6303cbb620293bf6,
        0xd93c6dabba8ba44f,
        0x4800013b93e07009
    ];

    /// COFACTOR^(-1) mod r =
    /// 54736692702506479266359545630135535008734683389654089840601071154286118289550
    const COFACTOR_INV: Fr =
        MontFp!("54736692702506479266359545630135535008734683389654089840601071154286118289550");
}

impl SWCurveConfig for Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = Fq::ZERO;

    /// COEFF_B = 33355508094144
    const COEFF_B: Fq = MontFp!("33355508094144");

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const GENERATOR: G1Affine = G1Affine::new_unchecked(G1_GENERATOR_X, G1_GENERATOR_Y);
}

/// G1_GENERATOR_X =
/// 21041479060334059994917916352561995599681496683378184338495044148554334865922285043275465528741964306980383475001455404610600675285973607295911295491046433550562665948850109
///
/// This is point (4, 9980063903480773747485374411020714173816276815606417861225439654554390745719005727274055955009294624578956029865160575967717579468338487582878433342359310352690204717350854)
/// removing the cofactor
pub const G1_GENERATOR_X: Fq = MontFp!("21041479060334059994917916352561995599681496683378184338495044148554334865922285043275465528741964306980383475001455404610600675285973607295911295491046433550562665948850109");

/// G1_GENERATOR_Y =
/// 1841897213774345215194080213544536275781413636928742996668240879239998760310806751411541978966337027054045888026250482482159097701195947482630556667627310741134774424209442
pub const G1_GENERATOR_Y: Fq = MontFp!("1841897213774345215194080213544536275781413636928742996668240879239998760310806751411541978966337027054045888026250482482159097701195947482630556667627310741134774424209442");
