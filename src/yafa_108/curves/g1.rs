use ark_ec::{
    models::{short_weierstrass::SWCurveConfig, CurveConfig},
    short_weierstrass::{Affine, Projective},
    AffineRepr, CurveGroup,
};
use ark_ff::{Field, MontFp};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;

use crate::yafa_108::{Fq, Fr};

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
    /// 600746357873917083993205832519317177466334726496174638989123538389804293458724181786500830516591
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0xed29225ad506d56f,
        0x4360c489b2a69d35,
        0xd497d8292a77e57a,
        0x2e5856438790676e,
        0x4800000815ea74e0
    ];

    /// COFACTOR^(-1) mod r =
    /// 29525709254955688747989059845541707767137581572857490857189440336502437830996
    const COFACTOR_INV: Fr =
        MontFp!("29525709254955688747989059845541707767137581572857490857189440336502437830996");
}

impl SWCurveConfig for Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = Fq::ZERO;

    /// COEFF_B = 5645376
    const COEFF_B: Fq = MontFp!("5645376");

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const GENERATOR: G1Affine = G1Affine::new_unchecked(G1_GENERATOR_X, G1_GENERATOR_Y);
}

/// G1_GENERATOR_X =
/// 34224882689487856897895479418745239985629366350609644482426457479519960220365127786522998842348723780309585364044838207298438594720820768669193775884940942820440641403063878
///
/// This is point (1, 2251806090139824558188356038789175761391042867695566584791619687372201424747844952708255342828857479046478251177789121867581508195508201416938853880984464178574376900833404)
/// removing the cofactor
pub const G1_GENERATOR_X: Fq = MontFp!("34224882689487856897895479418745239985629366350609644482426457479519960220365127786522998842348723780309585364044838207298438594720820768669193775884940942820440641403063878");

/// G1_GENERATOR_Y =
/// 531252248722162526044089347491731336050196735203224076668566468302402949866896014787473405742067724694482254210724441982798591548014930288574006468717544900780067911325970
pub const G1_GENERATOR_Y: Fq = MontFp!("531252248722162526044089347491731336050196735203224076668566468302402949866896014787473405742067724694482254210724441982798591548014930288574006468717544900780067911325970");
