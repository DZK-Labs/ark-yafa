use ark_ec::{
    models::CurveConfig,
    short_weierstrass::{Affine, Projective, SWCurveConfig},
    AffineCurve,
};
use ark_ff::{BitIteratorBE, Field, MontFp};
use ark_std::{vec, vec::Vec, One, Zero};

use crate::yafa_146::{Fq, Fq2, Fr, ATE_LOOP_COUNT};

pub type G2Affine = Affine<Parameters>;
pub type G2Projective = Projective<Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct G2Prepared {
    // Stores the coefficients of the line evaluations as calculated in
    // https://eprint.iacr.org/2013/722.pdf
    pub ell_coeffs: Vec<EllCoeff<Fq2>>,
    pub infinity: bool,
}

pub(crate) type EllCoeff<F> = (F, F, F);

#[derive(Clone, Copy, Debug)]
struct G2HomProjective {
    x: Fq2,
    y: Fq2,
    z: Fq2,
}

impl Default for G2Prepared {
    fn default() -> Self {
        Self::from(G2Affine::prime_subgroup_generator())
    }
}

impl From<G2Affine> for G2Prepared {
    fn from(q: G2Affine) -> Self {
        let two_inv = Fq::one().double().inverse().unwrap();
        match q.is_zero() {
            true => G2Prepared {
                ell_coeffs: vec![],
                infinity: true,
            },
            false => {
                let mut ell_coeffs = vec![];
                let mut r = G2HomProjective {
                    x: q.x,
                    y: q.y,
                    z: Fq2::one(),
                };

                for i in BitIteratorBE::new(ATE_LOOP_COUNT).skip(1) {
                    ell_coeffs.push(doubling_step(&mut r, &two_inv));

                    if i {
                        ell_coeffs.push(addition_step(&mut r, &q));
                    }
                }

                Self {
                    ell_coeffs,
                    infinity: false,
                }
            }
        }
    }
}

impl G2Prepared {
    pub fn is_zero(&self) -> bool {
        self.infinity
    }
}

fn doubling_step(r: &mut G2HomProjective, two_inv: &Fq) -> EllCoeff<Fq2> {
    // Formula for line function when working with
    // homogeneous projective coordinates.

    let mut a = r.x * &r.y;
    a.mul_assign_by_fp(two_inv);
    let b = r.y.square();
    let c = r.z.square();
    let e = Parameters::COEFF_B * &(c.double() + &c);
    let f = e.double() + &e;
    let mut g = b + &f;
    g.mul_assign_by_fp(two_inv);
    let h = (r.y + &r.z).square() - &(b + &c);
    let i = e - &b;
    let j = r.x.square();
    let e_square = e.square();

    r.x = a * &(b - &f);
    r.y = g.square() - &(e_square.double() + &e_square);
    r.z = b * &h;

    // This is a divisive twist
    (-h, j.double() + &j, i)
}

fn addition_step(r: &mut G2HomProjective, q: &G2Affine) -> EllCoeff<Fq2> {
    // Formula for line function when working with
    // homogeneous projective coordinates.

    let theta = r.y - &(q.y * &r.z);
    let lambda = r.x - &(q.x * &r.z);
    let c = theta.square();
    let d = lambda.square();
    let e = lambda * &d;
    let f = r.z * &c;
    let g = r.x * &d;
    let h = e + &f - &g.double();
    r.x = lambda * &h;
    r.y = theta * &(g - &h) - &(e * &r.y);
    r.z *= &e;
    let j = theta * &q.x - &(lambda * &q.y);

    // This is a divisive twist
    (lambda, -theta, j)
}

impl CurveConfig for Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;

    /// COFACTOR =
    /// 208944723538136275842934825487721719944357363218924936521565367644832496954484236054395518
    /// 029919319763436921041473784490674990167148701594903075840663632617329304570348691066252329
    /// 59381521895831602203960656224597442267949372368242654610238587983559434866653274386897396
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0x97eb95cb492a99f4,
        0xce1097b616174157,
        0x69ff3af87115ad37,
        0xaad0f9fa6c57adbb,
        0xe89b7dd49d7a672e,
        0x66d37be389df638a,
        0xf1d388886c9fcd6a,
        0x1e7cb685aa19089f,
        0x11c7fe3ecab2bfa4,
        0x28813796458eeeab,
        0x9a129319a58fe700,
        0x5f5465426db28467,
        0x64f89dae62b24f99,
        0xa200058c197e205
    ];

    /// COFACTOR^(-1) mod r =
    /// 52342251624758656962424635542600246432390708264947618736639330869844643473789
    const COFACTOR_INV: Fr =
        MontFp!("52342251624758656962424635542600246432390708264947618736639330869844643473789");
}

impl SWCurveConfig for Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq2 = Fq2::new(Fq::ZERO, Fq::ZERO);

    /// COEFF_B = G1::COEFF_B / u
    const COEFF_B: Fq2 = Fq2::new(Fq::ZERO, MontFp!("565347594816"));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const GENERATOR: G2Affine = G2Affine::new_unchecked(G2_GENERATOR_X, G2_GENERATOR_Y);
    // This is the point (u + 4, 56268234608410584509335952523906779591923024718516727382353769583
    // 9042556486544408093949957902074246729482291972632315061261331134959992728979959547596836023
    // 8639466912619219*u + 1549244716205334453943662290654166012055803483162652679442772269869164
    // 5581725493987465742793060036277134832165591892442449296288650854155894467334221898486129239
    // 512217529073)
    // multiplied by the cofactor
    //
    // The number of points on the curve
    // 1209707303679711869521701003500290607415757496777039484715548003544001706965971215131424169
    // 7637611826882570078699336838344331626954875886166837539349405006589346957661356360935844425
    // 3845876818888589332253391608769436909946181275415785184259524560205574601021756234710113763
    // 2689382898715418246998285567880568710254458669018456424437988969476952804
}

const G2_GENERATOR_X: Fq2 = Fq2::new(G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
const G2_GENERATOR_Y: Fq2 = Fq2::new(G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

/// G2_GENERATOR_X_C0 =
/// 15961858629678698558604263579817401801634171552903891025991809467507692432772682588292817173076251185109360722117296071793256403954089060793212318126143130609622773815531730
pub const G2_GENERATOR_X_C0: Fq = MontFp!("15961858629678698558604263579817401801634171552903891025991809467507692432772682588292817173076251185109360722117296071793256403954089060793212318126143130609622773815531730");

/// G2_GENERATOR_X_C1 =
/// 10832407179955914211847371507868395694831167875939483416974996320568813633110183337679716246669139630694946155994134164573070651417188923217059486369595504440405929405295070
pub const G2_GENERATOR_X_C1: Fq = MontFp!("10832407179955914211847371507868395694831167875939483416974996320568813633110183337679716246669139630694946155994134164573070651417188923217059486369595504440405929405295070");

/// G2_GENERATOR_Y_C0 =
/// 16822710565940580249223862576115209605427445899327931090428260033671198185206263769780729157475880870824879495366633381853655228190440466094750660128619185403553786706812293
pub const G2_GENERATOR_Y_C0: Fq = MontFp!("16822710565940580249223862576115209605427445899327931090428260033671198185206263769780729157475880870824879495366633381853655228190440466094750660128619185403553786706812293");

/// G2_GENERATOR_Y_C1 =
/// 1104233349431442808145740601052820001361225084204411121082176275621513935336905661254958546223874109479566034987625127812138372112294386337684996796070026833127481450661093
pub const G2_GENERATOR_Y_C1: Fq = MontFp!("1104233349431442808145740601052820001361225084204411121082176275621513935336905661254958546223874109479566034987625127812138372112294386337684996796070026833127481450661093");
