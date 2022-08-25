use ark_ff::{
    fields::fp3::{Fp3, Fp3Config},
    Field, MontFp,
};

use crate::yafa_108::Fq;

pub type Fq3 = Fp3<Fq3Config>;

pub struct Fq3Config;

impl Fp3Config for Fq3Config {
    type Fp = Fq;

    /// NONRESIDUE = 11
    /// It means that (x^3 - 11) is irreducible in F_{q^3}
    /// We also need "u" to be a nonresidue so it can be the twist
    const NONRESIDUE: Fq = MontFp!("11");

    /// t = (q^3 - 1) / pow(2, 3)
    const TWO_ADICITY: u32 = 3;


    #[rustfmt::skip]
    /// t_minus_one_div_2 = 2629663266467176870871002443124951184583329704009094068531710073760496
    /// 870167236050677745303113103585055358358466362691791465526226187396896633733264074557712305
    /// 436133938840477125775498584163824749636229813623010348686986267561874346621561015933614721
    /// 074794956143620396128239207898300512812188460816197391783882508336407354718487241526208132
    /// 900643660312817351203052784131581431362236654758101701422135319332524450434246494841700451
    /// 339645762359070996386180138129126352319421850061827361591505926342246237933107823656807
    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = &[
        0x9df37fa6adfa3f67,
        0x47e94b65dd587942,
        0xc80333b6e0d0ac78,
        0x2b6257722f37046f,
        0xfd092fd4d8b87b13,
        0x1675078c5bfca9dc,
        0x697423ba265a4597,
        0xbe63c85687f90125,
        0xfa5ea96d320f0c09,
        0x0cc15241ae6ff30e,
        0x51613fef090403f4,
        0x32757bd13a1882b3,
        0xeebbf697f70ab363,
        0x1f8381d48711040a,
        0x09db730c94d2c725,
        0x52898f07af37a87c,
        0x023f9943f709f256,
        0x5c6bd9ddfd3f3c70,
        0x8387b25e652d089b,
        0x59aedb9796302114,
        0xf7d178f3bdb92811,
        0x711bf5628c5ff289,
        0xf357a9ff2d4b9ce7,
        0xf612ccb8e8db8e91,
        0x57157e1a6e24f86d,
        0xe6b8e3c82c63af1a,
        0xb640003d666c6
    ];

    /// a quadratic nonresidue u
    /// u^T
    const QUADRATIC_NONRESIDUE_TO_T: Fq3 = Fq3::new(
        MontFp!("20352210931107105404090576606765520961773321791847852436891069667257662057768295150100850989487484739038066884340836714918897718329742678773708674671463699964752171757206793"),
        Fq::ZERO,
        Fq::ZERO,
    );

    /// NQR = 11
    /// NQR ^ ((q^0 - 1) / 3) = 1
    /// NQR ^ ((q^1 - 1) / 3) = 26889017910506495156320058129568772566653933834480888250851357865717474104627594400663718411564763129536634204743601058030506421687843844062271438195269332524424351101183828
    /// NQR ^ ((q^2 - 1) / 3) = 7891820029458153760463333779529065900060223375068625705890303802702606422769111189933903712579664581217514620699178559520642750983460643979247817048104789503567679423896212
    const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[
        Fq::ONE,
        MontFp!("26889017910506495156320058129568772566653933834480888250851357865717474104627594400663718411564763129536634204743601058030506421687843844062271438195269332524424351101183828"),
        MontFp!("7891820029458153760463333779529065900060223375068625705890303802702606422769111189933903712579664581217514620699178559520642750983460643979247817048104789503567679423896212"),
    ];

    /// NQR = 11
    /// NQR ^ (2 * (q^0 - 1) / 3) = 1
    /// NQR ^ (2 * (q^1 - 1) / 3) = 7891820029458153760463333779529065900060223375068625705890303802702606422769111189933903712579664581217514620699178559520642750983460643979247817048104789503567679423896212
    /// NQR ^ (2 * (q^2 - 1) / 3) = 26889017910506495156320058129568772566653933834480888250851357865717474104627594400663718411564763129536634204743601058030506421687843844062271438195269332524424351101183828
    const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[
        Fq::ONE,
        Self::FROBENIUS_COEFF_FP3_C1[2],
        Self::FROBENIUS_COEFF_FP3_C1[1],
    ];
}