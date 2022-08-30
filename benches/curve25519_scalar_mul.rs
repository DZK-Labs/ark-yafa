#[macro_use]
extern crate core;

use ark_bls12_381::Bls12_381;
use ark_bn254::Bn254;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ec::twisted_edwards::TECurveConfig;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{BigInteger, Field, One, PrimeField};
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::{AllocatedFp, FpVar};
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::non_zero_affine::NonZeroAffineVar;
use ark_r1cs_std::groups::curves::twisted_edwards::{AffineVar, MontgomeryAffineVar};
use ark_r1cs_std::prelude::CondSelectGadget;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use ark_std::marker::PhantomData;
use ark_std::ops::Add;
use ark_std::time::Instant;
use ark_std::{test_rng, UniformRand};

use ark_yafa::ed25519;

type Fr = ed25519::Fr;
type Fq = ed25519::Fq;

const NUM_SCALAR_MUL: usize = 1;

#[derive(Default, Copy, Clone)]
pub struct NativeCircuit<CF: PrimeField> {
    constraint_field_phantom: PhantomData<CF>,
}

impl ConstraintSynthesizer<Fq> for NativeCircuit<Fq> {
    fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
        let mut rng = test_rng();
        //
        // important: this circuit has constants that come from a pseudorandom generator,
        //      usually, this means that every time the circuit would be different, and
        //      proving will fail because it is a different circuit from the one used in
        //      the indexing phase.
        //
        //      fortunately, we use the test_rng which is fixed, so the two circuits would
        //      be the same. This should not be used in production.
        //

        let mut points = Vec::new();
        let mut point_vars = Vec::new();

        let mut scalars = Vec::new();
        let mut scalar_bit_vars = Vec::new();

        let num_bits = ed25519::Fr::MODULUS_BIT_SIZE;

        let generator = ed25519::EdwardsAffine::prime_subgroup_generator();

        for _ in 0..NUM_SCALAR_MUL {
            let scalar = Fr::rand(&mut rng);
            scalars.push(scalar);

            let point = generator.mul_bigint(scalar.into_bigint()).into_affine();
            points.push(point);

            let (x, y) = point.xy().unwrap();
            let x_var = FpVar::<Fq>::new_constant(cs.clone(), x.clone())?;
            let y_var = FpVar::<Fq>::new_constant(cs.clone(), y.clone())?;
            let point_var = AffineVar::<ed25519::EdwardsParameters, FpVar<Fq>>::new(x_var, y_var);
            point_vars.push(point_var);

            let mut cur = scalar.into_bigint();

            let mut scalar_bit_var = vec![];
            for _ in 0..num_bits {
                let bit = cur.is_odd();
                cur.div2();
                scalar_bit_var.push(Boolean::new_witness(cs.clone(), || Ok(bit))?);
            }
            scalar_bit_vars.push(scalar_bit_var);
        }

        let (placeholder_x_var, placeholder_y_var, placeholder_neg_y_var) = {
            let placeholder = ed25519::EdwardsAffine::rand(&mut rng);
            let (x, y) = placeholder.xy().unwrap();
            let placeholder_x_var = FpVar::<Fq>::new_constant(cs.clone(), x.clone())?;
            let placeholder_y_var = FpVar::<Fq>::new_constant(cs.clone(), y.clone())?;

            let placeholder_neg_y_var = placeholder_y_var.negate()?;

            (placeholder_x_var, placeholder_y_var, placeholder_neg_y_var)
        };

        let base_vars = {
            let mut base_vars = Vec::new();
            let mut cur = generator.into_projective();
            for _ in 0..num_bits {
                let point = cur.into_affine();
                let (x, y) = point.xy().unwrap();
                let x_var = FpVar::<Fq>::new_constant(cs.clone(), x.clone())?;
                let y_var = FpVar::<Fq>::new_constant(cs.clone(), y.clone())?;

                let point_var =
                    AffineVar::<ed25519::EdwardsParameters, FpVar<Fq>>::new(x_var, y_var);

                base_vars.push(point_var);
                cur.double_in_place();
            }
            base_vars
        };

        let placeholder_var = AffineVar::<ed25519::EdwardsParameters, FpVar<Fq>>::new(
            placeholder_x_var.clone(),
            placeholder_y_var,
        );
        let placeholder_neg_var = AffineVar::<ed25519::EdwardsParameters, FpVar<Fq>>::new(
            placeholder_x_var,
            placeholder_neg_y_var,
        );
        for i in 0..NUM_SCALAR_MUL {
            let mut res: AffineVar<ed25519::EdwardsParameters, FpVar<Fq>> = placeholder_var.clone();

            for j in 0..num_bits {
                let new_res = res.clone().add(&base_vars[i]);
                res = AffineVar::<ed25519::EdwardsParameters, FpVar<Fq>>::conditionally_select(
                    &scalar_bit_vars[i][j as usize],
                    &new_res,
                    &res,
                )?;
            }

            res = res.add(&placeholder_neg_var);
            res.enforce_equal(&point_vars[i])?;

            drop(res);
        }

        Ok(())
    }
}

#[derive(Clone)]
pub struct NonNativeAffineVar<P: TECurveConfig, CF: PrimeField>
where
    P::BaseField: PrimeField,
{
    x: NonNativeFieldVar<P::BaseField, CF>,
    y: NonNativeFieldVar<P::BaseField, CF>,
}

impl<P: TECurveConfig, CF: PrimeField> NonNativeAffineVar<P, CF>
where
    P::BaseField: PrimeField,
{
    pub fn new(
        x: NonNativeFieldVar<P::BaseField, CF>,
        y: NonNativeFieldVar<P::BaseField, CF>,
    ) -> Self {
        Self { x, y }
    }
}

impl<P: TECurveConfig, CF: PrimeField> Add<&Self> for NonNativeAffineVar<P, CF>
where
    P::BaseField: PrimeField,
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        let a = NonNativeFieldVar::<P::BaseField, CF>::Constant(P::COEFF_A);
        let d = NonNativeFieldVar::<P::BaseField, CF>::Constant(P::COEFF_D);

        // Compute U = (x1 + y1) * (x2 + y2)
        let u1 = (&self.x * &a.negate().unwrap()) + &self.y;
        let u2 = &other.x + &other.y;

        let u = u1 * &u2;

        // Compute v0 = x1 * y2
        let v0 = &other.y * &self.x;

        // Compute v1 = x2 * y1
        let v1 = &other.x * &self.y;

        // Compute C = d*v0*v1
        let v2 = &v0 * &v1 * d;

        // Compute x3 = (v0 + v1) / (1 + v2)
        let t0 = &v0 + &v1;
        let t1 = &v2 + &NonNativeFieldVar::<P::BaseField, CF>::Constant(P::BaseField::one());
        let t1_inv = t1.inverse().unwrap();
        let x3 = &t0 * &t1_inv;

        // Compute y3 = (U + a * v0 - v1) / (1 - v2)
        let t0 = &u + (&a * &v0) - &v1;
        let t1 = &NonNativeFieldVar::<P::BaseField, CF>::Constant(P::BaseField::one()) - &v2;
        let t1_inv = t1.inverse().unwrap();
        let y3 = &t0 * &t1_inv;

        Self { x: x3, y: y3 }
    }
}

impl<P: TECurveConfig, CF: PrimeField> EqGadget<CF> for NonNativeAffineVar<P, CF>
where
    P::BaseField: PrimeField,
{
    fn is_eq(&self, other: &Self) -> Result<Boolean<CF>, SynthesisError> {
        let is_x_eq = self.x.is_eq(&other.x)?;
        let is_y_eq = self.y.is_eq(&other.y)?;
        is_x_eq.and(&is_y_eq)
    }
}

impl<P: TECurveConfig, CF: PrimeField> CondSelectGadget<CF> for NonNativeAffineVar<P, CF>
where
    P::BaseField: PrimeField,
    P: Clone,
{
    fn conditionally_select(
        cond: &Boolean<CF>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        let x = NonNativeFieldVar::<P::BaseField, CF>::conditionally_select(
            cond,
            &true_value.x,
            &false_value.x,
        )?;
        let y = NonNativeFieldVar::<P::BaseField, CF>::conditionally_select(
            cond,
            &true_value.x,
            &false_value.x,
        )?;

        Ok(Self { x, y })
    }
}

#[derive(Default, Copy, Clone)]
pub struct NonNativeCircuit<F: PrimeField, CF: PrimeField> {
    field_phantom: PhantomData<F>,
    constraint_field_phantom: PhantomData<CF>,
}

impl<CF: PrimeField> ConstraintSynthesizer<CF> for NonNativeCircuit<Fq, CF> {
    fn generate_constraints(self, cs: ConstraintSystemRef<CF>) -> Result<(), SynthesisError> {
        let mut rng = test_rng();
        //
        // important: this circuit has constants that come from a pseudorandom generator,
        //      usually, this means that every time the circuit would be different, and
        //      proving will fail because it is a different circuit from the one used in
        //      the indexing phase.
        //
        //      fortunately, we use the test_rng which is fixed, so the two circuits would
        //      be the same. This should not be used in production.
        //

        let mut points = Vec::new();
        let mut point_vars = Vec::new();

        let mut scalars = Vec::new();
        let mut scalar_bit_vars = Vec::new();

        let num_bits = ed25519::Fr::MODULUS_BIT_SIZE;

        let generator = ed25519::EdwardsAffine::prime_subgroup_generator();

        for _ in 0..NUM_SCALAR_MUL {
            let scalar = Fr::rand(&mut rng);
            scalars.push(scalar);

            let point = generator.mul_bigint(scalar.into_bigint()).into_affine();
            points.push(point);

            let (x, y) = point.xy().unwrap();
            let x_var = NonNativeFieldVar::<Fq, CF>::new_constant(cs.clone(), x.clone())?;
            let y_var = NonNativeFieldVar::<Fq, CF>::new_constant(cs.clone(), y.clone())?;
            let point_var = NonNativeAffineVar::<ed25519::EdwardsParameters, CF>::new(x_var, y_var);
            point_vars.push(point_var);

            let mut cur = scalar.into_bigint();

            let mut scalar_bit_var = vec![];
            for _ in 0..num_bits {
                let bit = cur.is_odd();
                cur.div2();
                scalar_bit_var.push(Boolean::new_witness(cs.clone(), || Ok(bit))?);
            }
            scalar_bit_vars.push(scalar_bit_var);
        }

        let (placeholder_x_var, placeholder_y_var, placeholder_neg_y_var) = {
            let placeholder = ed25519::EdwardsAffine::rand(&mut rng);
            let (x, y) = placeholder.xy().unwrap();
            let placeholder_x_var =
                NonNativeFieldVar::<Fq, CF>::new_constant(cs.clone(), x.clone())?;
            let placeholder_y_var =
                NonNativeFieldVar::<Fq, CF>::new_constant(cs.clone(), y.clone())?;
            let placeholder_neg_y_var = placeholder_y_var.negate()?;

            (placeholder_x_var, placeholder_y_var, placeholder_neg_y_var)
        };

        let base_vars = {
            let mut base_vars = Vec::new();
            let mut cur = generator.into_projective();
            for _ in 0..num_bits {
                let point = cur.into_affine();
                let (x, y) = point.xy().unwrap();
                let x_var = NonNativeFieldVar::<Fq, CF>::new_constant(cs.clone(), x.clone())?;
                let y_var = NonNativeFieldVar::<Fq, CF>::new_constant(cs.clone(), y.clone())?;

                let point_var =
                    NonNativeAffineVar::<ed25519::EdwardsParameters, CF>::new(x_var, y_var);

                base_vars.push(point_var);
                cur.double_in_place();
            }
            base_vars
        };

        let placeholder_var = NonNativeAffineVar::<ed25519::EdwardsParameters, CF>::new(
            placeholder_x_var.clone(),
            placeholder_y_var,
        );
        let placeholder_neg_var = NonNativeAffineVar::<ed25519::EdwardsParameters, CF>::new(
            placeholder_x_var,
            placeholder_neg_y_var,
        );

        for i in 0..NUM_SCALAR_MUL {
            let mut res: NonNativeAffineVar<ed25519::EdwardsParameters, CF> =
                placeholder_var.clone();

            for j in 0..num_bits {
                let new_res = res.clone().add(&base_vars[i]);
                res = NonNativeAffineVar::<ed25519::EdwardsParameters, CF>::conditionally_select(
                    &scalar_bit_vars[i][j as usize],
                    &new_res,
                    &res,
                )?;
            }

            res = res.add(&placeholder_neg_var);
            res.enforce_equal(&point_vars[i])?;

            drop(res);
        }

        Ok(())
    }
}

fn test_gemini_yafa_108() {
    let native_circuit = NativeCircuit::default();
    let r1cs = ark_gemini::circuit::generate_relation::<Fq, NativeCircuit<Fq>>(native_circuit);

    println!("r1cs size: {}", r1cs.a.len());

    let num_constraints = 200000;
    let num_variables = 200000;
    let num_non_zero = 200000;

    let mut rng = test_rng();

    let timer = Instant::now();
    let ck = ark_gemini::kzg::CommitterKey::<ark_yafa::yafa_108::Yafa>::new(
        num_non_zero + num_variables + num_constraints,
        5,
        &mut rng,
    );
    println!("SRS generation: {}", timer.elapsed().as_secs_f64());

    let timer = Instant::now();
    for _ in 0..1 {
        let _ = ark_gemini::psnark::Proof::new_time(&r1cs, &ck);
    }
    println!(
        "Proof generation for one proof: {}",
        timer.elapsed().as_secs_f64()
    );
}

fn test_gemini_yafa_146() {
    let native_circuit = NativeCircuit::default();
    let r1cs = ark_gemini::circuit::generate_relation::<Fq, NativeCircuit<Fq>>(native_circuit);

    println!("r1cs size: {}", r1cs.a.len());

    let num_constraints = 200000;
    let num_variables = 200000;
    let num_non_zero = 200000;

    let mut rng = test_rng();

    let timer = Instant::now();
    let ck = ark_gemini::kzg::CommitterKey::<ark_yafa::yafa_146::Yafa>::new(
        num_non_zero + num_variables + num_constraints,
        5,
        &mut rng,
    );
    println!("SRS generation: {}", timer.elapsed().as_secs_f64());

    let timer = Instant::now();
    for _ in 0..1 {
        let _ = ark_gemini::psnark::Proof::new_time(&r1cs, &ck);
    }
    println!(
        "Proof generation for one proof: {}",
        timer.elapsed().as_secs_f64()
    );
}

fn test_groth16_bls12_381() {
    let mut rng = test_rng();

    let non_native_circuit = NonNativeCircuit::<Fq, ark_bls12_381::Fr>::default();
    let params = ark_groth16::generate_random_parameters::<Bls12_381, _, _>(
        non_native_circuit.clone(),
        &mut rng,
    )
    .unwrap();

    println!("constraints: {}", params.a_query.len());

    let timer = Instant::now();
    for _ in 0..1 {
        let _ = ark_groth16::create_random_proof(non_native_circuit.clone(), &params, &mut rng)
            .unwrap();
    }

    println!(
        "Proof generation for one proof: {}",
        timer.elapsed().as_secs_f64()
    );
}

fn test_gemini_bls12_381() {
    let non_native_circuit = NonNativeCircuit::<Fq, ark_bls12_381::Fr>::default();
    let r1cs = ark_gemini::circuit::generate_relation::<
        ark_bls12_381::Fr,
        NonNativeCircuit<Fq, ark_bls12_381::Fr>,
    >(non_native_circuit);

    println!("r1cs size: {}", r1cs.a.len());

    let num_constraints = 2000000;
    let num_variables = 2000000;
    let num_non_zero = 2000000;

    let mut rng = test_rng();

    let timer = Instant::now();
    let ck = ark_gemini::kzg::CommitterKey::<Bls12_381>::new(
        num_non_zero + num_variables + num_constraints,
        5,
        &mut rng,
    );
    println!("SRS generation: {}", timer.elapsed().as_secs_f64());

    let timer = Instant::now();
    for _ in 0..1 {
        let _ = ark_gemini::psnark::Proof::new_time(&r1cs, &ck);
    }
    println!(
        "Proof generation for one proof: {}",
        timer.elapsed().as_secs_f64()
    );
}

fn test_groth16_bn254() {
    let mut rng = test_rng();

    let non_native_circuit = NonNativeCircuit::<Fq, ark_bn254::Fr>::default();
    let params = ark_groth16::generate_random_parameters::<Bn254, _, _>(
        non_native_circuit.clone(),
        &mut rng,
    )
    .unwrap();

    println!("constraints: {}", params.a_query.len());

    let timer = Instant::now();
    for _ in 0..1 {
        let _ = ark_groth16::create_random_proof(non_native_circuit.clone(), &params, &mut rng)
            .unwrap();
    }

    println!(
        "Proof generation for one proof: {}",
        timer.elapsed().as_secs_f64()
    );
}

fn test_gemini_bn254() {
    let non_native_circuit = NonNativeCircuit::<Fq, ark_bn254::Fr>::default();
    let r1cs = ark_gemini::circuit::generate_relation::<
        ark_bn254::Fr,
        NonNativeCircuit<Fq, ark_bn254::Fr>,
    >(non_native_circuit);

    println!("r1cs size: {}", r1cs.a.len());

    let num_constraints = 2000000;
    let num_variables = 2000000;
    let num_non_zero = 2000000;

    let mut rng = test_rng();

    let timer = Instant::now();
    let ck = ark_gemini::kzg::CommitterKey::<Bn254>::new(
        num_non_zero + num_variables + num_constraints,
        5,
        &mut rng,
    );
    println!("SRS generation: {}", timer.elapsed().as_secs_f64());

    let timer = Instant::now();
    for _ in 0..1 {
        let _ = ark_gemini::psnark::Proof::new_time(&r1cs, &ck);
    }
    println!(
        "Proof generation for one proof: {}",
        timer.elapsed().as_secs_f64()
    );
}

fn main() {
    println!("Gemini on Yafa 108:");
    test_gemini_yafa_108();
    println!("Gemini on Yafa 146:");
    test_gemini_yafa_146();
    println!("Gemini on BN254:");
    test_gemini_bn254();
    println!("Gemini on BLS12-381:");
    test_gemini_bls12_381();
    println!("Groth16 on BN254:");
    test_groth16_bn254();
    println!("Groth16 on BLS12-381:");
    test_groth16_bls12_381();
}
