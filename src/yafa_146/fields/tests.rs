use ark_algebra_test_templates::{
    fields::*, generate_field_serialization_test, generate_field_test,
};
use ark_ff::{Field, One, PrimeField, UniformRand, Zero};
use ark_serialize::{buffer_bit_byte_size, CanonicalSerialize};
use ark_std::{rand::Rng, test_rng};
use core::ops::{AddAssign, MulAssign, SubAssign};

use crate::yafa_146::*;

generate_field_test!(yafa_146; fq2; fq6; fq12; mont(9, 4); );
generate_field_serialization_test!(yafa_146; fq2; fq6; fq12;);
