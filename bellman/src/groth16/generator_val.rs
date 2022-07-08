use rand_core::RngCore;
use std::ops::{AddAssign, MulAssign};
use std::sync::Arc;

use ff::{Field, PrimeField};
use group::{prime::PrimeCurveAffine, Curve, Group, Wnaf, WnafGroup};
use pairing::Engine;

use super::{Parameters, VerifyingKey};

use crate::{Circuit, ConstraintSystem, Index, LinearCombination, SynthesisError, Variable};

use crate::domain::{EvaluationDomain, Scalar};

use crate::multicore::Worker;

use bls12_381::Bls12;

use crate::groth16::mpc::{ParameterPair, TauParameterPair};

pub fn matrix_eval<E>(
    matrix: &Vec<Vec<(E::Fr, usize)>>,
    tau_g1: &Vec<E::G1Affine>,
    tau_g2: &Vec<E::G2Affine>,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>)
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    unimplemented!();
}

pub fn matrix_mul_tau<E>(matrix: &Vec<Vec<(E::Fr, usize)>>, tau: &Vec<E::G1Affine>) -> E::G1Affine
where
    E: Engine,
    E::G1: WnafGroup,
{
    assert_eq!(matrix[0].len(), tau.len());
    let mut result: Vec<E::G1> = Vec::new();
    for i in 0..matrix.len() {
        let mut result_part: Vec<E::G1> = Vec::new();
        for j in 0..matrix[i].len() {
            assert_eq!(matrix[i][j].1, j);
            result_part.push((tau[j] * matrix[i][j].0));
        }
        if result.len() == 0 {
            result = result_part;
        } else {
            result = tau_list_add::<E>(&result, &result_part);
        }
    }
    return tau_list_sum::<E>(&result);
}

pub fn tau_list_add<E>(v1: &Vec<E::G1>, v2: &Vec<E::G1>) -> Vec<E::G1>
where
    E: Engine,
    E::G1: WnafGroup,
{
    assert_eq!(v1.len(), v2.len());
    let mut result: Vec<E::G1> = Vec::new();
    for i in 0..v1.len() {
        result.push(v1[i] + v2[2])
    }
    result
}

pub fn tau_list_sum<E>(v1: &Vec<E::G1>) -> E::G1Affine
where
    E: Engine,
    E::G1: WnafGroup,
{
    let mut result = v1[0];
    if v1.len() > 1 {
        for i in 1..v1.len() {
            result = result + v1[i];
        }
    }
    return result.to_affine();
}
