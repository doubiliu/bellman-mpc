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
pub fn generator_val<E: Engine>() -> ParameterPair<E> {
    //多人运算验证，产生一个随机点，初始值是g1/g2，
    //[1,2,3,4]这几个人运算之后得到4!*g1或4!*g2
    unimplemented!();
}

pub fn generator_val_from<E: Engine>(p: &ParameterPair<E>) -> ParameterPair<E> {
    //多人运算验证，产生一个随机点，初始值是p，
    //[1,2,3,4]这几个人运算之后得到4!*p
    unimplemented!();
}

pub fn generator_tau<E: Engine>() -> TauParameterPair<E> {
    //多人运算验证，产生一系列指数增长的随机点，初始值是[g1,g1,g1,g1]，
    //[1,2,3,4]这几个人运算之后得到[[g1,g1,g1,g1],[g1,2*g1,3*g1,4*g1],[g1,4*g1,9*g1,16*g1].....]
    unimplemented!();
}

pub fn generator_tau_from<E: Engine>(tau: &ParameterPair<E>) -> TauParameterPair<E> {
    //多人运算验证，产生一系列指数增长的随机点，初始值是[tau0,tau1,tau2,tau3]，
    //[1,2,3,4]这几个人运算之后得到[[tau0,tau1,tau2,tau3],[tau0,2*tau1,3*tau2,4*tau3],[tau0,4*tau1,9*tau2,16*tau3].....]
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
