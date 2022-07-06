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
#[derive(Clone, Copy)]
pub struct ParameterPair<E: Engine> {
    pub g1_result: Option<E::G1Affine>,
    pub g2_result: Option<E::G2Affine>,
    pub g1_mine: Option<E::G1Affine>,
    pub g2_mine: Option<E::G2Affine>,
}

impl<E: Engine> ParameterPair<E> {
    fn get_g1(&self) -> E::G1Affine {
        self.g1_result.unwrap()
    }
}

pub fn initParameterList<E>() -> Vec<ParameterPair<E>>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let mut newp = ParameterPair {
        g1_result: Some(E::G1Affine::generator()),
        g2_result: Some(E::G2Affine::generator()),
        g1_mine: None,
        g2_mine: None,
    };
    return [newp].to_vec();
}

pub fn ParamterListExcute<E: Engine>(
    mut vec: Vec<ParameterPair<E>>,
    p: ParameterPair<E>,
) -> Vec<ParameterPair<E>>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let mut newp = ParameterPair {
        g1_result: None,
        g2_result: None,
        g1_mine: Some(E::G1Affine::generator()),
        g2_mine: Some(E::G2Affine::generator()),
    };
    let x = p.clone();
    let len = vec.len();
    if (len != 0) {
        let b = verify_mpc_g1(&p, &vec);
        assert_eq!(b, true);
        newp = x;
    }
    vec.push(newp);
    return vec;
}
pub fn mpc_common_paramters_custom_generator<E>(
    paramter_last: &ParameterPair<E>,
    my_alpha: E::Fr,
) -> Result<ParameterPair<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    mpc_common_paramters_custom::<E>(
        E::G1Affine::generator(),
        E::G2Affine::generator(),
        paramter_last,
        my_alpha,
    )
}

pub fn mpc_common_paramters_custom<E>(
    g1: E::G1Affine,
    g2: E::G2Affine,
    paramter_last: &ParameterPair<E>,
    my_alpha: E::Fr,
) -> Result<ParameterPair<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    //let g2 = E::G2::generator();
    //let g1 = E::G1::generator();
    let g1_before = paramter_last.g1_result.unwrap();
    let g1_after = (g1_before * my_alpha).to_affine();
    let g2_before = paramter_last.g2_result.unwrap();
    let g2_after = (g2_before * my_alpha).to_affine();
    let g1_mine = (g1 * my_alpha).to_affine();
    let g2_mine = (g2 * my_alpha).to_affine();
    let result = ParameterPair {
        g1_result: Some(g1_after),
        g2_result: Some(g2_after),
        g1_mine: Some(g1_mine),
        g2_mine: Some(g2_mine),
    };
    return Ok(result);
}

pub fn mpc_bad_paramters_custom<E>(
    paramter_last: &ParameterPair<E>,
    my_alpha: E::Fr,
) -> Result<ParameterPair<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let g2 = E::G2::generator();
    let g1 = E::G1::generator();
    let g1_before = paramter_last.g1_result.unwrap();
    let g1_after = (g1 * my_alpha).to_affine();
    let g2_before = paramter_last.g2_result.unwrap();
    let g2_after = (g2 * my_alpha).to_affine();
    let g1_mine = (g1).to_affine();
    let g2_mine = (g2 * my_alpha).to_affine();
    let result = ParameterPair {
        g1_result: Some(g1_after),
        g2_result: Some(g2_after),
        g1_mine: Some(g1_mine),
        g2_mine: Some(g2_mine),
    };
    return Ok(result);
}

pub fn verify_mpc_g1<E>(new_paramter: &ParameterPair<E>, paramters: &Vec<ParameterPair<E>>) -> bool
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let g1 = E::G1::generator();
    let g2 = E::G2::generator();
    let mut result = E::pairing(&new_paramter.g1_mine.unwrap(), &g2.to_affine())
        == E::pairing(&g1.to_affine(), &new_paramter.g2_mine.unwrap());
    let index = paramters.len();

    if (index >= 1) {
        let paramter_last = new_paramter.g1_result.unwrap();
        let paramter_part2 = new_paramter.g2_mine.unwrap();
        let paramter_part1 = paramters[index - 1].g1_result.unwrap();
        /*let paramter_last = paramters[index - 1].g1_result.unwrap();
        let paramter_part2 = paramters[index - 1].g2_mine.unwrap();
        let paramter_part1 = paramters[index - 2].g1_result.unwrap();*/
        result = result
            && E::pairing(&paramter_last, &g2.to_affine())
                == E::pairing(&paramter_part1, &paramter_part2);
    }
    result
}

#[derive(Clone)]
pub struct TauParameterPair<E: Engine> {
    pub list: Vec<ParameterPair<E>>,
}

impl<E: Engine> TauParameterPair<E> {
    pub fn get_g1(&self) -> Vec<E::G1Affine> {
        let mut result = Vec::new();
        for i in 0..self.list.len() {
            result.push(self.list[i].get_g1())
        }
        result
    }
}

pub fn initTauParameterList<E>(n: usize) -> Vec<TauParameterPair<E>>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let mut result = Vec::new();
    for i in 0..n {
        result.push(ParameterPair {
            g1_result: Some(E::G1Affine::generator()),
            g2_result: Some(E::G2Affine::generator()),
            g1_mine: None,
            g2_mine: None,
        });
    }

    return [TauParameterPair { list: result }].to_vec();
}
pub fn mpc_common_tauparamters_custom_generator<E>(
    tauparamter_last: &TauParameterPair<E>,
    my_x: Vec<E::Fr>,
) -> Result<TauParameterPair<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    mpc_common_tauparamters_custom(
        E::G1Affine::generator(),
        E::G2Affine::generator(),
        tauparamter_last,
        my_x,
    )
}

pub fn mpc_common_tauparamters_custom<E>(
    g1: E::G1Affine,
    g2: E::G2Affine,
    tauparamter_last: &TauParameterPair<E>,
    my_x: Vec<E::Fr>,
) -> Result<TauParameterPair<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    //let g2 = E::G2::generator();
    //let g1 = E::G1::generator();
    let mut result_list = Vec::new();
    let len = my_x.len();
    assert_eq!(len, tauparamter_last.list.len());
    for i in 0..len {
        let xg1_before = tauparamter_last.list[i].g1_result.unwrap();
        let xg1_after = (xg1_before * my_x[i]).to_affine();
        let xg2_before = tauparamter_last.list[i].g2_result.unwrap();
        let xg2_after = (xg2_before * my_x[i]).to_affine();
        let xg1_mine = (g1 * my_x[i]).to_affine();
        let xg2_mine = (g2 * my_x[i]).to_affine();

        let result = ParameterPair {
            g1_result: Some(xg1_after),
            g2_result: Some(xg2_after),
            g1_mine: Some(xg1_mine),
            g2_mine: Some(xg2_mine),
        };
        //let result = mpc_common_paramters_custom(&tauparamter_last.list[i], my_x[i]).unwrap();
        result_list.push(result);
    }
    return Ok(TauParameterPair { list: result_list });
}

pub fn TauParamterListExcute<E: Engine>(
    mut vec: Vec<TauParameterPair<E>>,
    p: TauParameterPair<E>,
) -> Vec<TauParameterPair<E>>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let b = verify_mpc_x(&p, &vec);
    assert_eq!(b, true);
    vec.push(p);
    return vec;
}

pub fn verify_x_pow<E>(new_xparamter: &TauParameterPair<E>) -> bool
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let mut result = true;
    let t_len = new_xparamter.list.len();
    let g2 = E::G2::generator().to_affine();
    if t_len > 1 {
        for i in 1..t_len {
            result = result
                && E::pairing(
                    &new_xparamter.list[i - 1].g1_result.unwrap(),
                    &new_xparamter.list[0].g2_result.unwrap(),
                ) == E::pairing(&new_xparamter.list[i].g1_result.unwrap(), &g2);
        }
    }
    return result;
}

pub fn verify_mpc_x<E>(
    new_xparamter: &TauParameterPair<E>,
    paramters: &Vec<TauParameterPair<E>>,
) -> bool
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    //let g2 = E::G2::generator().to_affine();
    let mut result = verify_x_pow(&new_xparamter);
    let index = paramters.len();
    let mut list = Vec::new();
    for i in 0..index {
        list.push(paramters[i].list[0].clone());
    }
    result = result && verify_mpc_g1(&new_xparamter.list[0], &list);
    result
}

use crate::groth16::generator::KeypairAssembly;
use crate::groth16::generator_val::*;
//公共计算产生的参数，要传给进行非通用计算的用户
pub struct CommonParamter<E: Engine> {
    pub alpha: ParameterPair<E>,
    pub beta: ParameterPair<E>,
    pub tau: TauParameterPair<E>,
    pub alpha_mul_tau: TauParameterPair<E>,
    pub beta_mul_tau: TauParameterPair<E>,
}

impl<E: Engine> CommonParamter<E> {
    pub fn matrix(
        &self,
        at_aux: &Vec<Vec<(E::Fr, usize)>>,
        bt_aux: &Vec<Vec<(E::Fr, usize)>>,
        ct_aux: &Vec<Vec<(E::Fr, usize)>>,
    ) -> CommonParamterMatrix<E> {
        //let tau_matrix = matrix_mul_tau::<E>(&ct_aux, &self.tau.get_g1());
        unimplemented!()
    }
}

pub fn mpc_common_paramters_custom_all<E>() -> CommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    //下边每一次生成的时候单独验证，但是有耦合的地方
    //eg:alpha_mul_tau 与 alpha
    //时候也应该考虑验证吧
    /*let alpha = generator_val();
    let beta = generator_val();
    let tau = generator_tau();
    let alpha_mul_tau = generator_tau_from(&alpha);
    let beta_mul_tau = generator_tau_from(&beta);
    CommonParamter {
        alpha,
        beta,
        tau,
        alpha_mul_tau,
        beta_mul_tau,
    }*/
    unimplemented!()
}

//非通用计算参数
pub struct UnCommonParamter<E: Engine> {
    pub delta: ParameterPair<E>,
    pub gamma: ParameterPair<E>,
    pub ic: TauParameterPair<E>,
    pub l: TauParameterPair<E>,
    pub h: TauParameterPair<E>,
}

pub struct CommonParamterMatrix<E: Engine> {
    pub tau_matrix: TauParameterPair<E>,
    pub alpha_mul_tau_matrix: TauParameterPair<E>,
    pub beta_mul_tau_matrix: TauParameterPair<E>,
}

//要传给进行非通用计算的用户使用公共计算产生的参数和电路产生的矩阵运算非通用计算参数
pub fn mpc_uncommon_paramters_custom_all<E: Engine>(
    common_paramter_matrix: &CommonParamterMatrix<E>,
) -> UnCommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    //(belat*ux + alpha*vx + wx)/gamma 为例，
    //u是我们拿的公开数据，beta*x就是传的beta_mul_tau这个数组
    //计算出 u * beta_mul_tau
    //再一次写出 v * alpha_mul_tau 和 w * tau
    //作为一个整体的起始点G，去和每一个人的1/gamma做运算
    //即 generator_val_from(G)。
    //而这要求在每一次上传时
    //另t(x)似乎也需要一些数组
    /*let gamma = generator_val();
    let delta = generator_val();
    let beta_mul_tau_div_gamma = generator_tau()
    let beta_mul_tau_u = matrix_mul_tau(&at_aux, &common_paramter.beta_mul_tau);
    let alpha_mul_tau_v = matrix_mul_tau(&bt_aux, &common_paramter.alpha_mul_tau);
    let tau_w = matrix_mul_tau(&ct_aux, &common_paramter.tau);
    let sum = tau_list_sum(beta_mul_tau_u, alpha_mul_tau_v, tau_w);*/

    unimplemented!()
}
