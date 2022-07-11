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
pub struct ParameterPair<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub g1_result: Option<E::G1>,
    pub g2_result: Option<E::G2>,
    pub g1_mine: Option<E::G1>,
    pub g2_mine: Option<E::G2>,
}

impl<E> ParameterPair<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    fn get_g1(&self) -> E::G1 {
        self.g1_result.unwrap()
    }
    fn get_g2(&self) -> E::G2 {
        self.g2_result.unwrap()
    }
}

pub fn initParameterList<E>() -> Vec<ParameterPair<E>>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let mut newp = ParameterPair {
        g1_result: Some(E::G1::generator()),
        g2_result: Some(E::G2::generator()),
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
        g1_mine: Some(E::G1::generator()),
        g2_mine: Some(E::G2::generator()),
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
        E::G1::generator(),
        E::G2::generator(),
        paramter_last,
        my_alpha,
    )
}

pub fn mpc_common_paramters_custom<E>(
    g1: E::G1,
    g2: E::G2,
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
    let g1_after = (g1_before * my_alpha);
    let g2_before = paramter_last.g2_result.unwrap();
    let g2_after = (g2_before * my_alpha);
    let g1_mine = (g1 * my_alpha);
    let g2_mine = (g2 * my_alpha);
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
    let g1_after = (g1 * my_alpha);
    let g2_before = paramter_last.g2_result.unwrap();
    let g2_after = (g2 * my_alpha);
    let g1_mine = (g1);
    let g2_mine = (g2 * my_alpha);
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
    let mut result = E::pairing(&new_paramter.g1_mine.unwrap().to_affine(), &g2.to_affine())
        == E::pairing(&g1.to_affine(), &new_paramter.g2_mine.unwrap().to_affine());
    let index = paramters.len();

    if (index >= 1) {
        let paramter_last = new_paramter.g1_result.unwrap();
        let paramter_part2 = new_paramter.g2_mine.unwrap();
        let paramter_part1 = paramters[index - 1].g1_result.unwrap();
        /*let paramter_last = paramters[index - 1].g1_result.unwrap();
        let paramter_part2 = paramters[index - 1].g2_mine.unwrap();
        let paramter_part1 = paramters[index - 2].g1_result.unwrap();*/
        result = result
            && E::pairing(&paramter_last.to_affine(), &g2.to_affine())
                == E::pairing(&paramter_part1.to_affine(), &paramter_part2.to_affine());
    }
    result
}

#[derive(Clone)]
pub struct TauParameterPair<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub list: Vec<ParameterPair<E>>,
}

impl<E> TauParameterPair<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub fn get_g1(&self) -> Vec<E::G1> {
        let mut result = Vec::new();
        for i in 0..self.list.len() {
            result.push(self.list[i].get_g1())
        }
        result
    }
    pub fn get_g1_affine(&self) -> Vec<E::G1Affine> {
        let mut result = Vec::new();
        for i in 0..self.list.len() {
            result.push(self.list[i].get_g1().to_affine())
        }
        result
    }

    pub fn get_g2(&self) -> Vec<E::G2> {
        let mut result = Vec::new();
        for i in 0..self.list.len() {
            result.push(self.list[i].get_g2())
        }
        result
    }

    pub fn get_g2_affine(&self) -> Vec<E::G2Affine> {
        let mut result = Vec::new();
        for i in 0..self.list.len() {
            result.push(self.list[i].get_g2().to_affine())
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
            g1_result: Some(E::G1::generator()),
            g2_result: Some(E::G2::generator()),
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
        E::G1::generator(),
        E::G2::generator(),
        tauparamter_last,
        my_x,
    )
}

pub fn mpc_common_tauparamters_custom<E>(
    g1: E::G1,
    g2: E::G2,
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
        let xg1_after = (xg1_before * my_x[i]);
        let xg2_before = tauparamter_last.list[i].g2_result.unwrap();
        let xg2_after = (xg2_before * my_x[i]);
        let xg1_mine = (g1 * my_x[i]);
        let xg2_mine = (g2 * my_x[i]);

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

pub fn TauParamterListExcute<E>(
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
                    &new_xparamter.list[i - 1].g1_result.unwrap().to_affine(),
                    &new_xparamter.list[0].g2_result.unwrap().to_affine(),
                ) == E::pairing(&new_xparamter.list[i].g1_result.unwrap().to_affine(), &g2);
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
//公共计算产生的参数，要传给进行非通用计算的用户
pub struct CommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub alpha: ParameterPair<E>,
    pub beta: ParameterPair<E>,
    pub tau: TauParameterPair<E>,
    pub alpha_mul_tau: TauParameterPair<E>,
    pub beta_mul_tau: TauParameterPair<E>,
}

impl<E> CommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub fn to_storage_format(&self) -> CommonParamterInStorage<E> {
        CommonParamterInStorage {
            alpha_g1: self.alpha.get_g1(),
            alpha_g2: self.alpha.get_g2(),
            beta_g1: self.beta.get_g1(),
            beta_g2: self.beta.get_g2(),
            tau_g1: self.tau.get_g1(),
            tau_g2: self.tau.get_g2(),
            alpha_mul_tau_g1: self.alpha_mul_tau.get_g1(),
            alpha_mul_tau_g2: self.alpha_mul_tau.get_g2(),
            beta_mul_tau_g1: self.beta_mul_tau.get_g1(),
            beta_mul_tau_g2: self.beta_mul_tau.get_g2(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct CommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub alpha_g1: E::G1,
    pub alpha_g2: E::G2,
    pub beta_g1: E::G1,
    pub beta_g2: E::G2,
    pub tau_g1: Vec<E::G1>,
    pub tau_g2: Vec<E::G2>,
    pub alpha_mul_tau_g1: Vec<E::G1>,
    pub alpha_mul_tau_g2: Vec<E::G2>,
    pub beta_mul_tau_g1: Vec<E::G1>,
    pub beta_mul_tau_g2: Vec<E::G2>,
}

fn list_mul_matrix<E>(
    list_g1: &Vec<E::G1>,
    list_g2: &Vec<E::G2>,
    matrix: &Vec<Vec<(E::Fr, usize)>>,
) -> (Vec<E::G1>, Vec<E::G2>)
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let len = matrix[0].len();
    assert!(len == list_g2.len());
    assert!(len == list_g1.len());
    let mut result_g1 = Vec::new();
    let mut result_g2 = Vec::new();
    for i in 0..matrix.len() {
        for j in 0..len {
            if i == 0 {
                result_g1.push(list_g1[j] * matrix[i][j].0);
                result_g2.push(list_g2[j] * matrix[i][j].0);
            } else {
                result_g1[j] = result_g1[j] + list_g1[j] * matrix[i][j].0;
                result_g2[j] = result_g2[j] + list_g2[j] * matrix[i][j].0;
            }
        }
    }
    (result_g1, result_g2)
}

impl<E> CommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub fn matrix(
        &self,
        at_aux: &Vec<Vec<(E::Fr, usize)>>,
        bt_aux: &Vec<Vec<(E::Fr, usize)>>,
        ct_aux: &Vec<Vec<(E::Fr, usize)>>,
        num_inputs: usize,
        num_aux: usize,
    ) -> CommonParamterMatrix<E>
    where
        E: Engine,
        E::G1: WnafGroup,
        E::G2: WnafGroup,
    {
        let (alpha_matrix_g1, alpha_matrix_g2) =
            list_mul_matrix::<E>(&self.alpha_mul_tau_g1, &self.alpha_mul_tau_g2, bt_aux);
        let (beta_matrix_g1, beta_matrix_g2) =
            list_mul_matrix::<E>(&self.beta_mul_tau_g1, &self.beta_mul_tau_g2, at_aux);
        let (tau_matrix_g1, tau_matrix_g2) =
            list_mul_matrix::<E>(&self.tau_g1, &self.tau_g2, ct_aux);
        let mut matrixed_g1_front = E::G1::identity();
        let mut matrixed_g2_front = E::G2::identity();
        let mut matrixed_g1_back = E::G1::identity();
        let mut matrixed_g2_back = E::G2::identity();

        for i in 0..num_aux + num_inputs {
            if i < num_aux {
                matrixed_g1_front =
                    matrixed_g1_front + alpha_matrix_g1[i] + beta_matrix_g1[i] + tau_matrix_g1[i];
                matrixed_g2_front =
                    matrixed_g2_front + alpha_matrix_g2[i] + beta_matrix_g2[i] + tau_matrix_g2[i];
            } else {
                matrixed_g1_back =
                    matrixed_g1_front + alpha_matrix_g1[i] + beta_matrix_g1[i] + tau_matrix_g1[i];
                matrixed_g2_back =
                    matrixed_g2_front + alpha_matrix_g2[i] + beta_matrix_g2[i] + tau_matrix_g2[i];
            }
        }

        CommonParamterMatrix {
            matrixed_g1_front: matrixed_g1_front - E::G1::identity(),
            matrixed_g2_front: matrixed_g2_front - E::G2::identity(),
            matrixed_g1_back: matrixed_g1_back - E::G1::identity(),
            matrixed_g2_back: matrixed_g2_back - E::G2::identity(),
        }
    }
}
/*
use std::collections::HashMap;
pub struct CommonParamtersMap<E: Engine> {
    len: i32,
    map: HashMap<i32, CommonParamterInStorage<E>>,
}

impl<E: Engine> CommonParamtersMap<E> {
    pub fn get_last(&self) -> CommonParamterInStorage<E> {
        let last_index = self.len;
        self.map[&last_index].clone()
    }
}*/

pub fn make_new_paramter<E>(
    x: &u64,
    pointg1: &E::G1,
    pointg2: &E::G2,
    baseg1: &E::G1,
    baseg2: &E::G2,
) -> ParameterPair<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    ParameterPair::<E> {
        g1_result: Some(*pointg1 * E::Fr::from(*x)),
        g2_result: Some(*pointg2 * E::Fr::from(*x)),
        g1_mine: Some(*baseg1 * E::Fr::from(*x)),
        g2_mine: Some(*baseg2 * E::Fr::from(*x)),
    }
}

pub fn make_new_tau_paramter<E>(
    x: &u64,
    pointg1_list: &Vec<E::G1>,
    pointg2_list: &Vec<E::G2>,
) -> TauParameterPair<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let base_g1 = E::G1::generator();
    let base_g2 = E::G2::generator();
    let mut list: Vec<ParameterPair<E>> = Vec::new();
    assert_eq!(pointg1_list.len(), pointg1_list.len());
    for i in 0..pointg1_list.len() {
        list.push(make_new_paramter(
            x,
            &pointg1_list[i],
            &pointg2_list[i],
            &base_g1,
            &base_g2,
        ));
    }
    TauParameterPair { list }
}

pub fn initial_common_paramters<E>() -> CommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let len = 0;
    let init_paramters = 0;
    unimplemented!()
}

pub fn mpc_common_paramters_generator<E>(
    storage: &CommonParamterInStorage<E>,
    (alpha, beta, tau): (u64, u64, u64),
) -> CommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let g1 = E::G1::generator();
    let g2 = E::G2::generator();
    //从common里拿出alpha数据，计算本次alpha
    let new_alpha = make_new_paramter::<E>(&alpha, &storage.alpha_g1, &storage.alpha_g2, &g1, &g2);

    //从common里拿出alpha数据，计算本次beta
    let new_beta = make_new_paramter::<E>(&beta, &storage.beta_g1, &storage.beta_g2, &g1, &g2);

    //从common里拿出x数据，计算本次x[]
    let new_tau = make_new_tau_paramter::<E>(&tau, &storage.tau_g1, &storage.tau_g2);

    //计算alpha*x[]
    let new_alpha_mul_tau = make_new_tau_paramter::<E>(
        &(alpha * tau),
        &storage.alpha_mul_tau_g1,
        &storage.alpha_mul_tau_g2,
    );
    //计算beta*x[]

    let new_beta_mul_tau = make_new_tau_paramter::<E>(
        &(beta * tau),
        &storage.beta_mul_tau_g1,
        &storage.beta_mul_tau_g2,
    );
    //return
    CommonParamter {
        alpha: new_alpha,
        beta: new_beta,
        tau: new_tau,
        alpha_mul_tau: new_alpha_mul_tau,
        beta_mul_tau: new_beta_mul_tau,
    }
}

pub fn verify_new_paramter<E>(paramter: &ParameterPair<E>, baseg1: &E::G1, baseg2: &E::G2) -> bool
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    E::pairing(
        &paramter.g1_result.unwrap().to_affine(),
        &E::G2Affine::generator(),
    ) == E::pairing(&baseg1.to_affine(), &paramter.g2_mine.unwrap().to_affine())
}

pub fn verify_common_paramter<E>(
    storage: &CommonParamterInStorage<E>,
    new_paramter: &CommonParamter<E>,
) -> CommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let len = new_paramter.tau.list.len();
    let mut result = (len == new_paramter.alpha_mul_tau.list.len())
        && (new_paramter.beta_mul_tau.list.len() == len);
    //从stoarge和新参数验证alpha
    let result_alpha =
        verify_new_paramter(&new_paramter.alpha, &storage.alpha_g1, &storage.alpha_g2);

    //从stoarge和新参数验证beta
    let result_beta = verify_new_paramter(&new_paramter.beta, &storage.beta_g1, &storage.beta_g2);

    let mut result_tau = true;
    let mut result_alpha_tau = true;
    let mut result_beta_tau = true;

    for i in 0..len {
        //x[]自身验证指数数组
        if i > 1 {
            result_tau = result_tau
                && (E::pairing(
                    &new_paramter.tau.list[i].g1_result.unwrap().to_affine(),
                    &E::G2Affine::generator(),
                ) == E::pairing(
                    &new_paramter.tau.list[i - 1].g1_result.unwrap().to_affine(),
                    &new_paramter.tau.list[0].g2_result.unwrap().to_affine(),
                ));
        }
        //从storage和新参数里验证 alpha*x
        result_alpha_tau = result_alpha_tau
            && verify_new_paramter(
                &new_paramter.alpha_mul_tau.list[i],
                &storage.alpha_mul_tau_g1[i],
                &storage.alpha_mul_tau_g2[i],
            );
        //从storage和新参数里验证 beta*x
        result_beta_tau = result_beta_tau
            && verify_new_paramter(
                &new_paramter.beta_mul_tau.list[i],
                &storage.beta_mul_tau_g1[i],
                &storage.beta_mul_tau_g2[i],
            );
    }

    assert_eq!(
        true,
        result && result_alpha && result_beta && result_tau && result_alpha_tau && result_beta_tau
    );
    new_paramter.to_storage_format()
}

pub fn mpc_common_paramters_custom_all<E>() -> CommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    //on-chain
    let mut paramter_in_storage = initial_common_paramters::<E>();
    //under-chain
    let player1_common = mpc_common_paramters_generator(&paramter_in_storage, (1, 2, 3));
    //on-chain
    paramter_in_storage = verify_common_paramter(&paramter_in_storage, &player1_common);
    //under-chain
    let player2_common = mpc_common_paramters_generator(&paramter_in_storage, (2, 3, 4));
    //on-chain
    paramter_in_storage = verify_common_paramter(&paramter_in_storage, &player2_common);
    //under-chain
    let player3_common = mpc_common_paramters_generator(&paramter_in_storage, (3, 4, 5));
    //on-chain
    paramter_in_storage = verify_common_paramter(&paramter_in_storage, &player3_common);

    //on-chain
    paramter_in_storage
}

//非通用计算参数
pub struct UnCommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub delta: ParameterPair<E>,
    pub gamma: ParameterPair<E>,
    pub ic: TauParameterPair<E>,
    pub l: TauParameterPair<E>,
    pub h: TauParameterPair<E>,
}
#[derive(Debug, Clone)]
pub struct UnCommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub alpha_g1: E::G1Affine,
    /*pub alpha_g2: E::G2Affine,
    pub beta_g1: E::G1Affine,
    pub beta_g2: E::G2Affine,
    pub tau_g1: Vec<E::G1Affine>,
    pub tau_g2: Vec<E::G2Affine>,
    pub alpha_mul_tau_g1: Vec<E::G1Affine>,
    pub alpha_mul_tau_g2: Vec<E::G2Affine>,
    pub beta_mul_tau_g1: Vec<E::G1Affine>,
    pub beta_mul_tau_g2: Vec<E::G2Affine>,*/
}
/*
pub struct UnCommonParamtersMap<E: Engine> {
    len: i32,
    map: HashMap<i32, UnCommonParamterInStorage<E>>,
}

impl<E: Engine> UnCommonParamtersMap<E> {
    pub fn get_last(&self) -> UnCommonParamterInStorage<E> {
        let last_index = self.len;
        self.map[&last_index].clone()
    }
}*/

pub struct CommonParamterMatrix<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    pub matrixed_g1_front: E::G1,
    pub matrixed_g2_front: E::G2,
    pub matrixed_g1_back: E::G1,
    pub matrixed_g2_back: E::G2,
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
    //on-chain
    let mut paramter_in_storage = initial_uncommon_paramters::<E>(common_paramter_matrix);
    //under-chain
    let player1_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (1, 2));
    //on-chain
    paramter_in_storage = verify_uncommon_paramter(&paramter_in_storage, &player1_common);
    //under-chain
    let player2_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (2, 3));
    //on-chain
    paramter_in_storage = verify_uncommon_paramter(&paramter_in_storage, &player2_common);
    //under-chain
    let player3_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (3, 4));
    //on-chain
    paramter_in_storage = verify_uncommon_paramter(&paramter_in_storage, &player3_common);

    //on-chain
    paramter_in_storage;

    unimplemented!()
}

pub fn initial_uncommon_paramters<E>(
    common_paramter_matrix: &CommonParamterMatrix<E>,
) -> UnCommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let len = 0;
    let init_paramters = 0;
    unimplemented!()
}

pub fn mpc_uncommon_paramters_generator<E>(
    storage: &UnCommonParamterInStorage<E>,
    (gamma, delta): (i32, i32),
) -> UnCommonParamter<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    unimplemented!();
}

pub fn verify_uncommon_paramter<E>(
    storage: &UnCommonParamterInStorage<E>,
    new_paramter: &UnCommonParamter<E>,
) -> UnCommonParamterInStorage<E>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    unimplemented!();
}
