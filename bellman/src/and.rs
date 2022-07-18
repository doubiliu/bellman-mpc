// For randomness (during paramgen and proof generation)
use rand::thread_rng;

// For benchmarking
use std::time::{Duration, Instant};

// Bring in some tools for using finite fiels
use ff::Field;

// We're going to use the BLS12-381 pairing-friendly elliptic curve.
use bls12_381::{Bls12, Scalar};

// We're going to use the Groth16 proving system.
use crate::and_mod::{AndDemo, RangeDemo};
use crate::groth16::mpc::{
    initParameterList, initTauParameterList, mpc_bad_paramters_custom, mpc_common_paramters_custom,
    mpc_common_paramters_custom_generator, mpc_common_tauparamters_custom,
    mpc_common_tauparamters_custom_generator, ParameterPair, ParamterListExcute,
    TauParamterListExcute,
};
use crate::groth16::{
    batch, create_proof, create_random_proof, generate_parameters, generate_random_parameters,
    prepare_verifying_key, verify_proof, Parameters, Proof,
};

#[test]
fn test_anddemo_bls12() {
    let mut rng = thread_rng();
    let params = {
        let c = AndDemo {
            a: None,
            b: None,
            tt: None,
        };

        //generate_parameters::<Bls12, _>(c, g1, g2, alpha, beta, gamma, delta, tau).unwrap()
        generate_random_parameters::<Bls12, _, _>(c, &mut rng).unwrap()
    };
    let pvk = prepare_verifying_key(&params.vk);

    println!("Creating proofs...");

    let r = Scalar::from(27134u64);
    let s = Scalar::from(17146u64);
    let c_except = Scalar::zero();
    let proof = {
        let c = AndDemo {
            a: Some(true),
            b: Some(false),
            tt: None,
        };
        //create_random_proof(c, &params, &mut rng).unwrap()
        create_proof(c, &params, r, s).unwrap()
    };

    assert!(verify_proof(&pvk, &proof, &[c_except]).is_ok());

    let proof_a = proof.a.to_uncompressed();

    let proof_b = proof.b.to_uncompressed();

    let proof_c = proof.c.to_uncompressed();

    //println!("A Point: {:?}", proof_a);

    //println!("B Point: {:?}", proof_b);

    //println!("C Point: {:?}", proof_c);

    let vk = params.vk;

    let alpha_g1 = vk.alpha_g1.to_uncompressed();

    let beta_g2 = vk.beta_g2.to_uncompressed();
    // gamma in g2 for verifying. Never the point at infinity.
    let gamma_g2 = vk.gamma_g2.to_uncompressed();

    let delta_g2 = vk.delta_g2.to_uncompressed();

    //println!("alpha_g1 Point: {:?}", alpha_g1);

    //println!("beta_g2 Point: {:?}", beta_g2);

    //println!("gamma_g2 Point: {:?}", gamma_g2);

    //println!("delta_g2 Point: {:?}", delta_g2);

    let mut acc = pvk.ic[0].to_uncompressed();

    let pvkiclen = pvk.ic.len();

    let ic1 = pvk.ic[1].to_uncompressed();

    //println!("ic0 Point: {:?}", acc);

    //println!("ic1 Point: {:?}", ic1);

    //println!("ic len: {:?}", pvkiclen);
}
#[test]
fn test_rangedemo_bls12() {
    let g1 = Scalar::one();
    let g2 = Scalar::one();
    let alpha = Scalar::from(48577);
    let beta = Scalar::from(22580);
    let gamma = Scalar::from(53332);
    let delta = Scalar::from(5481);
    let tau = Scalar::from(3673);
    let mut rng = thread_rng();
    let params = {
        let c = RangeDemo {
            a: Some(1u64),
            b: Some(2u64),
            n: Some(4u64),
            w: Some(9u64),
            wArray: Some([0u64, 0u64, 0u64, 0u64]),
            less_or_equal: Some(1u64),
            less: Some(1u64),
            not_all_zeros: Some(1u64),
            tt: None,
        };

        //generate_parameters::<Bls12, _>(c, g1, g2, alpha, beta, gamma, delta, tau).unwrap()
        generate_random_parameters::<Bls12, _, _>(c, &mut rng).unwrap()
    };
    let pvk = prepare_verifying_key(&params.vk);

    println!("Creating proofs...");

    let r = Scalar::from(27134);
    let s = Scalar::from(17146);
    let c_except = Scalar::from(2u64);
    let proof = {
        let c = RangeDemo {
            a: Some(1u64),
            b: Some(2u64),
            n: Some(4u64),
            w: Some(9u64),
            wArray: Some([1u64, 0u64, 0u64, 1u64]),
            less_or_equal: Some(1u64),
            less: Some(1u64),
            not_all_zeros: Some(1u64),

            tt: None,
        };
        //create_random_proof(c, &params, &mut rng).unwrap()
        create_proof(c, &params, r, s).unwrap()
    };

    assert!(verify_proof(&pvk, &proof, &[c_except]).is_ok());

    let proof_a = proof.a.to_uncompressed();

    let proof_b = proof.b.to_uncompressed();

    let proof_c = proof.c.to_uncompressed();

    println!("A Point: {:?}", proof_a);

    println!("B Point: {:?}", proof_b);

    println!("C Point: {:?}", proof_c);

    let vk = params.vk;

    let alpha_g1_beta_g2 = pvk.alpha_g1_beta_g2;

    let ic = pvk.ic;
}
/*
#[test]
pub fn test_mpc_alpha() {
    let player1 = Scalar::from(1);
    let player2 = Scalar::from(2);
    let player3 = Scalar::from(3);
    //链上：初始化参数列表，一个g1
    let mut list = initParameterList::<Bls12>();
    //链下：player1计算
    let player1_param = mpc_common_paramters_custom_generator::<Bls12>(&list[0], player1);
    //链上：验证player1，成功则添加
    list = ParamterListExcute::<Bls12>(list, player1_param.unwrap());
    //链下：player1计算
    let player2_param = mpc_common_paramters_custom_generator::<Bls12>(&list[1], player2);
    //链上：验证player2，成功则添加
    list = ParamterListExcute::<Bls12>(list, player2_param.unwrap());
    //链下：player3计算
    let player3_param = mpc_common_paramters_custom_generator::<Bls12>(&list[2], player3);
    //捣乱的
    //let badman3_params = mpc_bad_paramters_custom::<Bls12>(&list[2], player3);
    //链上：验证player3，成功则添加
    list = ParamterListExcute::<Bls12>(list, player3_param.unwrap());
    //list = ParamterListExcute::<Bls12>(list, badman3_params.unwrap());
    //验证list[3]与6*g1是否相等'
    let g1 = bls12_381::G1Affine::generator();
    let scalar = player1 * player2 * player3;
    let standard = bls12_381::G1Affine::from(g1 * scalar);
    assert_eq!(list[3].g1_result.unwrap(), standard);
}

#[test]
pub fn test_mpc_tau() {
    let player3 = [Scalar::from(2), Scalar::from(4), Scalar::from(8)].to_vec();
    let player2 = [Scalar::from(3), Scalar::from(9), Scalar::from(27)].to_vec();
    let player1 = [Scalar::from(1), Scalar::from(1), Scalar::from(1)].to_vec();
    let player0 = [Scalar::from(1), Scalar::from(1), Scalar::from(1)].to_vec();

    //链上：初始化参数列表
    let mut list = initTauParameterList::<Bls12>(player1.len());

    let player0_param = mpc_common_tauparamters_custom_generator::<Bls12>(&list[0], player0);

    println!("player0 create custom done");
    //链上：验证player1，成功则添加
    list = TauParamterListExcute::<Bls12>(list, player0_param.unwrap());
    println!("player0 verify done");

    //链下：player1计算
    let player1_param = mpc_common_tauparamters_custom_generator::<Bls12>(&list[1], player1);
    println!("player1 create custom done");
    //链上：验证player1，成功则添加
    list = TauParamterListExcute::<Bls12>(list, player1_param.unwrap());
    println!("player1 verify done");
    //链下：player1计算
    let player2_param = mpc_common_tauparamters_custom_generator::<Bls12>(&list[2], player2);
    println!("player2 create custom done");
    //链上：验证player2，成功则添加
    list = TauParamterListExcute::<Bls12>(list, player2_param.unwrap());
    println!("player2 verify done");
    //链下：player3计算
    let player3_param = mpc_common_tauparamters_custom_generator::<Bls12>(&list[3], player3);
    println!("player3 create custom done");
    //捣乱的
    //let badman3_params = mpc_bad_paramters_custom::<Bls12>(&list[2], player3);
    //链上：验证player3，成功则添加
    list = TauParamterListExcute::<Bls12>(list, player3_param.unwrap());
    println!("player3 verify done");
    //list = ParamterListExcute::<Bls12>(list, badman3_params.unwrap());
    let g1 = bls12_381::G1Affine::generator();
    let scalar = [6, 36, 216];
    let standard = [
        bls12_381::G1Affine::from(g1 * Scalar::from(scalar[0])),
        bls12_381::G1Affine::from(g1 * Scalar::from(scalar[1])),
        bls12_381::G1Affine::from(g1 * Scalar::from(scalar[2])),
    ];
    assert_eq!(list[4].list[0].g1_result.unwrap(), standard[0]);
    assert_eq!(list[4].list[1].g1_result.unwrap(), standard[1]);
    assert_eq!(list[4].list[2].g1_result.unwrap(), standard[2]);
}*/
