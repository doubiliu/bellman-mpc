// For randomness (during paramgen and proof generation)
use rand::thread_rng;

// For benchmarking
use std::time::{Duration, Instant};

// Bring in some tools for using finite fiels
use ff::Field;

// We're going to use the BLS12-381 pairing-friendly elliptic curve.
use bls12_381::{Bls12, Scalar};

// We're going to use the Groth16 proving system.
use crate::groth16::{
    batch, create_proof, create_random_proof, generate_parameters, generate_random_parameters,
    prepare_verifying_key, verify_proof, Parameters, Proof,
};

use crate::range_mod::*;

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
            w:Some(9u64),
            wArray:Some([0u64,0u64,0u64,0u64]),
            less_or_equal:Some(1u64),
            less:Some(1u64),
            not_all_zeros:Some(1u64),
            crArray:Some([0u64,0u64,0u64,0u64]),
            tt: None,
        };

        //generate_parameters::<Bls12, _>(c, g1, g2, alpha, beta, gamma, delta, tau).unwrap()
        generate_random_parameters::<Bls12, _, _>(c, &mut rng).unwrap()
    };
    let pvk = prepare_verifying_key(&params.vk);

    println!("Creating proofs...");

    let r = Scalar::from(27134);
    let s = Scalar::from(17146);
    let c_except = Scalar::from(4u64);
    let proof = {
        let c = RangeDemo {
            a: Some(1u64),
            b: Some(2u64),
            n: Some(4u64),

            w:Some(9u64),
            wArray:Some([1u64,0u64,0u64,1u64]),
            less_or_equal:Some(1u64),
            less:Some(1u64),
            not_all_zeros:Some(1u64),
            crArray:Some([1u64,1u64,1u64,1u64]),

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
