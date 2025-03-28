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

/// Generates a random common reference string for
/// a circuit.
pub fn generate_random_parameters<E, C, R>(
    circuit: C,
    mut rng: &mut R,
) -> Result<Parameters<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
    C: Circuit<E::Fr>,
    R: RngCore,
{
    let g1 = E::G1::generator();
    let g2 = E::G2::generator();
    let alpha = E::Fr::from(6);
    let beta = E::Fr::from(24);
    let gamma = E::Fr::from(6);
    let delta = E::Fr::from(24);
    let tau = E::Fr::from(2);
    generate_parameters::<E, C>(circuit, g1, g2, alpha, beta, gamma, delta, tau)
}

/// This is our assembly structure that we'll use to synthesize the
/// circuit into a QAP.
pub struct KeypairAssembly<Scalar: PrimeField> {
    num_inputs: usize,
    num_aux: usize,
    num_constraints: usize,
    at_inputs: Vec<Vec<(Scalar, usize)>>,
    bt_inputs: Vec<Vec<(Scalar, usize)>>,
    ct_inputs: Vec<Vec<(Scalar, usize)>>,
    at_aux: Vec<Vec<(Scalar, usize)>>,
    bt_aux: Vec<Vec<(Scalar, usize)>>,
    ct_aux: Vec<Vec<(Scalar, usize)>>,
}

impl<Scalar: PrimeField> ConstraintSystem<Scalar> for KeypairAssembly<Scalar> {
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<Scalar, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't even invoke the
        // function for obtaining one.

        let index = self.num_aux;
        self.num_aux += 1;
        self.at_aux.push(vec![]);
        self.bt_aux.push(vec![]);
        self.ct_aux.push(vec![]);

        Ok(Variable(Index::Aux(index)))
    }

    fn alloc_input<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<Scalar, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't even invoke the
        // function for obtaining one.

        let index = self.num_inputs;
        self.num_inputs += 1;

        self.at_inputs.push(vec![]);
        self.bt_inputs.push(vec![]);
        self.ct_inputs.push(vec![]);

        Ok(Variable(Index::Input(index)))
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<Scalar>) -> LinearCombination<Scalar>,
        LB: FnOnce(LinearCombination<Scalar>) -> LinearCombination<Scalar>,
        LC: FnOnce(LinearCombination<Scalar>) -> LinearCombination<Scalar>,
    {
        fn eval<Scalar: PrimeField>(
            l: LinearCombination<Scalar>,
            inputs: &mut [Vec<(Scalar, usize)>],
            aux: &mut [Vec<(Scalar, usize)>],
            this_constraint: usize,
        ) {
            for (index, coeff) in l.0 {
                println!("这次进入的参数是：{:?}", coeff);
                match index {
                    Variable(Index::Input(id)) => inputs[id].push((coeff, this_constraint)),
                    Variable(Index::Aux(id)) => aux[id].push((coeff, this_constraint)),
                }
            }
        }

        eval(
            a(LinearCombination::zero()),
            &mut self.at_inputs,
            &mut self.at_aux,
            self.num_constraints,
        );
        eval(
            b(LinearCombination::zero()),
            &mut self.bt_inputs,
            &mut self.bt_aux,
            self.num_constraints,
        );
        eval(
            c(LinearCombination::zero()),
            &mut self.ct_inputs,
            &mut self.ct_aux,
            self.num_constraints,
        );

        self.num_constraints += 1;
    }

    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn pop_namespace(&mut self) {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }
}

use crate::groth16::mpc::{
    mpc_common_paramters_custom_all, mpc_uncommon_paramters_custom_all, CommonParamter,
    UnCommonParamter,
};
#[allow(clippy::too_many_arguments)]
pub fn generate_parameters_mpc<E, C>(
    circuit: C,
    g1: E::G1,
    g2: E::G2,
) -> Result<Parameters<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
    C: Circuit<E::Fr>,
{
    let mut assembly = KeypairAssembly {
        num_inputs: 0,
        num_aux: 0,
        num_constraints: 0,
        at_inputs: vec![],
        bt_inputs: vec![],
        ct_inputs: vec![],
        at_aux: vec![],
        bt_aux: vec![],
        ct_aux: vec![],
    };

    // Allocate the "one" input variable
    assembly.alloc_input(|| "", || Ok(E::Fr::one()))?;

    circuit.synthesize(&mut assembly)?;
    for i in 0..assembly.num_inputs {
        assembly.enforce(|| "", |lc| lc + Variable(Index::Input(i)), |lc| lc, |lc| lc);
    }
    let mut a = vec![E::G1Affine::identity(); assembly.num_inputs + assembly.num_aux];
    let mut b_g1 = vec![E::G1Affine::identity(); assembly.num_inputs + assembly.num_aux];
    let mut b_g2 = vec![E::G2Affine::identity(); assembly.num_inputs + assembly.num_aux];

    let cp = mpc_common_paramters_custom_all::<E>();
    let cp_m = cp.matrix(
        &assembly.at_aux,
        &assembly.bt_aux,
        &assembly.ct_aux,
        assembly.num_inputs,
        assembly.num_aux,
        assembly.num_constraints,
    );
    let ucp = mpc_uncommon_paramters_custom_all::<E>(&cp_m);
    let vk = VerifyingKey::<E> {
        alpha_g1: cp.alpha_g1.to_affine(),
        beta_g1: cp.beta_g1.to_affine(),
        beta_g2: cp.beta_g2.to_affine(),
        gamma_g2: ucp.gamma_g2.to_affine(),
        delta_g1: ucp.delta_g1.to_affine(),
        delta_g2: ucp.delta_g2.to_affine(),
        ic: vec_to_list::<E>(&ucp.kin_g1),
    };
    Ok(Parameters {
        vk,
        h: Arc::new(vec_to_list::<E>(&ucp.h_g1)),
        l: Arc::new(vec_to_list::<E>(&ucp.kout_g1)),
        // Filter points at infinity away from A/B queries
        a: Arc::new(
            a.into_iter()
                .filter(|e| bool::from(!e.is_identity()))
                .collect(),
        ),
        b_g1: Arc::new(
            b_g1.into_iter()
                .filter(|e| bool::from(!e.is_identity()))
                .collect(),
        ),
        b_g2: Arc::new(
            b_g2.into_iter()
                .filter(|e| bool::from(!e.is_identity()))
                .collect(),
        ),
    })
}

/// Create parameters for a circuit, given some toxic waste.
#[allow(clippy::too_many_arguments)]
pub fn generate_parameters<E, C>(
    circuit: C,
    g1: E::G1,
    g2: E::G2,
    alpha: E::Fr,
    beta: E::Fr,
    gamma: E::Fr,
    delta: E::Fr,
    tau: E::Fr,
) -> Result<Parameters<E>, SynthesisError>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
    C: Circuit<E::Fr>,
{
    let mut assembly = KeypairAssembly {
        num_inputs: 0,
        num_aux: 0,
        num_constraints: 0,
        at_inputs: vec![],
        bt_inputs: vec![],
        ct_inputs: vec![],
        at_aux: vec![],
        bt_aux: vec![],
        ct_aux: vec![],
    };
    // Allocate the "one" input variable
    assembly.alloc_input(|| "", || Ok(E::Fr::one()))?;

    // Synthesize the circuit.
    circuit.synthesize(&mut assembly)?;
    println!("没声明时");
    println!("at_aux{:?}", assembly.at_aux);
    println!("bt_aux{:?}", assembly.bt_aux);
    println!("ct_aux{:?}", assembly.ct_aux);
    // Input constraints to ensure full density of IC query
    // x * 0 = 0
    for i in 0..assembly.num_inputs {
        assembly.enforce(|| "", |lc| lc + Variable(Index::Input(i)), |lc| lc, |lc| lc);
    }
    println!("声明后");
    println!("at_aux{:?}", assembly.at_aux);
    println!("bt_aux{:?}", assembly.bt_aux);
    println!("ct_aux{:?}", assembly.ct_aux);
    println!("num_aux :{}", assembly.num_aux);
    println!("num_input :{}", assembly.num_inputs);

    {
        let fuck = assembly.at_aux[0][0].0;
        let alpha_front = g1 * fuck;
    }

    // Create bases for blind evaluation of polynomials at tau
    let powers_of_tau = vec![Scalar::<E::Fr>(E::Fr::zero()); assembly.num_constraints];

    let mut powers_of_tau = EvaluationDomain::from_coeffs(powers_of_tau)?;
    let cp = mpc_common_paramters_custom_all::<E>();
    let cp_m = cp.matrix_test(
        &assembly.at_aux,
        &assembly.bt_aux,
        &assembly.ct_aux,
        assembly.num_inputs,
        assembly.num_aux,
        assembly.num_constraints,
    );

    let ucp = mpc_uncommon_paramters_custom_all::<E>(&cp_m);

    // Compute G1 window table
    let mut g1_wnaf = Wnaf::new();
    let g1_wnaf = g1_wnaf.base(g1, {
        // H query
        (powers_of_tau.as_ref().len() - 1)
        // IC/L queries
        + assembly.num_inputs + assembly.num_aux
        // A query
        + assembly.num_inputs + assembly.num_aux
        // B query
        + assembly.num_inputs + assembly.num_aux
    });

    // Compute G2 window table
    let mut g2_wnaf = Wnaf::new();
    let g2_wnaf = g2_wnaf.base(g2, {
        // B query
        assembly.num_inputs + assembly.num_aux
    });

    let gamma_inverse = {
        let inverse = gamma.invert();
        if bool::from(inverse.is_some()) {
            Ok(inverse.unwrap())
        } else {
            Err(SynthesisError::UnexpectedIdentity)
        }
    }?;
    let delta_inverse = {
        let inverse = delta.invert();
        if bool::from(inverse.is_some()) {
            Ok(inverse.unwrap())
        } else {
            Err(SynthesisError::UnexpectedIdentity)
        }
    }?;

    let worker = Worker::new();

    let mut h = vec![E::G1Affine::identity(); powers_of_tau.as_ref().len() - 1];
    {
        // Compute powers of tau
        {
            let powers_of_tau = powers_of_tau.as_mut();
            worker.scope(powers_of_tau.len(), |scope, chunk| {
                for (i, powers_of_tau) in powers_of_tau.chunks_mut(chunk).enumerate() {
                    scope.spawn(move |_scope| {
                        let mut current_tau_power = tau.pow_vartime(&[(i * chunk) as u64]);

                        for p in powers_of_tau {
                            p.0 = current_tau_power;
                            current_tau_power.mul_assign(&tau);
                        }
                    });
                }
            });
        }
        // coeff = t(x) / delta
        let mut coeff = powers_of_tau.z(&tau);
        coeff.mul_assign(&delta_inverse);

        // Compute the H query with multiple threads
        worker.scope(h.len(), |scope, chunk| {
            for (h, p) in h
                .chunks_mut(chunk)
                .zip(powers_of_tau.as_ref().chunks(chunk))
            {
                let mut g1_wnaf = g1_wnaf.shared();

                scope.spawn(move |_scope| {
                    // Set values of the H query to g1^{(tau^i * t(tau)) / delta}
                    let h_proj: Vec<_> = p[..h.len()]
                        .iter()
                        .map(|p| {
                            // Compute final exponent
                            let mut exp = p.0;
                            exp.mul_assign(&coeff);

                            // Exponentiate
                            g1_wnaf.scalar(&exp)
                        })
                        .collect();

                    // Batch normalize
                    E::G1::batch_normalize(&h_proj, h);
                });
            }
        });
    }

    // Use inverse FFT to convert powers of tau to Lagrange coefficients
    powers_of_tau.ifft(&worker);
    let powers_of_tau = powers_of_tau.into_coeffs();
    /*let p_l = powers_of_tau.len();
    let mut tau_list_pl_g1 = vec![E::G1Affine::identity(); p_l];
    let mut tau_list_pl_g2 = vec![E::G2Affine::identity(); p_l];
    for i in 0..p_l {
        tau_list_pl_g1[i] = (E::G1Affine::generator() * powers_of_tau[i].0).to_affine();

        tau_list_pl_g2[i] = (E::G2Affine::generator() * powers_of_tau[i].0).to_affine();
    };*/

    let mut a = vec![E::G1Affine::identity(); assembly.num_inputs + assembly.num_aux];
    let mut b_g1 = vec![E::G1Affine::identity(); assembly.num_inputs + assembly.num_aux];
    let mut b_g2 = vec![E::G2Affine::identity(); assembly.num_inputs + assembly.num_aux];
    let mut ic = vec![E::G1Affine::identity(); assembly.num_inputs];
    let mut l = vec![E::G1Affine::identity(); assembly.num_aux];

    #[allow(clippy::too_many_arguments)]
    fn eval<E: Engine>(
        // wNAF window tables
        g1_wnaf: &Wnaf<usize, &[E::G1], &mut Vec<i64>>,
        g2_wnaf: &Wnaf<usize, &[E::G2], &mut Vec<i64>>,

        // Lagrange coefficients for tau
        powers_of_tau: &[Scalar<E::Fr>],

        // QAP polynomials
        at: &[Vec<(E::Fr, usize)>],
        bt: &[Vec<(E::Fr, usize)>],
        ct: &[Vec<(E::Fr, usize)>],

        // Resulting evaluated QAP polynomials
        a: &mut [E::G1Affine],
        b_g1: &mut [E::G1Affine],
        b_g2: &mut [E::G2Affine],
        ext: &mut [E::G1Affine],

        // Inverse coefficient for ext elements
        inv: &E::Fr,

        // Trapdoors
        alpha: &E::Fr,
        beta: &E::Fr,

        // Worker
        worker: &Worker,
    ) {
        // Sanity check
        assert_eq!(a.len(), at.len());
        assert_eq!(a.len(), bt.len());
        assert_eq!(a.len(), ct.len());
        assert_eq!(a.len(), b_g1.len());
        assert_eq!(a.len(), b_g2.len());
        assert_eq!(a.len(), ext.len());

        // Evaluate polynomials in multiple threads
        worker.scope(a.len(), |scope, chunk| {
            for ((((((a, b_g1), b_g2), ext), at), bt), ct) in a
                .chunks_mut(chunk)
                .zip(b_g1.chunks_mut(chunk))
                .zip(b_g2.chunks_mut(chunk))
                .zip(ext.chunks_mut(chunk))
                .zip(at.chunks(chunk))
                .zip(bt.chunks(chunk))
                .zip(ct.chunks(chunk))
            {
                let mut g1_wnaf = g1_wnaf.shared();
                let mut g2_wnaf = g2_wnaf.shared();

                scope.spawn(move |_scope| {
                    let mut a_proj = vec![E::G1::identity(); a.len()];
                    let mut b_g1_proj = vec![E::G1::identity(); b_g1.len()];
                    let mut b_g2_proj = vec![E::G2::identity(); b_g2.len()];
                    let mut ext_proj = vec![E::G1::identity(); ext.len()];

                    for ((((((a, b_g1), b_g2), ext), at), bt), ct) in a_proj
                        .iter_mut()
                        .zip(b_g1_proj.iter_mut())
                        .zip(b_g2_proj.iter_mut())
                        .zip(ext_proj.iter_mut())
                        .zip(at.iter())
                        .zip(bt.iter())
                        .zip(ct.iter())
                    {
                        fn eval_at_tau<S: PrimeField>(
                            powers_of_tau: &[Scalar<S>],
                            p: &[(S, usize)],
                        ) -> S {
                            let mut acc = S::zero();

                            for &(ref coeff, index) in p {
                                let mut n = powers_of_tau[index].0;

                                n.mul_assign(coeff);
                                acc.add_assign(&n);
                            }

                            acc
                        }

                        // Evaluate QAP polynomials at tau
                        let mut at = eval_at_tau(powers_of_tau, at);
                        let mut bt = eval_at_tau(powers_of_tau, bt);
                        let ct = eval_at_tau(powers_of_tau, ct);

                        // Compute A query (in G1)
                        if !at.is_zero_vartime() {
                            *a = g1_wnaf.scalar(&at);
                        }

                        // Compute B query (in G1/G2)
                        if !bt.is_zero_vartime() {
                            *b_g1 = g1_wnaf.scalar(&bt);
                            *b_g2 = g2_wnaf.scalar(&bt);
                        }

                        at *= beta;
                        bt *= alpha;

                        let mut e = at;
                        e.add_assign(&bt);
                        e.add_assign(&ct);
                        e.mul_assign(inv);

                        *ext = g1_wnaf.scalar(&e);
                    }

                    // Batch normalize
                    E::G1::batch_normalize(&a_proj, a);
                    E::G1::batch_normalize(&b_g1_proj, b_g1);
                    E::G2::batch_normalize(&b_g2_proj, b_g2);
                    E::G1::batch_normalize(&ext_proj, ext);
                });
            }
        });
    }

    // Evaluate for inputs.
    eval::<E>(
        &g1_wnaf,
        &g2_wnaf,
        &powers_of_tau,
        &assembly.at_inputs,
        &assembly.bt_inputs,
        &assembly.ct_inputs,
        &mut a[0..assembly.num_inputs],
        &mut b_g1[0..assembly.num_inputs],
        &mut b_g2[0..assembly.num_inputs],
        &mut ic,
        &gamma_inverse,
        &alpha,
        &beta,
        &worker,
    );

    // Evaluate for auxiliary variables.
    eval::<E>(
        &g1_wnaf,
        &g2_wnaf,
        &powers_of_tau,
        &assembly.at_aux,
        &assembly.bt_aux,
        &assembly.ct_aux,
        &mut a[assembly.num_inputs..],
        &mut b_g1[assembly.num_inputs..],
        &mut b_g2[assembly.num_inputs..],
        &mut l,
        &delta_inverse,
        &alpha,
        &beta,
        &worker,
    );
    assert_eq!(cp.tau_g1[1], g1 * E::Fr::from(2));
    assert_eq!(cp.alpha_mul_tau_g1[0], g1 * E::Fr::from(6));
    assert_eq!(cp.alpha_mul_tau_g1[1], g1 * E::Fr::from(6 * 2));
    assert_eq!(cp.beta_mul_tau_g1[0], g1 * E::Fr::from(24));
    assert_eq!(cp.beta_mul_tau_g1[1], g1 * E::Fr::from(24 * 2));
    //assert_eq!(ic[1] * gamma_inverse, g1 * E::Fr::from(0)); //?
    /*for i in 0..200 {
        println!("第{}个是吗", i);
        assert_ne!(ic[1], (g1 * E::Fr::from(i)).to_affine());
    }
    assert_eq!(ic[1], (g1 * E::Fr::from(0)).to_affine());*/
    // Don't allow any elements be unconstrained, so that
    // the L query is always fully dense.
    for e in l.iter() {
        if e.is_identity().into() {
            return Err(SynthesisError::UnconstrainedVariable);
        }
    }

    assert_eq!(h[0], ucp.h_g1[0].to_affine());
    assert_eq!(h[1], ucp.h_g1[1].to_affine());
    let g1 = g1.to_affine();
    let g2 = g2.to_affine();
    let vk = VerifyingKey::<E> {
        alpha_g1: (g1 * alpha).to_affine(),
        beta_g1: (g1 * beta).to_affine(),
        beta_g2: (g2 * beta).to_affine(),
        gamma_g2: (g2 * gamma).to_affine(),
        delta_g1: (g1 * delta).to_affine(),
        delta_g2: (g2 * delta).to_affine(),
        ic,
    };
    //检查vk
    assert_eq!(vk.alpha_g1, cp.alpha_g1.to_affine());
    assert_eq!(vk.beta_g1, cp.beta_g1.to_affine());
    assert_eq!(vk.beta_g2, cp.beta_g2.to_affine());
    assert_eq!(vk.gamma_g2, ucp.gamma_g2.to_affine());
    assert_eq!(vk.delta_g1, ucp.delta_g1.to_affine());
    assert_eq!(vk.delta_g2, ucp.delta_g2.to_affine());
    Ok(Parameters {
        vk,
        h: Arc::new(h),
        l: Arc::new(l),

        // Filter points at infinity away from A/B queries
        a: Arc::new(
            a.into_iter()
                .filter(|e| bool::from(!e.is_identity()))
                .collect(),
        ),
        b_g1: Arc::new(
            b_g1.into_iter()
                .filter(|e| bool::from(!e.is_identity()))
                .collect(),
        ),
        b_g2: Arc::new(
            b_g2.into_iter()
                .filter(|e| bool::from(!e.is_identity()))
                .collect(),
        ),
    })
}

pub fn vec_to_list<E>(list: &Vec<E::G1>) -> Vec<E::G1Affine>
where
    E: Engine,
    E::G1: WnafGroup,
    E::G2: WnafGroup,
{
    let mut result = Vec::new();
    for i in 0..list.len() {
        result.push(list[i].to_affine());
    }
    result
}
