#[cfg(test)]
mod mpc_tests {
    use super::*;

    use crate::groth16::mpc::*;
    use bls12_381::{Bls12, Scalar};

    #[test]
    fn all_test() {
        let g1 = bls12_381::G1Affine::generator();
        let g2 = bls12_381::G2Affine::generator();
        let mut paramter_in_storage =
            initial_common_paramters::<Bls12>(8 /*这个数表示最高到x的几次幂*/);
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

        let at_aux = vec![
            vec![(Scalar::from(1), 0usize), (Scalar::from(2), 1usize)],
            vec![],
        ];
        let bt_aux = vec![
            vec![(Scalar::from(1), 0usize), (Scalar::from(2), 1usize)],
            vec![(Scalar::from(3), 0usize), (Scalar::from(4), 1usize)],
        ];
        let ct_aux = vec![vec![], vec![]];

        let common_paramter_matrix = paramter_in_storage.matrix(&at_aux, &bt_aux, &ct_aux, 2, 2, 4);
        let mut paramter_in_storage = initial_uncommon_paramters::<Bls12>(&common_paramter_matrix);
        //under-chain
        let player1_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (1, 2));
        paramter_in_storage = verify_uncommon_paramter(
            &common_paramter_matrix,
            &paramter_in_storage,
            &player1_common,
        );
        let player2_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (2, 3));
        //on-chain
        paramter_in_storage = verify_uncommon_paramter(
            &common_paramter_matrix,
            &paramter_in_storage,
            &player2_common,
        );
        //under-chain
        let player3_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (3, 4));
        //on-chain
        paramter_in_storage = verify_uncommon_paramter(
            &common_paramter_matrix,
            &paramter_in_storage,
            &player3_common,
        );
    }
    #[test]
    fn test_h() {
        let g1 = bls12_381::G1Affine::generator();
        println!("g1*-1:{:?}", bls12_381::G1Affine::from(-g1));
    }
    #[test]
    fn HexToScalar() {}

    #[test]
    #[test]
    fn common_works() {
        //验证通用系数
        let g1 = bls12_381::G1Affine::generator();
        let g2 = bls12_381::G2Affine::generator();

        let mut paramter_in_storage =
            initial_common_paramters::<Bls12>(8 /*这个数表示最高到x的几次幂*/);
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
        /*let alpha_g1 = paramter_in_storage.alpha_g1;
        let alpha_g2 = paramter_in_storage.alpha_g2;
        let beta_g1 = paramter_in_storage.beta_g1;
        let beta_g2 = paramter_in_storage.beta_g2;
        let tau_g1 = paramter_in_storage.tau_g1;
        let tau_g2 = paramter_in_storage.tau_g2;
        let alpha_mul_tau_g1 = paramter_in_storage.alpha_mul_tau_g1;
        let alpha_mul_tau_g2 = paramter_in_storage.alpha_mul_tau_g2;
        let beta_mul_tau_g1 = paramter_in_storage.beta_mul_tau_g1;
        let beta_mul_tau_g2 = paramter_in_storage.beta_mul_tau_g2;

        assert_eq!(alpha_g1, g1 * Scalar::from(6));
        assert_eq!(alpha_g2, g2 * Scalar::from(6));
        assert_eq!(beta_g1, g1 * Scalar::from(24));
        assert_eq!(beta_g2, g2 * Scalar::from(24));
        assert_eq!(tau_g1.len(), 8);
        assert_eq!(tau_g1[0], g1 * Scalar::from(1));
        assert_eq!(tau_g2[0], g2 * Scalar::from(1));
        assert_eq!(tau_g1[1], g1 * Scalar::from(60));
        assert_eq!(tau_g2[1], g2 * Scalar::from(60));
        assert_eq!(tau_g1[2], g1 * Scalar::from(3600));
        assert_eq!(tau_g2[2], g2 * Scalar::from(3600));
        assert_eq!(alpha_mul_tau_g1.len(), 8);
        assert_eq!(alpha_mul_tau_g1[0], g1 * Scalar::from(1 * 6));
        assert_eq!(alpha_mul_tau_g2[0], g2 * Scalar::from(1 * 6));
        assert_eq!(alpha_mul_tau_g1[1], g1 * Scalar::from(60 * 6));
        assert_eq!(alpha_mul_tau_g2[1], g2 * Scalar::from(60 * 6));
        assert_eq!(alpha_mul_tau_g1[2], g1 * Scalar::from(3600 * 6));
        assert_eq!(alpha_mul_tau_g2[2], g2 * Scalar::from(3600 * 6));
        assert_eq!(beta_mul_tau_g1.len(), 8);
        assert_eq!(beta_mul_tau_g1[0], g1 * Scalar::from(1 * 24));
        assert_eq!(beta_mul_tau_g2[0], g2 * Scalar::from(1 * 24));
        assert_eq!(beta_mul_tau_g1[1], g1 * Scalar::from(60 * 24));
        assert_eq!(beta_mul_tau_g2[1], g2 * Scalar::from(60 * 24));*/

        let matrix1 = vec![
            vec![(Scalar::from(1), 0usize), (Scalar::from(1), 1usize)],
            vec![],
        ];
        let matrix2 = vec![
            vec![(Scalar::from(1), 0usize)],
            vec![(Scalar::from(1), 1usize)],
        ];
        let matrix3: Vec<Vec<(Scalar, usize)>> = Vec::new();

        let cp_m = paramter_in_storage.matrix_test(&matrix1, &matrix2, &matrix3, 2, 2, 4);

        assert_eq!(cp_m.matrixed_g1_front[0], g1 * Scalar::from(24 * 61 + 6));
        assert_eq!(cp_m.matrixed_g1_front[1], g1 * Scalar::from(6 * 60));
    }

    #[test]
    fn matrix_works() {
        let g1 = bls12_381::G1Affine::generator();
        let g2 = bls12_381::G2Affine::generator();
        let x1 = g1 * Scalar::from(60);
        let x2 = g1 * Scalar::from(240);
        println!("{:?}", bls12_381::G1Affine::from(x1));
        println!("{:?}", bls12_381::G1Affine::from(x2));
        //let matrix3:Vec<Vec<(Scalar, usize)>> = Vec::new();
    }
    #[test]
    fn uncommonn_works() {
        let g1 = bls12_381::G1Affine::generator();
        let g2 = bls12_381::G2Affine::generator();

        let common_paramter_matrix = CommonParamterMatrix::<Bls12> {
            matrixed_g1_front: vec![g1 * Scalar::from(6), g1 * Scalar::from(12)],
            matrixed_g2_front: vec![g2 * Scalar::from(6), g2 * Scalar::from(12)],
            matrixed_g1_back: vec![g1 * Scalar::from(24), g1 * Scalar::from(48)],
            matrixed_g2_back: vec![g2 * Scalar::from(24), g2 * Scalar::from(48)],
            matrixed_h_g1: vec![
                g1 * Scalar::from(2),
                g1 * Scalar::from(4),
                g1 * Scalar::from(6),
                g1 * Scalar::from(8),
            ],
            matrixed_h_g2: vec![
                g2 * Scalar::from(2),
                g2 * Scalar::from(4),
                g2 * Scalar::from(6),
                g2 * Scalar::from(8),
            ],
        };
        let mut paramter_in_storage = initial_uncommon_paramters::<Bls12>(&common_paramter_matrix);
        //under-chain
        let player1_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (1, 2));
        //on-chain
        paramter_in_storage = verify_uncommon_paramter(
            &common_paramter_matrix,
            &paramter_in_storage,
            &player1_common,
        );
        /*assert_eq!(
            paramter_in_storage.kin_g1[0],
            common_paramter_matrix.matrixed_g1_front[0]
        );
        assert_eq!(
            paramter_in_storage.kin_g1[1],
            common_paramter_matrix.matrixed_g1_front[1]
        );*/
        assert_eq!(paramter_in_storage.gamma_g2, g2 * Scalar::from(1));
        //under-chain
        let player2_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (2, 3));
        /*assert_eq!(
            player2_common.ic.list[0].g1_result.unwrap(),
            common_paramter_matrix.matrixed_g1_front[0] * Scalar::from(2).invert().unwrap()
        );
        assert_eq!(
            player2_common.ic.list[1].g1_result.unwrap(),
            common_paramter_matrix.matrixed_g1_front[1] * Scalar::from(2).invert().unwrap()
        );*/
        assert_eq!(
            player2_common.delta.g2_result.unwrap(),
            g2 * Scalar::from(6)
        );
        //assert_eq!(paramter_in_storage.gamma_g2, g2 * Scalar::from(1));
        assert_eq!(player2_common.gamma.g2_mine.unwrap(), g2 * Scalar::from(2));

        //on-chain
        paramter_in_storage = verify_uncommon_paramter(
            &common_paramter_matrix,
            &paramter_in_storage,
            &player2_common,
        );
        //under-chain
        let player3_common = mpc_uncommon_paramters_generator(&paramter_in_storage, (3, 4));
        //on-chain
        paramter_in_storage = verify_uncommon_paramter(
            &common_paramter_matrix,
            &paramter_in_storage,
            &player3_common,
        );
        //under-chain

        assert_eq!(paramter_in_storage.gamma_g1, g1 * Scalar::from(6));
        assert_eq!(paramter_in_storage.gamma_g2, g2 * Scalar::from(6));
        assert_eq!(paramter_in_storage.delta_g1, g1 * Scalar::from(24));
        assert_eq!(paramter_in_storage.delta_g2, g2 * Scalar::from(24));
        let gamma_invert = Scalar::from(6).invert().unwrap();
        let delta_invert = Scalar::from(24).invert().unwrap();
        assert_eq!(
            paramter_in_storage.kin_g1[0],
            common_paramter_matrix.matrixed_g1_front[0] * gamma_invert
        );
        assert_eq!(
            paramter_in_storage.kin_g1[1],
            common_paramter_matrix.matrixed_g1_front[1] * gamma_invert
        );

        assert_eq!(
            paramter_in_storage.kout_g1[0],
            common_paramter_matrix.matrixed_g1_back[0] * delta_invert
        );
        assert_eq!(
            paramter_in_storage.kout_g1[1],
            common_paramter_matrix.matrixed_g1_back[1] * delta_invert
        );

        assert_eq!(
            paramter_in_storage.h_g1[0],
            common_paramter_matrix.matrixed_h_g1[0] * delta_invert
        );
        assert_eq!(
            paramter_in_storage.h_g1[1],
            common_paramter_matrix.matrixed_h_g1[1] * delta_invert
        );
        assert_eq!(
            paramter_in_storage.h_g1[2],
            common_paramter_matrix.matrixed_h_g1[2] * delta_invert
        );
        assert_eq!(
            paramter_in_storage.h_g1[3],
            common_paramter_matrix.matrixed_h_g1[3] * delta_invert
        );

        //on-chain
        //paramter_in_storage
    }
}
