use bls12_381::{
    multi_miller_loop, pairing, G1Affine, G1Projective, G2Affine, G2Prepared, G2Projective, Gt,
    MillerLoopResult, Scalar,
};
use std::mem::transmute;

use std::convert::AsMut;
fn clone_into_array<A, T>(slice: &[T]) -> A
where
    A: Default + AsMut<[T]>,
    T: Clone,
{
    let mut a = A::default();
    <A as AsMut<[T]>>::as_mut(&mut a).clone_from_slice(slice);
    a
}

use crate::gt_utils::{adc, mac, sbb};
/// p = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
const MODULUS: [u64; 6] = [
    0xb9fe_ffff_ffff_aaab,
    0x1eab_fffe_b153_ffff,
    0x6730_d2a0_f6b0_f624,
    0x6477_4b84_f385_12bf,
    0x4b1b_a7b6_434b_acd7,
    0x1a01_11ea_397f_e69a,
];

/// INV = -(p^{-1} mod 2^64) mod 2^64
const INV: u64 = 0x89f3_fffc_fffc_fffd;

pub fn gt_format(tmp: &[u64; 72]) -> [u8; 576] {
    let c0 = fp6_to_bytes_format(&tmp[0..36]);
    let c1 = fp6_to_bytes_format(&tmp[36..72]);
    let mut res = [0u8; 576];
    res[0..288].copy_from_slice(&c1);
    res[288..576].copy_from_slice(&c0);
    res
}

pub fn fp2_to_bytes_format(tmp: &[u64]) -> [u8; 96] {
    let c0 = fp_to_bytes_format(&tmp[0..6]);
    let c1 = fp_to_bytes_format(&tmp[6..12]);
    let mut res = [0u8; 96];
    res[0..48].copy_from_slice(&c1);
    res[48..96].copy_from_slice(&c0);
    res
}

pub fn fp6_to_bytes_format(tmp: &[u64]) -> [u8; 288] {
    let c0 = fp2_to_bytes_format(&tmp[0..12]);
    let c1 = fp2_to_bytes_format(&tmp[12..24]);
    let c2 = fp2_to_bytes_format(&tmp[24..36]);
    let mut res = [0u8; 288];
    res[0..96].copy_from_slice(&c2);
    res[96..192].copy_from_slice(&c1);
    res[192..288].copy_from_slice(&c0);
    res
}

pub fn fp_to_bytes_format(tmp: &[u64]) -> [u8; 48] {
    let tmp = montgomery_reduce(
        tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], 0, 0, 0, 0, 0, 0,
    );
    let mut res = [0; 48];
    res[0..8].copy_from_slice(&tmp[5].to_be_bytes());
    res[8..16].copy_from_slice(&tmp[4].to_be_bytes());
    res[16..24].copy_from_slice(&tmp[3].to_be_bytes());
    res[24..32].copy_from_slice(&tmp[2].to_be_bytes());
    res[32..40].copy_from_slice(&tmp[1].to_be_bytes());
    res[40..48].copy_from_slice(&tmp[0].to_be_bytes());

    res
}

pub(crate) const fn montgomery_reduce(
    t0: u64,
    t1: u64,
    t2: u64,
    t3: u64,
    t4: u64,
    t5: u64,
    t6: u64,
    t7: u64,
    t8: u64,
    t9: u64,
    t10: u64,
    t11: u64,
) -> [u64; 6] {
    // The Montgomery reduction here is based on Algorithm 14.32 in
    // Handbook of Applied Cryptography
    // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

    let k = t0.wrapping_mul(INV);
    let (_, carry) = mac(t0, k, MODULUS[0], 0);
    let (r1, carry) = mac(t1, k, MODULUS[1], carry);
    let (r2, carry) = mac(t2, k, MODULUS[2], carry);
    let (r3, carry) = mac(t3, k, MODULUS[3], carry);
    let (r4, carry) = mac(t4, k, MODULUS[4], carry);
    let (r5, carry) = mac(t5, k, MODULUS[5], carry);
    let (r6, r7) = adc(t6, 0, carry);

    let k = r1.wrapping_mul(INV);
    let (_, carry) = mac(r1, k, MODULUS[0], 0);
    let (r2, carry) = mac(r2, k, MODULUS[1], carry);
    let (r3, carry) = mac(r3, k, MODULUS[2], carry);
    let (r4, carry) = mac(r4, k, MODULUS[3], carry);
    let (r5, carry) = mac(r5, k, MODULUS[4], carry);
    let (r6, carry) = mac(r6, k, MODULUS[5], carry);
    let (r7, r8) = adc(t7, r7, carry);

    let k = r2.wrapping_mul(INV);
    let (_, carry) = mac(r2, k, MODULUS[0], 0);
    let (r3, carry) = mac(r3, k, MODULUS[1], carry);
    let (r4, carry) = mac(r4, k, MODULUS[2], carry);
    let (r5, carry) = mac(r5, k, MODULUS[3], carry);
    let (r6, carry) = mac(r6, k, MODULUS[4], carry);
    let (r7, carry) = mac(r7, k, MODULUS[5], carry);
    let (r8, r9) = adc(t8, r8, carry);

    let k = r3.wrapping_mul(INV);
    let (_, carry) = mac(r3, k, MODULUS[0], 0);
    let (r4, carry) = mac(r4, k, MODULUS[1], carry);
    let (r5, carry) = mac(r5, k, MODULUS[2], carry);
    let (r6, carry) = mac(r6, k, MODULUS[3], carry);
    let (r7, carry) = mac(r7, k, MODULUS[4], carry);
    let (r8, carry) = mac(r8, k, MODULUS[5], carry);
    let (r9, r10) = adc(t9, r9, carry);

    let k = r4.wrapping_mul(INV);
    let (_, carry) = mac(r4, k, MODULUS[0], 0);
    let (r5, carry) = mac(r5, k, MODULUS[1], carry);
    let (r6, carry) = mac(r6, k, MODULUS[2], carry);
    let (r7, carry) = mac(r7, k, MODULUS[3], carry);
    let (r8, carry) = mac(r8, k, MODULUS[4], carry);
    let (r9, carry) = mac(r9, k, MODULUS[5], carry);
    let (r10, r11) = adc(t10, r10, carry);

    let k = r5.wrapping_mul(INV);
    let (_, carry) = mac(r5, k, MODULUS[0], 0);
    let (r6, carry) = mac(r6, k, MODULUS[1], carry);
    let (r7, carry) = mac(r7, k, MODULUS[2], carry);
    let (r8, carry) = mac(r8, k, MODULUS[3], carry);
    let (r9, carry) = mac(r9, k, MODULUS[4], carry);
    let (r10, carry) = mac(r10, k, MODULUS[5], carry);
    let (r11, _) = adc(t11, r11, carry);

    // Attempt to subtract the modulus, to ensure the value
    // is smaller than the modulus.
    subtract_p(&[r6, r7, r8, r9, r10, r11])
}

#[inline]
const fn subtract_p(tmp: &[u64; 6]) -> [u64; 6] {
    let (r0, borrow) = sbb(tmp[0], MODULUS[0], 0);
    let (r1, borrow) = sbb(tmp[1], MODULUS[1], borrow);
    let (r2, borrow) = sbb(tmp[2], MODULUS[2], borrow);
    let (r3, borrow) = sbb(tmp[3], MODULUS[3], borrow);
    let (r4, borrow) = sbb(tmp[4], MODULUS[4], borrow);
    let (r5, borrow) = sbb(tmp[5], MODULUS[5], borrow);

    // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
    // borrow = 0x000...000. Thus, we use it as a mask!
    let r0 = (tmp[0] & borrow) | (r0 & !borrow);
    let r1 = (tmp[1] & borrow) | (r1 & !borrow);
    let r2 = (tmp[2] & borrow) | (r2 & !borrow);
    let r3 = (tmp[3] & borrow) | (r3 & !borrow);
    let r4 = (tmp[4] & borrow) | (r4 & !borrow);
    let r5 = (tmp[5] & borrow) | (r5 & !borrow);

    [r0, r1, r2, r3, r4, r5]
}

#[test]
fn test_format() {
    let bytes: [u8; 576] = [
        197, 133, 31, 160, 51, 228, 114, 25, 56, 37, 119, 253, 118, 43, 211, 151, 249, 205, 107,
        201, 111, 84, 206, 200, 20, 6, 212, 102, 115, 62, 246, 206, 128, 55, 132, 129, 39, 52, 17,
        166, 37, 216, 198, 63, 138, 68, 243, 19, 149, 105, 157, 46, 176, 49, 99, 210, 125, 126,
        121, 247, 130, 164, 104, 157, 146, 234, 57, 141, 36, 41, 155, 156, 170, 7, 49, 225, 162,
        28, 128, 244, 102, 176, 188, 189, 50, 7, 108, 161, 120, 4, 54, 186, 175, 164, 60, 8, 65,
        182, 22, 9, 219, 97, 226, 89, 13, 150, 62, 178, 244, 182, 22, 39, 69, 156, 189, 160, 16,
        91, 229, 200, 168, 237, 77, 156, 217, 11, 219, 11, 197, 170, 253, 87, 191, 158, 248, 140,
        94, 122, 119, 158, 146, 183, 214, 18, 53, 95, 225, 176, 136, 81, 200, 95, 101, 99, 9, 143,
        58, 110, 160, 52, 44, 214, 42, 224, 166, 38, 49, 219, 11, 153, 154, 125, 169, 90, 111, 252,
        16, 194, 137, 235, 245, 85, 47, 161, 137, 136, 111, 146, 58, 112, 35, 23, 120, 135, 130,
        113, 41, 143, 88, 147, 133, 117, 171, 17, 134, 91, 246, 67, 223, 159, 39, 236, 245, 170,
        131, 49, 246, 157, 201, 138, 225, 215, 115, 250, 176, 153, 76, 166, 166, 118, 225, 100, 31,
        143, 56, 88, 140, 167, 159, 23, 18, 239, 42, 202, 17, 10, 42, 103, 107, 241, 163, 42, 181,
        185, 17, 13, 110, 5, 157, 105, 208, 18, 68, 164, 165, 91, 26, 34, 119, 1, 29, 192, 41, 85,
        115, 108, 222, 206, 224, 102, 57, 195, 221, 159, 30, 167, 245, 5, 121, 198, 98, 176, 161,
        136, 10, 211, 4, 131, 252, 53, 93, 106, 197, 90, 13, 41, 31, 168, 166, 52, 200, 208, 199,
        7, 55, 218, 194, 48, 84, 205, 240, 10, 80, 128, 247, 127, 194, 240, 174, 46, 215, 226, 166,
        93, 36, 9, 86, 81, 27, 121, 118, 6, 46, 159, 19, 254, 24, 73, 35, 200, 209, 226, 244, 27,
        86, 60, 159, 69, 158, 76, 193, 227, 211, 185, 83, 94, 232, 163, 32, 0, 167, 33, 30, 18, 10,
        130, 204, 154, 197, 65, 131, 97, 175, 21, 177, 58, 153, 36, 140, 101, 149, 124, 185, 134,
        168, 28, 114, 56, 235, 115, 188, 52, 116, 71, 73, 215, 86, 82, 139, 74, 80, 234, 2, 25,
        164, 139, 109, 206, 134, 12, 248, 211, 163, 4, 170, 110, 104, 251, 135, 74, 166, 24, 38,
        207, 32, 185, 27, 231, 131, 187, 69, 57, 167, 146, 172, 119, 82, 42, 160, 70, 240, 148,
        159, 229, 14, 252, 247, 88, 96, 120, 243, 205, 88, 113, 246, 69, 249, 130, 27, 6, 193, 124,
        103, 229, 219, 159, 170, 71, 248, 3, 87, 230, 52, 97, 165, 219, 120, 128, 110, 138, 153,
        67, 154, 236, 215, 28, 102, 55, 153, 26, 154, 89, 170, 177, 68, 238, 66, 8, 47, 246, 160,
        201, 250, 223, 5, 182, 227, 155, 21, 142, 194, 63, 241, 74, 13, 186, 134, 12, 177, 255, 82,
        106, 160, 242, 15, 232, 108, 144, 26, 114, 72, 202, 148, 118, 20, 133, 176, 3, 62, 24, 131,
        117, 226, 228, 206, 64, 221, 175, 103, 245, 252, 165, 38, 229, 210, 150, 109, 154, 66, 34,
        31, 134, 73, 159, 126, 25,
    ];
    let g1 = G1Affine::generator();
    let g2 = G2Affine::generator();
    let gt = pairing(&g1, &g2);

    let tmp_gt = unsafe { transmute::<[u8; 576], [u64; 72]>(bytes) };
    let bytes = gt_format(&tmp_gt);
    println!("gt:");
    for i in bytes.iter() {
        print!("{},", i);
    }
    println!("");
    println!("gt+gt:");

    let gtplusgt = gt + gt;
    let tmp_gt = unsafe { transmute::<Gt, [u64; 72]>(gtplusgt) };
    let bytes = gt_format(&tmp_gt);
    for i in bytes.iter() {
        print!("{},", i);
    }
    println!("");

    println!("gt*3:");

    let gtmul3 = gt * Scalar::from(3);
    let tmp_gt = unsafe { transmute::<Gt, [u64; 72]>(gtmul3) };
    let bytes = gt_format(&tmp_gt);
    for i in bytes.iter() {
        print!("{},", i);
    }
    println!("");
    println!("gt*-3:");

    let tmp_gt = unsafe { transmute::<Gt, [u64; 72]>(-gtmul3) };
    let bytes = gt_format(&tmp_gt);
    for i in bytes.iter() {
        print!("{},", i);
    }
    println!("");
}
