use bls12_381::Scalar;
use ff::PrimeField;

use crate::{Circuit, ConstraintSystem, LinearCombination, SynthesisError};

pub struct RangeDemo<Scalar: PrimeField> {
    pub a: Option<u64>,
    pub b: Option<u64>,
    pub n: Option<u64>,
    pub w:Option<u64>,
    pub wArray:Option<[u64;4]>,
    pub less_or_equal:Option<u64>,
    pub less:Option<u64>,
    pub not_all_zeros:Option<u64>,
    pub crArray:Option<[u64;4]>,
    pub tt: Option<Scalar>,
}

impl<Scalar: PrimeField> Circuit<Scalar> for RangeDemo<Scalar> {
    fn synthesize<CS: ConstraintSystem<Scalar>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        let a_var = Scalar::from(self.a.unwrap());
        let b_var = Scalar::from(self.b.unwrap());
        let n_var = Scalar::from(self.n.unwrap());
        let w_var = Scalar::from(self.w.unwrap());
        let mut wArray_var= vec![];
        for i in 0 .. self.wArray.unwrap().len(){
            let wArray_v = cs.alloc(||"",||Ok(Scalar::from(*self.wArray.unwrap().get(i).unwrap())));
            wArray_var.push(wArray_v.unwrap());
        }
        let mut crArray_var= vec![];
        for i in 0 .. self.crArray.unwrap().len() {
            let cArray_v = cs.alloc(||"",||Ok(Scalar::from(*self.crArray.unwrap().get(i).unwrap())));
            crArray_var.push(cArray_v.unwrap());
        }

        let not_all_zeros_var = Scalar::from(self.not_all_zeros.unwrap());
        let less_or_equal_var = Scalar::from(self.less_or_equal.unwrap());
        let less_var = Scalar::from(self.less.unwrap());

        let mut a = cs.alloc(|| "a", ||  Ok(a_var))?;
        let mut b = cs.alloc(|| "b", || Ok(b_var))?;

        let mut w = cs.alloc(|| "w", || Ok(w_var))?;
        let mut not_all_zeros = cs.alloc(|| "not_all_zeros", || Ok(not_all_zeros_var))?;
        let mut less_or_equal = cs.alloc(|| "less_or_equal", || Ok(less_or_equal_var))?;
        let mut less = cs.alloc(|| "less", || Ok(less_var))?;

        let t=1<<(self.n.unwrap()-1u64);
        cs.enforce(
            || "w=2^n+b-a",
            |lc| lc + w,
            |lc| lc + CS::one(),
            |lc| lc+(Scalar::from(t),CS::one())+b-a,
        );

        let mut lc1 = LinearCombination::<Scalar>::zero();
        for i in 0..wArray_var.len(){
            lc1 = lc1 + (Scalar::from(1<<i),wArray_var[i]);
        }
        lc1=lc1-w;
        cs.enforce(
            || "2^0*w0+.......-w=0",
            |_| lc1,
            |lc| lc + CS::one(),
            |lc| lc,
        );

        for i in 0..wArray_var.len() {
            cs.enforce(
                || "w0(1-w0)=0",
                |lc| lc + wArray_var[i],
                |lc| lc + CS::one()-wArray_var[i],
                |lc| lc,
            );
        }

        cs.enforce(
            || "w0=cr0",
            |lc| lc + wArray_var[0],
            |lc| lc + CS::one(),
            |lc| lc+crArray_var[0],
        );

        for i in 1..crArray_var.len() {
            cs.enforce(
                || "(cr_(i-1)-1)(wi-1)=1-cr_i",
                |lc| lc + crArray_var[i-1]-CS::one(),
                |lc| lc + wArray_var[i]-CS::one(),
                |lc| lc+CS::one()- crArray_var[i],
            );
        }

        cs.enforce(
            || "not_all_zeros=cr_n",
            |lc| lc + not_all_zeros,
            |lc| lc + CS::one(),
            |lc| lc+crArray_var[crArray_var.len()-1],
        );

        cs.enforce(
            || "wn=less_or_equal*wn",
            |lc| lc + wArray_var[wArray_var.len()-1],
            |lc| lc + less_or_equal,
            |lc| lc+wArray_var[wArray_var.len()-1],
        );

        cs.enforce(
            || "wn*less_or_equal=less",
            |lc| lc + wArray_var[wArray_var.len()-1],
            |lc| lc + not_all_zeros,
            |lc| lc+less,
        );
        Ok(())
    }
}

/*impl<Scalar: PrimeField> Circuit<Scalar> for RangeDemo<Scalar> {
    fn synthesize<CS: ConstraintSystem<Scalar>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        let a_var = Scalar::from(self.a.unwrap());
        let b_var = Scalar::from(self.b.unwrap());
        let n_var = Scalar::from(self.n.unwrap());

        let mut a = cs.alloc(|| "a", ||  Ok(a_var))?;
        let mut b = cs.alloc(|| "b", || Ok(b_var))?;
        let mut n = cs.alloc_input(|| "n", || Ok(n_var))?;

        let mut lc = LinearCombination::<Scalar>::zero();
        lc = lc + (Scalar::from(3u64),a)+b;
        cs.enforce(
            || "n=3a+b",
            |_| lc,
            |lc| lc+CS::one(),
            |lc| lc+n,
        );
        Ok(())
    }*/