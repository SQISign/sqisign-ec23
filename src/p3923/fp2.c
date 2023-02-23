#include "fp2.h"

const fp2 fp2_0 = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };

/* Arithmetic modulo X^2 + 1 */

// 5 + 4i is a quadratic non-residue
fp2 fp2_non_residue() {
  fp2 res;
  fp_add3(&res.im, &fp_1, &fp_1);
  fp_add2(&res.im, &res.im);
  fp_add3(&res.re, &res.im, &fp_1);
  return res;
}

void fp2_mul3_old(fp2 *x, fp2 const *y, fp2 const *z) {
  fp xsum, xim;
  fp_add3(&xsum, &y->re, &y->im);
  fp_add3(&xim, &z->re, &z->im);
  fp_mul2(&xsum, &xim);
  fp_mul3(&xim, &y->im, &z->im);
  fp_mul3(&x->re, &y->re, &z->re);
  fp_sub3(&x->im, &xsum, &xim);
  fp_sub2(&x->im, &x->re);
  fp_sub2(&x->re, &xim);
}

void fp2_mul3(fp2 *x, fp2 const *y, fp2 const *z)
{
  fp t;

  fp2_mul_c0(&t, y, z);              // c0 = a0*b0 - a1*b1
  fp2_mul_c1(&x->im, y, z);          // c1 = a0*b1 + a1*b0 
  //x->re = t;
  x->re.x.c[0] = t.x.c[0]; x->re.x.c[1] = t.x.c[1]; x->re.x.c[2] = t.x.c[2]; x->re.x.c[3] = t.x.c[3];
}

void fp2_sq2_old(fp2 *x, fp2 const *y) {
  fp sum, diff;
  fp_add3(&sum, &y->re, &y->im);
  fp_sub3(&diff, &y->re, &y->im);
  fp_mul3(&x->im, &y->re, &y->im);
  fp_add2(&x->im, &x->im);
  fp_mul3(&x->re, &sum, &diff);
}

void fp2_sq2(fp2* x, fp2 const* y) {
    fp2 t;

    fp2_sq_c0(&t, y);               // c0 = (a0+a1)(a0-a1)
    fp2_sq_c1(&x->im, y);           // c1 = 2a0*a1
    //x->re = t;
    x->re.x.c[0] = t.re.x.c[0]; x->re.x.c[1] = t.re.x.c[1]; x->re.x.c[2] = t.re.x.c[2]; x->re.x.c[3] = t.re.x.c[3];
}

void fp2_inv(fp2 *x) {
  fp inorm, im2;
  fp_sq2(&inorm, &x->re);
  fp_sq2(&im2, &x->im);
  fp_add2(&inorm, &im2);
  fp_inv(&inorm);
  fp_mul2(&x->re, &inorm);
  fp_mul2(&x->im, &inorm);
  fp_neg1(&x->im);
}

bool fp2_issquare(const fp2 *x) {
  fp inorm, im2;
  fp_sq2(&inorm, &x->re);
  fp_sq2(&im2, &x->im);
  fp_add2(&inorm, &im2);
  return fp_issquare(&inorm);
}

void fp2_frob2(fp2 *x, const fp2 *y) {
  x->re = y->re;
  fp_neg2(&x->im, &y->im);
}

void fp2_exp(fp2 *res, fp2 const *x, uintbig const *k)
{
    if (fp2_iszero(x)) { *res = *x; return; }
    const fp2 xcopy = *x;
    *res = fp2_1;

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i));
    do {
        fp2_sq1(res);
        if (uintbig_bit(k, i)) {
            fp2_mul2(res, &xcopy);
        }
    } while (i--);
}



// dlp of h in basis g, which has order ell, naive implementation
bool fp2_dlp_naive(long *res, const fp2 *h, const fp2 *g, long ell) {
    long logarithm = 0;
    fp2 x = fp2_1;

    for (int i = 0; i < ell; ++i) {
        if (fp2_equal(h,&x)) { *res = logarithm; return true;}
        logarithm++;
        fp2_mul2(&x,g);
    }

    return false;
}


void fp2_sqrt(fp2 *x) {
    if (fp_iszero(&x->im)) {
        fp x_re_copy = x->re;

        if (fp_issquare(&x_re_copy)) {
            fp_sqrt(&x->re);
            return;
        }
        else {
            fp_neg2(&x->im, &x->re);
            fp_sqrt(&x->im);
            fp_set(&x->re, 0);
            return;
        }
    }

    fp sdelta, re, tmp1, tmp2, inv2, im;

    // sdelta = sqrt(re^2 + im^2)
    fp_sq2(&sdelta, &x->re);
    fp_sq2(&tmp1, &x->im);
    fp_add2(&sdelta, &tmp1);

    fp_sqrt(&sdelta);

    fp_set(&inv2,2);
    fp_inv(&inv2);

    fp_add3(&re,&x->re,&sdelta);
    fp_mul2(&re,&inv2);
    tmp2 = re;

    if (!fp_issquare(&tmp2)) {
        fp_sub3(&re,&x->re,&sdelta);
        fp_mul2(&re,&inv2);
    }

    fp_sqrt(&re);

    im = re;

    fp_inv(&im);
    fp_mul2(&im,&inv2);
    fp_mul2(&im,&x->im);

    x->re = re;
    x->im = im;
}

