#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "fp.h"
#include "fp2.h"

int main(){
  srand48(1);
  fp a, b, c, d, e, f;

  // Testing distributivity + commutativity
  for (int i = 0; i < 100; i++) {
    fp_random(&a);
    fp_random(&b);
    fp_random(&c);

    fp_add3(&d, &b, &c);
    fp_mul2(&d, &a);
    fp_mul3(&e, &b, &a);
    fp_mul3(&f, &c, &a);
    fp_add2(&e, &f);
    fp_sub2(&d, &e);

    assert(fp_iszero(&d));

    fp_neg2(&b, &a);
    fp_add2(&a, &b);
    assert(fp_iszero(&a));
  }

  // Testing inverse
  for (int i = 0; i < 100; i++) {
    fp_random(&a);
    b = a;
    fp_inv(&a);
    fp_mul3(&c, &a, &b);

    assert(c.x.c[0] == fp_1.x.c[0]);
    assert(c.x.c[1] == fp_1.x.c[1]);
    assert(c.x.c[2] == fp_1.x.c[2]);
    assert(c.x.c[3] == fp_1.x.c[3]);
    assert(fp_issquare(&a) == fp_issquare(&b));
  }

  /* Same for GF(p^2) */
  fp2 A, B, C, D, E, F;

  // Testing distributivity + commutativity
  for (int i = 0; i < 100; i++) {
    fp2_random(&A);
    fp2_random(&B);
    fp2_random(&C);

    fp2_add3(&D, &B, &C);
    fp2_mul2(&D, &A);
    fp2_mul3(&E, &B, &A);
    fp2_mul3(&F, &C, &A);
    fp2_add2(&E, &F);
    fp2_sub2(&D, &E);

    assert(fp2_iszero(&D));
    
    fp2_neg2(&B, &A);
    fp2_add2(&A, &B);
    assert(fp2_iszero(&A));
  }

  // Testing inverse
  for (int i = 0; i < 100; i++) {
    fp2_random(&A);
    B = A;
    fp2_inv(&A);
    fp2_mul3(&C, &A, &B);

    assert(C.re.x.c[0] == fp_1.x.c[0]);
    assert(C.re.x.c[1] == fp_1.x.c[1]);
    assert(C.re.x.c[2] == fp_1.x.c[2]);
    assert(C.re.x.c[3] == fp_1.x.c[3]);
    assert(fp_iszero(&C.im));
  }

  // Testing square root
  for (int i = 0; i < 100; i++) {
    fp2_random(&A);
    fp2_sq2(&B, &A);
    A = B;
    fp2_sqrt(&B);
    fp2_sq1(&B);

    assert(fp2_equal(&A, &B));
  }

  // Testing new GF(p) multiplication
  for (int i = 0; i < 100; i++) {
      fp_random(&a);
      fp_random(&b);
      fp_mul3(&c, &a, &b);
      fp_mul3_old(&d, &a, &b);
      fp_sub2(&d, &c);

      assert(fp_iszero(&d));
  }

  // Testing new GF(p^2) multiplication
  for (int i = 0; i < 100; i++) {
      fp2_random(&A);
      fp2_random(&B);
      fp2_mul3(&C, &A, &B);
      fp2_mul3_old(&D, &A, &B);
      fp2_sub2(&D, &C);

      assert(fp2_iszero(&D));
  }

  // Testing new GF(p^2) squaring
  for (int i = 0; i < 100; i++) {
      fp2_random(&A);
      fp2_sq2(&C, &A);
      fp2_sq2_old(&D, &A);
      fp2_sub2(&D, &C);

      assert(fp2_iszero(&D));
  }
  
  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}

