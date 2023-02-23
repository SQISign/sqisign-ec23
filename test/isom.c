#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "isomorphism.h"
#include "isogenies.h"
#include "constants.h"
#include "curve.h"

static isog_degree mont_order(const proj *P, const proj *A, bool plus) {
  
    const long *fact, *mult;
    long len;
      

    if(plus) {
        fact = p_plus_fact; mult = p_plus_mult;
        len = p_plus_len;
    }
    else {
        fact = p_minus_fact; mult = p_minus_mult;
        len = p_minus_len;
    }

  proj tmp;
  uintbig cof;
  isog_degree deg = degree_co((isog_degree){ 0 }, mult, len);
  for (int j = 0; j < len; j++) {
    degree_unset(&deg, j);
    degree_to_uint(&cof, deg, fact, len);
    xMUL(&tmp, A, P, &cof);
    uintbig_set(&cof, fact[j]);
    uint8_t v = 0;
    for ( ; !mont_iszero(&tmp); v++) {
      xMUL(&tmp, A, &tmp, &cof);
      assert(v <= mult[j]);
    }
    degree_set(&deg, j, v);
  }
  return deg;
}

int main() {
  srand48(1);
    
  isomorphism isom1, isom2;
  proj A; init_curve(&A);
  proj B, P, Q, j1, j2;

  for (int i = 0; i < 10; i++) {
    fp2_random(&P.x); P.z = fp2_1;

    bool plus = is_p_plus_one_side(&P, &A);

    if (plus) {
      xMUL(&P, &A, &P, &p_plus_odd_cofactor);
    } else {
      xMUL(&P, &A, &P, &p_minus_odd_cofactor);
    }
    
    isog_degree deg = mont_order(&P, &A, plus);
    B = A;
    rand_isom(&isom1, &B);
    jinv256(&j1, &A);
    jinv256(&j2, &B);
    assert(mont_equal(&j1,&j2));

    Q = P;
    mont_isom_apply(&isom1, &Q);
    assert(is_p_plus_one_side(&Q, &B) == plus);
    assert(deg.val == mont_order(&Q, &B, plus).val);
    
    mont_isom(&isom2, &A, &B);
    mont_isom_apply(&isom2, &P);
    assert(is_p_plus_one_side(&P, &B) == plus);
    assert(deg.val == mont_order(&P, &B, plus).val);
  }

  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
