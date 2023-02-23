#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"
#include "two_walks.h"
#include "curve.h"
#include "printing.h"
int main() {
  srand48(1);

  two_walk phi;
  init_curve(&phi.A);
  proj A, B, C, P, Pbak, Pbbak;

  for (int i = 0; i < 100; i++) {
    // sample point on curve
    while (true) {
      fp2_random(&phi.ker.x); phi.ker.z = fp2_1;
      if (!is_on_curve(&phi.ker, &phi.A))
	continue;
      // multiply by cofactor
      xMUL(&phi.ker, &phi.A, &phi.ker, &p_even_cofactor);
      // computer order
      P = phi.ker;
      for (phi.len = 0; !mont_iszero(&P) && !fp2_iszero(&P.x); phi.len++)
	xDBL(&P, &phi.A, &P);
      if (mont_iszero(&P) && phi.len > 0)
	break;
    }
    
    // Sample point
    fp2_random(&P.x); P.z = fp2_1;
    bool oncurve = is_on_curve(&P, &phi.A);
    
    // Evaluate isogeny
    A = phi.A;
    Pbak = P;
    eval_walk(&phi, &B, &P);
    Pbbak = P;
    assert(is_on_curve(&P, &B) == oncurve);
    
    // Evaluate dual
    eval_dual(&phi, &B, &P);
    // Check dual isogeny equation
    for (int i = 0; i < phi.len; i++)
      xDBL(&Pbak, &phi.A, &Pbak);
    assert(mont_equal(&Pbak, &P));

    // Double check dual
    isomorphism isom;
    proj j1, j2;
    // construct dual
    dual_walk(&phi);
    jinv256(&j1, &B);
    jinv256(&j2, &phi.A);
    assert(mont_equal(&j1, &j2));
    // move Pbbak to new curve
    mont_isom(&isom, &B, &phi.A);
    mont_isom_apply(&isom, &Pbbak);
    assert(is_on_curve(&Pbbak, &phi.A) == oncurve);
    // evaluate dual
    eval_walk(&phi, &C, &Pbbak);
    assert(is_on_curve(&Pbbak, &C) == oncurve);
    // check consistency
    jinv256(&j1, &C);
    jinv256(&j2, &A);    
    assert(mont_equal(&j1, &j2));
    mont_isom(&isom, &C, &A);
    mont_isom_apply(&isom, &Pbbak);
    if (fp2_iszero(&j1.x) && !mont_equal(&Pbbak, &P))
      mont0_dist(&Pbbak, &Pbbak);
    if (fp2_iszero(&j1.x) && !mont_equal(&Pbbak, &P))
      mont0_dist(&Pbbak, &Pbbak);
    if (fp2_iszero(&A.x) && !mont_equal(&Pbbak, &P))
      mont0_dist(&Pbbak, &Pbbak);
    assert(mont_equal(&Pbbak, &P));
  }
  
  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
