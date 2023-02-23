#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"
#include "curve.h"
#include "two_walks.h"

int main() {
  srand48(1);

  proj A; init_curve(&A);

  proj P, Q, PQ, T, TT;
  uintbig a, b;
  find_basis_2e(&P, &Q, &PQ, &A);
  for (int i = 0; i < 50; i++) {
    uintbig_random(&a);
    uintbig_random(&b);
    xBIDIM(&T, &A, &P, &a, &Q, &b, &PQ);
    assert(bidim_log_2e(&a, &b, &A, &T, &P, &Q, &PQ, two_tors_height));
    xBIDIM(&TT, &A, &P, &a, &Q, &b, &PQ);
    assert(mont_equal(&TT, &T));
  }

  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
