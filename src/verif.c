#include "verif.h"
#include <assert.h>



void decompress_old(two_walk *walk, proj *A, const uintbig *zip, long len) {
  // long mask = (1 << two_tors_height) - 1;
  // long hint_mask = (0xf << two_tors_height);
  // uintbig a;
  proj P, Q, PQ;
  for (int i = 0; i < len; i++) {
    // uintbig_set(&a, zip[i] & mask);
    // long hint = (zip[i] & hint_mask) >> two_tors_height;
    // get the next kernel
    normalize_proj(A);
    find_basis_2e(&P, &Q, &PQ, A);  // TODO: use point Q from previous step + hint
    xBIDIM(&walk[i].ker, A, &P, &zip[i], &Q, &uintbig_1, &PQ);
    walk[i].A = *A;
    walk[i].len = two_tors_height;
    // take the next step
    eval_walk(walk+i, A, &Q);
  }
}

void challenge_alt(proj *A, const uintbig *m) {
  proj H, K, tmp;
  uintbig_set(&H.x.re.x, m->c[0]);
  uintbig_set(&H.x.im.x, m->c[1]);
  uintbig_set(&H.z.re.x, m->c[2]);
  uintbig_set(&H.z.im.x, m->c[3]);

  isog_degree deg = { 0 };
  degree_set(&deg, 0, p_plus_mult[0]);
  odd_isogeny phi;

  while(true) {
    fp_add2(&H.x.re, &fp_1);
    if (!is_on_curve(&H, A))
      continue;
    xMUL(&K, A, &H, &p_plus_odd_cofactor);
    isog_degree cof = degree_co(deg, p_plus_mult, p_plus_len);
    uintbig a;
    degree_to_uint(&a, cof, p_plus_fact, p_plus_len);
    xMUL(&K, A, &K, &a);
    tmp = K;
    uintbig_set(&a, p_plus_fact[0]);
    for (int i = 1; i < p_plus_mult[0]; i++) {
      xMUL(&tmp, A, &tmp, &a);
    }
    if (!mont_iszero(&tmp))
      break;
  }

  phi.kernel_plus = K;
  phi.kernel_minus.x = fp2_1;
  phi.kernel_minus.z = fp2_0;
  phi.deg_plus = deg;
  phi.deg_minus.val = 0;
  eval(A, &phi, &H);

  //
  degree_set(&deg, 0, p_minus_mult[0]);
  while(true) {
    fp_add2(&H.x.re, &fp_1);
    if (is_on_curve(&H, A))
      continue;
    xMUL(&K, A, &H, &p_minus_odd_cofactor);
    isog_degree cof = degree_co(deg, p_minus_mult, p_minus_len);
    uintbig a;
    degree_to_uint(&a, cof, p_minus_fact, p_minus_len);
    xMUL(&K, A, &K, &a);
    tmp = K;
    uintbig_set(&a, p_minus_fact[0]);
    for (int i = 1; i < p_minus_mult[0]; i++) {
      xMUL(&tmp, A, &tmp, &a);
    }
    if (!mont_iszero(&tmp))
      break;
  }

  phi.kernel_minus = K;
  phi.kernel_plus.x = fp2_1;
  phi.kernel_plus.z = fp2_0;
  phi.deg_minus = deg;
  phi.deg_plus.val = 0;
  eval(A, &phi, &H);
}
