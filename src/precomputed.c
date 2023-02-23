#include <assert.h>
#include "precomputed.h"
#include "curve.h"
#include "uintbig.h"

struct precomp_struct global_setup;

void init_precomputations() {
  init_curve(&global_setup.E0);
  init_precomputations_generated();

  long var = fetch_var();
  GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC);
  GEN p = strtoi(p_str);
  /*   GEN B = alg_hilbert(nf, stoi(-1), negi(p), var, 0); */

  long q = q_norm;
  // long c = -1;
  GEN multtable = mkvec4(
    mkmat4(mkcol4s(1,0,0,0), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), mkcol4s(0,0,0,1)),
    mkmat4(mkcol4s(0,1,0,0), mkcol4s(-q,0,0,0), mkcol4s(0,0,0,-1), mkcol4s(0,0,q,0)),
    mkmat4(mkcol4s(0,0,1,0), mkcol4s(0,0,0,1), mkcol4(gneg(p),gen_0,gen_0,gen_0), mkcol4(gen_0,gneg(p),gen_0,gen_0)),
    mkmat4(mkcol4s(0,0,0,1), mkcol4s(0,0,-q,0), mkcol4(gen_0,p,gen_0,gen_0), mkcol4(gneg(gmul(stoi(q),p)),gen_0,gen_0,gen_0))
    );

  GEN B = alg_csa_table(nf, multtable, var,0);

  GEN B_1 = mkcol4s(1,0,0,0);
  GEN B_i = mkcol4s(0,1,0,0);
  GEN B_j = mkcol4s(0,0,1,0);
  GEN B_ji = mkcol4s(0,0,0,1);


  /*   GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2 */
  /*   GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2 */

  global_setup.p = p;
  global_setup.q = q;
  global_setup.B = B; // the quaternion algebra
  global_setup.qf = mkmat4(mkcol4s(1,0,0,0),
			   mkcol4s(0,q,0,0),
			   mkcol4(gen_0,gen_0,p,gen_0),
			   mkcol4(gen_0,gen_0,gen_0,gmul(p,stoi(q)))); // quadratic form defined by the reduced norm

  // global_setup.torsion_fm = Z_factor_limit(strtoi(
  //    all_the_torsion_str
  //    ), 70000000);

  global_setup.gen_p_plus_fact = cgetg(3, t_MAT);
  gel(global_setup.gen_p_plus_fact,1) = cgetg(p_plus_len+1, t_COL);
  gel(global_setup.gen_p_plus_fact,2) = cgetg(p_plus_len+1, t_COL);
  global_setup.gen_p_minus_fact = cgetg(3, t_MAT);
  gel(global_setup.gen_p_minus_fact,1) = cgetg(p_minus_len+1, t_COL);
  gel(global_setup.gen_p_minus_fact,2) = cgetg(p_minus_len+1, t_COL);

  global_setup.gen_p_plus_primary = cgetg(p_plus_len+1, t_VEC);
  global_setup.gen_p_minus_primary = cgetg(p_minus_len+1, t_VEC);

  for (int i = 0; i < p_plus_len; ++i) {
    gel(gel(global_setup.gen_p_plus_fact,1),i+1) = stoi(p_plus_fact[i]);
    gel(gel(global_setup.gen_p_plus_fact,2),i+1) = stoi(p_plus_mult[i]);
    gel(global_setup.gen_p_plus_primary,i+1) = powuu(p_plus_fact[i],p_plus_mult[i]);
  }

  for (int i = 0; i < p_minus_len; ++i) {
    gel(gel(global_setup.gen_p_minus_fact,1),i+1) = stoi(p_minus_fact[i]);
    gel(gel(global_setup.gen_p_minus_fact,2),i+1) = stoi(p_minus_mult[i]);
    gel(global_setup.gen_p_minus_primary,i+1) = powuu(p_minus_fact[i],p_minus_mult[i]);
  }

  global_setup.gen_odd_torsion = gmul(ZV_prod(global_setup.gen_p_plus_primary),
				      ZV_prod(global_setup.gen_p_minus_primary));

  // GEN B1 = mkcol4(gen_1,gen_0,gen_0,gen_0);
  // GEN B2 = mkcol4(ghalf,ghalf,gen_0,gen_0);
  // GEN B3 = mkcol4(gen_0,gen_0,ghalf,gneg(ghalf));
  // GEN B4 = mkcol4(gen_0,gdiv(gen_1,stoi(q)),gen_0,gneg(gdiv(stoi(c),stoi(q))));

  GEN B1 = mkcol4(gen_1,gen_0,gen_0,gen_0);
  GEN B2 = mkcol4(gen_0,gen_1,gen_0,gen_0);
  GEN B3 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
  GEN B4 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2

  global_setup.O0 = alglathnf(B,mkmat4(B1,B2,B3,B4), gen_0); // the cannonical maximal order
  /*   global_setup.O0 = alglathnf(B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0); */


  global_setup.one = B_1;
  global_setup.i = B_i;
  global_setup.j = B_j;
  global_setup.ji = B_ji;

  global_setup.O0_b1 = B1;
  global_setup.O0_b2 = B2;
  global_setup.O0_b3 = B3;
  global_setup.O0_b4 = B4;
  global_setup.O0_to_standard = mkmat4(B1, B2, B3, B4);
  global_setup.standard_to_O0 = RgM_inv(global_setup.O0_to_standard);

  global_setup.smooth_famat_for_klpt =Z_factor_limit(global_setup.gen_odd_torsion, 7045009+1);

  for (int i = 0; i < EXTREMAL_ORDERS_N; ++i) {
    global_setup.orders[i] = NULL;
  }

  long r = 3;
  while (r < EXTREMAL_ORDERS_N) {
      if (kronecker(gneg(stoi(r)),global_setup.p) == -1)  {
          global_setup.orders[r] = malloc(sizeof(p_extremal_maximal_order));
          *(global_setup.orders[r]) = get_p_extremal_maximal_order_from_scratch(r);
      }
      r = itos(nextprime(stoi(r+1)));
  }

  if (q < EXTREMAL_ORDERS_N) {
      free(global_setup.orders[q]);
      global_setup.orders[q] = malloc(sizeof(p_extremal_maximal_order));
      (*(global_setup.orders[q])).order = global_setup.O0;
      (*(global_setup.orders[q])).i = global_setup.i;
      (*(global_setup.orders[q])).j = global_setup.j;
      (*(global_setup.orders[q])).q = global_setup.q;
  }
  // printf("DONE\n");

  // output(global_setup.smooth_famat_for_klpt);

  /*   // output(alg_O0_to_standard(mkcol4s(1,0,0,0))); */
  /*   // output(alg_O0_to_standard(mkcol4s(0,1,0,0))); */
  /*   // output(alg_O0_to_standard(mkcol4s(0,0,1,0))); */
  /*   // output(alg_O0_to_standard(mkcol4s(0,0,0,1))); */

}
