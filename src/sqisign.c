#define _XOPEN_SOURCE
#define _unused(x) ((void)(x))

#include <stdint.h>
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "tedwards.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "sqisign.h"
#include "printing.h"
// #include "isomorphism.h"
// #include "uintbig.h"


void init_compressed_sig(compressed_signature *comp_sigma) {
  comp_sigma->zip=malloc(4*sizeof(uint64_t)*signing_length_two_tors_height_step);
}

void free_compressed_sig(compressed_signature *comp_sigma) {
  free(comp_sigma->zip);
}

//if two tors_height is too small this can fail if there is not enough torsion (more than usual...)
void keygen(public_key *pk, secret_key *sk) {

    unsigned int length_NI = security_level/2;
    //it is important that this integer is a multiple of two_tors_height
    unsigned int e_tau = (double)  ( ((security_level*2) / two_tors_height)+1)*two_tors_height;
    GEN NI = NULL;
    GEN NJ = powiu(gen_2, e_tau);

    do {
        NI = randomprime(powiu(gen_2, length_NI));
    } while (Fp_issquare(gen_2,NI));

    GEN gamma = NULL;
    while (!gamma) {
        // parity option is 1 so the 2-walk is not backtracking
        //we use this function because we need gamma to have integer coefficients
        gamma = norm_equation_special(global_setup.p, gmul(NI,NJ), 1, true);

    }

    gamma = gtrans(gamma);

    sk->I_large_prime = lideal_create(global_setup.B, global_setup.O0, gamma, NI);
    sk->I_two = lideal_create(global_setup.B, global_setup.O0, alg_conj(global_setup.B, gamma), NJ);

    init_trivial_two_walk_long(&sk->phi_two);
    pk-> E =global_setup.E0;
    proj Q;
    ideal_to_isogeny_1728(&sk->phi_two,sk->I_two,&Q);
    eval_walk(&sk->phi_two.phi[sk->phi_two.len-1],&pk->E,&Q);


    //old method keygen below, might be useful
    // ideal_to_isogeny_O0_two_long(&sk->phi_two, &sk->I_T, &sk->phi_T, sk->I_two, true);



    // if (!sk->phi_T.phi2_set) {
    //     assert(sk->phi_T.phi2_dual_set);
    //     sk->phi_T.phi2 = sk->phi_T.phi2_dual;
    //     sk->phi_T.middle = sk->phi_T.target;
    //     dual(&sk->phi_T.middle, &sk->phi_T.phi2);
    //     sk->phi_T.phi2_set = true;
    // }

    // if (!sk->phi_T.phi2_dual_set) {
    //     assert(sk->phi_T.phi2_set);
    //     sk->phi_T.phi2_dual = sk->phi_T.phi2;
    //     sk->phi_T.target = sk->phi_T.middle;
    //     dual(&sk->phi_T.target, &sk->phi_T.phi2_dual);
    //     sk->phi_T.phi2_dual_set = true;
    // }

    // pk->E = sk->phi_T.target;

}

void commitment(GEN *coeff, GEN *I, GEN *I_even, odd_isogeny *phi_com,two_walk_long *phi_com_even){
    pari_sp ltop = avma;

    GEN coeff_plus_1 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_plus_2 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_minus_1 = cgetg(p_minus_len+1, t_VEC);
    GEN coeff_minus_2 = cgetg(p_minus_len+1, t_VEC);

    for (int i = 0; i < p_plus_len; ++i) {
        gel(coeff_plus_1,i+1) = powuu(p_plus_fact[i], p_plus_mult[i] - p_plus_mult_com[i]);
        gel(coeff_plus_2,i+1) = randomi(powuu(p_plus_fact[i],p_plus_mult_com[i]));
        gel(coeff_plus_2,i+1) = gmul(gel(coeff_plus_2,i+1), gel(coeff_plus_1,i+1));

        gel(coeff_plus_1,i+1) = gmod(gel(coeff_plus_1,i+1), powuu(p_plus_fact[i],p_plus_mult[i]));
        gel(coeff_plus_2,i+1) = gmod(gel(coeff_plus_2,i+1), powuu(p_plus_fact[i],p_plus_mult[i]));
    }

    for (int i = 0; i < p_minus_len; ++i) {
        gel(coeff_minus_1,i+1) = powuu(p_minus_fact[i], p_minus_mult[i] - p_minus_mult_com[i]);;
        gel(coeff_minus_2,i+1) = randomi(powuu(p_minus_fact[i],p_minus_mult_com[i]));
        gel(coeff_minus_2,i+1) = gmul(gel(coeff_minus_2,i+1), gel(coeff_minus_1,i+1));

        gel(coeff_minus_1,i+1) = gmod(gel(coeff_minus_1,i+1), powuu(p_minus_fact[i],p_minus_mult[i]));
        gel(coeff_minus_2,i+1) = gmod(gel(coeff_minus_2,i+1), powuu(p_minus_fact[i],p_minus_mult[i]));
    }

    GEN coeff_plus = mkvec2(coeff_plus_1, coeff_plus_2);
    GEN coeff_minus = mkvec2(coeff_minus_1, coeff_minus_2);
    *coeff = mkvec2(coeff_plus, coeff_minus);

    *I = kernel_to_ideal_O0_T(*coeff);
    // GEN I_even_int;

    //this if we want to have 2^* isogenies in the commitment
    if (need_even_commit) {
      *I_even = lideal_random_2e(global_setup.B, global_setup.O0,two_tors_height);
      //this loop is necessary to ensure that ideal_to_isogeny_1728 works
      while (!lideal_equals( lideal_create(global_setup.B, global_setup.O0,lideal_generator(*I_even),gen_2),lideal_create(global_setup.B, global_setup.O0,mkcol4(gen_1,gen_1,gen_0,gen_0),gen_2)  )) {
        *I_even = lideal_random_2e(global_setup.B, global_setup.O0,two_tors_height);
      }
      proj A = global_setup.E0;
      ideal_to_isogeny_1728(phi_com_even,*I_even,&A);
    }
    else {
      *I_even = lideal_create( global_setup.B, global_setup.O0,global_setup.one,gen_1 );
    }



    proj ker = coeff_to_E0(gel(*coeff,1), false);
    proj ker_twist = coeff_to_E0(gel(*coeff,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(*I)));


    phi_com->kernel_plus = ker;
    phi_com->kernel_minus = ker_twist;
    phi_com->deg_plus = deg;
    phi_com->deg_minus = deg_twist;
    gerepileall(ltop, 3, I, I_even, coeff);
}



void deterministic_point(proj *P, proj const *A, long ell, long e, bool twist, GEN seed) {
    pari_sp ltop = avma;
    GEN rand_state = getrand();
    setrand(seed);

    unsigned short newseed[3] = {1,2,3};
    unsigned short *oldptr = seed48(newseed);

    uintbig cofactor;
    uintbig_add3(&cofactor, &p, &uintbig_1);
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell);
    }
    proj Z;

    while (1) {
        fp2_random(&P->x); fp2_random(&P->z);
        if (twist == is_on_curve(P, A)) continue;
        xMUL(P, A, P, &cofactor);
        Z = *P;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&Z, A, &Z, &ell_big);
        }
        if (!fp2_iszero(&Z.z)) {
            // a final test
            xMUL(&Z, A, &Z, &ell_big);
            assert(fp2_iszero(&Z.z));
            return;
        }
    }

    setrand(rand_state);
    seed48(oldptr);
    avma = ltop;
}


void deterministic_basis(point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
    proj P1, P2;
    point P1_mul, P2_mul, tmp;
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    fp2 weil;

    proj E;
    mont_to_ted(&E, A, twist);

    GEN seed = gen_1;

    deterministic_point(&P1, A, ell, e, twist, seed);

    seed = gadd(seed, gen_1);

    mont_to_ted_point(P1_ted, A, &P1);
    P1_mul = *P1_ted;
    for (int i = 0; i < e-1; ++i) {
        ted_mul(&P1_mul, &P1_mul, &E, &ell_big);
    }

    assert(ted_is_on_curve(&P1_mul,&E));
    assert(!ted_iszero(&P1_mul));
    ted_mul(&tmp, &P1_mul, &E, &ell_big);
    assert(ted_iszero(&tmp));

    do {
        deterministic_point(&P2, A, ell, e, twist, seed);
        mont_to_ted_point(P2_ted, A, &P2);
        P2_mul = *P2_ted;
        for (int i = 0; i < e-1; ++i) {
            ted_mul(&P2_mul, &P2_mul, &E, &ell_big);
        }
        assert(ted_is_on_curve(&P2_mul,&E));
        assert(!ted_iszero(&P2_mul));
        ted_mul(&tmp, &P2_mul, &E, &ell_big);
        assert(ted_iszero(&tmp));

        ted_weil(&weil, &E, &P1_mul, &P2_mul, &ell_big);
        fp2_sub2(&weil, &fp2_1);
        seed = gadd(seed, gen_1);
    } while (fp2_iszero(&weil));
}




bool mont_dlp(GEN *a, GEN *b, const proj *A, const proj *P, const proj *P1,
    const proj *P2, const proj *P1plusP2, long ell, long e, bool twist) {
    proj Q;
    proj E;
    point basis_ted[3], P_ted;
    mont_to_ted(&E, A, twist);

    mont_to_ted_point(&P_ted, A, P);
    mont_to_ted_point(&basis_ted[0], A, P1);
    mont_to_ted_point(&basis_ted[1], A, P2);

    ted_add(&basis_ted[2], &E, &basis_ted[0], &basis_ted[1]);
    ted_to_mont_point(&Q, &basis_ted[2]);
    if (!mont_equal(&Q, P1plusP2)) { ted_neg(&basis_ted[1], &basis_ted[1]); }
    ted_add(&basis_ted[2], &E, &basis_ted[0], &basis_ted[1]);
    ted_to_mont_point(&Q, &basis_ted[2]);
    assert(mont_equal(&Q, P1plusP2));

    bool found = ted_bidim_log(a, b, &E, &P_ted, &basis_ted[0], &basis_ted[1], ell, e);

    return found;
}


//in fact only work for ell=2
bool mont_power_dlp(uintbig *a,const proj *A, const proj *Q, const proj *P,const proj *PQ,long e) {
  return mont_two_DLP(a,A,Q,P,PQ,e);
}



//renormalize the two walk into a walk of step having length two_tors_height
//zip contains a compressed representation of this walk
// return the length of the last step
long normalized_walk(two_walk *w,uintbig *zip, long *n) {
  long tot = 0;
  for (int i=0;i<*n;i++) {
    tot+=w[i].len;
  }
  long step= tot / two_tors_height;
  long add_step=1;
  if (tot == step*two_tors_height) {
    add_step=0;
  }
  long index_w = 0;
  two_walk norm_walk[step+add_step];
  uintbig big_res;
  // proj Pc[step+add_step],Qc[step+add_step],PQc[step+add_step],Ac[step+add_step];
  while ( w[index_w].len == two_tors_height ) {
    norm_walk[index_w]=w[index_w];
    //this computes a basis <P,Q> and the value log such that the kernel is generated by
    // P + log*Q
    // currently the algorithm is pushing the basis through and computing a DLP
    // for this one we already have the kernel but doing a bidimensionnal DLP appears difficult
    // as we would need some extra point (the difference between the kernel and the points of the basis ?)
    proj P1,P2,P3,dummy_A;
    normalize_proj(&w[index_w].A);
    find_basis_2e(&P1,&P2,&P3,&w[index_w].A);
    // Pc[0]=P1;
    // Qc[0]=P2;
    // PQc[0]=P3;
    // Ac[0]=w[index_w].A;
    proj push_points[3];
    push_points[0]=P2;
    push_points[1]=P1;
    push_points[2]=P3;
    eval_walk_mult(&w[index_w],&dummy_A,push_points,3);
    // what we do when P is the kernel generator ?

    bool dlp=mont_power_dlp(&zip[index_w],&dummy_A,&(push_points[0]),&(push_points[1]),&(push_points[2]),two_tors_height);

    // uintbig_set(&zip[index_w],res);
    uintbig_set(&big_res,pow(2,two_tors_height));
    // zip[index_w]=pow(2,two_tors_height) - zip[index_w];
    uintbig_sub3(&zip[index_w],&big_res,&zip[index_w]);

    assert(dlp);
    _unused(dlp);

    index_w++;
    isomorphism isom;
    mont_isom(&isom,&w[index_w].A,&dummy_A);
    mont_isom_apply(&isom,&w[index_w].ker);
    // w[index_w].A= dummy_A;
  }

  proj A = w[index_w].A;
  proj P = w[index_w].ker;
  long order_P = w[index_w].len;
  assert(mont_equal(&A,&w[index_w].A));


  for (int index=index_w; index < step; index++ ){
    norm_walk[index].A=A;
    norm_walk[index].len=two_tors_height;
    proj P1,P2,P3;
    proj dummy_A;
    normalize_proj(&A);
    find_basis_2e(&P1,&P2,&P3,&A);
    // Pc[index]=P1;
    // Qc[index]=P2;
    // PQc[index]=P3;
    // Ac[index]=A;
    #ifndef NDEBUG
    proj test_order;
    test_order=P1;
    for (int i=1;i<two_tors_height; i++) {
      xDBL(&test_order,&A,&test_order);
    }
    assert(!mont_iszero(&test_order));
    xDBL(&test_order,&A,&test_order);
    assert(mont_iszero(&test_order));
    uintbig a;
    assert(!mont_power_dlp(&a,&A,&P2,&P1,&P3,two_tors_height));
    uintbig log_test;
    long mult_fac=55;
    uintbig x_test;
    uintbig_set(&log_test,mult_fac);
    proj test_point,test_pt2;
    xMUL(&test_point,&A,&P1,&log_test);
    uintbig_set(&log_test,mult_fac+1);
    xMUL(&test_pt2,&A,&P1,&log_test);

    assert(mont_power_dlp(&x_test,&A,&test_point,&P1,&test_pt2,two_tors_height));
    #endif

    long w_step=1;
    long remain_step = order_P;
    while (remain_step+w[index_w + w_step].len <two_tors_height) {
      remain_step+=w[index_w+w_step].len;
      w_step++;
    }
    two_walk loc_phi[w_step+1];
    loc_phi[0].A = A;
    loc_phi[0].ker=P;
    loc_phi[0].len =order_P;


    if (w_step > 1) {
      for (int i=1; i < w_step ;i++){
        loc_phi[i] = w[index_w+i];
      }
      // remain_step+=w[index_w+1].len;
    }
    remain_step = two_tors_height - remain_step;

    if (remain_step == w[index_w+w_step].len) {
        loc_phi[w_step] = w[index_w+w_step];
        index_w += w_step+1;
        order_P = w[index_w].len;
        P = w[index_w].ker;
    }
    else {
      loc_phi[w_step].A = w[index_w+w_step].A;
      loc_phi[w_step].len=remain_step;
      xDBL(&loc_phi[w_step].ker,&w[index_w+w_step].A,&w[index_w+w_step].ker);
      for (int i=0; i < w[index_w+w_step].len - remain_step -1 ; i++ ) {
        xDBL(&loc_phi[w_step].ker,&w[index_w+w_step].A,&loc_phi[w_step].ker);
      }
      index_w += w_step;
      order_P = w[index_w].len - remain_step;
      P=w[index_w].ker;
      isomorphism isom;
      eval_walk_isom_mult(&isom,&(loc_phi[w_step]), &dummy_A,&(loc_phi[w_step]),&P,1);
    }
    proj push_points[3];
    push_points[0]=P2;push_points[1]=P1;push_points[2]=P3;

    for (int i=0; i < w_step +1 ;i ++) {
      // long log;
      isomorphism isom;
      eval_walk_isom_mult(&isom,&(loc_phi[i]),&A,&(loc_phi[i]),push_points,3);
      if (i != w_step) {
        isomorphism isom;
        mont_isom(&isom,&loc_phi[i+1].A,&A);
        mont_isom_apply(&isom,&loc_phi[i+1].ker);
        loc_phi[i+1].A=A;
      }
      if (i == w_step) {
        isomorphism isom;
        mont_isom(&isom,&dummy_A,&A);
        mont_isom_apply(&isom,&P);
      }
    }
    uintbig x;

    if (mont_iszero(&(push_points[0])) ) {
      uintbig_set(&x,1);
      uintbig_set(&zip[index],0);
    }
    if (mont_iszero(&(push_points[1])) ){
      uintbig_set(&x,0);
      uintbig_set(&zip[index],1);
    }
    else {
      uintbig_set(&x,1);
      bool dlp =mont_power_dlp(&zip[index],&A,&(push_points[0]),&(push_points[1]),&(push_points[2]),two_tors_height );
      // zip[index]=pow(2,two_tors_height) - zip[index];
      // uintbig_set(&y,zip[index]);
      // uintbig_set(&zip[index],res);
      uintbig_set(&big_res,pow(2,two_tors_height));
      uintbig_sub3(&zip[index],&big_res,&zip[index]);
      assert(dlp);
      _unused(dlp);
      #ifndef NDEBUG
      proj verif_dlp;
      verif_dlp=push_points[1];
      xMUL(&verif_dlp,&A,&verif_dlp,&zip[index]);
      assert(mont_equal(&verif_dlp,&(push_points[0])));
      #endif
    }

    xBIDIM(&norm_walk[index].ker,&(loc_phi[0].A),&P2,&x,&P1,&zip[index],&P3);
    isomorphism isom;
    proj A_target;
    eval_walk_isom(&isom,&norm_walk[index],&A_target,NULL,&norm_walk[index],NULL);
    #ifndef NDEBUG
    proj verif=P;
    for (int i=1; i < order_P; i++) {
      xDBL(&verif,&A,&verif);
    }
    assert(!mont_iszero(&verif));
    xDBL(&verif,&A,&verif);
    assert(mont_iszero(&verif));
    // proj proj_tmp = {fp2_1, fp2_0};
    proj j1,j2;
    jinv256(&j1,&A);
    // eval_walk(&norm_walk[index],&A_target,&proj_tmp);;
    jinv256(&j2,&A_target);
    assert(mont_equal(&j1, &j2));
    #endif

  }
  if (add_step == 1){
    norm_walk[step].A=A;
    norm_walk[step].ker=P;
    norm_walk[step].len=order_P;

    //compute the compression coeff for the last one
    isomorphism isom;

    proj P1,P2,P3,A_target;
    // eval_walk_isom(&isom,&norm_walk[step],&A_target,NULL,&norm_walk[step],NULL);
    normalize_proj(&A);
    find_basis_2e(&P1,&P2,&P3,&A);
    two_walk phi_test;
    phi_test.A=A;
    phi_test.len=order_P;
    proj push_points[3];
    push_points[0]=P2;
    push_points[1]=P1;
    push_points[2]=P3;
    for (int i=0;i < two_tors_height - order_P; i++) {
      for (int j=0; j<3; j++){
        xDBL(&(push_points[j]),&A,&(push_points[j]));
      }
    }
    P1=push_points[1];
    P2=push_points[0];
    P3=push_points[2];
    eval_walk_isom_mult(&isom,&norm_walk[step],&A_target,&norm_walk[step],push_points,3);
    //what we do when P1 is the kernel generator ?
    bool dlp=mont_power_dlp(&zip[step],&A_target,&(push_points[0]),&(push_points[1]),&(push_points[2]),order_P);

    // uintbig_set(&zip[step],res);
    uintbig_set(&big_res,pow(2,two_tors_height));
    uintbig_sub3(&zip[step],&big_res,&zip[step]);
    // zip[step]= pow(2,order_P) - zip[step];
    assert(dlp);
    _unused(dlp);
    #ifndef NDEBUG
    proj A_test;
    // uintbig_set(&a,zip[step]);
    xBIDIM(&phi_test.ker,&phi_test.A,&P1,&zip[step],&P2,&uintbig_1,&P3);
    isomorphism isom2;
    eval_walk_isom(&isom2,&phi_test,&A_test,NULL,&phi_test,NULL);
    proj j1,j2;
    jinv256(&j1,&A_target);
    jinv256(&j2,&A_test);
    assert(mont_equal(&j1,&j2));
    #endif

  }
  for ( int i=0; i<step+add_step; i++){
    w[i]=norm_walk[i];
  }
  *n=step+add_step;
  #ifndef NDEBUG
  proj Q, PQ;
  A=w[0].A;
  two_walk walk[step+add_step];
  for (int i = 0; i < step; i++) {
    // uintbig_set(&a, zip[i] );
    // long hint = (zip[i] & hint_mask) >> two_tors_height;
    // get the next kernel
    normalize_proj(&A);
    find_basis_2e(&P, &Q, &PQ, &A);  // TODO: use point Q from previous step + hint

    // assert(mont_equal(&A,&Ac[i]));
    // assert(mont_equal(&P,&Pc[i]));
    // assert(mont_equal(&Q,&Qc[i]));
    // assert(mont_equal(&PQ,&PQc[i]));
    xBIDIM(&(walk[i].ker), &A, &P, &zip[i], &Q, &uintbig_1, &PQ);
    // assert(mont_equal(&walk[i].ker,&w[i].ker));
    walk[i].A = A;
    walk[i].len = two_tors_height;
    // take the next step
    isomorphism isom;
    eval_walk_isom(&isom,&walk[i], &A, NULL,&walk[i],NULL);
    // A=walk[i].A;
    proj j1,j2;
    jinv256(&j1,&A);
    jinv256(&j2,&w[i+1].A);
    assert(mont_equal(&j1,&j2));
  }
  #endif

  return order_P;

}

//this function is used to hash a coefficient for the even challenge from the message
//right now this is the identity
static void hash_to_even_chall(uintbig *even_chall,const uintbig *m) {
  even_chall->c[0] = 3*m->c[0]+2;
  even_chall->c[1] = 3*m->c[1]+7;
  even_chall->c[2] =5*m->c[2]+5;
  even_chall->c[3] = 7*m->c[3]+2;
}


void challenge(proj *E_cha, two_walk_long* phi_chall_even, const uintbig *m, const proj *E_com, const proj *basis_plus, const proj *basis_minus, GEN *dlog, proj *basis_two){
    pari_sp ltop = avma;
    //unsigned short newseed[3] = {1,2,3};
    //unsigned short *oldptr = seed48(newseed);

    odd_isogeny phi;
    proj A = *E_com;
    long index;
    bool twist;
    uintbig k;
    long ell;
    proj H, P, Z,X;
    _unused(X);
    uintbig_set(&H.x.re.x, m->c[0]);
    uintbig_set(&H.x.im.x, m->c[1]);
    uintbig_set(&H.z.re.x, m->c[2]);
    uintbig_set(&H.z.im.x, m->c[3]);


    uintbig cofactor_plus = {1,0,0,0}, cofactor_minus = {1,0,0,0};
    uintbig order_plus = {1,0,0,0}, order_minus = {1,0,0,0};
    isog_degree deg;
    degree_one(&deg);


    // compute cofactor and degree of the 'plus' part
    for (int i = 0; i < p_plus_len; ++i) {
        ell = p_plus_fact[i];
        for (int j = 0; j < p_plus_mult[i] - p_plus_mult_cha[i]; ++j){
            uintbig_mul3_64(&cofactor_plus, &cofactor_plus, ell);
        }
        for (int j = 0; j < p_plus_mult_cha[i]; ++j){
            uintbig_mul3_64(&order_plus, &order_plus, ell);
            index = ell_to_index(p_plus_fact[i], &twist);
            degree_set(&deg, index, p_plus_mult_cha[i]);
        }

    }


    bool bad;

    // find the 'plus' part of the kernel
    while (1) {
        //fp2_random(&P->x); fp2_random(&P->z);
        fp_add2(&H.x.re, &fp_1);
        if (!is_on_curve(&H, E_com)) continue;
        xMUL(&P, E_com, &H, &cofactor_plus);
        xMUL(&P, E_com, &P, &p_plus_odd_cofactor);

        bad = false;
        for (int i = 0; i < p_plus_len; ++i) {
            if (p_plus_mult_cha[i] > 0) {
                ell = p_plus_fact[i];
                Z = P;
                uintbig_div3_64(&k, &order_plus, ell);

                xMUL(&Z, E_com, &Z, &k);
                if (mont_iszero(&Z)) { bad = true; break; }

                #ifndef NDEBUG
                uintbig ell_big;
                uintbig_set(&ell_big, ell);
                xMUL(&Z, E_com, &Z, &ell_big);
                assert(mont_iszero(&Z));
                #endif
            }
        }
        if (bad) continue;
        else break;

    }

    // At this point, P generates the 'plus' part of the kernel
    GEN dlog_plus = NULL;
    if (basis_plus) {
        long len = p_plus_len;
        GEN coeff_1 = cgetg(len+1, t_VEC);
        GEN coeff_2 = cgetg(len+1, t_VEC);
        proj Q;
        _unused(Q);
        for (int i = 0; i < len; ++i) {

            gel(coeff_1,i+1) = gen_0;
            gel(coeff_2,i+1) = gen_0;

            if (0 < p_plus_mult_cha[i]) {
                GEN a,b;

                long ell = p_plus_fact[i];

                // find basis of the basis of the ell-torsion
                proj basis_i[3], P_i;

                uintbig_div3_64(&k, &order_plus, ell);
                for (int j = 1; j < p_plus_mult_cha[i]; ++j){
                  uintbig_div3_64(&k, &k, ell);
                }
                xMUL(&basis_i[0], E_com, &basis_plus[0], &k);
                xMUL(&basis_i[1], E_com, &basis_plus[1], &k);
                xMUL(&basis_i[2], E_com, &basis_plus[2], &k);
                xMUL(&P_i, E_com, &P, &k);




                bool dlp = mont_dlp(&a, &b, &A, &P_i, &basis_i[0], &basis_i[1], &basis_i[2],
                    ell, p_plus_mult_cha[i], false);
                assert(dlp);
                _unused(dlp);
                if (!dlp ){
                  printf("no dlps !!!\n");
                }


                if (dlp) {
                  gel(coeff_1,i+1) = a;
                  gel(coeff_2,i+1) = b;
                  #ifndef NDEBUG
                  uintbig a_big,b_big;
                  gentobig(&a_big, gel(coeff_1,i+1));
                  gentobig(&b_big, gel(coeff_2,i+1));
                  // test correctness dlp
                  xBIDIM(&Z, &A, &basis_i[0], &a_big, &basis_i[1], &b_big, &basis_i[2]);
                  assert(mont_equal(&Z,&P_i));
                  #endif
                }

            }
        }

        //dlog_plus = torsion_crt_compose(mkvec2(coeff_1, coeff_2), false);
        dlog_plus = mkvec2(coeff_1, coeff_2);

        GEN dlog_plus_composed = torsion_crt_compose(dlog_plus, false);
        _unused(dlog_plus_composed);
        #ifndef NDEBUG
        uintbig a_big, b_big;

        gentobig(&a_big, gel(dlog_plus_composed, 1));
        gentobig(&b_big, gel(dlog_plus_composed, 2));
        xBIDIM(&Q, &A, &basis_plus[0], &a_big, &basis_plus[1], &b_big, &basis_plus[2]);


        // // final test
        // assert(mont_equal(&Q,&P));
        // this test doesn't pass because the mont points are only up to sign,
        // so P = ±Q at every prime, but the choice of sign may defer between primes

        // so instead, we test correctness at each prime
        for (int i = 0; i < len; ++i) {
            long ell = p_plus_fact[i];
            k = order_plus;
            for (int j = 0; j < p_plus_mult_cha[i]; ++j){
              uintbig_div3_64(&k, &k, ell);
            }
            xMUL(&Z, E_com, &P, &k);
            xMUL(&X, E_com, &Q, &k);
            assert(mont_equal(&Z,&X));
        }

        // test that the order is correct
        xMUL(&Z, E_com, &P, &order_plus);
        assert(mont_iszero(&Z));
        xMUL(&Z, E_com, &Q, &order_plus);
        assert(mont_iszero(&Z));

        #endif

    }


    phi.kernel_plus = P;
    phi.kernel_minus.x = fp2_1;
    phi.kernel_minus.z = fp2_0;
    phi.deg_plus = deg;
    phi.deg_minus.val = 0;


    long len_points = (basis_plus) ? 3 : 1;
    if (basis_two) len_points += 3;

    proj points[len_points];
    points[0] = P;

    if (basis_minus) {
        points[0] = basis_minus[0];
        points[1] = basis_minus[1];
        points[2] = basis_minus[2];
    }
    if (basis_two) {
        points[3] = basis_two[0];
        points[4] = basis_two[1];
        points[5] = basis_two[2];
    }
    eval_mult(&A, &phi, points, len_points);

    degree_one(&deg);

    // compute cofactor and degree of the 'minus' part
    for (int i = 0; i < p_minus_len; ++i) {
        ell = p_minus_fact[i];
        for (int j = 0; j < p_minus_mult[i] - p_minus_mult_cha[i]; ++j){
            uintbig_mul3_64(&cofactor_minus, &cofactor_minus, ell);
        }
        for (int j = 0; j < p_minus_mult_cha[i]; ++j){
            uintbig_mul3_64(&order_minus, &order_minus, ell);
            index = ell_to_index(p_minus_fact[i], &twist);
            degree_set(&deg, index, p_minus_mult_cha[i]);
        }

    }

    // find the 'minus' part of the kernel
    while (1) {
        fp_add2(&H.x.re, &fp_1);
        //fp2_random(&P->x); fp2_random(&P->z);
        if (is_on_curve(&H, &A)) continue;
        xMUL(&P, &A, &H, &cofactor_minus);
        xMUL(&P, &A, &P, &p_minus_odd_cofactor);


        bad = false;
        for (int i = 0; i < p_minus_len; ++i) {
            if (p_minus_mult_cha[i] > 0) {
                ell = p_minus_fact[i];
                Z = P;
                uintbig_div3_64(&k, &order_minus, ell);

                xMUL(&Z, &A, &Z, &k);
                if (mont_iszero(&Z)) { bad = true; break; }

                #ifndef NDEBUG
		            uintbig ell_big;
                uintbig_set(&ell_big, ell);
                xMUL(&Z, &A, &Z, &ell_big);
                assert(mont_iszero(&Z));
                #endif
            }
        }
        if (bad) continue;
        else break;
    }

    // At this point, P generates the 'minus' part of the kernel

    GEN dlog_minus = NULL;
    if (basis_minus) {
        long len = p_minus_len;
        GEN coeff_1 = cgetg(len+1, t_VEC);
        GEN coeff_2 = cgetg(len+1, t_VEC);
        proj Q;
        _unused(Q);
        for (int i = 0; i < len; ++i) {

            gel(coeff_1,i+1) = gen_0;
            gel(coeff_2,i+1) = gen_0;

            // if (0 < p_minus_mult_cha[i]) {
            //     GEN a,b;


            //     bool dlp = mont_dlp(&a, &b, &A, &P, &points[0], &points[1], &points[2],
            //         p_minus_fact[i], p_minus_mult_cha[i], true);
            //     assert(dlp);
            //     _unused(dlp);
            //     gel(coeff_1,i+1) = a;
            //     gel(coeff_2,i+1) = b;
            // }



            if (0 < p_minus_mult_cha[i]) {
                GEN a,b;

                long ell = p_minus_fact[i];

                // find basis of the basis of the ell-torsion
                proj basis_i[3], P_i;

                uintbig_div3_64(&k, &order_minus, ell);
                for (int j = 1; j < p_minus_mult_cha[i]; ++j){
                  uintbig_div3_64(&k, &k, ell);
                }
                xMUL(&basis_i[0], &A, &points[0], &k);
                xMUL(&basis_i[1], &A, &points[1], &k);
                xMUL(&basis_i[2], &A, &points[2], &k);
                xMUL(&P_i, &A, &P, &k);

                xMUL(&Z, &A, &basis_i[0], &order_minus);
                assert(mont_iszero(&Z));
                xMUL(&Z, &A, &basis_i[1], &order_minus);
                assert(mont_iszero(&Z));
                xMUL(&Z, &A, &basis_i[2], &order_minus);
                assert(mont_iszero(&Z));
                xMUL(&Z, &A, &P_i, &order_minus);
                assert(mont_iszero(&Z));


                bool dlp = mont_dlp(&a, &b, &A, &P_i, &basis_i[0], &basis_i[1], &basis_i[2],
                    ell, p_minus_mult_cha[i], true);



                assert(dlp);
                _unused(dlp);

                gel(coeff_1,i+1) = a;
                gel(coeff_2,i+1) = b;

                #ifndef NDEBUG
                uintbig a_big,b_big;
                gentobig(&a_big, gel(coeff_1,i+1));
                gentobig(&b_big, gel(coeff_2,i+1));
                // test correctness dlp
                xBIDIM(&Z, &A, &basis_i[0], &a_big, &basis_i[1], &b_big, &basis_i[2]);
                assert(mont_equal(&Z,&P_i));
                #endif
            }


        }

        //dlog_minus = torsion_crt_compose(mkvec2(coeff_1, coeff_2), true);

        dlog_minus = mkvec2(coeff_1, coeff_2);

        GEN dlog_minus_composed = torsion_crt_compose(dlog_minus, true);
        _unused(dlog_minus_composed);

        #ifndef NDEBUG
        uintbig a_big, b_big;

        gentobig(&a_big, gel(dlog_minus_composed, 1));
        gentobig(&b_big, gel(dlog_minus_composed, 2));
        xBIDIM(&Q, &A, &points[0], &a_big, &points[1], &b_big, &points[2]);


        // // final test
        // assert(mont_equal(&Q,&P));
        // this test doesn't pass because the mont points are only up to sign,
        // so P = ±Q at every prime, but the choice of sign may defer between primes

        // so instead, we test correctness at each prime
        for (int i = 0; i < len; ++i) {
            long ell = p_minus_fact[i];
            k = order_minus;
            for (int j = 0; j < p_minus_mult_cha[i]; ++j){
              uintbig_div3_64(&k, &k, ell);
            }
            xMUL(&Z, &A, &P, &k);
            xMUL(&X, &A, &Q, &k);
            assert(mont_equal(&Z,&X));
        }

        // test that the order is correct
        xMUL(&Z, &A, &P, &order_minus);
        assert(mont_iszero(&Z));
        xMUL(&Z, &A, &Q, &order_minus);
        assert(mont_iszero(&Z));

        #endif

    }

    len_points = (basis_two) ? 3 : 1;

    if (basis_two) {
        points[0] = points[3];
        points[1] = points[4];
        points[2] = points[5];
    }

    phi.kernel_minus = P;
    phi.kernel_plus.x = fp2_1;
    phi.kernel_plus.z = fp2_0;
    phi.deg_minus = deg;
    phi.deg_plus.val = 0;

    eval_mult(&A, &phi, points, len_points);

    if (basis_two) {
        basis_two[0] = points[0];
        basis_two[1] = points[1];
        basis_two[2] = points[2];
    }

    //add a 2^two_tors_height isogeny to the challenge isogeny
    GEN dlog_a = NULL; GEN dlog_b = NULL;
    if (need_even_chall) {
      normalize_proj(&A);
      proj P1,P2,P12;
      find_basis_2e(&P1, &P2, &P12, &A);
      //reducing the order by 2
      xDBL(&P1,&A,&P1);
      xDBL(&P2,&A,&P2);
      xDBL(&P12,&A,&P12);



      uintbig even_chall;


      //this is not properly hashing the message to a challenge
      hash_to_even_chall(&even_chall,m);


      two_walk isog_chall;
      // the -1 is really important
      isog_chall.len = two_tors_height-1;
      isog_chall.A = A;


      uintbig one;
      uintbig_set(&one,1);
      xBIDIM(&isog_chall.ker,&A,&P2,&one,&P1,&even_chall,&P12);
      #ifndef NDEBUG

      #endif


      if (basis_two) {
        //computing the coefficients that we need to find the ideal
        uintbig x1,x2;
        proj2 xyP1,xyP2,xyQ,xyP12;
        xtoxy(&xyP1,&A,&points[0]);
        xtoxy(&xyP2,&A,&points[1]);
        xtoxy(&xyQ,&A,&isog_chall.ker);
        proj2 inter;
        xyADD(&inter,&A,&xyP1,&xyP2);
        proj test_sign;
        xytox(&test_sign,&inter);
        if(!mont_equal(&test_sign,&points[2])) {
          xyNEG(&xyP2,&xyP2);
        }
        xyADD(&xyP12,&A,&xyP1,&xyP2);
        proj2 xyPQ1,xyPQ2,xyPQ12;
        xyADD(&xyPQ1,&A,&xyP1,&xyQ);
        xyADD(&xyPQ2,&A,&xyP2,&xyQ);
        xyADD(&xyPQ12,&A,&xyP12,&xyQ);
        proj PQ1,PQ2,PQ12;
        xytox(&PQ1,&xyPQ1);
        xytox(&PQ2,&xyPQ2);
        xytox(&PQ12,&xyPQ12);
        bool dlpt1 = mont_bidim_two_DLP(&x1,&x2,&A,&isog_chall.ker,&points[0],&points[1],&points[2],&PQ1,&PQ2,&PQ12,two_tors_height-1);
        assert(dlpt1);
        _unused(dlpt1);
        dlog_a = bigtogen(&x1);
        dlog_b = bigtogen(&x2);


      }
      #ifndef NDEBUG
      proj R= isog_chall.ker;

      for (int i =0;i< two_tors_height -2;i++) {
        xDBL(&R,&A,&R);
      }
      assert(!mont_iszero(&R));
      xDBL(&R,&A,&R);
      assert(mont_iszero(&R));

      #endif

      two_walk_stol(phi_chall_even, &isog_chall);
      eval_walk(&isog_chall,&A,&P12);
    }




    *E_cha = A;

    if (dlog && !need_even_chall) { *dlog = gerepilecopy(ltop, mkvec2(dlog_plus, dlog_minus)); }
    else if (dlog && need_even_chall) { *dlog = gerepilecopy(ltop, mkvec4(dlog_plus, dlog_minus,dlog_a,dlog_b)); }

    //seed48(oldptr);
    //avma = ltop;
}

//function to test the decompression
void decompress(two_walk *walk, proj *A, const uintbig *zip, long len,long last_step) {

  proj P, Q, PQ;

  for (int i = 0; i < len-1; i++) {
    // printf("zip %ld ",zip[i]);

    // uintbig_set(&a, zip[i]);
    // long hint = (zip[i] & hint_mask) >> two_tors_height;
    // get the next kernel
    normalize_proj(A);
    find_basis_2e(&P, &Q, &PQ, A);  // TODO: use point Q from previous step + hint
    xBIDIM(&(walk[i].ker), A, &P, &zip[i], &Q, &uintbig_1, &PQ);
    // printf(" k %ld ",walk[i].ker.x.re.x.c[0]);
    walk[i].A = *A;
    walk[i].len = two_tors_height;
    // take the next step
    isomorphism isom;
    eval_walk_isom(&isom,&walk[i], A, NULL,&walk[i],NULL);
    // printf("\n");
  }
  //last step of smaller size
  // uintbig_set(&a, zip[len-1]);
  // long hint = (zip[len-1] & hint_mask) >> two_tors_height;
  // get the next kernel
  normalize_proj(A);
  find_basis_2e(&P, &Q, &PQ, A);
  // printf(" A %ld \n",A->x.re.x.c[0]);
  for (int i=0; i < two_tors_height - last_step ; i++){
    xDBL(&P,A,&P);
    xDBL(&Q,A,&Q);
    xDBL(&PQ,A,&PQ);
  }
  xBIDIM(&walk[len-1].ker, A, &P, &zip[len-1], &Q, &uintbig_1, &PQ);
  walk[len-1].A = *A;
  walk[len-1].len = last_step;
  isomorphism isom;
  eval_walk_isom(&isom,&walk[len-1],A,NULL,&walk[len-1], NULL);
  // printf(" A %ld \n",A->x.re.x.c[0]);
}



void decompress_new(two_walk *walk,proj *A, const uintbig *zip, long len,long last_step) {
  proj P, Q, PQ,temp;
  //first_step
  normalize_proj(A);
  find_basis_2e(&P, &Q, &PQ, A);
  xBIDIM(&(walk[0].ker), A, &P, &zip[0], &Q, &uintbig_1, &PQ);
  walk[0].A = *A;
  walk[0].len = two_tors_height;
  eval_walk(&walk[0],A,&P);
  //next steps
  for (int i =1 ; i< len-1; i++) {
    deterministic_second_point(&Q,&P,A,two_tors_height);
    xBILIFT(&PQ,&temp,&P,&Q,A);
    xBIDIM(&(walk[i].ker), A, &P, &zip[i], &Q, &uintbig_1, &PQ);
    walk[i].A = *A;
    walk[i].len = two_tors_height;
    eval_walk(&walk[i],A,&P);
  }

  //last step of smaller size

  for (int i=0; i < two_tors_height - last_step ; i++){
    xDBL(&P,A,&P);
  }
  // normalize_proj(A);normalize_proj(&P);
  deterministic_second_point(&Q,&P,A,last_step);
  for (int i=0; i < two_tors_height - last_step ; i++){
    xDBL(&Q,A,&Q);
  }
  xBILIFT(&PQ,&temp,&P,&Q,A);
  // printf("A = "); print_proj_hash(A); printf("\n");
  // printf("P = "); print_proj_hash(&P); printf("\n");
  // printf("Q = "); print_proj_hash(&Q); printf("\n");
  // printf("PQ "); print_proj_hash(&PQ); printf("\n");
  xBIDIM(&walk[len-1].ker, A, &P, &zip[len-1], &Q, &uintbig_1, &PQ);
  walk[len-1].A = *A;
  walk[len-1].len = last_step;
  eval_walk(&walk[len-1],A,&P);
}

// basis_two is the image of the basis of the two torsion through phi_com and phi_cha
void response(two_walk_long *sigma,uintbig *zip,  GEN coeff_ker_challenge_commitment, const secret_key *sk, const proj *basis_two, const proj *E_cha){
    pari_sp ltop = avma;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;

    GEN I = sk->I_large_prime;
    GEN I_two = sk->I_two;

    // STEP 1: compute the ideal of phi_challenge_commitment
    GEN I_phipsi = kernel_to_ideal_O0_T(coeff_ker_challenge_commitment);
    // STEP 2: find J of norm a power of two such that J inter I is equivalent to I_phipsi
    GEN beta = lideal_isom(I_two, I); // I_two*beta = I
    GEN alpha = gmul(beta, lideal_norm(I_two)); // n(alpha) = n(I)
    GEN gamma = lideal_generator_coprime(I_phipsi, gmul(gen_2, lideal_norm(I)));
    GEN generator = algmul(A, gamma, alpha);
    GEN norm = gmul(lideal_norm(I_two), lideal_norm(I_phipsi));
    GEN K = lideal_create(A, order, generator, norm);

    assert(lideal_isom(I_phipsi, lideal_inter(K,I))); // K inter I is equivalent to I_phipsi


    GEN J;
    GEN n;
    do{
        J = klpt_general_power(I, K, gen_2); // J inter I is equivalent to I_phipsi
        // printf("path length %ld\n", Z_lval(lideal_norm(J),2));
        alg_primitive(&n, A, order, algmul(A, lideal_generator(J), lideal_generator(I_two)));

        // backtracking?
    } while(gcmp(n,gen_1) != 0);
    // printf("path length %ld\n", Z_lval(lideal_norm(J),2));

    // STEP 3: compute L of norm n(J) such that L inter sk->I_T is equivalent to I_phipsi

    GEN L = lideal_inter(J,I);
    assert(lideal_isom(I_phipsi, L));

    beta = lideal_isom(I, sk->I_T); // I*beta = I_T;

    L = lideal_mul(L,beta);

    assert(lideal_isom(I_phipsi, L));
    assert(gcmp(lideal_norm(L), gmul(lideal_norm(sk->I_T), lideal_norm(J))) == 0);


    GEN dummy_ideal;
    special_isogeny dummy_isogeny;

    // STEP 4: convert to an isogeny


    #ifndef NDEBUG
    alg_primitive(&n, A, order, algmul(A, lideal_generator(L), lideal_generator(sk->I_two)));
    assert(gcmp(n,gen_1) == 0);
    #endif

    long delta = 14;
    // long len_tail = two_tors_height + delta;
    GEN L_cut = lideal_create(A, order, lideal_generator(L), shifti(lideal_norm(L), -len_tail));
    // clock_t t = tic();
    ideal_to_isogeny_two(sigma, &dummy_ideal, &dummy_isogeny, L_cut, sk->I_T, sk->I_two, &sk->phi_T, &sk->phi_two, false);
    // printf("ideal_to_isog %f ms \n",toc(t));

    if (len_tail!=0) {
      GEN gen_tail = lideal_isom(L, I_phipsi); // L*gen_tail = I_phipsi;
      gen_tail = gmul(gen_tail, lideal_norm(L));
      GEN L_tail = lideal_create(A, order, gen_tail, powuu(2,two_tors_height));

      two_walk phi_tail_dual;
      phi_tail_dual.A = *E_cha;
      phi_tail_dual.len = two_tors_height;

      GEN v_tail = ideal_to_kernel_O0_ell(L_tail, 2);
      uintbig x, y;
      gentobig(&x, gel(v_tail,1));
      gentobig(&y, gel(v_tail,2));

      xBIDIM(&phi_tail_dual.ker, &phi_tail_dual.A, &basis_two[0], &x, &basis_two[1], &y, &basis_two[2]);


      isomorphism isom;
      proj L_cut_target, phi_tail_dual_target, proj_tmp = {fp2_1, fp2_0};
      eval_walk(&sigma->phi[sigma->len-1], &L_cut_target, &proj_tmp);

      eval_walk_isom(&isom, &phi_tail_dual, &phi_tail_dual_target, NULL, &phi_tail_dual, NULL);

      two_walk eta;

      bool done;

      proj from = L_cut_target;
      proj to = phi_tail_dual_target; // phi_2 source

      done = MITM(&eta, &from, &to, delta);
      assert(done);

      two_walk phi_tail = phi_tail_dual;
      dual_walk(&phi_tail);

      two_walk_composition_sl(sigma, &eta, sigma);
      two_walk_composition_sl(sigma, &phi_tail, sigma);
    }


    normalized_walk(sigma->phi,zip,&(sigma->len));
    // sigma->len = normalized_n;
    // assert(mont_equal(&sigma->phi[0].A)    //

    avma = ltop;
}

// basis_two is the image of the basis of the two torsion through phi_com and phi_cha
void response_new(two_walk_long *sigma,uintbig *zip,  GEN coeff_ker_challenge_commitment, GEN I_even , const secret_key *sk, const proj *basis_two, const proj *E_cha, bool *test){
    pari_sp ltop = avma;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;

    GEN I = sk->I_large_prime;
    GEN I_two = sk->I_two;

    // STEP 1: compute the ideal of phi_challenge_commitment
    GEN I_phipsi = kernel_to_ideal_O0_T(coeff_ker_challenge_commitment);
    I_phipsi = lideal_inter(I_phipsi,I_even);
    // STEP 2: find J of norm a power of two such that J inter I is equivalent to I_phipsi
    GEN beta = lideal_isom(I_two, I); // I_two*beta = I
    GEN alpha = gmul(beta, lideal_norm(I_two)); // n(alpha) = n(I)
    GEN gamma = lideal_generator_coprime(I_phipsi, gmul(gen_2, lideal_norm(I)));
    GEN generator = algmul(A, gamma, alpha);
    GEN norm = gmul(lideal_norm(I_two), lideal_norm(I_phipsi));
    GEN K = lideal_create(A, order, generator, norm);


    gamma =lideal_isom(I_phipsi, lideal_inter(K,I)); // K inter I is equivalent to I_phipsi
    assert(gamma);
    // GEN gam_t = gmul( lideal_norm(lideal_inter(K,I)),gamma);

    GEN J;
    GEN n;
    do{
        J = klpt_general_power(I, K, gen_2); //
        // printf("path length %ld\n", Z_lval(lideal_norm(J),2));
        alg_primitive(&n, A, order, algmul(A, lideal_generator(J), lideal_generator(I_two)));

        // backtracking?
    } while(gcmp(n,gen_1) != 0);


    // STEP 3: compute L of norm n(J)2^* such that L inter sk->I_T is equivalent to I_phipsi
    GEN L = lideal_inter(J,I);

    //Compute last_I which is the end of I_two
    beta = lideal_isom(I, I_two); // I*beta = I_two;
    assert(beta);    //option 1
    GEN last_I = lideal_create(A,lideal_right_order(I_two),alg_conj(A,lideal_generator_coprime(I_two,gen_1)),gpowgs(gen_2,two_tors_height));
    GEN alpha_last =lideal_generator_coprime(last_I,gen_1);
    assert(gequal(lideal_norm(last_I),gpowgs(gen_2,two_tors_height)) );
    alpha_last= algmul(A,alpha_last,alg_conj(A,beta));
    alpha_last= algmul(A,beta,alpha_last);
    alpha_last = algmul(A,alg_scalar(A,gdiv(gen_1,algnorm(A,beta,0))),alpha_last);

    last_I = lideal_create(A,lideal_right_order(I),alpha_last,gpowgs(gen_2,two_tors_height));
    assert(gequal(lideal_norm(last_I),gpowgs(gen_2,two_tors_height)) );



    // STEP 4: convert to an isogeny
    long delta = 14;

    GEN L_cut = lideal_create(A, lideal_order(L), lideal_generator(L), shifti(lideal_norm(L), -len_tail));


    //the value test first indicates if the beginning curve has a special property
    //after the computation it will indicate if the computation as failed somehow
    *test = false;
    //endpoint will indicate the endpoint of the translated isogeny
    proj endpoint,A_fin;
    ideal_to_isogeny_two_new(sigma,sk->I_two,sk->I_large_prime,L_cut,last_I,&sk->phi_two,zip,test,&endpoint);
    A_fin = sigma->phi[sigma->len-1].A;
    // deterministic_second_point(&Q,&sigma->phi[sigma->len-1].ker,&sigma->phi[sigma->len-1].A,two_tors_height);
    eval_walk(&sigma->phi[sigma->len-1],&A_fin,&endpoint);
    if (*test) {
      //len_tail should either be 0, delta or two_tors_height + delta and indicates if some of the translation can be done from the end of the chain, in the event of a yes, then the translation is done on an ideal of norm len_tail
      // fow now this is always 0 or delta because the precomputation on the two_torsion is not initialized.
      //it is set to delta when it allows us to do one less step
       if (len_tail != 0) {


         if (len_tail == delta) {

           //this is to decide if we will need an iso for the mitm

           gamma = lideal_isom(I_phipsi, L);

           assert(gamma);
           gamma= gmul(lideal_norm(I_phipsi),gamma);
           GEN X;
           alglatcontains(A,order,gamma,&X);
           assert(X);
           long need_iso = Z_lval(content(X),2);


           proj Q;
           Q = endpoint;
           for (int i =0 ; i<two_tors_height -  len_tail; i++) {
             xDBL(&Q,&A_fin,&Q);
           }
           proj R,QR,temp;
           proj points[3];
           // normalize_proj(&Q);normalize_proj(&A_fin);
           deterministic_second_point(&R,&Q,&A_fin,len_tail);
           for (int i =0 ; i<two_tors_height -  len_tail; i++) {
             xDBL(&R,&A_fin,&R);
           }
           xBILIFT(&QR,&temp,&Q,&R,&A_fin);
           points[0]= R;points[1]=Q;points[2]=QR;
           two_walk eta;
           //
           bool done;
           proj from = A_fin;
           proj to = *E_cha;
           isomorphism isomo;
           if (need_iso) {
            rand_isom(&isomo,&to);
           }
           done = MITM(&eta, &from, &to, delta);
           //just in case but shound not happen
           while (!done) {
             rand_isom(&isomo,&to);
             done = MITM(&eta, &from, &to, delta);
           }
           assert(done);
           assert(mont_equal(&eta.A,&A_fin));
           two_walk_composition_sl(sigma, &eta, sigma);
           // sigma->phi[sigma->len-1] = eta;

           assert(mont_equal(&eta.A,&A_fin));
           eval_walk_mult(&eta,&A_fin,points,3);
           #ifndef  NDEBUG
           proj j1,j2;
           jinv256(&j1,E_cha);
           jinv256(&j2,&A_fin);
           assert(mont_equal(&j1,&j2));
           #endif
           uintbig big_res;
           done  = mont_power_dlp(&zip[sigma->len-1],&A_fin,&points[0],&points[1],&points[2],len_tail);
           assert(done);
           _unused(done);

           uintbig_set(&big_res,pow(2,len_tail));
           uintbig_sub3(&zip[sigma->len-1],&big_res,&zip[sigma->len-1]);


         }
         else {
           GEN gen_tail = lideal_isom(L, I_phipsi); // L*gen_tail = I_phipsi;
           assert(gen_tail);
           gen_tail = algmul(A, alg_scalar(A,lideal_norm(L)),gen_tail );
           GEN L_tail = lideal_create(A, order, gen_tail, powuu(2,two_tors_height));
           two_walk phi_tail_dual;
           phi_tail_dual.A = *E_cha;
           phi_tail_dual.len = two_tors_height;

           GEN v_tail = ideal_to_kernel_O0_ell(L_tail, 2);
           uintbig x, y;
           gentobig(&x, gel(v_tail,1));
           gentobig(&y, gel(v_tail,2));

           xBIDIM(&phi_tail_dual.ker, &phi_tail_dual.A, &basis_two[0], &x, &basis_two[1], &y, &basis_two[2]);


           isomorphism isom;
           proj L_cut_target, phi_tail_dual_target, proj_tmp = {fp2_1, fp2_0};
           eval_walk(&sigma->phi[sigma->len-1], &L_cut_target, &proj_tmp);

           eval_walk_isom(&isom, &phi_tail_dual, &phi_tail_dual_target, NULL, &phi_tail_dual, NULL);

           two_walk eta;
           //
           bool done;
           //
           proj from = L_cut_target;
           proj to = phi_tail_dual_target; // phi_2 source
           //
           done = MITM(&eta, &from, &to, delta);
           assert(done);
           //
           two_walk phi_tail = phi_tail_dual;
           dual_walk(&phi_tail);
           //
           two_walk_composition_sl(sigma, &eta, sigma);
           two_walk_composition_sl(sigma, &phi_tail, sigma);
           // custom size
           uintbig zip_test2[3];
           long len_end = 3;
           two_walk phi_test[3];
           phi_test[0] = sigma->phi[sigma->len-3];
           phi_test[1] = sigma->phi[sigma->len-2];
           phi_test[2] = sigma->phi[sigma->len-1];
           // phi_test.len = 3;


           normalized_walk(phi_test,zip_test2,&len_end);
           zip[signing_length_two_tors_height_step -3] = zip_test2[0];
           zip[signing_length_two_tors_height_step -2] = zip_test2[1];
           zip[signing_length_two_tors_height_step -1] = zip_test2[2];

           sigma->phi[sigma->len-3] = phi_test[0];
           sigma->phi[sigma->len-2] = phi_test[1];
           sigma->phi[sigma->len-1] = phi_test[2];
         }


      }


      #ifndef NDEBUG
      proj A_check = (sigma->phi)[0].A;;
      proj j1,j2,j3;
      jinv256(&j1,E_cha);
      jinv256(&j3,&A_fin);
      two_walk check[signing_length_two_tors_height_step];
      long last_step= sigma->phi[sigma->len-1].len;
      decompress_new(check,&A_check,zip,(sigma->len),last_step);
      jinv256(&j2,&A_check);
      assert(mont_equal(&j1,&j3));
      assert(mont_equal(&j1,&j2));

      #endif
    }






    avma = ltop;
}

bool simple_check_signature(const two_walk_long *sigma, const uintbig *zip, const public_key *pk, proj *E_cha) {
    proj j1,j2,j3;
    jinv256(&j1, &pk->E);
    jinv256(&j2, &sigma->phi[0].A);
    assert(mont_equal(&j1, &j2));

    proj sigma_target = sigma->phi[0].A, proj_tmp = {fp2_1, fp2_0};
    for (int i = 0; i < sigma->len; ++i) {
        jinv256(&j1, &sigma_target);
        jinv256(&j2, &sigma->phi[i].A);
        if (!mont_equal(&j1, &j2)) return false;
        sigma_target = sigma->phi[i].A;
        eval_walk(&sigma->phi[i], &sigma_target, &proj_tmp);
    }

    jinv256(&j1, &sigma_target);
    normalize_proj(&j1);
    jinv256(&j2, E_cha);
    normalize_proj(&j2);
    assert(mont_equal(&j1,&j2));
    proj A_target=(sigma->phi)[0].A;
    two_walk check[sigma->len];
    long last_step= sigma->phi[sigma->len-1].len;
    decompress(check,&A_target,zip,(sigma->len),last_step);
    jinv256(&j3, &A_target);
    // printf("%ld %ld \n",j2.x.re.x.c[0],j2.x.re.x.c[1]);
    return (mont_equal(&j1, &j2) && mont_equal(&j1,&j3));
}


bool simple_check_signature_new(const two_walk_long *sigma, const uintbig *zip, const public_key *pk, proj *E_cha) {
    proj j1,j2,j3;
    jinv256(&j1, &pk->E);
    jinv256(&j2, &sigma->phi[0].A);
    assert(mont_equal(&j1, &j2));

    proj sigma_target = sigma->phi[0].A, proj_tmp = {fp2_1, fp2_0};
    for (int i = 0; i < sigma->len; ++i) {
        jinv256(&j1, &sigma_target);
        jinv256(&j2, &sigma->phi[i].A);
        if (!mont_equal(&j1, &j2)) return false;
        sigma_target = sigma->phi[i].A;
        eval_walk(&sigma->phi[i], &sigma_target, &proj_tmp);
    }

    jinv256(&j1, &sigma_target);
    normalize_proj(&j1);
    jinv256(&j2, E_cha);
    normalize_proj(&j2);
    assert(mont_equal(&j1,&j2));
    proj A_target=(sigma->phi)[0].A;
    two_walk check[sigma->len];
    long last_step= sigma->phi[sigma->len-1].len;
    decompress_new(check,&A_target,zip,(sigma->len),last_step);
    jinv256(&j3, &A_target);
    // printf("%ld %ld \n",j2.x.re.x.c[0],j2.x.re.x.c[1]);
    return (mont_equal(&j1, &j2) && mont_equal(&j1,&j3));
}


void sign(compressed_signature *comp_sigma, const secret_key *sk, const public_key *pk, const uintbig *m) {
    pari_sp ltop = avma;

    GEN coeff_com, I_com,I_com_even;
    odd_isogeny phi_com;
    two_walk_long phi_com_even;
    init_trivial_two_walk_long(&phi_com_even);
    proj E_cha;
    _unused(pk);
    //printf("Commitment\n");
    commitment(&coeff_com, &I_com, &I_com_even, &phi_com,&phi_com_even);

    comp_sigma->E_com = global_setup.E0;


    // compute the image of a basis of the torsion used for the challenge

    proj points[9];
    points[0] = torsion_basis_sum[0];
    points[1] = torsion_basis_sum[1];
    points[2] = torsion_basis_sum[2];
    points[3] = torsion_basis_twist_sum[0];
    points[4] = torsion_basis_twist_sum[1];
    points[5] = torsion_basis_twist_sum[2];
    points[6] = torsion_basis_two[0];
    points[7] = torsion_basis_two[1];
    points[8] = torsion_basis_two[2];

    eval_mult(&comp_sigma->E_com, &phi_com, points, 9);

    proj basis_plus[3], basis_minus[3], basis_two[3];
    basis_plus[0] = points[0];
    basis_plus[1] = points[1];
    basis_plus[2] = points[2];
    basis_minus[0] = points[3];
    basis_minus[1] = points[4];
    basis_minus[2] = points[5];
    basis_two[0] = points[6];
    basis_two[1] = points[7];
    basis_two[2] = points[8];

    uintbig ell_big;
    for (int i = 0; i < p_plus_len; ++i) {
        uintbig_set(&ell_big, p_plus_fact[i]);
        for (int j = 0; j < p_plus_mult[i] - p_plus_mult_cha[i]; ++j){
            for (int l = 0; l < 3; ++l) {
                xMUL(&basis_plus[l], &comp_sigma->E_com, &basis_plus[l], &ell_big);
            }
        }
    }

    for (int i = 0; i < p_minus_len; ++i) {
        uintbig_set(&ell_big, p_minus_fact[i]);
        for (int j = 0; j < p_minus_mult[i] - p_minus_mult_cha[i]; ++j){
            for (int l = 0; l < 3; ++l) {
                xMUL(&basis_minus[l], &comp_sigma->E_com, &basis_minus[l], &ell_big);
            }
        }
    }

    //printf("Challenge\n");

    // comp_sigma->E_com= Sigma->E_com;

    GEN dlog;
    two_walk_long phi_chall_even;
    challenge(&E_cha, &phi_chall_even ,m, &comp_sigma->E_com, basis_plus, basis_minus, &dlog, basis_two);

    GEN coeff_ker_challenge_commitment = gadd(coeff_com,dlog);


    #ifndef NDEBUG

    odd_isogeny phi_test;
    GEN I = kernel_to_ideal_O0_T(coeff_ker_challenge_commitment);
    proj ker = coeff_to_E0(gel(coeff_ker_challenge_commitment,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff_ker_challenge_commitment,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(I)));

    phi_test.kernel_plus = ker;
    phi_test.kernel_minus = ker_twist;
    phi_test.deg_plus = deg;
    phi_test.deg_minus = deg_twist;

    proj A = global_setup.E0, H = {fp2_1, fp2_0};
    eval(&A, &phi_test, &H);

    assert(mont_equal(&A, &E_cha));

    #endif


    // uint64_t zip[signing_length_two_tors_height_step];
    two_walk_long sigma;
    init_trivial_two_walk_long(&sigma);
    response(&sigma, comp_sigma->zip, coeff_ker_challenge_commitment, sk, basis_two, &E_cha);
    // response_new(&sigma, comp_sigma->zip, coeff_ker_challenge_commitment, sk, basis_two, &E_cha);
    // Sigma->sigma = sigma;
    // zip_copy(comp_sigma,zip,signing_length_two_tors_height_step);


    // for (int i=0;i<signing_length_two_tors_height_step ;i ++ ){
    //     comp_sigma->zip[i]=zip[i];
    // }


    assert(simple_check_signature(&sigma,comp_sigma->zip, pk, &E_cha));

    free_two_walk_long(&sigma);

    avma = ltop;
}

void sign_new(compressed_signature *comp_sigma, const secret_key *sk, const public_key *pk, const uintbig *m) {
    pari_sp ltop = avma;

    GEN coeff_com, I_com,I_even;
    two_walk_long phi_com_even,phi_chall_even;
    init_trivial_two_walk_long(&phi_com_even);
    init_trivial_two_walk_long(&phi_chall_even);
    odd_isogeny phi_com;
    proj E_cha;
    _unused(pk);



    // printf("Commitment\n");
    commitment(&coeff_com, &I_com, &I_even, &phi_com,&phi_com_even);

    comp_sigma->E_com = global_setup.E0;


    // compute the image of a basis of the torsion used for the challenge

    proj points[11];
    points[0] = torsion_basis_sum[0];
    points[1] = torsion_basis_sum[1];
    points[2] = torsion_basis_sum[2];


    points[3] = torsion_basis_twist_sum[0];
    points[4] = torsion_basis_twist_sum[1];
    points[5] = torsion_basis_twist_sum[2];

    if (need_even_chall) {
        proj P;
        P.z = fp2_1;
        random_Fp_point_two_f(&P,&comp_sigma->E_com);
        proj2 Pt1,Pt2,Pt3;
        xtoxy(&Pt1,&comp_sigma->E_com,&P);

        proj2 P_test_corr;
        montxy0_frob(&P_test_corr,&Pt1);
        if (!xy_equal(&P_test_corr,&Pt1)) {
          mont0_dist(&P,&P);
          montxy0_dist(&Pt1,&Pt1);
        }
        #ifndef NDEBUG
        montxy0_frob(&P_test_corr,&Pt1);
        assert(xy_equal(&P_test_corr,&Pt1));
        #endif
        montxy0_dist(&Pt2,&Pt1);
        xyADD(&Pt3,&comp_sigma->E_com,&Pt1,&Pt2);
        xyNEG(&Pt1,&Pt1);
        xyADD(&Pt2,&comp_sigma->E_com,&Pt1,&Pt2);
        xyADD(&Pt1,&comp_sigma->E_com,&Pt3,&Pt2);


        xytox(&points[6],&Pt3);
        xytox(&points[7],&Pt2);
        xytox(&points[8],&Pt1);

        #ifndef NDEBUG
        proj R= points[6];proj S= points[7];

        for (int i =0;i< two_tors_height -2;i++) {
          xDBL(&R,&comp_sigma->E_com,&R);
          xDBL(&S,&comp_sigma->E_com,&S);
        }
        assert(!mont_iszero(&R));
        assert(!mont_iszero(&S));
        assert(!mont_equal(&R,&S));
        xDBL(&R,&comp_sigma->E_com,&R);
        xDBL(&S,&comp_sigma->E_com,&S);
        assert(mont_iszero(&R));
        assert(mont_iszero(&S));

        #endif
    }
    else {
      points[6] = torsion_basis_two[0];
      points[7] = torsion_basis_two[1];
      points[8] = torsion_basis_two[2];
    }




    points[9] = phi_com.kernel_plus;
    points[10] = phi_com.kernel_minus;
    if (need_even_commit) {
      eval_walk_long_mult(&phi_com_even,&comp_sigma->E_com,points,11);
      phi_com.kernel_plus=points[9];
      phi_com.kernel_minus=points[10];
    }
    if (len_tail == 0 && !need_even_chall){
      eval_mult(&comp_sigma->E_com, &phi_com, points, 6);
    }
    else {
      eval_mult(&comp_sigma->E_com, &phi_com, points, 9);
    }


    proj basis_plus[3], basis_minus[3], basis_two[3];
    basis_plus[0] = points[0];
    basis_plus[1] = points[1];
    basis_plus[2] = points[2];
    basis_minus[0] = points[3];
    basis_minus[1] = points[4];
    basis_minus[2] = points[5];
    basis_two[0] = points[6];
    basis_two[1] = points[7];
    basis_two[2] = points[8];

    uintbig ell_big;
    for (int i = 0; i < p_plus_len; ++i) {
        uintbig_set(&ell_big, p_plus_fact[i]);
        for (int j = 0; j < p_plus_mult[i] - p_plus_mult_cha[i]; ++j){
            for (int l = 0; l < 3; ++l) {
                xMUL(&basis_plus[l], &comp_sigma->E_com, &basis_plus[l], &ell_big);
            }
        }
    }

    for (int i = 0; i < p_minus_len; ++i) {
        uintbig_set(&ell_big, p_minus_fact[i]);
        for (int j = 0; j < p_minus_mult[i] - p_minus_mult_cha[i]; ++j){
            for (int l = 0; l < 3; ++l) {
                xMUL(&basis_minus[l], &comp_sigma->E_com, &basis_minus[l], &ell_big);
            }
        }
    }

    // printf("Challenge\n");

    // comp_sigma->E_com= Sigma->E_com;

    GEN dlog;
    if (len_tail == 0 && !need_even_chall){
      challenge(&E_cha, &phi_chall_even, m, &comp_sigma->E_com, basis_plus, basis_minus, &dlog, NULL);
    }
    else {
      challenge(&E_cha, &phi_chall_even ,m, &comp_sigma->E_com, basis_plus, basis_minus, &dlog, basis_two);
    }


    GEN coeff_ker_challenge_commitment = gadd(coeff_com,mkvec2(gel(dlog,1),gel(dlog,2)));
    //compute the even ideal for the challenge when needed
    if (need_even_chall) {
      kernel_to_ideal_two_f(gel(dlog,3),gel(dlog,4),two_tors_height-1,&I_even);
    }


    #ifndef NDEBUG

    odd_isogeny phi_test;
    GEN I = kernel_to_ideal_O0_T(coeff_com);
    proj ker = coeff_to_E0(gel(coeff_com,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff_com,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(I)));

    phi_test.kernel_plus = ker;
    phi_test.kernel_minus = ker_twist;
    phi_test.deg_plus = deg;
    phi_test.deg_minus = deg_twist;

    proj A = global_setup.E0, H = {fp2_1, fp2_0};

    if (need_even_commit ) {
      two_walk_long phi_two_test;
      init_trivial_two_walk_long(&phi_two_test);
      push_two_walk_long_through_odd_isogeny(&phi_two_test, &phi_com_even, &phi_test, &A);
      A = phi_two_test.phi[0].A;
      eval_walk_long_mult(&phi_two_test, &A, &H, 1);
    }
    else {
      eval_mult(&A,&phi_test,&H,1);
    }


    assert(mont_equal(&A, &comp_sigma->E_com));


    // // why doesnt the following work?
    // A = global_setup.E0; H = (proj){fp2_1, fp2_0};
    // phi_test = push_odd_isogeny_through_two_walk_long(&phi_test, &A, &phi_com_even);
    // eval(&A, &phi_test, &H);
    // assert(mont_equal(&A, &comp_sigma->E_com));


    #endif


    #ifndef NDEBUG
    {
    odd_isogeny phi_test;
    GEN I = kernel_to_ideal_O0_T(coeff_ker_challenge_commitment);
    proj ker = coeff_to_E0(gel(coeff_ker_challenge_commitment,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff_ker_challenge_commitment,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(I)));

    phi_test.kernel_plus = ker;
    phi_test.kernel_minus = ker_twist;
    phi_test.deg_plus = deg;
    phi_test.deg_minus = deg_twist;

    proj A = global_setup.E0, H = {fp2_1, fp2_0};
    if (need_even_commit) {
      two_walk_long phi_two_test;
      init_trivial_two_walk_long(&phi_two_test);
      push_two_walk_long_through_odd_isogeny(&phi_two_test, &phi_com_even, &phi_test, &A);
      A = phi_two_test.phi[0].A;
      eval_walk_long_mult(&phi_two_test, &A, &H, 1);
    }
    else {
      eval_mult(&A,&phi_test,&H,1);
      eval_walk_long_mult(&phi_chall_even, &A, &H, 1);
    }

    if (need_even_chall) {

    }

    assert(mont_equal(&A, &E_cha));
    }
    #endif
    // TOC(tcom,"com + chall");

    // uint64_t zip[signing_length_two_tors_height_step];
    two_walk_long sigma;
    init_trivial_two_walk_long(&sigma);
    // response(&sigma, comp_sigma->zip, coeff_ker_challenge_commitment, sk, basis_two, &E_cha);
    bool test;

    response_new(&sigma, comp_sigma->zip, coeff_ker_challenge_commitment, I_even, sk, basis_two, &E_cha,&test);

    if (test) {
      assert(simple_check_signature_new(&sigma,comp_sigma->zip, pk, &E_cha));
    }
    else {
      printf("computation failed due to norm equation failure \n");
    }

    free_two_walk_long(&sigma);

    avma = ltop;
}





bool verif(compressed_signature *comp_sigma, const public_key *pk,const uintbig *m){
    proj A_chall = comp_sigma->E_com;
    two_walk_long phi_chall_even;
    init_trivial_two_walk_long(&phi_chall_even);
    challenge(&A_chall,&phi_chall_even, m, &A_chall, NULL, NULL, NULL, NULL);
    proj A_check=pk->E;
    two_walk walk_check[signing_length_two_tors_height_step];
    // for (int i=0;i< 30;i++){
    //   printf("%ld ",comp_sigma->zip[i]);
    // }
    // printf("\n");

    decompress(walk_check,&A_check,comp_sigma->zip,signing_length_two_tors_height_step,last_step_length);

    proj j1,j2;
    jinv256(&j1,&A_chall);
    normalize_proj(&j1);
    // printf("%ld %ld \n",j1.x.re.x.c[0],j1.x.re.x.c[1]);
    jinv256(&j2,&A_check);
    normalize_proj(&j2);
    // printf("%ld %ld \n",j2.x.re.x.c[0],j2.x.re.x.c[1]);
    return mont_equal(&j1,&j2);
}

bool verif_new(compressed_signature *comp_sigma, const public_key *pk,const uintbig *m){
    proj A_chall = comp_sigma->E_com;
    two_walk_long phi_chall_even;
    init_trivial_two_walk_long(&phi_chall_even);
    challenge(&A_chall,&phi_chall_even, m, &A_chall, NULL, NULL, NULL, NULL);
    proj A_check=pk->E;
    two_walk walk_check[signing_length_two_tors_height_step];
    // for (int i=0;i< 30;i++){
    //   printf("%ld ",comp_sigma->zip[i]);
    // }
    // printf("\n");

    decompress_new(walk_check,&A_check,comp_sigma->zip,signing_length_two_tors_height_step,last_step_length);

    proj j1,j2;
    jinv256(&j1,&A_chall);
    normalize_proj(&j1);
    // printf("%ld %ld \n",j1.x.re.x.c[0],j1.x.re.x.c[1]);
    jinv256(&j2,&A_check);
    normalize_proj(&j2);
    // printf("%ld %ld \n",j2.x.re.x.c[0],j2.x.re.x.c[1]);
    return mont_equal(&j1,&j2);
}
