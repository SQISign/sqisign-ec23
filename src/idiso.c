#define _unused(x) ((void)(x))

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "printing.h"
#include "sqisign.h"


odd_isogeny trivial_odd_isogeny() {
    odd_isogeny phi;

    degree_one(&phi.deg_plus);
    degree_one(&phi.deg_minus);
    phi.kernel_plus.x = fp2_1;
    phi.kernel_plus.z = fp2_0;
    phi.kernel_minus.x = fp2_1;
    phi.kernel_minus.z = fp2_0;

    return phi;
}

special_isogeny trivial_special_isogeny() {
    special_isogeny phi;

    phi.source = global_setup.E0;
    phi.target = global_setup.E0;
    phi.middle = global_setup.E0;

    phi.phi2_set = true;
    phi.phi2_dual_set = true;

    phi.phi1 = trivial_odd_isogeny();
    phi.phi2 = trivial_odd_isogeny();
    phi.phi2_dual = trivial_odd_isogeny();

    return phi;
}



void random_two_walk(two_walk *phi){
    pari_sp ltop = avma;

    phi->A = global_setup.E0;

    phi->len = two_tors_height;

    const proj *basis = torsion_basis_two;
    uintbig x, y;
    proj P;

    // do {
        gentobig(&x, gen_1);
        gentobig(&y, randomi(powuu(2,two_tors_height)));

        xBIDIM(&phi->ker, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);

        gentobig(&x, powuu(2,two_tors_height-1));
        xMUL(&P, &global_setup.E0, &phi->ker, &x);
        assert(!mont_iszero(&P));
    // } while (!fp2_iszero(&P.x)); // checks that it is not the degree 2 endomorphism, when j = 1728

    gentobig(&x, gen_2);
    xMUL(&P, &global_setup.E0, &P, &x);
    assert(mont_iszero(&P));

    avma = ltop;
}


void free_two_walk_long(two_walk_long *phi) {
    if (phi->len) { free(phi->phi); phi->len = 0; }
}

void init_trivial_two_walk_long(two_walk_long *phi) {
    phi->phi = NULL;
    phi->len = 0;
}

static GEN alg_standard_to_O0(GEN elt) {
    return RgM_RgC_mul(global_setup.standard_to_O0, elt);
}


GEN kernel_to_ideal_gen_action(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {

    if (e < 1) {
        return mkcol4s(1,0,0,0);
    }

    pari_sp ltop = avma;

    GEN gelle = powuu(ell,e);

    GEN v1 = gmul(m1, v);
    GEN v2 = gmul(m2, v);
    GEN v3 = gmul(m3, v);
    GEN v4 = gmul(m4, v);

    GEN matsys = mkmat4(v1,v2,v3,v4);
    GEN ker = matkermod(matsys, gelle, NULL);

    GEN sol, sol_reduced;

    int dim_ker = lg(ker)-1;


    if (dim_ker == 0) {
        avma = ltop;
        return mkcol4s(1,0,0,0);
    }
    for (int i = 1; i <= dim_ker; ++i){
        sol = gel(ker,i);
        sol_reduced = gmodgs(sol,ell);
        if (!isexactzero(sol_reduced)) break;
    }

    return gerepilecopy(ltop, sol);
}


GEN kernel_to_ideal_action_O0(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {

    pari_sp ltop = avma;

    GEN gelle = powuu(ell,e);
    GEN gen_v = kernel_to_ideal_gen_action(v, m1, m2, m3, m4, ell, e);



    GEN generator = gmul(gel(gen_v,1),global_setup.O0_b1);
    generator = gadd(generator, gmul(gel(gen_v,2),global_setup.O0_b2));
    generator = gadd(generator, gmul(gel(gen_v,3),global_setup.O0_b3));
    generator = gadd(generator, gmul(gel(gen_v,4),global_setup.O0_b4));

    GEN ideal = lideal_create(global_setup.B, global_setup.O0, generator, gelle);

    return gerepilecopy(ltop, ideal);
}

// endo is an endomorphism expressed in the basis of the order whose action on the torsion correspond to m1 m2 m3 m4
// i.e. endo acts on the torsion as endo[0]*m1 + endo[1]*m2 + endo[2]*m3 + endo[3]*m4
GEN endo_to_kernel_action(GEN endo, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {
    pari_sp ltop = avma;



    GEN gelle = gpowgs(stoi(ell),e);

    GEN endo_m = gmul(gel(endo,1),m1);
    endo_m = gadd(endo_m, gmul(gel(endo,2),m2));
    endo_m = gadd(endo_m, gmul(gel(endo,3),m3));
    endo_m = gadd(endo_m, gmul(gel(endo,4),m4));



    GEN ker = matkermod(endo_m, gelle, NULL);

    GEN sol, sol_reduced;

    int dim_ker = lg(ker)-1;

    if (dim_ker == 0) {
        return mkcol2(gen_0,gen_0);
    }

    for (int i = 1; i <= dim_ker; ++i){
        sol = gel(ker,i);
        sol_reduced = gmodgs(sol,ell);
        if (!isexactzero(sol_reduced)) break;
        assert(i < dim_ker);
    }


    // remains to compute sol[1]*P1 + sol[2]*P2 where P1,P2 is a basis of the torsion

    return gerepilecopy(ltop, sol);
}

GEN ideal_to_kernel_action_O0(GEN I, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {
    pari_sp ltop = avma;

    GEN generator = lideal_generator(I);
    GEN endo = alg_standard_to_O0(generator);
    GEN sol = endo_to_kernel_action(endo, m1, m2, m3, m4, ell, e);
    // remains to compute sol[1]*P1 + sol[2]*P2 where P1,P2 is a basis of the torsion

    return gerepilecopy(ltop, sol);
}

void action_from_elle(GEN *m1, GEN *m2, GEN *m3, GEN *m4, long ell, long e) {
    pari_sp ltop = avma;

    GEN gelle = powuu(ell,e);

    *m1 = mkmat2(mkcol2s(1,0),mkcol2s(0,1));
    if (ell == 2) {
        *m2 = gmod(global_setup.action_two_2, gelle);
        *m3 = gmod(global_setup.action_two_3, gelle);
        *m4 = gmod(global_setup.action_two_4, gelle);

    }
    else {
        bool twist;
        unsigned long index = ell_to_index(ell, &twist);
        if (!twist) {
            *m2 = gmod(global_setup.action_2[index], gelle);
            *m3 = gmod(global_setup.action_3[index], gelle);
            *m4 = gmod(global_setup.action_4[index], gelle);
        }
        else {
            *m2 = gmod(global_setup.action_twist_2[index], gelle);
            *m3 = gmod(global_setup.action_twist_3[index], gelle);
            *m4 = gmod(global_setup.action_twist_4[index], gelle);
        }
    }

    gerepileall(ltop, 4, m1,m2,m3,m4);
}

GEN kernel_to_ideal_gen_O0_ell(GEN v, long ell, long *e) {
    pari_sp ltop = avma;
    long e_max = ell_to_e(ell);
    GEN gcd = ggcd(gel(v,1),gel(v,2));
    long e_diff = (isexactzero(gcd)) ? e_max : Z_lval(gcd, ell);

    if (e_diff >= e_max) {
        avma = ltop;
        *e = 0;
        return mkcol4s(1,0,0,0);
    }

    GEN generator;
    GEN m1,m2,m3,m4;
    *e = (e_diff < e_max) ? e_max - e_diff : 0;
    action_from_elle(&m1, &m2, &m3, &m4, ell, *e);
    v = gdiv(v, powuu(ell,e_diff)); // from ell^e_max torsion to ell^e torsion
    generator = kernel_to_ideal_gen_action(v, m1, m2, m3, m4, ell, *e);
    return gerepilecopy(ltop, generator);
}

GEN kernel_to_ideal_O0_ell(GEN v, long ell) {
    pari_sp ltop = avma;
    long e_max = ell_to_e(ell);
    GEN gcd = ggcd(gel(v,1),gel(v,2));
    long e_diff = (isexactzero(gcd)) ? e_max : Z_lval(gcd, ell);

    if (e_diff >= e_max) {
        avma = ltop;
        return lideal_create(global_setup.B, global_setup.O0, alg_scalar(global_setup.B,gen_1), gen_1);
    }

    GEN I;
    GEN m1,m2,m3,m4;
    long e = (e_diff < e_max) ? e_max - e_diff : 0;
    action_from_elle(&m1, &m2, &m3, &m4, ell, e);
    v = gdiv(v, powuu(ell,e_diff)); // from ell^e_max torsion to ell^e torsion
    I = kernel_to_ideal_action_O0(v, m1, m2, m3, m4, ell, e);
    return gerepilecopy(ltop, I);
}

// assume I is a cyclic ideal
GEN ideal_to_kernel_O0_ell(GEN I, long ell) {
    pari_sp ltop = avma;
    long e = Z_lval(lideal_norm(I), ell);
    long e_max = ell_to_e(ell);
    long e_diff = e_max - e;
    GEN v;
    GEN m1,m2,m3,m4;
    action_from_elle(&m1, &m2, &m3, &m4, ell, e);
    v = ideal_to_kernel_action_O0(I, m1, m2, m3, m4, ell, e);
    v = gmul(v, powuu(ell,e_diff)); // from ell^e torsion to ell^e_max torsion
    return gerepilecopy(ltop, v);
}

// assume I is a cyclic ideal
GEN endo_to_kernel_O0_ell(GEN alpha, long ell, long e) {
    pari_sp ltop = avma;
    long e_max = ell_to_e(ell);
    long e_diff = e_max - e;
    GEN v;
    GEN m1,m2,m3,m4;
    action_from_elle(&m1, &m2, &m3, &m4, ell, e);

    GEN endo = alg_standard_to_O0(alpha);

    v = endo_to_kernel_action(endo, m1, m2, m3, m4, ell, e);
    v = gmul(v, powuu(ell,e_diff)); // from ell^e torsion to ell^e_max torsion
    return gerepilecopy(ltop, v);
}



GEN ideal_to_kernel_O0_T(GEN I, GEN fact_norm) {
    pari_sp ltop = avma;
    long len = lg(gel(fact_norm,1)), ell, e_max, e;
    GEN coeff_1 = zerovec(on_curve_len), coeff_2 = zerovec(on_curve_len);
    GEN coeff_twist_1 = zerovec(on_twist_len), coeff_twist_2 = zerovec(on_twist_len);
    GEN alpha = lideal_generator(I);
    long index;
    bool twist;
    GEN v;

    for (int i = 1; i < len; ++i) {
        ell = itos(gel(gel(fact_norm,1),i));
        e_max = itos(gel(gel(fact_norm,2),i));
        index = ell_to_index(ell, &twist);
        e = Z_lval(lideal_norm(I), ell);
        v = endo_to_kernel_O0_ell(alpha, ell, e);
        if (!twist) {
            gel(coeff_1,index+1) = gel(v,1);
            gel(coeff_2,index+1) = gel(v,2);
        }
        else {
            gel(coeff_twist_1,index+1) = gel(v,1);
            gel(coeff_twist_2,index+1) = gel(v,2);
        }
    }

    return gerepilecopy(ltop, mkvec2(mkvec2(coeff_1,coeff_2),mkvec2(coeff_twist_1,coeff_twist_2)));
}

GEN torsion_crt_compose (GEN coeff, bool twist) {
    pari_sp ltop = avma;
    GEN p_primary = ((twist && curve_order_is_p_plus_one) || (!twist && !curve_order_is_p_plus_one)) ?
                        global_setup.gen_p_minus_primary :
                        global_setup.gen_p_plus_primary;

    GEN a = ZV_chinese(gel(coeff,1), p_primary, NULL);
    GEN b = ZV_chinese(gel(coeff,2), p_primary, NULL);

    return gerepilecopy(ltop, mkcol2(a,b));
}

GEN torsion_crt_decompose (GEN v, bool twist) {
    pari_sp ltop = avma;
    long len;
    GEN p_primary;
    if ((twist && curve_order_is_p_plus_one) || (!twist && !curve_order_is_p_plus_one)) {
        len = p_minus_len;
        p_primary = global_setup.gen_p_minus_primary;
    }
    else {
        len = p_plus_len;
        p_primary = global_setup.gen_p_plus_primary;

    }

    GEN coeff_1 = cgetg(len+1, t_VEC);
    GEN coeff_2 = cgetg(len+1, t_VEC);

    for (int i = 0; i < len; ++i) {
        gel(coeff_1,i+1) = gmod(gel(v,1), gel(p_primary,i+1));
        gel(coeff_2,i+1) = gmod(gel(v,2), gel(p_primary,i+1));
    }

    return gerepilecopy(ltop, mkvec2(coeff_1,coeff_2));
}

proj vec_to_E0(GEN v, bool twist) {
    proj pt;
    const proj *basis = (twist) ? torsion_basis_twist_sum : torsion_basis_sum;
    uintbig x, y;

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));

    // TODO: replace the following with a small multiplication and a combination of the form P + yQ (i.e., x = 1)
    xBIDIM(&pt, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);

    return pt;
}


void famat_to_degree(isog_degree *deg_curve, isog_degree *deg_twist, GEN f) {
    degree_one(deg_curve);
    degree_one(deg_twist);

    if (lg(f) == 1) return;

    long m = lg(gel(f,1)) - 1;
    long ell, e, index;
    bool twist;

    for (long i = 1; i <= m; ++i) {
        ell = itos_or_0(gel(gel(f,1),i));
        e = itos_or_0(gel(gel(f,2),i));
        index = ell_to_index(ell, &twist);
        if (!twist) {
            degree_set(deg_curve, index, e);
        }
        else {
            degree_set(deg_twist, index, e);
        }
    }
}

proj coeff_to_E0(GEN coeff, bool twist) {
    GEN v = torsion_crt_compose (coeff, twist);
    return vec_to_E0(v, twist);
}

odd_isogeny ideal_to_isogeny_O0_T(GEN I, GEN factorisation_norm) {
    odd_isogeny phi;
    GEN coeff = ideal_to_kernel_O0_T(I,factorisation_norm);
    proj ker = coeff_to_E0(gel(coeff,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, factorisation_norm);


    if (curve_order_is_p_plus_one) {
        phi.kernel_plus = ker;
        phi.kernel_minus = ker_twist;
        phi.deg_plus = deg;
        phi.deg_minus = deg_twist;
    }
    else {
        phi.kernel_plus = ker_twist;
        phi.kernel_minus = ker;
        phi.deg_plus = deg_twist;
        phi.deg_minus = deg;

    }

    return phi;
}



two_walk ideal_to_isogeny_O0_two(GEN I) {
    pari_sp ltop = avma;
    two_walk phi;
    phi.A = global_setup.E0;

    GEN v = ideal_to_kernel_O0_ell(I, 2);

    const proj *basis = torsion_basis_two;
    uintbig x, y;

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));

    // TODO: replace the following with a small multiplication and a combination of the form P + yQ (i.e., x = 1)
    xBIDIM(&phi.ker, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);



    phi.len = Z_lval(lideal_norm(I), 2);

    avma = ltop;
    return phi;
}




GEN kernel_to_ideal_O0_T(GEN coeff) {
    pari_sp ltop = avma;

    GEN I, v, N = gen_1;
    long ell, e, e_max;

    GEN x, ideal_generator_decomposed = mkvec4(
        cgetg(on_curve_len + on_twist_len + 1, t_VEC),
        cgetg(on_curve_len + on_twist_len + 1, t_VEC),
        cgetg(on_curve_len + on_twist_len + 1, t_VEC),
        cgetg(on_curve_len + on_twist_len + 1, t_VEC));

    for (int i = 0; i < on_curve_len; ++i) {

        ell = on_curve_fact[i];
        e_max = on_curve_mult[i];
        v = mkcol2(gel(gel(gel(coeff,1),1),i+1), gel(gel(gel(coeff,1),2),i+1));
        x = kernel_to_ideal_gen_O0_ell(v, ell, &e);
        N = gmul(N, powuu(ell,e));
        gel(gel(ideal_generator_decomposed,1),i+1) = gel(x,1);
        gel(gel(ideal_generator_decomposed,2),i+1) = gel(x,2);
        gel(gel(ideal_generator_decomposed,3),i+1) = gel(x,3);
        gel(gel(ideal_generator_decomposed,4),i+1) = gel(x,4);
    }

    for (int i = 0; i < on_twist_len; ++i) {
        ell = on_twist_fact[i];
        e_max = on_twist_mult[i];
        v = mkcol2(gel(gel(gel(coeff,2),1),i+1), gel(gel(gel(coeff,2),2),i+1));
        x = kernel_to_ideal_gen_O0_ell(v, ell, &e);
        N = gmul(N, powuu(ell,e));
        gel(gel(ideal_generator_decomposed,1),on_curve_len+i+1) = gel(x,1);
        gel(gel(ideal_generator_decomposed,2),on_curve_len+i+1) = gel(x,2);
        gel(gel(ideal_generator_decomposed,3),on_curve_len+i+1) = gel(x,3);
        gel(gel(ideal_generator_decomposed,4),on_curve_len+i+1) = gel(x,4);
    }

    GEN p_primary = shallowconcat(global_setup.gen_p_plus_primary, global_setup.gen_p_minus_primary);
    if (!curve_order_is_p_plus_one) {
        p_primary = shallowconcat(global_setup.gen_p_minus_primary, global_setup.gen_p_plus_primary);
    }

    GEN a = ZV_chinese(gel(ideal_generator_decomposed,1), p_primary, NULL);
    GEN b = ZV_chinese(gel(ideal_generator_decomposed,2), p_primary, NULL);
    GEN c = ZV_chinese(gel(ideal_generator_decomposed,3), p_primary, NULL);
    GEN d = ZV_chinese(gel(ideal_generator_decomposed,4), p_primary, NULL);

    GEN generator = gmul(a,global_setup.O0_b1);
    generator = gadd(generator, gmul(b,global_setup.O0_b2));
    generator = gadd(generator, gmul(c,global_setup.O0_b3));
    generator = gadd(generator, gmul(d,global_setup.O0_b4));

    I = lideal_create(global_setup.B, global_setup.O0, generator, N);

    return gerepilecopy(ltop, I);
}



// Evaluate special isogeny phi : source -> target at point P.
// sets P to the image point
proj eval_special(proj *range, special_isogeny *phi, const proj *P) {
    proj A = phi->source;
    proj Q = *P;
    isomorphism isom;

    eval(&A, &phi->phi1, &Q);

    if (!phi->phi2_set) {
        assert(phi->phi2_dual_set);
        phi->phi2 = phi->phi2_dual;
        phi->middle = phi->target;
        dual(&phi->middle, &phi->phi2);
        phi->phi2_set = true;
    }

    #ifndef NDEBUG
    proj j1,j2;
    jinv256(&j1, &A);
    jinv256(&j2, &phi->middle);
    // if (!mont_equal(&j1,&j2)) fprintf(stderr,"ERROR in eval_special: the two isogenies do not meet in the middle\n");
    assert(mont_equal(&j1,&j2));
    #endif

    mont_isom(&isom, &A, &phi->middle);
    mont_isom_apply(&isom, &Q);
    A = phi->middle;
    eval(&A, &phi->phi2, &Q);

    *range = A;

    return Q;
}









// Evaluate special isogeny phi : source -> target at point P.
// sets P to the images points
void eval_special_mult(proj *range, special_isogeny *phi, proj *P, long len) {
    proj A = phi->source;
    isomorphism isom;

    eval_mult(&A, &phi->phi1, P, len);

    if (!phi->phi2_set) {
        assert(phi->phi2_dual_set);
        phi->phi2 = phi->phi2_dual;
        phi->middle = phi->target;
        dual(&phi->middle, &phi->phi2);
        phi->phi2_set = true;
    }

    #ifndef NDEBUG
    proj j1,j2;
    jinv256(&j1, &A);
    jinv256(&j2, &phi->middle);
    // if (!mont_equal(&j1,&j2)) fprintf(stderr,"ERROR in eval_special: the two isogenies do not meet in the middle\n");
    assert(mont_equal(&j1,&j2));
    #endif

    mont_isom(&isom, &A, &phi->middle);

    for (int i = 0; i < len; ++i) {
        mont_isom_apply(&isom, P+i);
    }

    A = phi->middle;
    eval_mult(&A, &phi->phi2, P,len);

    *range = A;
}



void two_walk_stol(two_walk_long *phil, const two_walk *phi) {
    if (phi->len  == 0) { free_two_walk_long(phil); init_trivial_two_walk_long(phil); return; }
    two_walk_long res;
    res.len = 1;
    res.phi = malloc(sizeof(two_walk)*(1));
    res.phi[0] = *phi;
    free_two_walk_long(phil);
    *phil = res;
}

void two_walk_composition_ll(two_walk_long *phi, const two_walk_long *phi2, const two_walk_long *phi1){
    two_walk_long res;
    res.len = phi1->len + phi2->len;
    res.phi = malloc(sizeof(two_walk)*(res.len));
    for (int i = 0; i < phi1->len; ++i) {
        res.phi[i] = phi1->phi[i];
    }
    for (int i = 0; i < phi2->len; ++i) {
        res.phi[i+phi1->len] = phi2->phi[i];
    }
    free_two_walk_long(phi);
    *phi = res;
}

void copy_two_walk_long(two_walk_long *copy, const two_walk_long *phi) {
    two_walk_long triv;
    init_trivial_two_walk_long(&triv);
    two_walk_composition_ll(copy, phi, &triv);
}

void two_walk_composition_ls(two_walk_long *phi, const two_walk_long *phi2, const two_walk *phi1) {
    if (phi1->len == 0) { copy_two_walk_long(phi, phi2); return; }
    two_walk_long phi1l;
    init_trivial_two_walk_long(&phi1l);
    two_walk_stol(&phi1l, phi1);
    two_walk_composition_ll(phi, phi2, &phi1l);
    free_two_walk_long(&phi1l);
}

void two_walk_composition_sl(two_walk_long *phi, const two_walk *phi2, const two_walk_long *phi1) {
    if (phi2->len == 0) { copy_two_walk_long(phi, phi1); return; }
    two_walk_long phi2l;
    init_trivial_two_walk_long(&phi2l);
    two_walk_stol(&phi2l, phi2);
    two_walk_composition_ll(phi, &phi2l, phi1);
    free_two_walk_long(&phi2l);
}

void two_walk_composition_ss(two_walk_long *phi, const two_walk *phi2, const two_walk *phi1) {
    if (phi2->len == 0) {
        if (phi1->len == 0) { free_two_walk_long(phi); init_trivial_two_walk_long(phi); }
        else { two_walk_stol(phi, phi1); return; }
    }
    if (phi1->len == 0) { two_walk_stol(phi, phi2); return; }
    two_walk_long phi1l;
    init_trivial_two_walk_long(&phi1l);
    two_walk_stol(&phi1l, phi1);
    two_walk_long phi2l;
    init_trivial_two_walk_long(&phi2l);
    two_walk_stol(&phi2l, phi2);
    two_walk_composition_ll(phi, &phi2l, &phi1l);
    free_two_walk_long(&phi1l);
    free_two_walk_long(&phi2l);
}

odd_isogeny push_odd_isogeny_through_two_walk(const odd_isogeny *phi_odd, proj *phi_odd_source, const two_walk *phi_two) {
    odd_isogeny phi_odd_isom = *phi_odd;

    if (phi_two->len == 0) { return phi_odd_isom;}
    two_walk phi_two_isom = *phi_two;
    isomorphism isom1, isom2;

    #ifndef NDEBUG
    proj j1,j2;
    jinv256(&j1, phi_odd_source);
    jinv256(&j2, &phi_two->A);
    assert(mont_equal(&j1,&j2));
    #endif

    if (!mont_equal(phi_odd_source,&phi_two_isom.A)) {

        #ifndef NDEBUG
        // Check the j-invariant is not 0 or 1728 (otherwise the isomorphism is ambiguous)
        proj j1728;
        fp2_set(&j1728.x,1728);
        fp2_set(&j1728.z,256);
        assert(!mont_equal(&j1,&j1728) && !mont_iszero(&j1));
        #endif

        mont_isom(&isom2, phi_odd_source, &phi_two_isom.A);
        mont_isom_apply(&isom2, &phi_odd_isom.kernel_plus);
        mont_isom_apply(&isom2, &phi_odd_isom.kernel_minus);

    }


    eval_walk_isom(&isom1, &phi_two_isom, NULL, NULL, &phi_two_isom, NULL);

    mont_isom_apply(&isom1, &phi_odd_isom.kernel_plus);
    mont_isom_apply(&isom1, &phi_odd_isom.kernel_minus);


    #ifndef NDEBUG
    uintbig k;
    proj P = phi_two_isom.ker;
    gentobig(&k, powuu(2,phi_two_isom.len-1));
    xMUL(&P, &phi_two_isom.A, &P, &k);
    assert(!mont_iszero(&P));
    assert(!fp2_iszero(&P.x));
    gentobig(&k, powuu(2,1));
    xMUL(&P, &phi_two_isom.A, &P, &k);
    assert(mont_iszero(&P));

    jinv256(&j1, phi_odd_source);
    jinv256(&j2, &phi_two_isom.A);
    assert(mont_equal(&j1,&j2));
    #endif


    odd_isogeny res = phi_odd_isom;
    eval_walk(&phi_two_isom, phi_odd_source, &res.kernel_plus);
    eval_walk(&phi_two_isom, phi_odd_source, &res.kernel_minus);
    return res;
}

odd_isogeny push_odd_isogeny_through_two_walk_long(const odd_isogeny *phi_odd, proj *phi_odd_source, const two_walk_long *phi_two) {
    odd_isogeny res = *phi_odd;
    for (int i = 0; i < phi_two->len; ++i) {
        res = push_odd_isogeny_through_two_walk(&res, phi_odd_source, &phi_two->phi[i]);
    }
    return res;
}

two_walk push_two_walk_through_odd_isogeny(const two_walk *phi_two, const odd_isogeny *phi_odd, const proj *phi_odd_source) {
    two_walk res = *phi_two;
    isomorphism isom;

    mont_isom(&isom, &phi_two->A, phi_odd_source);
    mont_isom_apply(&isom, &res.ker);

    res.A = *phi_odd_source;

    eval(&res.A, phi_odd, &res.ker);

    return res;
}

void push_two_walk_long_through_odd_isogeny(two_walk_long *phi, const two_walk_long *phi_two, const odd_isogeny *phi_odd, const proj *phi_odd_source) {
    two_walk_long res;

    init_trivial_two_walk_long(&res);
    copy_two_walk_long(&res, phi_two);
    odd_isogeny phi_odd_i = *phi_odd;
    proj phi_odd_source_i = *phi_odd_source;


    for (int i = 0; i < phi_two->len; ++i) {
        res.phi[i] = push_two_walk_through_odd_isogeny(&res.phi[i], &phi_odd_i, &phi_odd_source_i);
        if (i+1 < phi_two->len) {
            #ifndef NDEBUG
            proj j1,j2;
            jinv256(&j1, &phi_odd_source_i);
            jinv256(&j2, &phi_two->phi[i].A);
            assert(mont_equal(&j1,&j2));
            #endif
            phi_odd_i = push_odd_isogeny_through_two_walk(&phi_odd_i, &phi_odd_source_i, &phi_two->phi[i]);
        }
    }

    free_two_walk_long(phi);
    copy_two_walk_long(phi, &res);
}

two_walk push_two_walk_through_special_isogeny(const two_walk *phi_two, special_isogeny *phi_special) {
    two_walk res = *phi_two;
    proj E0 = phi_special->source;
    isomorphism isom;

    mont_isom(&isom, &phi_two->A, &E0);
    mont_isom_apply(&isom, &res.ker);

    res.ker = eval_special(&res.A, phi_special, &res.ker);

    return res;
}


special_isogeny special_ideal_to_isogeny(GEN J, GEN I, const two_walk_long *phi_I) {
    pari_sp ltop = avma;
    special_isogeny phi;
    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);
    GEN H1 = lideal_create(A, order, lideal_generator(J), ggcd(global_setup.gen_odd_torsion,lideal_norm(J)));
    GEN beta = lideal_isom(I, J); // I*beta = J
    GEN alpha = gmul(alg_conj(A,beta), lideal_norm(I)); // J = chi_I(alpha) = I conj(alpha)/n(I)
    GEN H2 = lideal_create(A, order, alpha, gdiv(lideal_norm(J),lideal_norm(H1)));

    GEN fm = famat_mul(global_setup.gen_p_plus_fact,global_setup.gen_p_minus_fact);
    GEN fm_norm1 = famat_Z_gcd(fm, lideal_norm(H1));
    phi.phi1 = ideal_to_isogeny_O0_T(H1, fm_norm1);
    phi.source = global_setup.E0;

    GEN fm_norm2 = famat_Z_gcd(fm, lideal_norm(H2));
    phi.phi2_dual = ideal_to_isogeny_O0_T(H2, fm_norm2);
    phi.target = global_setup.E0;

    // push phi->phi2 though phi_I
    phi.phi2_dual = push_odd_isogeny_through_two_walk_long(&phi.phi2_dual, &phi.target, phi_I);

    #ifndef NDEBUG
    odd_isogeny X = phi.phi1,Y = phi.phi2_dual;
    proj E1 = phi.source, E2 = phi.target;
    dual(&E1, &X);
    dual(&E2, &Y);
    assert(mont_equal(&E1,&E2));
    #endif


    phi.phi2_dual_set = true;
    phi.phi2_set = false;
    avma = ltop;
    return phi;
}

void eval_walk_long_mult(const two_walk_long *phi, proj *B, proj *P, long cardinality) {
    isomorphism isom;

    two_walk dummy;
    for (int i = 0; i < phi->len; ++i) {
        mont_isom(&isom, B, &phi->phi[i].A);
        for (int i = 0; i < cardinality; ++i) {
            mont_isom_apply(&isom, &P[i]);
        }
        eval_walk_isom_mult(&isom, &dummy, B, &phi->phi[i], P, cardinality);
    }
}


GEN dual_coeff(GEN coeff, isog_degree deg_plus, isog_degree deg_minus) {
    pari_sp ltop = avma;
    GEN coeff_dual = coeff;

    isog_degree deg_curve, deg_twist;
    if (curve_order_is_p_plus_one){
        deg_curve = deg_plus;
        deg_twist = deg_minus;
    }
    else {
        deg_curve = deg_minus;
        deg_twist = deg_plus;
    }

    for (int i = 0; i < on_curve_len; ++i) {
        long ell = on_curve_fact[i];
        long e = on_curve_mult[i];
        GEN c1 = gel(gel(gel(coeff,1), 1), i+1);
        GEN c2 = gel(gel(gel(coeff,1), 2), i+1);
        long d = degree_get(deg_curve, i);
        if (Z_lval(c1, ell) == e-d) {
            gel(gel(gel(coeff_dual,1), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,1), 2), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
        }
        else if (Z_lval(c2, ell) == e-d) {
            gel(gel(gel(coeff_dual,1), 1), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
            gel(gel(gel(coeff_dual,1), 2), i+1) = gen_0;
        }
        else {
            assert(isexactzero(c1) && isexactzero(c2));
            gel(gel(gel(coeff_dual,1), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,1), 2), i+1) = gen_0;
        }

    }


    for (int i = 0; i < on_twist_len; ++i) {
        long ell = on_twist_fact[i];
        long e = on_twist_mult[i];
        GEN c1 = gel(gel(gel(coeff,2), 1), i+1);
        GEN c2 = gel(gel(gel(coeff,2), 2), i+1);
        long d = degree_get(deg_twist, i);
        if (Z_lval(c1, ell) == e-d) {
            gel(gel(gel(coeff_dual,2), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,2), 2), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
        }
        else if (Z_lval(c2, ell) == e-d) {
            gel(gel(gel(coeff_dual,2), 1), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
            gel(gel(gel(coeff_dual,2), 2), i+1) = gen_0;
        }
        else {
            assert(isexactzero(c1) && isexactzero(c2));
            gel(gel(gel(coeff_dual,2), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,2), 2), i+1) = gen_0;
        }
    }
    return gerepilecopy(ltop, coeff_dual);
}



// PRF to generate points
static void hash(proj *P, int i) {
  uintbig_set(&P->x.re.x, 3 * i + 13);
  uintbig_set(&P->z.re.x, 5 * i * i + 17);
  uintbig_set(&P->x.im.x, 7 * i * i * i + 19);
  uintbig_set(&P->z.im.x, 11 * i * i * i + 23);
}

//compute deterministically a second point Q of order 2^two_tors_height
//so that P,[2^(two_tors_height-f)] Q is a basis of the 2^f torsion.

void deterministic_second_point(proj *Q,const proj* P, proj* A,long f) {
  bool oncurve = class_mod_4 == 3;
  long cnt = 1;
  proj P2 = *P;
  // normalize_proj(A);
  // normalize_proj(P);
  for (int i = 1; i < f; i++)
    xDBL(&P2, A, &P2);
  assert(!mont_iszero(&P2));
  while (true) {
    hash(Q, cnt++);
    if (is_on_curve(Q, A) != oncurve)
      continue;
    // multiply by cofactor
    xMUL(Q, A, Q, &p_even_cofactor);
    // check it has maximal order
    proj Q2 = *Q;
    for (int i = 1; i < two_tors_height; i++)
      xDBL(&Q2, A, &Q2);
    if (!mont_iszero(&Q2) && !mont_equal(&Q2, &P2))
      break;
  }

}

//alternative method to compute the ideal_to_isogeny algorithm
//phi is the output, an isogeny of degree 2f
//I_two is an idea of norm 2^*. Phi_I is implicitly defined as the correpsonding isogeny.
//I is equivalent to I_two of smaller norm
//phi_I_basis  is the image through phi_I of a basis of the odd_torsion
//L_kernel is the kernel of the dual of the last 2f isogeny composing phi_I.
// K is the ideal that needs to be translated
//the result is the compressed value that will be used to compress the representation of the result isogeny
void ideal_to_isogeny_small_two_f(two_walk *phi, uintbig *res, eichler_package *eich ,  GEN I, GEN I_two,
   GEN K, GEN L , proj* phi_I_basis, proj *base_curve, long f,bool special_fixed, bool special_bad_norm){

  proj L_kernel = phi_I_basis[6];
  GEN A = lideal_algebra(I);
  GEN twof = gpowgs(gen_2,f);

  GEN order = lideal_order(I);
  assert(!lideal_equals(K,L));

  GEN fm = famat_sqr(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact));
  GEN gamma;
  //this is the first step
  if (special_fixed) {
    eichler_norm_equation_special_smooth_small_curve_fixed(eich,I, L, fm, gen_2);
  }
  //in case we need to  use  the special  extremal order
  else if (special_bad_norm) {
    eichler_norm_equation_special_smooth_small_curve(eich,I, L, fm, gen_2);
  }
  //
  else {
    //this is VERY ad hoc
    GEN remove_fact = mkmat2(mkcol4s(6983,4019,3691,4283),mkcol4s(2,2,2,2));
    fm = famat_div(fm,remove_fact);
    remove_fact = mkmat2(mkcols(2713),mkcols(2));
    fm = famat_div(fm,remove_fact);
    eichler_norm_equation_special_smooth(eich,I, L, fm, gen_2);
  }
  GEN gamma_wit;
  if (eich->delta) {
     gamma_wit= eich->gamma;
    assert(lideal_equals(eich->J,lideal_mul(eich->I,alg_conj(A,eich->delta) )));

    // this should only happen when we used the other extremal orders
    if (!eich->gamma ){
        // old code
        // gamma = lideal_isom(eich->J,I_two); // I_two = eich.J * gamma;

        // new code
        gamma = lideal_isom(I,I_two); // OR(I_two) = bar(gamma) * OR(I) * gamma / norm = bar(gamma) * OR(eich->I) * gamma / norm;
        gamma = algmul(A,eich->delta,gamma); // OR(eich->J) = bar(eich->delta) * OR(eich->I) * eich->delta, so
                                            // we get OR(I_two) = bar(gamma) * OR(eich->J) * gamma / norm;


    }
    //this should happen at the first step only when gamma has not been initialized already and so it init at 1
    else if (gequal( algnorm(A,eich->gamma,0),gen_1)) {
       gamma = lideal_isom(eich->J,I_two);
    }
    else {
      gamma = algmul(A,eich->delta,eich->gamma);
    }
    assert(lideal_isom(I,I_two));
    eich->gamma = gamma;
    assert(gamma);
    GEN beta = algmul(A,eich->beta,gamma);
    beta = algmul(A,alg_conj(A,gamma),beta);
    GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,gamma,0)),gen_0,gen_0,gen_0);
    beta = algmul(A,beta,overnorm);
    assert(alglatcontains(A,lideal_right_order(I_two),beta,NULL));

    GEN gen_two =lideal_generator_coprime(I_two,algnorm(A,beta,0));
    GEN H1 = lideal_create(A,order,algmul(A, gen_two ,beta), ggcd(global_setup.gen_odd_torsion,algnorm(A,beta,0)) );
    GEN H2 = lideal_create(A,order,algmul(A, gen_two ,alg_conj(A,beta) ), gdiv(algnorm(A,beta,0),lideal_norm(H1))  );
    assert(lideal_isom( lideal_inter(I_two,H1),lideal_inter(I_two,H2) ));
    //computing psi1 corresponding to H1
    odd_isogeny psi_1; //= ideal_to_isogeny_O0_T(H1_odd, famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H1_odd)));

    GEN factorisation_norm1 = famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H1));
    GEN coeff1 = ideal_to_kernel_O0_T(H1,factorisation_norm1);

    GEN v1_plus = torsion_crt_compose(gel(coeff1,1), false);
    GEN v1_minus = torsion_crt_compose(gel(coeff1,2), true);

    proj ker1_plus, ker1_minus;

    uintbig x1, y1;
    gentobig(&x1, gel(v1_plus,1));
    gentobig(&y1, gel(v1_plus,2));
    xBIDIM(&ker1_plus, base_curve, &phi_I_basis[0], &x1, &phi_I_basis[1], &y1, &phi_I_basis[2]);
    gentobig(&x1, gel(v1_minus,1));
    gentobig(&y1, gel(v1_minus,2));
    xBIDIM(&ker1_minus, base_curve, &phi_I_basis[3], &x1, &phi_I_basis[4], &y1, &phi_I_basis[5]);

    isog_degree deg1_plus, deg1_minus;
    famat_to_degree(&deg1_plus, &deg1_minus, factorisation_norm1);
    psi_1.kernel_plus = ker1_plus;
    psi_1.kernel_minus = ker1_minus;
    psi_1.deg_plus = deg1_plus;
    psi_1.deg_minus = deg1_minus;

    //computing psi2 corresponding to H2
    odd_isogeny psi_2; //= ideal_to_isogeny_O0_T(H1_odd, famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H1_odd)));

    GEN factorisation_norm2 = famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H2));
    GEN coeff2 = ideal_to_kernel_O0_T(H2,factorisation_norm2);

    GEN v2_plus = torsion_crt_compose(gel(coeff2,1), false);
    GEN v2_minus = torsion_crt_compose(gel(coeff2,2), true);

    proj ker2_plus, ker2_minus;
    uintbig x2, y2;
    gentobig(&x2, gel(v2_plus,1));
    gentobig(&y2, gel(v2_plus,2));
    xBIDIM(&ker2_plus, base_curve, &phi_I_basis[0], &x2, &phi_I_basis[1], &y2, &phi_I_basis[2]);
    gentobig(&x2, gel(v2_minus,1));
    gentobig(&y2, gel(v2_minus,2));
    xBIDIM(&ker2_minus, base_curve, &phi_I_basis[3], &x2, &phi_I_basis[4], &y2, &phi_I_basis[5]);

    isog_degree deg2_plus, deg2_minus;
    famat_to_degree(&deg2_plus, &deg2_minus, factorisation_norm2);

    psi_2.kernel_plus = ker2_plus;
    psi_2.kernel_minus = ker2_minus;
    psi_2.deg_plus = deg2_plus;
    psi_2.deg_minus = deg2_minus;

    //computing the required torsion points
    proj Kd;
    deterministic_second_point(&Kd,&L_kernel,base_curve,f);

    for (int i = 0; i <two_tors_height - f ; i++) {
      xDBL(&Kd,base_curve,&Kd);
    }


    proj P3,P4;
    xBILIFT(&P3,&P4,&L_kernel,&Kd,base_curve);

    //pushing the basis through the two isogenies of the  endomorphism
    proj P1[2];
    P1[0] = L_kernel;
    P1[1] = Kd;
    proj A1 = *base_curve;
    eval_mult(&A1,&psi_1,P1,2);


    proj pt[3];
    pt[0] = L_kernel;
    pt[1] = Kd;
    pt[2] = P3;
    proj A2 = *base_curve;
    assert(!mont_equal(&A2,&A1));

    eval_mult(&A2,&psi_2,pt,3);
    #ifndef NDEBUG
    proj j1,j2;
    jinv256(&j1,&A1);
    normalize_proj(&j1);
    jinv256(&j2,&A2);
    normalize_proj(&j2);
    assert(mont_equal(&j1,&j2));
    #endif
    isomorphism isomo;
    mont_isom(&isomo,&A1,&A2);
    mont_isom_apply(&isomo,&P1[0]);
    mont_isom_apply(&isomo,&P1[1]);

    //compute the DLPs
    uintbig x3,x4;
    proj2 xyP1,xyP2,xyQ,xyR,xyP12;
    xtoxy(&xyP1,&A2,&pt[0]);
    xtoxy(&xyP2,&A2,&pt[1]);
    xtoxy(&xyQ,&A2,&P1[0]);
    xtoxy(&xyR,&A2,&P1[1]);
    xyNEG(&xyR,&xyR);
    xyNEG(&xyQ,&xyQ);
    proj2 inter;
    xyADD(&inter,&A2,&xyP1,&xyP2);
    proj test_sign;
    xytox(&test_sign,&inter);
    if(!mont_equal(&test_sign,&pt[2])) {
      xyNEG(&xyP2,&xyP2);
    }
    xyADD(&xyP12,&A2,&xyP1,&xyP2);
    proj2 xyPQ1,xyPQ2,xyPQ12,xyPR1,xyPR2,xyPR12;
    xyADD(&xyPQ1,&A2,&xyP1,&xyQ);
    xyADD(&xyPQ2,&A2,&xyP2,&xyQ);
    xyADD(&xyPQ12,&A2,&xyP12,&xyQ);
    xyADD(&xyPR1,&A2,&xyP1,&xyR);
    xyADD(&xyPR2,&A2,&xyP2,&xyR);
    xyADD(&xyPR12,&A2,&xyP12,&xyR);
    proj PQ1,PQ2,PQ12,PR1,PR2,PR12;
    xytox(&PQ1,&xyPQ1);
    xytox(&PQ2,&xyPQ2);
    xytox(&PQ12,&xyPQ12);
    xytox(&PR1,&xyPR1);
    xytox(&PR2,&xyPR2);
    xytox(&PR12,&xyPR12);
    bool dlpt1 = mont_bidim_two_DLP(&x1,&x2,&A2,&P1[0],&pt[0],&pt[1],&pt[2],&PQ1,&PQ2,&PQ12,f);
    assert(dlpt1);
    dlpt1 = mont_bidim_two_DLP(&x3,&x4,&A2,&P1[1],&pt[0],&pt[1],&pt[2],&PR1,&PR2,&PR12,f);
    assert(dlpt1);
    _unused(dlpt1);

    #ifndef NDEBUG


    proj Pt1,Pt2;
    xBIDIM(&Pt1,&A2,&pt[0],&x1,&pt[1],&x2,&pt[2]);
    assert(mont_equal(&Pt1,&P1[0]));
    xBIDIM(&Pt2,&A2,&pt[0],&x3,&pt[1],&x4,&pt[2]);
    assert(mont_equal(&Pt2,&P1[1]));

    #endif

    //compute the coefficients x1_gen,x2_gen of the kernel in the basis L_kernel,Kd
    GEN x1_gen = bigtogen(&x1);
    GEN x2_gen = bigtogen(&x2);
    GEN x3_gen = bigtogen(&x3);
    GEN x4_gen = bigtogen(&x4);

    GEN H2_norm = gmod(lideal_norm(H2),twof);

    x1_gen = gmod(gmul(x1_gen,H2_norm),twof);
    x2_gen = gmod(gmul(x2_gen,H2_norm),twof);
    x3_gen = gmod(gmul(x3_gen,H2_norm),twof);
    x4_gen = gmod(gmul(x4_gen,H2_norm),twof);
    //use the trace to verify the signs
    GEN trace = gmod(gmul(gen_2,gel(beta,1)),twof);
    if (! gequal(trace,gmod( gadd(x1_gen,x4_gen),twof)  )) {
      x1_gen =gmod(gneg(x1_gen),twof);
      x2_gen =gmod(gneg(x2_gen),twof);

      if (! gequal(trace,gmod( gadd(x1_gen,x4_gen),twof)  )) {
        x4_gen =gmod(gneg(x4_gen),twof);
        x3_gen =gmod(gneg(x3_gen),twof);
        if (! gequal(trace,gmod( gadd(x1_gen,x4_gen),twof)  )) {
          x1_gen =gmod(gneg(x1_gen),twof);
          x2_gen =gmod(gneg(x2_gen),twof);
        }
      }
    }
    assert(gequal(trace,gmod( gadd(x1_gen,x4_gen),twof)));

    GEN coeffs = find_coeff(K,L,eich->beta,eich->delta,lideal_right_order(eich->J));
    assert(coeffs);

    x1_gen = gmod(gmul(x1_gen,gel(coeffs,2)) ,twof);
    x2_gen = gmod(gmul(x2_gen,gel(coeffs,2)) ,twof);
    x1_gen = gmod(gadd(x1_gen,gel(coeffs,1)) ,twof);
    gentobig(&x1,x1_gen);
    gentobig(&x2,x2_gen);


    phi->A = *base_curve;
    phi->len = f;

    uintbig res_temp;
    //this is in the event were this is the first step of the signature, we do a special procedure for compression
    if (special_fixed) {
      xBIDIM(&phi->ker,base_curve,&L_kernel,&x1,&Kd,&x2,&P3);
      proj P,Q,PQ,dummy_A;
      find_basis_2e(&P,&Q,&PQ,&phi->A);
      dummy_A = phi->A;
      proj push_points[3];
      push_points[0]=Q;
      push_points[1]=P;
      push_points[2]=PQ;
      phi_I_basis[6] = P;
      eval_walk_mult(phi,&dummy_A,push_points,3);

      bool dlp=mont_power_dlp(&res_temp,&dummy_A,&(push_points[0]),&(push_points[1]),&(push_points[2]),two_tors_height);
      assert(dlp);
      _unused(dlp);


      #ifndef NDEBUG
      Pt1 = PQ;

      for (int i=1;i<f;i++) {
        xDBL(&Pt1,&phi->A,&Pt1);
        // normalize_proj(&P_ker_test);
      }
      assert(!mont_iszero(&Pt1));
      xDBL(&Pt1,&phi->A,&Pt1);
      assert(mont_iszero(&Pt1));
      Pt1 = push_points[1];
      for (int i=1;i<f;i++) {
        xDBL(&Pt1,&dummy_A,&Pt1);
        // normalize_proj(&P_ker_test);
      }
      assert(!mont_iszero(&Pt1));
      xDBL(&Pt1,&dummy_A,&Pt1);
      assert(mont_iszero(&Pt1));
      #endif


      #ifndef NDEBUG
      xMUL(&push_points[1],&dummy_A,&push_points[1], &res_temp);
      assert(mont_equal(&push_points[0],&push_points[1]));

      #endif

      uintbig big_res;
      uintbig_set(&big_res,pow(2,two_tors_height));
      uintbig_sub3(&res_temp,&big_res,&res_temp);
      *res = res_temp;
      xBIDIM(&(phi->ker), &phi->A, &P, &res_temp, &Q, &uintbig_1, &PQ);

      #ifndef NDEBUG
      proj P_test;
      xBIDIM(&(P_test), &phi->A, &P, &res_temp, &Q, &uintbig_1, &PQ);
      Pt1 = P_test;
      for (int i=1;i<f;i++) {
        xDBL(&Pt1,&phi->A,&Pt1);
      }
      assert(!mont_iszero(&Pt1));
      xDBL(&Pt1,&phi->A,&Pt1);
      assert(mont_iszero(&Pt1));
      two_walk phi_test;
      phi_test.A  = phi->A;
      phi_test.ker= P_test;
      phi_test.len=f;
      proj A_test = phi->A;
      eval_walk(&phi_test, &A_test ,&PQ);

      #endif

    }
    else {
        x1_gen = gmod(gmul(x1_gen,ginvmod(x2_gen,twof)),twof);
      gentobig(&res_temp,x1_gen);
      *res = res_temp;
      xBIDIM(&phi->ker,base_curve,&L_kernel,&res_temp,&Kd,&uintbig_1,&P3);
      #ifndef NDEBUG
      proj P_ker_test = phi->ker;
      // normalize_proj(&P_ker_test);
      // normalize_proj(base_curve);
      for (int i=1;i<f;i++) {
        xDBL(&P_ker_test,base_curve,&P_ker_test);
        normalize_proj(&P_ker_test);
      }
      assert(!mont_iszero(&P_ker_test));
      xDBL(&P_ker_test,base_curve,&P_ker_test);
      assert(mont_iszero(&P_ker_test));
      #endif
    }


    // avma = ltop;

    //happens when we went through the extremal order
    if (!gamma_wit) {
      eich->J = I;
      eich->delta = mkcol4(gen_1,gen_0,gen_0,gen_0);
      eich->gamma = lideal_isom(I,I_two);
    }

    // return res;
  }


}


// T = global_setup.gen_odd_torsion
// f = two_tors_height
// I is a left O_0-ideal of norm dividing T^2 \ell^{2f+delta}
// J is a left O_0-ideal containing I of norm gcd(T^2,n(I))
// K is a left O_0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi * phi_J
// Finds L equivalent to I of norm dividing T^2
// phi_K_basis is the image through phi_K of a basis of the odd torsion, then
// phi is applied to it (a basis = 6 points: A, B and A+B for the curve and its twist)
void ideal_to_isogeny_two_2f_delta(two_walk_long *phi, GEN *L,
    special_isogeny *phi_L, GEN I, GEN J, GEN K,
    proj *phi_K_basis, proj *phi_K_target, int delta, GEN I_long) {
    pari_sp ltop = avma;
    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);
    isomorphism isom;
    clock_t t;
    #ifndef NDEBUG
    proj j1,j2;
    #endif
    long len_step = Z_lval(lideal_norm(I), 2);
    long e1 = (len_step < two_tors_height) ? len_step : two_tors_height;
    if (e1 > len_step) e1 = len_step;
    long e2 = len_step - delta - e1;
    if (e2 < 0) e2 = 0;

    long dist = len_step - e1 - e2;
    if (dist < 0) dist = 0;

    if (dist > 22) fprintf(stderr,"Warning: MITM distance is %ld\n", dist);



    GEN L_;
    printf("klpt...\n");
    if (len_step + Z_lval(lideal_norm(K), 2) < dbllog2r(itor(global_setup.p,10))/2. + 10 ) {
        GEN M;
        GEN alpha = lideal_isom(J, K); // J*alpha = K

        M = lideal_mul(I, alpha); // I*alpha, equivalent to I, but norm a power of 2

        L_ = klpt_special_smooth_small_2e_input(M, famat_sqr(global_setup.smooth_famat_for_klpt));
    }
    else {
        L_ = klpt_special_smooth(I, famat_sqr(global_setup.smooth_famat_for_klpt));
    }
    printf("klpt done\n");

    assert(lideal_isom(I, L_));

    GEN a = lideal_isom(J, K); // J*a = K
    if (gcmp(lideal_norm(K), gen_1) == 0) { a = alg_scalar(A,gen_1); /* make sure we don't apply a distorsion */ }

    GEN M = lideal_mul(I, a);
    assert(lideal_isom(L_,M));
    GEN b = lideal_isom(L_,M); // L_*gamma = M
    GEN gamma = gmul(b, lideal_norm(L_));

    assert(alglatcontains(A, lideal_lattice(K), gamma, NULL));
    assert(alglatcontains(A, lideal_lattice(L_), alg_conj(A, gamma), NULL));
    assert(gcmp(algnorm(A,gamma,0), gmul(powuu(2,len_step),gmul(lideal_norm(K),lideal_norm(L_)))) == 0);

    GEN n;
    alg_primitive(&n, A, order, gamma);
    assert(gcmp(n,gen_1) == 0);



    GEN H1_odd = lideal_create(A, order, gamma, ggcd(global_setup.gen_odd_torsion, lideal_norm(L_)));



    odd_isogeny psi_1;

    GEN factorisation_norm = famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H1_odd));
    GEN coeff = ideal_to_kernel_O0_T(H1_odd,factorisation_norm);

    GEN v_curve = torsion_crt_compose(gel(coeff,1), false);
    GEN v_twist = torsion_crt_compose(gel(coeff,2), true);

    proj ker_curve, ker_twist;

    uintbig x, y;
    gentobig(&x, gel(v_curve,1));
    gentobig(&y, gel(v_curve,2));
    xBIDIM(&ker_curve, phi_K_target, &phi_K_basis[0], &x, &phi_K_basis[1], &y, &phi_K_basis[2]);
    gentobig(&x, gel(v_twist,1));
    gentobig(&y, gel(v_twist,2));
    xBIDIM(&ker_twist, phi_K_target, &phi_K_basis[3], &x, &phi_K_basis[4], &y, &phi_K_basis[5]);


    isog_degree deg_plus, deg_minus;
    if (curve_order_is_p_plus_one)
        famat_to_degree(&deg_plus, &deg_minus, factorisation_norm);
    else
        famat_to_degree(&deg_minus, &deg_plus, factorisation_norm);

    psi_1.deg_plus = deg_plus;
    psi_1.deg_minus = deg_minus;

    if (curve_order_is_p_plus_one) {
        psi_1.kernel_plus = ker_curve;
        psi_1.kernel_minus = ker_twist;
    }
    else {
        psi_1.kernel_plus = ker_twist;
        psi_1.kernel_minus = ker_curve;
    }

    proj psi_1_source = *phi_K_target;


    odd_isogeny psi_1_dual = psi_1;
    proj psi_1_dual_source = psi_1_source;


    if ((psi_1_dual.deg_plus.val != 0) || (psi_1_dual.deg_minus.val != 0)) {

        proj ker_dual[2];




        GEN coeff_dual = dual_coeff(coeff, psi_1.deg_plus, psi_1.deg_minus);

        GEN v_curve_dual = torsion_crt_compose(gel(coeff_dual,1), false);
        gentobig(&x, gel(v_curve_dual,1));
        gentobig(&y, gel(v_curve_dual,2));
        xBIDIM(&ker_dual[0], &psi_1_source, &phi_K_basis[0], &x, &phi_K_basis[1], &y, &phi_K_basis[2]);


        GEN v_twist_dual = torsion_crt_compose(gel(coeff_dual,2), true);
        gentobig(&x, gel(v_twist_dual,1));
        gentobig(&y, gel(v_twist_dual,2));
        xBIDIM(&ker_dual[1], &psi_1_source, &phi_K_basis[3], &x, &phi_K_basis[4], &y, &phi_K_basis[5]);


        printf("eval_mult 1...\n");
        t = tic();
        eval_mult(&psi_1_dual_source, &psi_1, ker_dual, 2);
        TOC(t,"eval_mult");

        if (curve_order_is_p_plus_one) {
            psi_1_dual.kernel_plus = ker_dual[0];
            psi_1_dual.kernel_minus = ker_dual[1];
        }
        else {
            psi_1_dual.kernel_plus = ker_dual[1];
            psi_1_dual.kernel_minus = ker_dual[0];
        }

    }


    GEN gamma_conj = alg_conj(A, gamma);
    GEN H2_odd = lideal_create(A, order, gamma_conj, gdiv(lideal_norm(L_), lideal_norm(H1_odd)));
    GEN H2_two = lideal_create(A, order, gamma_conj, powuu(2, e2));

    odd_isogeny psi_2 = ideal_to_isogeny_O0_T(H2_odd, famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H2_odd)));




    // TODO: redundant computation
    GEN beta = lideal_isom(I, L_); // I*beta = L_
    GEN I_long_next = lideal_mul(I_long, beta); // I_i_long*beta
    GEN I1_next = lideal_create(A, order, lideal_generator(I_long_next), powuu(2, two_tors_height));

    two_walk phi_1_next = ideal_to_isogeny_O0_two(I1_next);



    two_walk phi_2 = ideal_to_isogeny_O0_two(H2_two);

    proj pt;


    proj pt2[2];
    pt2[0] = phi_1_next.ker;
    pt2[1] = phi_2.ker;

        printf("eval_mult 2...\n");
        t = tic();

    eval_mult(&phi_1_next.A, &psi_2, pt2, 2);
        TOC(t,"eval_mult");

    phi_1_next.ker = pt2[0];
    phi_2.ker = pt2[1];
    phi_2.A = phi_1_next.A;


    two_walk phi_2_adjusted; // push kernel away from (0:1)
    proj phi_2_adjusted_target;

    pt = phi_2.ker;
    eval_walk_isom(&isom, &phi_2_adjusted, &phi_2_adjusted_target, &pt, &phi_2, &pt);


    two_walk eta;
    two_walk phi_2_dual = phi_2_adjusted;

    dual_walk(&phi_2_dual);


    // Check dual
    #ifndef NDEBUG
    proj pt_save;
    proj B = phi_2_adjusted.A;
    pt.z = fp2_1;
    do {
        fp2_random(&pt.x);
    } while(!is_on_curve(&pt, &phi_2_adjusted.A));

    pt_save = pt;


    eval_walk(&phi_2_adjusted, &B, &pt);

    mont_isom(&isom, &B, &phi_2_dual.A);
    mont_isom_apply(&isom, &pt);

    B = phi_2_dual.A;
    eval_walk(&phi_2_dual, &B, &pt);

    mont_isom(&isom, &B, &phi_2_adjusted.A);
    mont_isom_apply(&isom, &pt);

    jinv256(&j1, &phi_2_adjusted_target);
    jinv256(&j2, &phi_2_dual.A);
    assert(mont_equal(&j1,&j2));
    #endif


    proj from = psi_1_dual_source, from0 = psi_1_dual_source;
    proj to = phi_2_adjusted_target; // phi_2 source
    two_walk_long phi_2_dual_eta;
    init_trivial_two_walk_long(&phi_2_dual_eta);

    if (dist > 0) {
        bool done;
        assert(dist != 1); // for now, the case dist == 1 crashes
        done = MITM2(&eta, &from, &to, dist);
        assert(done);
        two_walk_composition_ss(&phi_2_dual_eta, &phi_2_dual, &eta);
    }
    else {
        two_walk_stol(&phi_2_dual_eta, &phi_2_dual);
    }





    phi_L->source = global_setup.E0;
    phi_L->phi1 = psi_2;

    phi_L->phi2_dual_set = false;

    // since psi_1_dual has already been computed...
    phi_L->middle = psi_1_dual_source; // = psi_1_dual_source
    phi_L->phi2 = push_odd_isogeny_through_two_walk_long(&psi_1_dual, &phi_L->middle, &phi_2_dual_eta);

    phi_L->phi2_set = true;



    // push phi_2 and phi_1_next simultaneously through phi_L->phi2
    two_walk phi2_pushed = phi_2;


    isomorphism isom3;
    mont_isom(&isom3, &phi_1_next.A, &phi_L->middle);
    mont_isom_apply(&isom3, &phi_1_next.ker);
    phi_1_next.A = phi_L->middle;



    mont_isom_apply(&isom3, &phi2_pushed.ker);
    phi2_pushed.A = phi_L->middle;

    pt2[0] = phi2_pushed.ker;
    pt2[1] = phi_1_next.ker;
        printf("eval_mult 3...\n");
        t = tic();

    eval_mult(&phi2_pushed.A, &phi_L->phi2, pt2, 2);
        TOC(t,"eval_mult");

    phi2_pushed.ker = pt2[0];
    phi_1_next.ker = pt2[1];
    phi_1_next.A = phi2_pushed.A;


    rand_isom(&isom, &phi2_pushed.A);
    mont_isom_apply(&isom, &phi2_pushed.ker);

    proj phi2_pushed_target, proj_tmp;
    eval_walk(&phi2_pushed, &phi2_pushed_target, &proj_tmp);

    two_walk phi2_pushed_dual = phi2_pushed;
    dual_walk(&phi2_pushed_dual);




    from = psi_1_source, from0 = psi_1_source;
    to = phi2_pushed_target; // phi_2 source
    two_walk mitm_up;
    if (dist > 0) {
        bool done;
        done = MITM2(&mitm_up, &from, &to, dist);
        assert(done);
        two_walk_composition_ss(phi, &phi2_pushed_dual, &mitm_up);
    }
    else {
        two_walk_stol(phi, &phi2_pushed_dual);
    }


    two_walk_composition_sl(phi, &phi_1_next, phi);

    eval_walk_long_mult(phi, phi_K_target, phi_K_basis, 6);

    *L = gerepilecopy(ltop, L_);
    assert(lideal_isom(I, *L));

    free_two_walk_long(&phi_2_dual_eta);
}








// T = global_setup.gen_odd_torsion
// I is a left O0-ideal of norm dividing T^2 2^e for some positive integer e
// J = I + O0*T^2
// K is a left O0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi * phi_J
// Finds L equivalent to I of norm dividing T^2
void ideal_to_isogeny_two(two_walk_long *phi_res, GEN *L, special_isogeny *phi_L,
    GEN I, GEN J, GEN K, const special_isogeny *phi_J, const two_walk_long *phi_K,
    bool endpoint_close_to_E0){
    pari_sp ltop = avma;
    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);

    assert(lideal_isom(J, K));


    #ifndef NDEBUG
    GEN X = lideal_create(A, order, lideal_generator(I), lideal_norm(J));
    assert(gcmp(algnorm(A,lideal_isom(J, X),0),gen_1) == 0);
    assert(lideal_isom(J, K));

    if (phi_K->len > 0) {
        proj j1,j2;
        proj P = phi_K->phi[phi_K->len-1].ker;
        proj E = phi_K->phi[phi_K->len-1].A;
        jinv256(&j1, &phi_J->target);
        eval_walk(&phi_K->phi[phi_K->len-1], &E, &P);
        jinv256(&j2, &E);
        assert(mont_equal(&j1,&j2));
        assert(phi_J->phi2_set);
        assert(phi_J->phi2_dual_set);
    }

    long len = 0;
    for (int i = 0; i < phi_K->len; ++i) {
        len += phi_K->phi[i].len;
    }
    #endif

    // TODO : adjust delta
    long delta = (two_tors_height < 14) ? two_tors_height : 14;
    long len_phi = Z_lval(lideal_norm(I), 2);
    long len_step = 2*two_tors_height + delta;

    if (len_phi == 0) {
        avma = ltop; *L = J; *phi_L = *phi_J;
        free_two_walk_long(phi_res);
        init_trivial_two_walk_long(phi_res);
        return;
    }

    long steps = len_phi / len_step;
    if (steps*len_step < len_phi) ++steps;


    proj phi_K_basis[6], phi_K_target = global_setup.E0;
    phi_K_basis[0] = torsion_basis_sum[0];
    phi_K_basis[1] = torsion_basis_sum[1];
    phi_K_basis[2] = torsion_basis_sum[2];
    phi_K_basis[3] = torsion_basis_twist_sum[0];
    phi_K_basis[4] = torsion_basis_twist_sum[1];
    phi_K_basis[5] = torsion_basis_twist_sum[2];

    special_isogeny phi_J_i;
    two_walk_long phi_K_i, phi[steps];
    init_trivial_two_walk_long(&phi_K_i);

    long len_first_step = (endpoint_close_to_E0) ? (len_phi % len_step) : len_step;

    GEN I_i_long = I;
    GEN I_i_short = lideal_create(A, order, lideal_generator(I), gmul(lideal_norm(J),powuu(2, len_first_step)));
    GEN J_i = J;
    GEN K_i = K;
    phi_J_i = *phi_J;
    copy_two_walk_long(&phi_K_i, phi_K);


    GEN alpha, beta;


    int len_phi1 = (len_first_step < two_tors_height) ? len_first_step : two_tors_height;
    GEN I1 = lideal_create(A, order, lideal_generator(I_i_short), powuu(2, len_phi1));
    // printf("")
    two_walk phi_1_0 = ideal_to_isogeny_O0_two(I1);
    two_walk phi_1 = push_two_walk_through_special_isogeny(&phi_1_0, &phi_J_i);

    two_walk_composition_sl(&phi_K_i, &phi_1, &phi_K_i);


    eval_walk_long_mult(&phi_K_i, &phi_K_target, phi_K_basis, 6);



    for (int i = 0; i < steps; ++i) {
        printf("STEP %d/%ld\n",i+1,steps);
        alpha = lideal_isom(J_i, K_i); // J_i*alpha = K_i
        if (gcmp(lideal_norm(K_i), gen_1) == 0) {
            alpha = alg_scalar(A, gen_1); /* make sure we don't apply a distorsion */

        }

        init_trivial_two_walk_long(&phi[i]);

        ideal_to_isogeny_two_2f_delta(&phi[i], &J_i, &phi_J_i,
                                      I_i_short, J_i, K_i,
                                      phi_K_basis, &phi_K_target, delta,
                                      I_i_long);

        assert(lideal_isom(J_i, I_i_short));

        two_walk_composition_ll(&phi_K_i, &phi[i], &phi_K_i);
        // update K, phi_K and I_i_*

        K_i = lideal_mul(I_i_short, alpha); // I_i_short*alpha
        beta = lideal_isom(I_i_short, J_i); // I_i*beta = J_i
        I_i_long = lideal_mul(I_i_long, beta); // I_i_long*beta
        I_i_short = lideal_create(A, order, lideal_generator(I_i_long), gmul(lideal_norm(J_i),powuu(2, len_step)));

    }

    *L = gerepilecopy(ltop, J_i);
    assert(lideal_isom(I, *L));
    *phi_L = phi_J_i;



    // reconstruct phi
    two_walk_long phi_full;
    init_trivial_two_walk_long(&phi_full);
    two_walk_stol(&phi_full, &phi_1);
    for (int i = 0; i < steps; ++i) {
        two_walk_composition_ll(&phi_full, &phi[i], &phi_full);
        free_two_walk_long(&phi[i]);
    }


    copy_two_walk_long(phi_res, &phi_full);
    free_two_walk_long(&phi_full);
    free_two_walk_long(&phi_K_i);

}

// I is a left O0-ideal of norm 2^*, phi_I is the isogeny corresponding to I
// J is a left O0-ideal equivalent to I of rigth order O
// K is equal to K0 inter J where K0 is the left O-ideal of norm 2 to be translated
// L is a left O-ideal of norm 2^f corresponding to the end of phi_I
// Finds phi_K0
//only works when phi_I has at least degree 2^f for now.
//the value endpoint_close_to_E0 first indicates if the beginning curve has a special property
//after the computation it will indicate if the computation as failed somehow
void ideal_to_isogeny_two_new(two_walk_long *phi_res,
    GEN I,GEN J, GEN K, GEN L, const two_walk_long *phi_I, uintbig *zip,
    bool *endpoint_close_to_E0,proj *endpoint){

    pari_sp ltop = avma;
    GEN A = lideal_algebra(J);

    long len_phi = Z_lval(lideal_norm(K), 2);
    long len_step = two_tors_height;
    if (len_phi == 0) {
        avma = ltop;
        free_two_walk_long(phi_res);
        init_trivial_two_walk_long(phi_res);
        return;
    }

    long steps = len_phi / len_step;
    if (steps*len_step < len_phi) ++steps;

    *endpoint= global_setup.E0;
    proj phi_I_basis[7];
    phi_I_basis[0] = torsion_basis_sum[0];
    phi_I_basis[1] = torsion_basis_sum[1];
    phi_I_basis[2] = torsion_basis_sum[2];
    phi_I_basis[3] = torsion_basis_twist_sum[0];
    phi_I_basis[4] = torsion_basis_twist_sum[1];
    phi_I_basis[5] = torsion_basis_twist_sum[2];

    proj P = phi_I_basis[2];
    xMUL(&P,endpoint,&P,&p_even_cofactor);
    assert(mont_iszero(&P));

    two_walk phi[steps];
    long len_first_step = len_step;


    assert(len_first_step == two_tors_height);
    GEN I_two_i = I;
    GEN I_i = J;
    assert(lideal_isom(I_i,I_two_i));
    GEN J_i = J;
    GEN L_i = L;
    GEN K_i_long = K;
    GEN al = lideal_generator_coprime(K_i_long,gen_1);



    GEN K_i_short = lideal_create(A,lideal_order(J_i),al,gmul(lideal_norm(J_i),gpowgs(gen_2,len_first_step)));
    K_i_short=lideal_create(A,lideal_right_order(K_i_short),alg_conj(A,lideal_generator_coprime(K_i_short,gen_1)),gpowgs(gen_2,len_first_step));
    K_i_short=lideal_create(A,lideal_right_order(J_i),alg_conj(A,lideal_generator_coprime(K_i_short,gen_1)),gpowgs(gen_2,len_first_step));
    assert(gequal(lideal_norm(K_i_short),gpowgs(gen_2,len_first_step)) );
    two_walk phi_dual;
    phi_dual = phi_I->phi[phi_I->len-1];
    if (*endpoint_close_to_E0) {
      isomorphism isom;
      // rand_isom(&isom,&phi_dual.A);
      mont_isom(&isom,endpoint,&phi_dual.A);
      mont_isom_apply(&isom,&phi_I_basis[0]);
      mont_isom_apply(&isom,&phi_I_basis[1]);
      mont_isom_apply(&isom,&phi_I_basis[2]);
      mont_isom_apply(&isom,&phi_I_basis[3]);
      mont_isom_apply(&isom,&phi_I_basis[4]);
      mont_isom_apply(&isom,&phi_I_basis[5]);
      // mont_isom_apply(&isom,&phi_dual.ker);
      *endpoint = phi_dual.A;
    }

    //pushing the points through the isogeny phi_I
    eval_walk_long_mult(phi_I, endpoint, phi_I_basis, 6);
    if (phi_I-> len > 1) {
      assert(phi_I->phi[phi_I->len-1].len == two_tors_height);
    }

    dual_walk_no_isom(&phi_dual);

    assert(mont_equal(&phi_dual.A,endpoint));
    normalize_proj(&phi_dual.ker);
    normalize_proj(&phi_dual.A);
    phi_I_basis[6] = phi_dual.ker;
    #ifndef NDEBUG
    proj P6;
    P6 = phi_dual.ker;
    assert(!mont_iszero(&P6));
    for (int i =1;i< len_step;i++) {
      xDBL(&P6,endpoint,&P6);
    }
    assert(!mont_iszero(&P6));
    xDBL(&P6,endpoint,&P6);
    assert(mont_iszero(&P6));
    #endif
    long f = len_step;
    GEN gamma,delta;
    eichler_package eich;
    eich.beta=NULL; eich.J = NULL; eich.delta=NULL; eich.I = NULL; eich.L=NULL;
    eich.gamma = NULL;
    for (int i = 0; i < steps-1; ++i) {
        assert(lideal_isom(I_i,I_two_i));

        bool special_fixed =  (i == 0 && !*endpoint_close_to_E0 ) ;
        bool other_special_extremal = ( (i==1 && two_tors_height < 60 && ! *endpoint_close_to_E0) || (*endpoint_close_to_E0 && (i+1)*len_step < 128 ) );

        //this  is the main function of the loop
        assert(alglatcontains(A,lideal_right_order(I_i),lideal_generator(L_i),NULL));
        // clock_t t=tic();
        ideal_to_isogeny_small_two_f(&phi[i],&zip[i],&eich,I_i,I_two_i,K_i_short,L_i,phi_I_basis,&(phi_dual.A),f,special_fixed,other_special_extremal);
        // printf("ideal_to_isogeny_small_two_f n° %d ",i);TOC(t,"");

        if (eich.delta) {

          delta = alg_conj(A,eich.delta);
          assert(lideal_isom(I_i,I_two_i));

          //compute the final step size if necessary
          if (i == steps-2) {

            f = len_phi - (i+1)*f;
          }
          //actualize all the ideals used
          J_i =eich.J;
          GEN twof  = gpowgs(gen_2,f);
          assert(lideal_isom(I_i,I_two_i));
          assert(lideal_isom(I_i,J_i));
          assert(lideal_equals(J_i,lideal_mul(I_i,delta)));

          K_i_long = lideal_mul(K_i_long,delta);
          GEN alpha_long = lideal_generator_coprime(K_i_long,gen_1);
          I_i = lideal_create(A,lideal_order(K_i_long),alpha_long,gmul(lideal_norm(J_i),gpowgs(gen_2,two_tors_height)));
          GEN Ori = lideal_right_order(I_i);
          L_i = lideal_create(A,Ori,alg_conj(A,lideal_generator_coprime(I_i,gen_1)),lideal_norm(K_i_short));
          assert(gequal(lideal_norm(L_i),lideal_norm(K_i_short)));
          assert(gequal(lideal_norm(I_i),gmul(lideal_norm(J_i),gpowgs(gen_2,two_tors_height))));
          K_i_short = lideal_create(A,lideal_order(K_i_long),alpha_long,gmul(lideal_norm(J_i),gpowgs(gen_2,two_tors_height+f)));
          K_i_short=lideal_create(A,lideal_right_order(K_i_short),alg_conj(A,lideal_generator_coprime(K_i_short,gen_1)),twof);
          K_i_short=lideal_create(A,Ori,alg_conj(A,lideal_generator_coprime(K_i_short,gen_1)),twof);
          assert(gequal(lideal_norm(K_i_short),twof) );
          gamma = eich.gamma;
          assert(gamma);
          I_two_i = lideal_mul(I_i,gamma);
          assert(lideal_isom(I_i,I_two_i));
          #ifndef NDEBUG
          proj P_test3 = phi[i].ker;
          proj A_test = phi[i].A;
          // normalize_proj(&phi[i].A);
          normalize_proj(&P_test3);
          // assert(!mont_iszero(&P_test3));
          for (int i =1;i<two_tors_height;i++){
            // printf("phi[i].A = "); print_proj_hash(&phi[i].A); printf("\n");
            xDBL(&P_test3,&A_test,&P_test3);
            normalize_proj(&P_test3);
          }
          assert(!mont_iszero(&P_test3));
          xDBL(&P_test3,&A_test,&P_test3);
          assert(mont_iszero(&P_test3));
          #endif
          //push  the points through the neew isogeny
          eval_walk_mult(&phi[i], endpoint, phi_I_basis, 7);
          phi_dual.A = *endpoint;
          phi_dual.ker=phi_I_basis[6];
          phi_dual.len = phi[i].len;
          #ifndef NDEBUG
          proj P_test = phi_dual.ker;
          for (int i =1;i<two_tors_height;i++){
            xDBL(&P_test,&phi_dual.A,&P_test);
          }
          assert(!mont_iszero(&P_test));
          xDBL(&P_test,&phi_dual.A,&P_test);
          assert(mont_iszero(&P_test));


          #endif
        }
        else {
          break;
        }


    }
    if (eich.delta) {
      #ifndef NDEBUG
      proj P_test = phi_dual.ker;
      for (int i =0;i<two_tors_height;i++){
        xDBL(&P_test,&phi_dual.A,&P_test);
      }
      assert(mont_iszero(&P_test));
      #endif
      //this is the final step
      //if the walk is too short, there might be a problem because special_fixed and special_step_1 are both false
      if (f != 0) {
        bool special_fixed =  false;
        //during keygen we need to use the extremal order again when two_tors_height is too ms
        bool special_step_1 = (*endpoint_close_to_E0 && two_tors_height < 59)  ;
        //compute the last kernel
        for (int i= 0; i< two_tors_height-f;i++){
          xDBL(&phi_I_basis[6],&phi_dual.A,&phi_I_basis[6]);
        }
        #ifndef NDEBUG
        proj P_test = phi_I_basis[6];
        for (int i =0;i<f;i++){
          xDBL(&P_test,&phi_dual.A,&P_test);
        }
        assert(mont_iszero(&P_test));
        #endif
        // normalize_proj(&phi_dual.A);
        // normalize_proj(&phi_I_basis[6]);
        ideal_to_isogeny_small_two_f(&phi[steps-1],&zip[steps-1],&eich,I_i,I_two_i,K_i_short,L_i,phi_I_basis,&phi_dual.A,f,special_fixed,special_step_1);
      }
      *endpoint = phi_I_basis[6];
      // eval_walk(&phi[steps-1],endpoint,&phi_I_basis[0]);

      // reconstruct phi
      two_walk_long phi_full;
      init_trivial_two_walk_long(&phi_full);
      two_walk_stol(&phi_full, &phi[0]);
      for (int i = 1; i < steps; ++i) {
          two_walk_composition_sl(&phi_full, &phi[i], &phi_full);
          // free(&phi[i]);
      }
      copy_two_walk_long(phi_res, &phi_full);
      free_two_walk_long(&phi_full);

      //tell if the computation succeeded
      *endpoint_close_to_E0 = true;
    }
    else {
      *endpoint_close_to_E0 = false;
    }


}




void ideal_to_isogeny_O0_two_long(two_walk_long *phi, GEN *L, special_isogeny *phi_L, GEN I,bool endpoint_close_to_E0) {
    pari_sp ltop = avma;
    GEN trivial_ideal = lideal_create(global_setup.B, global_setup.O0, mkcol4s(1,0,0,0), gen_1);
    special_isogeny triv_special = trivial_special_isogeny();
    two_walk_long triv_two;
    init_trivial_two_walk_long(&triv_two);

    ideal_to_isogeny_two(phi, L, phi_L, I, trivial_ideal, trivial_ideal, &triv_special, &triv_two,endpoint_close_to_E0);
    *L = gerepilecopy(ltop, *L);
}


void random_Fp_point_two_f(proj *P, proj *A) {
  bool oncurve = class_mod_4 == 3;
  while (true) {
    fp2_random_fp(&P->x);
    // fp2_random_fp(&P->z);

    if (is_on_curve(P, A) != oncurve)
      continue;
    // multiply by cofactor
    xMUL(P, A, P, &p_even_cofactor);
    // check it has maximal order
    proj P2 = *P;
    for (int i = 1; i < two_tors_height; i++)
      xDBL(&P2, A, &P2);

    if (!mont_iszero(&P2)){
      #ifndef NDEBUG
      xDBL(&P2, A, &P2);
      assert(mont_iszero(&P2));
      #endif
        break;
    }

  }

}



//this function will perform the translation from ideal to isogeny when the starting curve is j=j1728
//warning this probably does not work for any other curve !!
// The ideal I has norm 2^*
//will fail when I is not contained in O0< i+1,2>
void ideal_to_isogeny_1728(two_walk_long* phi ,GEN I,proj *pk) {
  pari_sp ltop = avma;
  GEN ell=gen_2;
  _unused(ell);
  long len_phi = Z_lval(lideal_norm(I), 2);
  long len_step = two_tors_height;
  if (len_phi == 0) {
      avma = ltop;
      free_two_walk_long(phi);
      init_trivial_two_walk_long(phi);
      return;
  }

  long steps = len_phi / len_step;
  if (steps*len_step < len_phi) ++steps;
  two_walk phi_I_one;
  two_walk_long phi_I_start;

  GEN order= lideal_order(I);
  GEN A = lideal_algebra(I);
  proj phi_I_target = global_setup.E0;

  //the first step is the tricky one
  //we need a setup to apply the algorithm in ideal_to_isogeny_small_two_f
  //a point P of order 2^f, the corresponding ideal and a "good" endomorphism of O0
  // for that we will use the point of order 2^f defined over Fp, the ideal generated by frob-1 and the endomorphis sqrt{-q}.
  proj P;
  P.z = fp2_1;
  random_Fp_point_two_f(&P,&phi_I_target);
  #ifndef NDEBUG
  proj P2;
  mont0_frob(&P2, &P);
  assert(mont_equal(&P2,&P));
  for (int i =1;i< len_step;i++) {
    xDBL(&P2,&phi_I_target,&P2);
  }
  assert(!mont_iszero(&P2));

  proj Poi;
  deterministic_second_point(&Poi,&P,&phi_I_target,len_step);
  proj poin2 = Poi;
  for (int i =1;i< len_step;i++) {
    xDBL(&poin2,&phi_I_target,&poin2);
  }
  assert(!mont_iszero(&poin2));
  assert(!mont_equal(&poin2,&P2));
  proj point2frob;
  mont0_frob(&point2frob, &poin2);
  assert(!mont_equal(&point2frob,&poin2));
  assert(!mont_equal(&point2frob,&P2));

  xDBL(&poin2,&phi_I_target,&poin2);
  assert(mont_iszero(&poin2));
  xDBL(&point2frob,&phi_I_target,&point2frob);
  assert(mont_iszero(&point2frob));
  xDBL(&P2,&phi_I_target,&P2);
  assert(mont_iszero(&P2));
  #endif

  GEN quat_j_minus_1 = mkcol4(gneg(gen_1),gen_0,gen_1,gen_0);
  GEN ellf = gpowgs(gen_2,two_tors_height);
  GEN ideal_P = lideal_create(A,order,quat_j_minus_1,ellf);
  assert(gequal(lideal_norm(ideal_P),ellf));
  GEN small_I = lideal_create(A,order,lideal_generator(I),ellf);
  assert(gequal(lideal_norm(small_I),ellf));

  proj2 Pt1,Pt2;
  xtoxy(&Pt1,&phi_I_target,&P);

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

  #ifndef NDEBUG
    proj P1,P3;
    xytox(&P3,&Pt2);
    xytox(&P1,&Pt1);
    assert(!xy_is_zero(&Pt2));
    assert(!mont_iszero(&P3));
    for (int i =1;i< two_tors_height;i++) {
      xDBL(&P1,&phi_I_target,&P1);
      xDBL(&P3,&phi_I_target,&P3);
    }
    assert(!mont_iszero(&P1));
    assert(!mont_iszero(&P3));
    xDBL(&P1,&phi_I_target,&P1);
    xDBL(&P3,&phi_I_target,&P3);
    assert(mont_iszero(&P3));
    assert(mont_iszero(&P1));
  #endif

  proj P_ker1;
  GEN coeffs = find_coeff(small_I,ideal_P,mkcol4(gen_0,gen_1,gen_0,gen_0),mkcol4(gen_1,gen_0,gen_0,gen_0),order);
  uintbig C,D;
  gentobig(&C,gel(coeffs,1));
  gentobig(&D,gel(coeffs,2));

  xyMUL(&Pt1,&phi_I_target,&Pt1,&C);
  xyMUL(&Pt2,&phi_I_target,&Pt2,&D);
  xyADD(&Pt1,&phi_I_target,&Pt1,&Pt2);
  xytox(&P_ker1,&Pt1);

  isomorphism isom;
  rand_isom(&isom,&phi_I_target);
  mont_isom_apply(&isom,&P_ker1);
  #ifndef NDEBUG
  proj P4;
  P4 = P_ker1;
  // assert(!xy_is_zero(&Pt2));
  assert(!mont_iszero(&P4));
  for (int i =1;i< two_tors_height;i++) {
    xDBL(&P4,&phi_I_target,&P4);
  }
  assert(!mont_iszero(&P4));
  xDBL(&P4,&phi_I_target,&P4);
  assert(mont_iszero(&P4));
  #endif

  phi_I_one.A = phi_I_target;
  phi_I_one.ker = P_ker1;
  phi_I_one.len = len_step;

  init_trivial_two_walk_long(&phi_I_start);
  two_walk_stol(&phi_I_start, &phi_I_one);





  if (steps > 1) {
    long bound_coeff_J  = 10;
    GEN delta;
    GEN J = lideal_equiv_prime_random(small_I,&delta,stoi(bound_coeff_J));
    assert(lideal_equals(small_I,lideal_mul(J,delta))); // small_I = J * delta
    GEN deltabis = lideal_isom(J,small_I);
    deltabis = alg_conj(A,delta);
    assert(lideal_equals(J,lideal_mul(small_I,deltabis)));

    assert(alglatcontains(A,lideal_right_order(small_I),lideal_generator_coprime(I,ell),NULL));
    GEN K;
    GEN alpha_L = algmul(A,alg_conj(A,lideal_generator_coprime(small_I,lideal_norm(J))),alg_conj(A,delta));
    alpha_L = gdiv( algmul(A,delta,alpha_L),algnorm(A,delta,0));

    assert(alglatcontains(A,lideal_right_order(J),alpha_L,NULL));
    K = lideal_mul(I,deltabis);
    GEN L = lideal_create(A,lideal_right_order(J),alpha_L,ellf);
    assert(gequal(lideal_norm(L),ellf));

    assert(gequal( lideal_norm(K),gmul(lideal_norm(J),gpowgs(ell,len_phi-len_step))));
    assert(lideal_isom(I,K));


    two_walk_long phi_I_end;
    init_trivial_two_walk_long(&phi_I_end);
    uintbig zip[steps++];
    bool endpoint_close_to_E0 = true;

    ideal_to_isogeny_two_new(&phi_I_end, small_I,J, K, L, &phi_I_start,zip, &endpoint_close_to_E0,pk);
    assert(endpoint_close_to_E0);
    two_walk_composition_ll(phi,&phi_I_end,&phi_I_start);
    free_two_walk_long(&phi_I_start);
    free_two_walk_long(&phi_I_end);
  }
  else {
    copy_two_walk_long(phi, &phi_I_start);
    free_two_walk_long(&phi_I_start);
  }


  avma = ltop;
}


//This takes two coefficients a,b and outputs the O0-ideal of norm 2^f whose kernel is (a+ib)P
// where P is the point of order 2^f defined over Fp.
//only works when f < two_tors_height for now (because we have not precomputed a basis of the 2^two_tors_height torsion)

void kernel_to_ideal_two_f(GEN a, GEN b,long f,GEN *I) {

  pari_sp ltop = avma;
  GEN A = global_setup.B;
  GEN order = global_setup.O0;
  GEN quat_j_minus_1 = mkcol4(gneg(gen_1),gen_0,gen_1,gen_0);
  GEN ellf = gpowgs(gen_2,f);
  GEN one_minus_i = mkcol4(gen_1,gneg(gen_1),gen_0,gen_0);
  GEN a_minus_ib = mkcol4(a,gneg(b),gen_0,gen_0);
  GEN ideal_P_2f = lideal_create(A,order,gdiv(algmul(A,quat_j_minus_1,one_minus_i),gen_2) , ellf );

  assert(gequal(lideal_norm(ideal_P_2f),ellf));

  GEN J = lideal_create(A,order,algmul(A,lideal_generator_coprime(ideal_P_2f,gen_2),a_minus_ib),ellf);
  assert(gequal(lideal_norm(J),ellf));

  *I = gerepilecopy(ltop,J);
}
