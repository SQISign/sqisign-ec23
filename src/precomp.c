#define _XOPEN_SOURCE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <math.h>
#include <assert.h>

#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"

#include "mont.h"
#include "tedwards.h"
#include "constants.h"
#include "precomputed.h"
#include "curve.h"
#include "two_walks.h"

struct quaternion_setup_t {
    GEN p; // the prime
    uintbig p_uint; // the prime
    GEN B; // the quaternion algebra
    GEN qf; // the quaternion algebra
    GEN O0; // the cannonical maximal order
    GEN one;
    GEN i;
    GEN j;
    GEN ji;
    GEN torsion_fm; // factorisation matrix of the available torsion

    GEN O0_b1;
    GEN O0_b2;
    GEN O0_b3;
    GEN O0_b4;
    GEN O0_to_standard;
    GEN standard_to_O0;

    proj E0;
};

struct quaternion_setup_t precomp_setup;


uintbig stobig(long long x) {
    uintbig x_big;
    uintbig_set(&x_big, x);
    return x_big;
}

char* pari_int_code(GEN i) {
    if (is_bigint(i)) return pari_sprintf("strtoi(\"%Ps\")", i);
    else return pari_sprintf("stoi(%PsULL)", i);
}

char* pari_2x2_matrix_code(GEN M) {
    return  pari_sprintf("mkmat2(mkcol2(%s,%s),mkcol2(%s,%s))", pari_int_code(gcoeff(M,1,1)),
            pari_int_code(gcoeff(M,2,1)),
            pari_int_code(gcoeff(M,1,2)),
            pari_int_code(gcoeff(M,2,2)));
}


char* fp_code(const fp *x) {
    uintbig p;
    fp_dec(&p, x);
    return pari_sprintf("{ 0x%lxULL, 0x%lxULL, 0x%lxULL, 0x%lxULL }", p.c[0], p.c[1], p.c[2], p.c[3]);
}

char* fp2_code(const fp2 *x) {
    return pari_sprintf("{ %s,\n  %s }", fp_code(&x->re), fp_code(&x->im));
}

char* proj_code(const proj *P) {
    if (fp2_iszero(&P->z))
        return pari_sprintf("{ %s,\n  %s }", fp2_code(&P->x), fp2_code(&P->z));
    else {
        fp2 x,z,tmp = P->z;
        fp2_inv(&tmp);
        fp2_mul3(&x,&P->x,&tmp);
        fp2_mul3(&z,&P->z,&tmp);
        return pari_sprintf("{ %s,\n  %s }", fp2_code(&x), fp2_code(&z));
    }
}

char* proj2_code(const proj2 *P) {
    if (fp2_iszero(&P->z))
        return pari_sprintf("{ %s,\n  %s,\n  %s }", fp2_code(&P->x), fp2_code(&P->y), fp2_code(&P->z));
    else {
        fp2 x,y,z,tmp = P->z;
        fp2_inv(&tmp);
        fp2_mul3(&x,&P->x,&tmp);
        fp2_mul3(&y,&P->y,&tmp);
        fp2_mul3(&z,&P->z,&tmp);
        return pari_sprintf("{ %s,\n  %s,\n  %s }", fp2_code(&x), fp2_code(&y), fp2_code(&z));
    }
}

char* ted_code(const point *P) {
    return pari_sprintf("{ %s,\n  %s,\n  %s,\n  %s }", fp2_code(&P->x), fp2_code(&P->y), fp2_code(&P->z), fp2_code(&P->t));
}

void print_loop(char* prefix, char* name1, char* name2) {
    printf("%sfor (int i = 0; i < 3; ++i) {\n",prefix);
    printf("%s\tfp_enc( &(&(&%s[i])->x)->re, &(%s[i][0][0]) );\n",prefix,name1,name2);
    printf("%s\tfp_enc( &(&(&%s[i])->x)->im, &(%s[i][0][1]) );\n",prefix,name1,name2);
    printf("%s\tfp_enc( &(&(&%s[i])->z)->re, &(%s[i][1][0]) );\n",prefix,name1,name2);
    printf("%s\tfp_enc( &(&(&%s[i])->z)->im, &(%s[i][1][1]) );\n",prefix,name1,name2);
    printf("%s}\n",prefix);
}

GEN norm0(GEN x) {
    return algnorm(precomp_setup.B, x,0);
}


void random_point(proj *P, proj const *A, long ell, long e, bool twist) {
    uintbig cofactor;
    if ((!twist && curve_order_is_p_plus_one) || (twist && !curve_order_is_p_plus_one))
        uintbig_add3(&cofactor, &p, &uintbig_1);
    else { uintbig_sub3(&cofactor, &p, &uintbig_1); }


    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell);
    }
    proj Z;




    while (1) {
        fp2_random(&P->x); P->z = fp2_1;
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
}

void random_basis(proj *P1, proj *P2, point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
    point P1_mul, P2_mul, tmp;
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    fp2 weil;




    proj E;
    mont_to_ted(&E, A, twist);

    random_point(P1, A, ell, e, twist);

    mont_to_ted_point(P1_ted, A, P1);
    P1_mul = *P1_ted;
    for (int i = 0; i < e-1; ++i) {
        ted_mul(&P1_mul, &P1_mul, &E, &ell_big);
    }

    assert(ted_is_on_curve(&P1_mul,&E));
    assert(!ted_iszero(&P1_mul));
    ted_mul(&tmp, &P1_mul, &E, &ell_big);
    assert(ted_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        mont_to_ted_point(P2_ted, A, P2);
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
    } while (fp2_iszero(&weil));

}

void random_basis_2e(proj *P1, proj *P2, proj const *A, long e, bool twist) {
    proj P1_mul, P2_mul, tmp;
    uintbig ell_big;
    long ell = 2;
    uintbig_set(&ell_big, ell);

    random_point(P1, A, ell, e, twist);

    P1_mul = *P1;
    for (int i = 0; i < e-1; ++i) {
        xMUL(&P1_mul, A, &P1_mul, &ell_big);
    }

    assert(is_on_curve(&P1_mul,A));
    assert(!mont_iszero(&P1_mul));
    xMUL(&tmp, A, &P1_mul, &ell_big);
    assert(mont_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        P2_mul = *P2;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&P2_mul, A, &P2_mul, &ell_big);
        }
        assert(is_on_curve(&P2_mul,A));
        assert(!mont_iszero(&P2_mul));
        xMUL(&tmp, A, &P2_mul, &ell_big);
        assert(mont_iszero(&tmp));

    } while (mont_equal(&P1_mul,&P2_mul));

}

bool fp2_ord_trivial(long *res, const fp2 *g, long bound) {
    long order = 1;
    fp2 x = *g;

    for (int i = 0; i < bound; ++i) {
        if (fp2_equal(&fp2_1,&x)) { *res = order; return true; }
        order++;
        fp2_mul2(&x,g);
    }

    return false;
}

GEN action_two_3_4(GEN m_i, GEN m_j, long e) {
    GEN m_ij2,m_1ji2;
    GEN gelle = stoi(1LL<<e);
    GEN gelle1 = stoi(1LL<<(e-1));
    GEN iMi, MM, test1, test2;
    GEN pplus14 = gdiv(gadd(precomp_setup.p,gen_1),stoi(4));
    GEN id = mkmat2(mkcol2s(1,0),mkcol2s(0,1));


    for (int i11 = 0; i11 < 2; ++i11){
        for (int i12 = 0; i12 < 2; ++i12){
            for (int i21 = 0; i21 < 2; ++i21){
                for (int i22 = 0; i22 < 2; ++i22){
                    m_ij2 = gdiv(gadd(m_i, m_j),gen_2);
                    if (i11) gcoeff(m_ij2,1,1) = gadd(gcoeff(m_ij2,1,1),gelle1);
                    if (i12) gcoeff(m_ij2,1,2) = gadd(gcoeff(m_ij2,1,2),gelle1);
                    if (i21) gcoeff(m_ij2,2,1) = gadd(gcoeff(m_ij2,2,1),gelle1);
                    if (i22) gcoeff(m_ij2,2,2) = gadd(gcoeff(m_ij2,2,2),gelle1);
                    printf("B\n");
                    output(gmul(m_i,m_ij2));
                    printf("C\n");
                    output(gmul(gmul(m_i,m_ij2),m_i));
                    printf("D\n");
                    output(gmod(gmul(gmul(m_i,m_ij2),m_i),gelle));
                    printf("E\n");

                    iMi = gmod(gmul(gmul(m_i,m_ij2),m_i),gelle);


                    printf("done\n");

                    test1 = gmod(gsub(m_ij2,m_i),gelle);

                    MM = gmod(gmul(m_ij2,m_ij2),gelle);
                    test2 = gmod(gneg(gmul(id,pplus14)),gelle);



                    if (gequal(iMi,test1) && gequal(MM,test2)) { // a candidate
                        m_1ji2 = gmod(gneg(gmul(m_ij2,m_i)),gelle);



                        printf("\t/* candidate %d%d%d%d */\n", i11,i12,i21,i22);
                        printf("\taction_two_3 = %s;\n", pari_2x2_matrix_code(m_1ji2));
                        printf("\taction_two_4 = %s;\n", pari_2x2_matrix_code(m_ij2));
                    }
                }
            }
        }
    }


    return NULL;
}


void compute_action(GEN *M, void (*endo)(point*, const point*), const point *P1, const point *P2, const proj *E, long ell, long e) {
    point Q;
    GEN a,b,c,d;
    assert(ted_is_on_curve(P1,E));
    endo(&Q,P1);
    // proj E1;
    // printf("E0 = %s\n", proj_code(E));
    // fp2_frob2(&E1.x, &E->x);
    // fp2_frob2(&E1.z, &E->z);
    assert(ted_is_on_curve(&Q,E));

    bool ok = ted_bidim_log(&a, &c, E, &Q, P1, P2, ell, e);
    assert(ok);

    endo(&Q,P2);
    ok = ted_bidim_log(&b, &d, E, &Q, P1, P2, ell, e);
    assert(ok);

    *M = mkmat2(mkcol2(a,c),mkcol2(b,d));
}

void check_action(GEN M, void (*endo)(point*, const point*), void (*mont_endo)(proj*, const proj*), const point *P1, const point *P2, const proj *E, const proj *E_mont, long ell, long e) {
    GEN x,y,u,v;
    uintbig X,Y,U,V;
    point Q, endoQ, endoQ_combination, tmp;

    GEN gelle = powuu(ell,e);

    point P12;
    proj M1,M2,M12,MQ,MendoQ,MendoQ_combination;

    ted_to_mont_point(&M1, P1);
    ted_to_mont_point(&M2, P2);
    ted_add(&P12, E, P1, P2);
    ted_to_mont_point(&M12, &P12);


    for (int i = 0; i < 10; ++i) {
        x = randomi(gelle);
        y = randomi(gelle);
        gentobig(&X, x);
        gentobig(&Y, y);

        ted_mul(&Q, P1, E, &X);
        ted_mul(&tmp, P2, E, &Y);
        ted_add(&Q, E, &Q, &tmp);

        endo(&endoQ,&Q);

        GEN A,B,C,D;
        A = gcoeff(M,1,1);
        B = gcoeff(M,1,2);
        C = gcoeff(M,2,1);
        D = gcoeff(M,2,2);

        u = gadd(gmul(A,x),gmul(B,y));
        v = gadd(gmul(C,x),gmul(D,y));
        gentobig(&U, u);
        gentobig(&V, v);

        ted_mul(&endoQ_combination, P1, E, &U);
        ted_mul(&tmp, P2, E, &V);
        ted_add(&endoQ_combination, E, &endoQ_combination, &tmp);

        assert(ted_equal(&endoQ,&endoQ_combination));


        // same check in montgomery form
        ted_to_mont_point(&MQ, &Q);
        mont_endo(&MendoQ,&MQ);

        xBIDIM(&MendoQ_combination, E_mont, &M1, &U, &M2, &V, &M12);

        assert(mont_equal(&MendoQ,&MendoQ_combination));
    }
}

// P12 = P1 + P2
void compute_action_2e(GEN *M, void (*endo)(proj*, const proj*), const proj *P1, const proj *P2, const proj *P12, const proj *A, long e) {
    proj Q, R;
    uintbig a,b,c,d,x,y;

    endo(&Q,P1);
    assert(is_on_curve(&Q,A));

    bidim_log_2e(&a, &c, A, &Q, P1, P2, P12, e);

    xBIDIM(&R, A, P1, &a, P2, &c, P12);
    assert(mont_equal(&Q, &R));

    endo(&Q,P2);
    bidim_log_2e(&b, &d, A, &Q, P1, P2, P12, e);

    xBIDIM(&R, A, P1, &b, P2, &d, P12);
    assert(mont_equal(&Q, &R));

    endo(&Q,P12);
    bidim_log_2e(&x, &y, A, &Q, P1, P2, P12, e);

    xBIDIM(&R, A, P1, &x, P2, &y, P12);
    assert(mont_equal(&Q, &R));

    uintbig ab, cd, ord = { 0 };
    ord.c[e / 64] = (uint64_t)1 << (e % 64);
    uintbig_add3(&ab, &a, &b);
    uintbig_add3(&cd, &c, &d);

    // Reduce mod 1 << e
    if (uintbig_bit(&ab, e))
      uintbig_sub3(&ab, &ab, &ord);
    if (uintbig_bit(&cd, e))
      uintbig_sub3(&cd, &cd, &ord);

    // Check that x = ±(a+b) and y = ±(c+d)
    if (!uintbig_equal(&x, &ab))
      uintbig_sub3(&ab, &ord, &ab);
    if (!uintbig_equal(&y, &cd))
      uintbig_sub3(&cd, &ord, &cd);
    // If not, change sign of b and d
    if (!uintbig_equal(&x, &ab) || !uintbig_equal(&y, &cd)) {
      uintbig_sub3(&b, &ord, &b);
      uintbig_sub3(&d, &ord, &d);
    }

    // Again
    uintbig_add3(&ab, &a, &b);
    uintbig_add3(&cd, &c, &d);

    // Reduce mod 1 << e
    if (uintbig_bit(&ab, e))
      uintbig_sub3(&ab, &ab, &ord);
    if (uintbig_bit(&cd, e))
      uintbig_sub3(&cd, &cd, &ord);

    // Check that x = ±(a+b) and y = ±(c+d), with the same sign !
    if (!uintbig_equal(&x, &ab)) {
      uintbig_sub3(&ab, &ord, &ab);
      uintbig_sub3(&cd, &ord, &cd);
    }
    assert(uintbig_equal(&x, &ab) && uintbig_equal(&y, &cd));

    *M = mkmat2(mkcol2(bigtogen(&a),bigtogen(&c)),mkcol2(bigtogen(&b),bigtogen(&d)));
}

bool check_action_2e(GEN M, void (*endo)(proj*, const proj*), void (*endoxy)(proj2*, const proj2*), const proj2 *P1xy, const proj2 *P2xy, const proj *E, long e) {
    uintbig X,Y,U,V;
    proj2 Qxy,Rxy,exy1,exy2,eQxy_combination,eQxy;

    GEN t;

    for (int i = 0; i < 20; ++i) {
        uintbig_random(&X);
        uintbig_random(&Y);

        endoxy(&exy1,P1xy);
        endoxy(&exy2,P2xy);

        xyMUL(&Qxy, E, P1xy, &X);
        // printf("xP1xy = %s\n",proj2_code(&Qxy));
        xyMUL(&Rxy, E, P2xy, &Y);
        // output(stoi(y));
        // printf("P2xy = %s\n",proj2_code(P2xy));
        // printf("yP2xy = %s\n",proj2_code(&Rxy));
        // assert(!xy_is_zero(&Rxy));
        xyADD(&Qxy, E, &Qxy, &Rxy);

        endoxy(&eQxy,&Qxy);



        uintbig a,b;
        proj pt1, pt2, pt12,eQ;
        proj2 P12xy;
        xytox(&pt1, P1xy);
        xytox(&pt2, P2xy);
        xyADD(&P12xy, E, P1xy, P2xy);
        xytox(&pt12, &P12xy);
        xytox(&eQ, &eQxy);
        assert(bidim_log_2e(&a, &b, E, &eQ, &pt1, &pt2, &pt12, e));

        // xMUL(&tmp, E, &pt1, &X);
        // printf("xP1 = %s\n",proj_code(&tmp));
        // xMUL(&tmp, E, &pt2, &Y);
        // printf("yP2 = %s\n",proj_code(&tmp));


        t = gmod(RgM_RgC_mul(M, mkcol2(bigtogen(&X),bigtogen(&Y))), mpshift(gen_1, e));
        //printf("e = %ld\n", e);
        // printf("bidim [%ld,%ld] or [%ld,%ld] \n",(1ULL<<e)-a,(1ULL<<e)-b,a,b);
        // printf("coeff ");
        // output(t);
        // printf("mat ");
        // output(M);




        gentobig(&U, gel(t,1));
        gentobig(&V, gel(t,2));

        xyMUL(&eQxy_combination, E, P1xy, &U);
        xyMUL(&Rxy, E, P2xy, &V);
        xyADD(&eQxy_combination, E, &eQxy_combination, &Rxy);
        proj2 neg;
        xyNEG(&neg, &eQxy);


        // printf("eQ = %s\n",proj_code(&eQ));
        // xBIDIM(&eQ_combination, E, &pt1, &U, &pt2, &V, &pt12);
        // printf("eQ_combination1 = %s\n",proj_code(&eQ_combination));
        // uintbig_set(&U, a);
        // uintbig_set(&V, b);
        // xBIDIM(&eQ_combination, E, &pt1, &U, &pt2, &V, &pt12);
        // printf("eQ_combination2 = %s\n",proj_code(&eQ_combination));
        // printf("eQxy = %s\n eQxy_combination = %s\n",proj2_code(&eQxy),proj2_code(&eQxy_combination));

        assert(xy_is_on_curve(E, &eQxy_combination));
        assert(xy_equal(&eQxy,&eQxy_combination) || xy_equal(&neg,&eQxy_combination));

        if (xy_equal(&neg,&eQxy_combination)) return false; // flip sign!
    }
    return true;
}


// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(80000000, 1<<18);

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }


    long var = fetch_var();
    GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC);

    GEN m1 = mkmat2(mkcol2s(1,0),mkcol2s(0,1));

    // Starting curve
    init_curve(&precomp_setup.E0);

    GEN p_pari = strtoi(p_str);

    // The generators of the quaternion algebra
    GEN B_1 = mkcol4s(1,0,0,0);
    GEN B_i = mkcol4s(0,1,0,0);
    GEN B_j = mkcol4s(0,0,1,0);
    GEN B_ji = mkcol4s(0,0,0,1);
    // //GEN B_ij = mkcol4s(0,0,0,-1);

    // The quaternion algebra (defined by mutliplication table)
    GEN multtable = mkvec4(
    mkmat4(mkcol4s(1,0,0,0), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), mkcol4s(0,0,0,1)),
    mkmat4(mkcol4s(0,1,0,0), mkcol4s(-q_norm,0,0,0), mkcol4s(0,0,0,-1), mkcol4s(0,0,q_norm,0)),
    mkmat4(mkcol4s(0,0,1,0), mkcol4s(0,0,0,1), mkcol4(gneg(p_pari),gen_0,gen_0,gen_0), mkcol4(gen_0,gneg(p_pari),gen_0,gen_0)),
    mkmat4(mkcol4s(0,0,0,1), mkcol4s(0,0,-q_norm,0), mkcol4(gen_0,p_pari,gen_0,gen_0), mkcol4(gneg(gmul(stoi(q_norm),p_pari)),gen_0,gen_0,gen_0))

    );
    GEN B = alg_csa_table(nf, multtable, var,0);

    // The canonical maximal order (see KLPT, Lemmas 2-4 and https://ia.cr/2018/371, Prop. 1)
    // WARNING: Does not apply to p = 1 mod 12
    GEN B1 = B_1, B2, B3, B4;
    if (class_mod_4 == 3) {
      B2 = B_i;                                   // i
      B3 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
      B4 = mkcol4(gen_0,ghalf,ghalf,gen_0);       // (i+j)/2
    } else {
      B2 = mkcol4(ghalf,ghalf,gen_0,gen_0);                 // (1+i)/2
      B3 = mkcol4(gen_0,gen_0,ghalf,gneg(ghalf));           // j(1-i)/2
      int c = -1;
      B4 = mkcol4(gen_0,gdiv(gen_1,stoi(q_norm)),
		  gen_0,gneg(gdiv(stoi(c),stoi(q_norm))));  // (1-cj)i/q   (c = √-p mod q)
    }
    GEN B_O0 = alglathnf(B,mkmat4(B1,B2,B3,B4), gen_0);

    precomp_setup.p = p_pari;
    precomp_setup.B = B; // the quaternion algebra
    precomp_setup.qf = mkmat4(mkcol4s(1,0,0,0),
                             mkcol4s(0,q_norm,0,0),
                             mkcol4(gen_0,gen_0,p_pari,gen_0),
                             mkcol4(gen_0,gen_0,gen_0,gmul(p_pari,stoi(q_norm)))); // quadratic form defined by the reduced norm

    precomp_setup.torsion_fm = Z_factor_limit(strtoi(
        all_the_torsion_str
        ), 70000000);

    precomp_setup.O0 = B_O0; // the canonical maximal order
    precomp_setup.one = B_1;
    precomp_setup.i = B_i;
    precomp_setup.j = B_j;
    precomp_setup.ji = B_ji;

    precomp_setup.O0_b1 = B1;
    precomp_setup.O0_b2 = B2;
    precomp_setup.O0_b3 = B3;
    precomp_setup.O0_b4 = B4;
    precomp_setup.O0_to_standard = mkmat4(B1, B2, B3, B4);
    precomp_setup.standard_to_O0 = RgM_inv(precomp_setup.O0_to_standard);







    proj E, E_twist;
    mont_to_ted(&E, &precomp_setup.E0, false);
    mont_to_ted(&E_twist, &precomp_setup.E0, true);








    GEN m_frob_two, m_dist_two;
    GEN gelle, inv2, invq;

    point basis_ted[on_curve_len][3];
    proj basis[on_curve_len][3];
















    // two-torsion

    proj basis_two[3];
    find_basis_2e(basis_two, basis_two + 1, basis_two + 2, &precomp_setup.E0);

    proj2 basisxy_two[3];

    xtoxy(basisxy_two + 0, &precomp_setup.E0, basis_two + 0);
    xtoxy(basisxy_two + 1, &precomp_setup.E0, basis_two + 1);

    xyADD(basisxy_two + 2, &precomp_setup.E0, basisxy_two + 0, basisxy_two + 1);

    proj pt;
    xytox(&pt, &basisxy_two[2]);

    if (!mont_equal(&pt, &basis_two[2])) {
      xyNEG(basisxy_two, basisxy_two);
      xyADD(basisxy_two + 2, &precomp_setup.E0, basisxy_two + 0, basisxy_two + 1);
      xytox(&pt, &basisxy_two[2]);
    }
    assert (mont_equal(&pt, &basis_two[2]));



    bool correct_sign;
    compute_action_2e(&m_dist_two, mont0_dist, &basis_two[0], &basis_two[1], &basis_two[2], &precomp_setup.E0, two_tors_height);

    do {
        // printf("check_action_2e call\n");
        correct_sign = check_action_2e(m_dist_two, mont0_dist, montxy0_dist, &basisxy_two[0], &basisxy_two[1], &precomp_setup.E0, two_tors_height);
        // printf("check_action_2e done\n");
        if (!correct_sign) { m_dist_two = gmod(gneg(m_dist_two), powuu(2,two_tors_height));}
    } while (!correct_sign);


    compute_action_2e(&m_frob_two, mont0_frob, &basis_two[0], &basis_two[1], &basis_two[2], &precomp_setup.E0, two_tors_height);
    do {
        correct_sign = check_action_2e(m_frob_two, mont0_frob, montxy0_frob, &basisxy_two[0], &basisxy_two[1], &precomp_setup.E0, two_tors_height);
        if (!correct_sign) { m_frob_two = gmod(gneg(m_frob_two), powuu(2,two_tors_height));}
    } while (!correct_sign);

    gelle = powuu(2,two_tors_height);

    // WARNING: Does not apply to p = 1 mod 12
    GEN action_two_2, action_two_3, action_two_4;
    if (class_mod_4 == 3) {
      action_two_2 = m_dist_two;
      // TODO: how to handle division by 2?
      //action_two_3 = gmod(gmul(gsub(m1,gmul(m_frob_two,m_dist_two)),inv2),gelle); // (1-ji)/2
      //action_two_4 = gmod(gmul(gadd(m_dist_two,m_frob_two),inv2),gelle); //(i+j)/2
      action_two_4 = action_two_3 = action_two_2;
    } else {
      int c = -1;
      action_two_2 = gmod(gneg(m_dist_two),gelle); //(1+i)/2
      action_two_3 = gmod(gmul(action_two_2,m_frob_two),gelle); //(j-ji)/2
      action_two_4 = gmod(gadd(action_two_2,gmul(action_two_3,stoi(c))),gelle); //(1+i+cj-cji)/2
      action_two_4 = gmod(gmul(action_two_4,gen_2),gelle); // (1+i+cj-cji)
      action_two_4 = gmod(gsub(gsub(action_two_4,gmul(m_frob_two,stoi(c))),m1),gelle); // (i-ji)
      action_two_4 = gmod(gmul(action_two_4,Fp_inv(stoi(q_norm), gelle)),gelle); //(i-ji)/q
    }


    GEN m_frob, m_dist;


    GEN action_2[on_curve_len],action_3[on_curve_len],action_4[on_curve_len];


    proj basis_twist[on_twist_len][3];
    point basis_twist_ted[on_twist_len][3];

    GEN action_twist_2[on_twist_len],action_twist_3[on_twist_len],action_twist_4[on_twist_len];

    for (int i = 0; i < on_twist_len; i++) {
        long ell = on_twist_fact[i], e = on_twist_mult[i];

	// TODO: find a uniform way to handle this:
        //if (ell == q) e++; // a correction for the fact that constants.h only counts the q^(e-1)-torsion

        gelle = powuu(ell,e);
        inv2 = Fp_inv(gen_2, gelle);
        random_basis(&basis_twist[i][0], &basis_twist[i][1], &basis_twist_ted[i][0], &basis_twist_ted[i][1], &precomp_setup.E0, ell, e, true);


        point test;
        proj test2;
        uintbig ell_big;
        gentobig(&ell_big, powuu(ell,e));
        ted_mul(&test, &basis_twist_ted[i][0], &E_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_mul(&test, &basis_twist_ted[i][1], &E_twist, &ell_big);
        assert(ted_iszero(&test));

        ted_to_mont_point(&test2, &basis_twist_ted[i][0]);
        assert(!is_on_curve(&test2, &precomp_setup.E0));
        xMUL(&test2, &precomp_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));

        ted_to_mont_point(&test2, &basis_twist_ted[i][1]);
        assert(!is_on_curve(&test2, &precomp_setup.E0));
        xMUL(&test2, &precomp_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));


        ted_add(&basis_twist_ted[i][2], &E_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1]);
        ted_to_mont_point(&basis_twist[i][2], &basis_twist_ted[i][2]);

        ted_mul(&test, &basis_twist_ted[i][2], &E_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_to_mont_point(&test2, &basis_twist_ted[i][2]);
        assert(!is_on_curve(&test2, &precomp_setup.E0));
        xMUL(&test2, &precomp_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));


        fprintf(stderr, "\\\\ computing m_dist %ld^%ld\n",ell,e);
        compute_action(&m_dist, ted0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, ell, e);
        // printf("dist check\n");
        check_action(m_dist, ted0_dist, mont0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e);
        // printf("dist ok\n");

        fprintf(stderr, "\\\\ computing m_frob %ld^%ld\n",ell,e);
        compute_action(&m_frob, ted0_frob_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, ell, e);
        // printf("frob check\n");
        check_action(m_frob,       ted0_frob_twist, mont0_frob, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e);
        // printf("frob ok\n");

	// WARNING: Does not apply to p = 1 mod 12
	if (class_mod_4 == 3) {
	  action_twist_2[i] = m_dist;
	  action_twist_3[i] = gmod(gmul(gsub(m1,gmul(m_frob,m_dist)),inv2),gelle); // (1-ji)/2
	  action_twist_4[i] = gmod(gmul(gadd(m_dist,m_frob),inv2),gelle); //(i+j)/2
	} else {
	  int c = -1;
	  action_twist_2[i] = gmod(gneg(m_dist),gelle); //(1+i)/2
	  action_twist_3[i] = gmod(gmul(action_twist_2[i],m_frob),gelle); //(j-ji)/2
	  action_twist_4[i] = gmod(gadd(action_twist_2[i],gmul(action_twist_3[i],stoi(c))),gelle); //(1+i+cj-cji)/2
	  action_twist_4[i] = gmod(gmul(action_twist_4[i],gen_2),gelle); // (1+i+cj-cji)
	  action_twist_4[i] = gmod(gsub(gsub(action_twist_4[i],gmul(m_frob,stoi(c))),m1),gelle); // (i-cji)
	  if (ell != q_norm) {
            action_twist_4[i] = gmod(gmul(action_twist_4[i], Fp_inv(stoi(q_norm), gelle) ),gelle); //(i-cji)/q
	  }
	  else { // only compute the action on ell^(e-1)
	    // TODO: needs some fixes
            action_twist_2[i] = gmod(action_twist_2[i],powuu(ell,e-1));
            action_twist_3[i] = gmod(action_twist_3[i],powuu(ell,e-1));
            //output(gmod(action_twist_4[i],stoi(ell)));
            assert(gisexactzero(gmod(action_twist_4[i],stoi(ell))));
            action_twist_4[i] = gdiv(action_twist_4[i],stoi(ell));
            action_twist_4[i] = gmod(action_twist_4[i],powuu(ell,e-1));

            uintbig ell_big;
            uintbig_set(&ell_big,ell);
            ted_mul(&basis_twist_ted[i][0], &basis_twist_ted[i][0], &E_twist, &ell_big);
            ted_mul(&basis_twist_ted[i][1], &basis_twist_ted[i][1], &E_twist, &ell_big);
            ted_mul(&basis_twist_ted[i][2], &basis_twist_ted[i][2], &E_twist, &ell_big);
            xMUL(&basis_twist[i][0], &precomp_setup.E0, &basis_twist[i][0], &ell_big);
            xMUL(&basis_twist[i][1], &precomp_setup.E0, &basis_twist[i][1], &ell_big);
            xMUL(&basis_twist[i][2], &precomp_setup.E0, &basis_twist[i][2], &ell_big);

            check_action(gmod(m_frob,powuu(ell,e-1)), ted0_frob_twist, mont0_frob, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e-1);
            check_action(gmod(m_dist,powuu(ell,e-1)), ted0_dist, mont0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e-1);

	  }
	}
    }









    for (int i = 0; i < on_curve_len; i++) {
        long ell = on_curve_fact[i], e = on_curve_mult[i];

        gelle = powuu(ell,e);
        inv2 = Fp_inv(gen_2, gelle);
        if (ell != q_norm) invq = Fp_inv(stoi(q_norm), gelle);









        random_basis(&basis[i][0], &basis[i][1], &basis_ted[i][0], &basis_ted[i][1], &precomp_setup.E0, ell, e, false);
        ted_add(&basis_ted[i][2], &E, &basis_ted[i][0], &basis_ted[i][1]);
        ted_to_mont_point(&basis[i][2], &basis_ted[i][2]);


        point test;
        proj test2, test3;

        // TEST that mont_add works properly
        mont_add(&test2, &precomp_setup.E0, &basis[i][0], &basis[i][1]);
        ted_neg(&test, &basis_ted[i][1]);
        ted_add(&test, &E, &basis_ted[i][0], &test);
        ted_to_mont_point(&test3, &test);
        assert(mont_equal(&test2,&basis[i][2]) || mont_equal(&test2,&test3));
        // END TEST



        fprintf(stderr, "\\\\ computing m_frob %ld^%ld\n",ell,e);

        compute_action(&m_frob, ted0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, ell, e);



        check_action(m_frob, ted0_frob, mont0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e);



        fprintf(stderr, "\\\\ computing m_dist %ld^%ld\n",ell,e);
        compute_action(&m_dist, ted0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, ell, e);
        check_action(m_dist, ted0_dist, mont0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e);

	// WARNING: Does not apply to p = 1 mod 12
	if (class_mod_4 == 3) {
	  action_2[i] = m_dist;
	  action_3[i] = gmod(gmul(gsub(m1,gmul(m_frob,m_dist)),inv2),gelle); // (1-ji)/2
	  action_4[i] = gmod(gmul(gadd(m_dist,m_frob),inv2),gelle); //(i+j)/2
	} else {
	  int c = -1;
	  action_2[i] = gmod(gneg(m_dist),gelle); //(1+i)/2
	  action_3[i] = gmod(gmul(action_2[i],m_frob),gelle); //(j-ji)/2
	  action_4[i] = gmod(gadd(action_2[i],gmul(action_3[i],stoi(c))),gelle); //(1+i+cj-cji)/2
	  action_4[i] = gmod(gmul(action_4[i],gen_2),gelle); // (1+i+cj-cji)
	  action_4[i] = gmod(gsub(gsub(action_4[i],gmul(m_frob,stoi(c))),m1),gelle); // (i-ji)
	  if (ell != q_norm)
            action_4[i] = gmod(gmul(action_4[i],invq),gelle); //(i-ji)/q
	  else { // only compute the action on ell^(e-1)
            action_2[i] = gmod(action_2[i],powuu(ell,e-1));
            action_3[i] = gmod(action_3[i],powuu(ell,e-1));
            assert(gisexactzero(gmod(action_4[i],stoi(ell))));
            action_4[i] = gdiv(action_4[i],stoi(ell));
            action_4[i] = gmod(action_4[i],powuu(ell,e-1));

            uintbig ell_big;
            uintbig_set(&ell_big,ell);
            ted_mul(&basis_ted[i][0], &basis_ted[i][0], &E, &ell_big);
            ted_mul(&basis_ted[i][1], &basis_ted[i][1], &E, &ell_big);
            ted_mul(&basis_ted[i][2], &basis_ted[i][2], &E, &ell_big);
            xMUL(&basis[i][0], &precomp_setup.E0, &basis[i][0], &ell_big);
            xMUL(&basis[i][1], &precomp_setup.E0, &basis[i][1], &ell_big);
            xMUL(&basis[i][2], &precomp_setup.E0, &basis[i][2], &ell_big);


            check_action(gmod(m_frob,powuu(ell,e-1)), ted0_frob, mont0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e-1);
            check_action(gmod(m_dist,powuu(ell,e-1)), ted0_dist, mont0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e-1);

	  }
	}
    }
























    point basis_ted_sum[3], basis_twist_ted_sum[3];
    proj basis_sum[3], basis_twist_sum[3];

    ted_add(&basis_ted_sum[0], &E, &basis_ted[0][0], &basis_ted[1][0]);
    ted_add(&basis_ted_sum[1], &E, &basis_ted[0][1], &basis_ted[1][1]);
    for (int i = 2; i < on_curve_len; i++) {
        ted_add(&basis_ted_sum[0], &E, &basis_ted_sum[0], &basis_ted[i][0]);
        ted_add(&basis_ted_sum[1], &E, &basis_ted_sum[1], &basis_ted[i][1]);
    }
    ted_add(&basis_ted_sum[2], &E, &basis_ted_sum[0], &basis_ted_sum[1]);
    ted_to_mont_point(&basis_sum[0], &basis_ted_sum[0]);
    ted_to_mont_point(&basis_sum[1], &basis_ted_sum[1]);
    ted_to_mont_point(&basis_sum[2], &basis_ted_sum[2]);

    ted_add(&basis_twist_ted_sum[0], &E_twist, &basis_twist_ted[0][0], &basis_twist_ted[1][0]);
    ted_add(&basis_twist_ted_sum[1], &E_twist, &basis_twist_ted[0][1], &basis_twist_ted[1][1]);
    for (int i = 2; i < on_twist_len; i++) {
        ted_add(&basis_twist_ted_sum[0], &E_twist, &basis_twist_ted_sum[0], &basis_twist_ted[i][0]);
        ted_add(&basis_twist_ted_sum[1], &E_twist, &basis_twist_ted_sum[1], &basis_twist_ted[i][1]);
    }
    ted_add(&basis_twist_ted_sum[2], &E_twist, &basis_twist_ted_sum[0], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[0], &basis_twist_ted_sum[0]);
    ted_to_mont_point(&basis_twist_sum[1], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[2], &basis_twist_ted_sum[2]);


    // printing stuff
    printf("#include \"precomputed.h\"\n\n");

    printf("// each basis entry is a triple of the form P,Q,P+Q\n");
    printf("// this is initializing the point using the classical representation {0,...,p-1} for elements in GF(p).\n");
    printf("// We don't use this representation for actual computation but rather the montgomery representation (the conversion is made in init_precomputations using the fp_enc function)\n");
    printf("// hence the ***_uintbig[] defined below should not be used in any actual piece of code.\n\n");

    printf("const uintbig torsion_basis_sum_uintbig[3][2][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", proj_code(&basis_sum[0]), proj_code(&basis_sum[1]), proj_code(&basis_sum[2]));

    printf("const uintbig torsion_basis_twist_sum_uintbig[3][2][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", proj_code(&basis_twist_sum[0]), proj_code(&basis_twist_sum[1]), proj_code(&basis_twist_sum[2]));

    printf("\n");


    printf("const uintbig torsion_basis_ted_sum_uintbig[3][4][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", ted_code(&basis_ted_sum[0]), ted_code(&basis_ted_sum[1]), ted_code(&basis_ted_sum[2]));

    printf("const uintbig torsion_basis_twist_ted_sum_uintbig[3][4][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", ted_code(&basis_twist_ted_sum[0]), ted_code(&basis_twist_ted_sum[1]), ted_code(&basis_twist_ted_sum[2]));

    printf("\n");


    printf("const uintbig torsion_basis_uintbig[%ld][3][2][2] = {\n", on_curve_len);
    for (int i = 0; i < on_curve_len; i++) {
        printf("{ %s,\n %s,\n %s },\n", proj_code(&basis[i][0]), proj_code(&basis[i][1]), proj_code(&basis[i][2]));
    }
    printf("};\n\n");


    printf("const uintbig torsion_basis_twist_uintbig[%ld][3][2][2] = {\n", on_twist_len);
    for (int i = 0; i < on_twist_len; i++) {
        printf("{ %s,\n %s,\n %s },\n", proj_code(&basis_twist[i][0]), proj_code(&basis_twist[i][1]), proj_code(&basis_twist[i][2]));
    }
    printf("};\n\n");


    printf("const uintbig torsion_basis_two_uintbig[3][2][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", proj_code(&basis_two[0]), proj_code(&basis_two[1]), proj_code(&basis_two[2]));




    printf("\n\n");
    printf("proj torsion_basis[%ld][3];\n", on_curve_len);
    printf("proj torsion_basis_sum[3];\n");
    printf("point torsion_basis_ted_sum[3];\n");
    printf("proj torsion_basis_twist[%ld][3];\n", on_twist_len);
    printf("proj torsion_basis_twist_sum[3];\n");
    printf("point torsion_basis_twist_ted_sum[3];\n");
    printf("proj torsion_basis_two[3];\n");



    printf("\n\n");
    printf("void init_precomputations_generated() {\n");


    printf("\tglobal_setup.action_2 = malloc(%ld*sizeof(GEN));\n", on_curve_len);
    printf("\tglobal_setup.action_3 = malloc(%ld*sizeof(GEN));\n", on_curve_len);
    printf("\tglobal_setup.action_4 = malloc(%ld*sizeof(GEN));\n", on_curve_len);
    printf("\tglobal_setup.action_twist_2 = malloc(%ld*sizeof(GEN));\n", on_twist_len);
    printf("\tglobal_setup.action_twist_3 = malloc(%ld*sizeof(GEN));\n", on_twist_len);
    printf("\tglobal_setup.action_twist_4 = malloc(%ld*sizeof(GEN));\n", on_twist_len);
    printf("\n");


    for (int i = 0; i < on_curve_len; i++) {
        printf("\tglobal_setup.action_2[%d] = %s;\n", i, pari_2x2_matrix_code(action_2[i]));
        printf("\tglobal_setup.action_3[%d] = %s;\n", i, pari_2x2_matrix_code(action_3[i]));
        printf("\tglobal_setup.action_4[%d] = %s;\n", i, pari_2x2_matrix_code(action_4[i]));
    }
    for (int i = 0; i < on_twist_len; i++) {
        printf("\tglobal_setup.action_twist_2[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_2[i]));
        printf("\tglobal_setup.action_twist_3[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_3[i]));
        printf("\tglobal_setup.action_twist_4[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_4[i]));
    }

    printf("\tglobal_setup.action_two_2 = %s;\n", pari_2x2_matrix_code(action_two_2));
    printf("\tglobal_setup.action_two_3 = %s;\n", pari_2x2_matrix_code(action_two_3));
    printf("\tglobal_setup.action_two_4 = %s;\n", pari_2x2_matrix_code(action_two_4));
    // action_two_3_4(m_dist_two, m_frob_two, two_tors_height);





    print_loop("\t","torsion_basis_two","torsion_basis_two_uintbig");
    print_loop("\t","torsion_basis_sum","torsion_basis_sum_uintbig");
    print_loop("\t","torsion_basis_twist_sum","torsion_basis_twist_sum_uintbig");
    printf("\tfor (int j=0;j<%ld;j++){\n", on_curve_len);
    print_loop("\t\t","torsion_basis[j]","torsion_basis_uintbig[j]");
    printf("\t}\n");
    printf("\tfor (int j=0;j<%ld;j++){\n", on_twist_len);
    print_loop("\t\t","torsion_basis_twist[j]","torsion_basis_twist_uintbig[j]");
    printf("\t}\n");

  printf("\
    \tfor (int i=0;i<3;i++){\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->x)->re, &(torsion_basis_ted_sum_uintbig[i][0][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->x)->im, &(torsion_basis_ted_sum_uintbig[i][0][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->y)->re, &(torsion_basis_ted_sum_uintbig[i][1][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->y)->im, &(torsion_basis_ted_sum_uintbig[i][1][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->z)->re, &(torsion_basis_ted_sum_uintbig[i][2][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->z)->im, &(torsion_basis_ted_sum_uintbig[i][2][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->t)->re, &(torsion_basis_ted_sum_uintbig[i][3][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->t)->im, &(torsion_basis_ted_sum_uintbig[i][3][1]) );\n\
  \t}\n\
  \tfor (int i=0;i<3;i++){\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->re, &(torsion_basis_twist_ted_sum_uintbig[i][0][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->im, &(torsion_basis_twist_ted_sum_uintbig[i][0][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->re, &(torsion_basis_twist_ted_sum_uintbig[i][1][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->im, &(torsion_basis_twist_ted_sum_uintbig[i][1][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->re, &(torsion_basis_twist_ted_sum_uintbig[i][2][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->im, &(torsion_basis_twist_ted_sum_uintbig[i][2][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->re, &(torsion_basis_twist_ted_sum_uintbig[i][3][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->im, &(torsion_basis_twist_ted_sum_uintbig[i][3][1]) );\n\
  \t}\n\n");




  //printf("/***** INSERT HERE INITIALISATION OF PARI VALUES *****/\n\n");





    printf("}\n\n");



    return 0;
}
