#define _unused(x) ((void)(x))

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pari/pari.h>
#include <assert.h>


#include "klpt.h"
#include "ideal.h"
#include "toolbox.h"
#include "constants.h"

//todo: remove this
//#include "precomputed.h"


GEN klpt_strong_approximation(GEN A, GEN order, GEN p, GEN N, GEN beta, GEN gamma, GEN L, GEN lambda) {
    pari_sp ltop = avma;

    long q = itos_or_0(algnorm(A,mkcol4s(0,1,0,0),0));
    assert(q != 0);
    GEN C0 = gel(beta,3);
    GEN D0 = gel(beta,4);

    // the following should not happen: make sure in klpt that both C0 and D0 are invertible
    int swap = 0;
    if (gcmp(ggcd(D0,N),gen_1) != 0) {
        if (q != 1) {
            avma = ltop;
            return NULL;
        }
        swap = 1;
        C0 = gel(beta,4);
        D0 = gel(beta,3);
    }

    //register the constant denominator of gamma
    GEN Xg;
    alglatcontains(A,order,gamma,&Xg);
    GEN ng = content(Xg);


    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda_2 = gmul(gmul(p,lambda),gen_2);
    GEN coeff_c = gmul(p_lambda_2,C0);
    GEN coeff_d = gmul(gmul(p_lambda_2,D0), stoi(q));
    GEN cst_term = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), N); // (L-lambda^2*beta_norm)/N

    GEN coeff_d_inv = Fp_inv(coeff_d,N);

    GEN cp0 = gen_0;
    GEN dp0 = Fp_mul(cst_term,coeff_d_inv,N);

    GEN lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c,coeff_d_inv, N))), mkcol2(gen_0,N)),N);
    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));

    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(4000), 2, 10); // returns at most 4000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v);
        }

        diff = gsub(target, close_v);
        // GEN norm = gmul(RgV_dotproduct(diff,diff), p);
        // TODO: this is a quick fix for q ≠ 1, it can be done better by addapting the quadratic form
        GEN norm = gmul(gadd(gsqr(gel(diff,1)), gmul(gsqr(gel(diff,2)), stoi(q))),p);
        if (gcmp(norm, L) <= 0) {
            GEN rhs = diviiexact(gsub(L, norm), gsqr(N));

            GEN ap, bp;
            if ((q != 1 && ispseudoprime(rhs,0) && cornacchia(stoi(q), rhs, &ap, &bp))
                || (q == 1 && cornacchia_extended(rhs, &ap, &bp) )) {
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);
                if (swap) {
                    c = gel(diff,2);
                    d = gel(diff,1);
                }
                GEN betap = mkcol4(a,b,c,d);

                if ( Z_lval(L,2) ==0  ) {
                  return gerepilecopy(ltop, betap);
                }
                else {
                  GEN X;
                  alglatcontains(A,order,algmul(A,gamma,betap),&X);
                  GEN n =content(X);
                  //checking that the constant factor is of the desired size to ensure fixed length
                  if (  Z_lval(n,2) == 1 + Z_lval(ng,2) ) {
                    return gerepilecopy(ltop, betap);
                  }
                  else {
                    //try to swap the values to improve the rate of success
                    //warning: it only works if p=3 mod 4 and we are not in the extremal case
                    GEN betap = mkcol4(b,a,c,d);
                    GEN X;
                    alglatcontains(A,order,algmul(A,gamma,betap),&X);
                    GEN n =content(X);
                    //checking that the constant factor is of the desired size to ensure fixed length
                    if (  Z_lval(n,2) == 1 + Z_lval(ng,2) ) {
                      return gerepilecopy(ltop, betap);
                    }

                  }

                }


            }
        }
    }

    avma = ltop;
    return NULL;
}


GEN klpt_strong_approximation_extremal(GEN A, p_extremal_maximal_order extremal, GEN order, GEN p, GEN N, GEN beta, GEN C0, GEN D0, GEN gamma, GEN L, GEN lambda) {
    pari_sp ltop = avma;

    // we made sure in klpt that both C0 and D0 are invertible
    assert(gequal1(ggcd(D0,N)));

    GEN q = stoi(extremal.q);

    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda_2 = gmul(gmul(p,lambda),gen_2);
    GEN coeff_c = gmul(p_lambda_2,C0);
    GEN coeff_d = gmul(gmul(p_lambda_2,D0), q);
    GEN cst_term = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), N); // (L-lambda^2*beta_norm)/N

    GEN coeff_d_inv = Fp_inv(coeff_d,N);

    GEN cp0 = gen_0;
    GEN dp0 = Fp_mul(cst_term,coeff_d_inv,N);

    GEN lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c,coeff_d_inv, N))), mkcol2(gen_0,N)),N);
    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));

    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(2000), 2, 10); // returns at most 2000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v);
        }

        diff = gsub(target, close_v);
        // norm = gmul(RgV_dotproduct(diff,diff), p);
        // TODO: this is a quick fix for q ≠ 1, it can be done better by addapting the quadratic form
        GEN norm = gmul(gadd(gsqr(gel(diff,1)), gmul(gsqr(gel(diff,2)), q)),p);

        if (gcmp(norm, L) <= 0) {
            GEN rhs = diviiexact(gsub(L, norm), gsqr(N));

            GEN ap, bp;
            if ((extremal.q != 1 && ispseudoprime(rhs,0) && cornacchia(q, rhs, &ap, &bp))
                || (extremal.q == 1 && cornacchia_extended(rhs, &ap, &bp) )) {
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);

                // GEN betap = mkcol4(a,b,c,d);

                GEN betap = mkcol4(a,gen_0,gen_0,gen_0);
                betap = gadd(betap,gmul(extremal.i,b));
                betap = gadd(betap,gmul(extremal.j,c));
                betap = gadd(betap,gmul(algmul(A,extremal.j,extremal.i),d));

                if ( Z_lval(L,2) ==0  ) {
                  return gerepilecopy(ltop, betap);
                }
                else {
                  GEN X;
                  GEN Xg;
                  alglatcontains(A,order,gamma,&Xg);
                  GEN ng = content(Xg);

                  alglatcontains(A,order,algmul(A,gamma,betap),&X);
                  GEN n =content(X);
                  //checking that the constant factor is of the desired length to ensure fixed length
                  if (  Z_lval(n,2) == 2 + Z_lval(ng,2) ) {
                    return gerepilecopy(ltop, betap);
                  }
                }


            }
        }
    }

    avma = ltop;
    return NULL;
}



//same as above but looks for solution where the coefficients might be half integers
GEN klpt_strong_approximation_divided_by_2(GEN A, GEN p, GEN N, GEN beta, GEN L, GEN lambda) {
    pari_sp ltop = avma;


    GEN C0 = gel(beta,3);
    GEN D0 = gel(beta,4);

    // the following should not happen: make sure in klpt that both C0 and D0 are invertible
    int swap = 0;
    if (gcmp(ggcd(D0,N),gen_1) != 0) {
        swap = 1;
        C0 = gel(beta,4);
        D0 = gel(beta,3);
    }

    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda_2 = gmul(gmul(p,lambda),gen_2);
    GEN coeff_c = gmul(p_lambda_2,C0);
    GEN coeff_d = gmul(p_lambda_2,D0);
    GEN cst_term = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), N); // (L-lambda^2*beta_norm)/N

    GEN coeff_d_inv = Fp_inv(coeff_d,N);

    GEN cp0 = gen_0;
    GEN dp0 = Fp_mul(cst_term,coeff_d_inv,N);

    GEN lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c,coeff_d_inv, N))), mkcol2(gen_0,N)),N);
    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));
    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(2000), 2, 10); // returns at most 2000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v);
        }

        diff = gsub(target, close_v);
        distance = gmul(RgV_dotproduct(diff,diff), p);
        if (gcmp(distance, L) <= 0 && (!gequal( gmod(gel(diff,1),gen_2),gmod(gel(diff,2),gen_2)))) {
            GEN rhs = diviiexact(gsub(L, distance), gsqr(N));

            GEN ap, bp;
            if (cornacchia_extended(rhs, &ap, &bp) ){
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);
                if (swap) {
                    c = gel(diff,2);
                    d = gel(diff,1);
                }

                // this condition should ensure that (a,b,c,d)/2 is inside O0 but not inside the special eichler order of level 2 that we want to avoid
                // if ( (!gequal( gmod(a,gen_2),gmod(b,gen_2)))) {

                  if (gequal( gmod(a,gen_2),gmod(d,gen_2))) {

                    GEN betap = mkcol4(gdiv(a,gen_2),gdiv(b,gen_2),gdiv(c,gen_2),gdiv(d,gen_2));
                    return gerepilecopy(ltop, betap);
                  }
                  else {
                    GEN a_temp = a;
                    a = b;
                    b = a_temp;
                    GEN betap = mkcol4(gdiv(a,gen_2),gdiv(b,gen_2),gdiv(c,gen_2),gdiv(d,gen_2));
                    return gerepilecopy(ltop, betap);
                  }

                // }


                //  if ( (!gequal( gmod(a,gen_2),gmod(b,gen_2)))  && (gequal( gmod(a,gen_2),gmod(d,gen_2))) && (gequal( gmod(b,gen_2),gmod(c,gen_2)))  )
                // {
                //   GEN betap = mkcol4(gdiv(a,gen_2),gdiv(b,gen_2),gdiv(c,gen_2),gdiv(d,gen_2));
                //   return gerepilecopy(ltop, betap);
                // }
                // else {
                //   // printf("found values but bad mod 2 \n");
                //   // output(gmod(a,gen_2));
                //   // output(gmod(b,gen_2));
                //   // output(gmod(c,gen_2));
                //   // output(gmod(d,gen_2));
                // }
                // if ( Z_lval(L,2) ==0  ) {
                //   return gerepilecopy(ltop, betap);
                // }
                // else {
                //   GEN X;
                //   GEN Xg;
                //   alglatcontains(A,order,gamma,&Xg);
                //   GEN ng = content(Xg);
                //   alglatcontains(A,order,algmul(A,gamma,betap),&X);
                //   GEN n =content(X);
                //   //checking that the constant factor is of the desired length to ensure fixed length
                //   if (  Z_lval(n,2) == 2 + Z_lval(ng,2) ) {
                //     return gerepilecopy(ltop, betap);
                //   }
                // }


            }
        }
    }

    avma = ltop;
    return NULL;
}

//same as above but looks for solution where the coefficients might be half integers
GEN klpt_strong_approximation_divided_by_2_extremal(GEN A, p_extremal_maximal_order extremal, GEN p, GEN N, GEN beta, GEN C0, GEN D0, GEN L, GEN lambda) {
    pari_sp ltop = avma;


    // we made sure in klpt that both C0 and D0 are invertible
    assert(gequal1(ggcd(D0,N)));

    GEN q = stoi(extremal.q);

    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda_2 = gmul(gmul(p,lambda),gen_2);
    GEN coeff_c = gmul(p_lambda_2,C0);
    GEN coeff_d = gmul(gmul(p_lambda_2,D0), q);
    GEN cst_term = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), N); // (L-lambda^2*beta_norm)/N

    GEN coeff_d_inv = Fp_inv(coeff_d,N);

    GEN cp0 = gen_0;
    GEN dp0 = Fp_mul(cst_term,coeff_d_inv,N);

    GEN lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c,coeff_d_inv, N))), mkcol2(gen_0,N)),N);
    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));
    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(2000), 2, 10); // returns at most 2000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v);
        }

        diff = gsub(target, close_v);
        // norm = gmul(RgV_dotproduct(diff,diff), p);
        // TODO: this is a quick fix for q ≠ 1, it can be done better by addapting the quadratic form
        GEN norm = gmul(gadd(gsqr(gel(diff,1)), gmul(gsqr(gel(diff,2)), q)),p);
        if (gcmp(norm, L) <= 0 && (!gequal( gmod(gel(diff,1),gen_2),gmod(gel(diff,2),gen_2)))) {
            GEN rhs = diviiexact(gsub(L, norm), gsqr(N));

            GEN ap, bp;
            if ((extremal.q != 1 && ispseudoprime(rhs,0) && cornacchia(q, rhs, &ap, &bp))
                || (extremal.q == 1 && cornacchia_extended(rhs, &ap, &bp) )) {
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);


                GEN betap = mkcol4(a,gen_0,gen_0,gen_0);
                betap = gadd(betap,gmul(extremal.i,b));
                betap = gadd(betap,gmul(extremal.j,c));
                betap = gadd(betap,gmul(algmul(A,extremal.j,extremal.i),d));

                betap = gdiv(betap,gen_2);

                if(alglatcontains(A,extremal.order,betap,NULL)){
                    return gerepilecopy(ltop, betap);
                }

            }
        }
    }

    avma = ltop;
    return NULL;
}





// as above but N is a power of 2
GEN klpt_strong_approximation_2e(GEN A, GEN p, GEN N, GEN beta, GEN L, GEN lambda) {
    pari_sp ltop = avma;

    GEN C0 = gel(beta,3);
    GEN D0 = gel(beta,4);


    long q = itos_or_0(algnorm(A,mkcol4s(0,1,0,0),0));
    assert(q != 0);

    int swap = 0;
    // if (gcmp(ggcd(D0,N),gen_1) != 0) {
    //     if (q > 2) {
    //         return NULL;
    //     }
    //     swap = 1;
    //     C0 = gel(beta,4);
    //     D0 = gel(beta,3);
    // }

    // if (q == 2 && gcmp(ggcd(C0,N),gen_1) != 0) {
    //     return NULL;
    // }

    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda = gmul(p,lambda);
    GEN coeff_c_half = gmul(p_lambda,C0);
    GEN coeff_d_half = gmul(gmul(p_lambda,D0),stoi(q));

    GEN cp0, dp0,lattice_basis;
    GEN cst_term_half = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), gmul(N,gen_2)); // (L-lambda^2*beta_norm)/(2*N)

    if (gcmp(ggcd(C0,N),gen_1) == 0) {
        GEN coeff_c_half_inv = Fp_inv(coeff_c_half,N);
        cp0 = Fp_mul(cst_term_half,coeff_c_half_inv,N);
        dp0 = gen_0;
        lattice_basis = gmul(mkmat2(mkcol2(Fp_mul(coeff_d_half,coeff_c_half_inv, N), gneg(gen_1)), mkcol2(N,gen_0)),N);

    }
    else if (q != 2 && gcmp(ggcd(D0,N),gen_1) == 0) {
        GEN coeff_d_half_inv = Fp_inv(coeff_d_half,N);
        dp0 = Fp_mul(cst_term_half,coeff_d_half_inv,N);
        cp0 = gen_0;
        lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c_half,coeff_d_half_inv, N))), mkcol2(gen_0,N)),N);
    }
    else if (q != 2 && gcmp(ggcd(C0,N),ghalf) == 0 && gcmp(ggcd(D0,N),ghalf) == 0) {
        GEN coeff_d_inv = Fp_inv(gmul(coeff_d_half,gen_2),N);
        dp0 = Fp_mul(gmul(cst_term_half,gen_2),coeff_d_inv,N);
        cp0 = gen_0;
        lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(gmul(coeff_c_half,gen_2),coeff_d_inv, N))), mkcol2(gen_0,N)),N);
    }
    else {
        return NULL;
    }

    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));

    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(2000), 2, 10); // returns at most 2000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v );
        }

        diff = gsub(target, close_v);
        // GEN norm = gmul(RgV_dotproduct(diff,diff), p);
        GEN norm = gmul(gadd(gsqr(gel(diff,1)), gmul(gsqr(gel(diff,2)), stoi(q))),p);



        if (gcmp(norm, L) <= 0) {
            GEN rhs = diviiexact(gsub(L, norm), gsqr(N));

            GEN ap, bp;
            if ((q != 1 && ispseudoprime(rhs,0) && cornacchia(stoi(q), rhs, &ap, &bp))
                || (q == 1 && cornacchia_extended(rhs, &ap, &bp) )) {
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);
                if (swap) {
                    c = gel(diff,2);
                    d = gel(diff,1);
                }

                GEN betap = mkcol4(a,b,c,d);

                return gerepilecopy(ltop, betap);
            }
        }
    }

    avma = ltop;
    return NULL;
}

// compute beta in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
GEN klpt_solve_beta(GEN A, GEN gamma, GEN J, GEN N) {
    pari_sp ltop = avma;
    GEN A_i = mkcol4s(0,1,0,0);
    GEN A_j = mkcol4s(0,0,1,0);

    GEN J_basis = lideal_basis(J);
    GEN J_scalar = lideal_scalar(J);
    GEN gamma_j = algmul(A,gamma,A_j);

    GEN matsys = gmul(mkmat2(gamma_j, algmul(A,gamma_j, A_i)), Q_denom(J_scalar));
    matsys = gconcat(matsys, gmul(J_basis, Q_remove_denom(J_scalar, NULL)));

    GEN ker;

    if (mpodd(N)) // prime case
        ker = matkermod(matsys, N, NULL); // flag = 1 because integral entries
    else // power of 2 case
        ker = matkermod(matsys, gmul(N,Q_denom(J_scalar)), NULL);

    unsigned long i = lg(ker)-1;
    GEN C0,D0;
    do {
        C0 = gel(gel(ker, i), 1);
        D0 = gel(gel(ker, i), 2);

        i--;
        if (i < 1) break;
    } while ((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0));

    if ((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0)) {
        avma = ltop;
        return NULL;
    }

    GEN beta = mkcol4(gen_0, gen_0, C0, D0); // beta in span(j,j*i) and gamma*beta in J
    return gerepilecopy(ltop, beta);
}


// compute beta in span(j,j*i) such that there is a zeta \in (Z/N(J)Z)^* with
// gamma*beta*delta - zeta in J*N
// N and N(J) are coprime
GEN klpt_solve_beta_special(GEN A, GEN gamma, GEN delta, GEN N, GEN J) {
    pari_sp ltop = avma;
    GEN A_1 = mkcol4s(1,0,0,0);
    GEN A_i = mkcol4s(0,1,0,0);
    GEN A_j = mkcol4s(0,0,1,0);
    GEN NJ = lideal_norm(J);

    GEN J_basis = lideal_basis(J);
    GEN J_scalar = lideal_scalar(J);

    GEN gamma_j = algmul(A,gamma,A_j);

    GEN gamma_j_delta = algmul(A,gamma_j,delta);
    GEN gamma_j_i_delta = algmul(A,algmul(A,gamma_j,A_i),delta);


    GEN matsys = gmul(mkmat3(gamma_j_delta, gamma_j_i_delta, A_1), Q_denom(J_scalar));

    matsys = gconcat(matsys, gmul(J_basis, gmul(Q_remove_denom(J_scalar, NULL),N)));


    GEN ker = matkermod(matsys, NJ, NULL); // flag = 1 because integral entries

    unsigned long i = lg(ker)-1;
    GEN C0,D0, zeta;
    do {
        C0 = gel(gel(ker, i), 1);
        D0 = gel(gel(ker, i), 2);
        zeta = gel(gel(ker, i), 3);
        i--;
        if (i < 1) break;
    } while (gcmp(zeta,gen_0) == 0);

    if (gcmp(zeta,gen_0) == 0) {
        avma = ltop;
        return NULL;
    }

    GEN beta = mkcol4(gen_0, gen_0, C0, D0); // beta in span(j,j*i) and gamma*beta in J
    return gerepilecopy(ltop, beta);
}


// same as above but using the input p-extremal order
// compute beta in span(j,j*i) such that there is a zeta \in (Z/N(J)Z)^* with
// gamma*beta*delta + zeta in J*N
// N and N(J) are coprime
// return beta, C0 and D0 with beta = C0*j + D0*j*i
GEN klpt_solve_beta_extremal(GEN A, p_extremal_maximal_order extremal, GEN gamma, GEN delta, GEN N, GEN J) {
    pari_sp ltop = avma;
    GEN A_1 = mkcol4s(1,0,0,0);
    GEN A_i = extremal.i;
    GEN A_j = extremal.j;
    GEN A_ji = algmul(A,A_j,A_i);
    GEN NJ = lideal_norm(J);

    GEN J_basis = lideal_basis(J);
    GEN J_scalar = lideal_scalar(J);

    GEN gamma_j = algmul(A,gamma,A_j);

    GEN gamma_j_delta = algmul(A,gamma_j,delta);
    GEN gamma_j_i_delta = algmul(A,algmul(A,gamma_j,A_i),delta);


    GEN matsys = gmul(mkmat3(gamma_j_delta, gamma_j_i_delta, A_1), Q_denom(J_scalar));

    matsys = gconcat(matsys, gmul(J_basis, gmul(Q_remove_denom(J_scalar, NULL),N)));


    GEN ker = matkermod(matsys, NJ, NULL); // flag = 1 because integral entries

    unsigned long i = lg(ker)-1;
    GEN C0,D0, zeta;
    do {
        C0 = gel(gel(ker, i), 1);
        D0 = gel(gel(ker, i), 2);
        zeta = gel(gel(ker, i), 3);
        i--;
        if (i < 1) break;
    } while (gcmp(zeta,gen_0) == 0);

    if (gcmp(zeta,gen_0) == 0) {
        // printf("No beta found (special)\n");
        avma = ltop;
        return NULL;
    }

    GEN beta = gadd(gmul(C0,A_j), gmul(D0,A_ji));//   mkcol4(gen_0, gen_0, C0, D0); // beta in span(j,j*i) and gamma*beta in J
    // output(beta);
    //tests
    GEN zeta_alg=mkcol4(zeta,gen_0,gen_0,gen_0);
    zeta=Fp_mul(zeta,gen_1,NJ);
    GEN beta_test=algmul(lideal_algebra(J),beta,delta);
    beta_test=algmul(lideal_algebra(J),gamma,beta_test);

    GEN delta_test=algadd(lideal_algebra(J),beta_test,zeta_alg);
    if (!alglatcontains(lideal_algebra(J),gel(J,1),delta_test,NULL) ) {
      printf("the found solution is not contained in the eichler order \n");
    }

    GEN res = mkvec3(beta,C0,D0);

    return gerepilecopy(ltop, res);
}

// runs KLPT for the left ideal I in the special order of the quaternion algebra A
// the result is an ideal equivalent to I of norm dividing the integer whose factorisation matrix is fm
// Assumes the basis of A is 1, i, j, j*i, where i^2 = -1 and j^2 = -p
GEN klpt_special_smooth_given_nearprime_J(GEN J0, GEN L0, GEN N, GEN fm) {
    pari_sp ltop = avma;

    GEN A = lideal_algebra(J0);
    GEN order = lideal_order(J0);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);
    long q = itos_or_0(algnorm(A,mkcol4s(0,1,0,0),0));

    GEN J = lideal_create(A, order, lideal_generator(J0), N);
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    // try 100 times then abandon
    for (int n = 1; n <= 100; ++n) {

        // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
        GEN fm_1, L1;

        GEN coeff = NULL;
        unsigned long ctr = 1;

        while(!coeff) {
            fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/10) + 3),4);
            L1 = famat_prod(fm_1);
            if (q == 1) coeff = norm_equation_special(p, gmul(L1,N), 0, false);
            else coeff = norm_equation_special_q(p,q,gmul(L1,N));
            ctr += 1;
        }

        GEN gamma = gtrans(coeff); // norm N*L1
        GEN fm_rem = famat_div(fm,fm_1); // remaining factors

        // compute beta in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        GEN beta = klpt_solve_beta(A, gamma, J, N);
        if (!beta) continue;
        GEN beta_norm = algnorm(A,beta,0);

        if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
            /// A problem with this beta!
            continue;
        }

        // generate L2 such that L2/beta_norm is a square modulo N
        int beta_norm_is_square = Fp_issquare(beta_norm,N), L2_is_square;
        GEN bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(80)), fm_2, L2; // (p*N^3)*MARGIN

        int safety_ctr = 0;
        GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
        GEN betap = NULL;

        do {

          do {
              fm_2 = famat_smallest(fm_rem, bound_L2);
              L2 = famat_prod(fm_2);
              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;
          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 1));
          do {
              fm_2 = famat_random(fm_rem, bound_L2,3);
              L2 = famat_prod(fm_2);
              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;
          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 5));

            do {
                fm_2 = famat_random(fm_rem, bound_L2,0);
                L2 = famat_prod(fm_2);
                L2_is_square = Fp_issquare(L2,N);
                safety_ctr++;
            } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 40));
            if (safety_ctr >= 40) { continue; } // not enough entropy in fm_rem to ensure a solution

            list_L2 = vec_append(list_L2,L2);

            // strong approximation
            GEN lambda = Fp_sqrt(Fp_div(L2,beta_norm,N), N);

            betap = klpt_strong_approximation(A, order, p, N, beta,gamma, L2, lambda);

            if (betap) {

                // TODO: implement and use lideal_mul_elt instead
                GEN norm = gmul(gmul(L1,L2),L0);
                GEN alpha = algmul(A,gamma,betap);
                GEN generator = algmul(A,lideal_generator_coprime(J0,norm),alg_conj(A,alpha));
                GEN newideal = lideal_create(A, order, generator, norm);

                return gerepilecopy(ltop, newideal);
            }
        } while (!betap && (safety_ctr < 40));

    }


    avma = ltop;
    return NULL;
}

GEN klpt_special_smooth(GEN I, GEN fm) {
    pari_sp ltop = avma;


    // find an equivalent nearprime ideal
    GEN J0 = lideal_equiv_nearprime(I,fm,0);

    GEN fm_0 = famat_Z_gcd(fm, lideal_norm(J0));
    GEN L0 = famat_prod(fm_0);
    GEN N = diviiexact(lideal_norm(J0), L0);

    GEN fm_without_L0 = famat_div(fm,fm_0);

    GEN klpt_sol = klpt_special_smooth_given_nearprime_J(J0,L0,N,fm_without_L0);

    long ctr = 0;
    while (!klpt_sol && ctr < 50) { // fallback strategy if first J0 failed (negligible probability)
        do {
            ctr++;
            // TODO: implement and use a randomised lideal_equiv_nearprime instead
            J0 = lideal_equiv_prime_random(I,NULL,stoi(5+ctr));
        } while (!J0);

        N = lideal_norm(J0);
        klpt_sol = klpt_special_smooth_given_nearprime_J(J0,gen_1,N,fm);
    }

    if (!klpt_sol) { avma = ltop; fprintf(stderr, "klpt_special_smooth did not find a solution\n"); return NULL; }
    else return gerepilecopy(ltop, klpt_sol);
}


GEN klpt_special_smooth_given_2e_J(GEN J, GEN fm) {
    pari_sp ltop = avma;

    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);

    GEN N = lideal_norm(J);
    GEN N_2 = gmul(N,gen_2);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);
    GEN p_div_N = gadd(truedivii(p, N), gen_1);
    long q = itos_or_0(algnorm(A,mkcol4s(0,1,0,0),0));


    // try 100 times then abandon
    for (int n = 1; n <= 100; ++n) {

        // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
        GEN fm_1, L1;

        GEN coeff = NULL;
        unsigned long ctr = 1;

        GEN O2 = alglatmul(A, order, alg_scalar(A, gen_2));
        GEN O4 = alglatmul(A, order, alg_scalar(A, stoi(4)));
        GEN K = lideal_create(A, order, mkcol4s(0,1,0,0), stoi(2));
        K = alglatmul(A, lideal_lattice(K), alg_scalar(A, stoi(2)));
        while(!coeff) {
            fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/20) + 1),0);
            L1 = famat_prod(fm_1);
            // coeff = norm_equation_special(p, gmul(L1,N), 1, false);
            if (q == 1) coeff = norm_equation_special(p, gmul(L1,N), 1, false);
            else if (q == 2) {
                // TODO: this is some unoptimised gymnastics to find a primitive solution when q ≠ 1
                coeff = norm_equation_special_q(p,q,gmul(gmul(L1,N),stoi(4)));
                if (coeff && alglatcontains(A, O2, gtrans(coeff), NULL)) {
                    if (q == 2 && !alglatcontains(A, K, gtrans(coeff), NULL)) {
                        //output(coeff);
                        coeff = NULL;
                    }
                    else if (alglatcontains(A, O4, gtrans(coeff), NULL)) {
                        //output(coeff);
                        coeff = NULL;
                    }
                    else {
                        coeff = gdiv(coeff, gen_2);
                    }
                }
                else { coeff = NULL; }
            }
            else {
                coeff = norm_equation_special_q(p,q,gmul(L1,N));
                if (coeff && alglatcontains(A, O2, gtrans(coeff), NULL)) {
                    coeff = NULL;
                }
            }
            ctr += 1;
        }
        GEN gamma = gtrans(coeff); // norm N*L1
        assert(gcmp(algnorm(A,gamma,0), gmul(L1,N)) == 0);
        assert(alglatcontains(A, order, gamma, NULL));
        assert(!alglatcontains(A, O2, gamma, NULL));
        GEN fm_rem = famat_div(fm,fm_1); // remaining factors
        // compute beta in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        GEN beta = klpt_solve_beta(A, gamma, J, N);
        if (!beta) continue;
        int beta_is_good = true;
        GEN beta_half = gdiv(beta, gen_2);
        while (alglatcontains(A, order, beta_half, NULL)) { // beta is divisible by 2
            // printf("HALF\n");
            beta = beta_half;
            beta_half = gdiv(beta, gen_2);
            if (!alglatcontains(A, lideal_lattice(J), algmul(A,gamma,beta), NULL)) {
                beta_is_good = false;
                break;
            }
        }
        if (!beta_is_good) continue;
        alglatcontains(A, order, beta, NULL);
        assert(alglatcontains(A, lideal_lattice(J), algmul(A,gamma,beta), NULL));

        GEN beta_norm = algnorm(A,beta,0);
        if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
            // output(beta);
            // output(ggcd(beta_norm, N));
            printf("Problem with beta!\n");
            continue;
        }
        // beta_norm is 1 mod 4 (exactly one coefficient of beta is odd)

        // generate L2 such that L2/beta_norm is a square modulo N
        GEN lambda = NULL;
        GEN bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(80)), fm_2, L2; // (p*N^3)*MARGIN

        int safety_ctr = 0;
        GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
        GEN betap = NULL;

        do {
            do {
                fm_2 = famat_random(fm_rem, bound_L2,0);
                L2 = famat_prod(fm_2);
                lambda = Zn_sqrt(Fp_div(L2,beta_norm,N_2), N_2);  // lambda is found mod 2*N instead of N

                safety_ctr++;
            } while (((!lambda) || RgV_isin(list_L2, L2)) && (safety_ctr < 40));
            if (safety_ctr >= 40) { continue; } // not enough entropy in fm_rem to ensure a solution

            list_L2 = vec_append(list_L2,L2);
            // strong approximation
            GEN betap = klpt_strong_approximation_2e(A, p, N, beta, L2, lambda);
            if (betap) {
                // TODO: implement and use lideal_mul_elt instead
                GEN norm = gmul(L1,L2);
                GEN alpha = algmul(A,gamma,betap);
                GEN generator = algmul(A,lideal_generator_coprime(J,norm),alg_conj(A,alpha));
                GEN newideal = lideal_create(A, order, generator, norm);

                return gerepilecopy(ltop, newideal);
            }
        }
        while ((!betap) && (safety_ctr < 40));
    }


    avma = ltop;
    return NULL;
}

GEN klpt_special_smooth_small_2e_input(GEN I, GEN fm){
    pari_sp ltop = avma;

    GEN A = lideal_algebra(I);
    GEN order = lideal_order(I);
    GEN J = I;

    long q = itos_or_0(algnorm(A,mkcol4s(0,1,0,0),0));

    if (q == 1 || q == 2) {
        GEN end_of_norm_2 = mkcol4s(1,1,0,0);
        if (q == 2) end_of_norm_2 = mkcol4s(0,1,0,0);

        GEN I_red = lideal_create(A, order, lideal_generator(I), gen_2);
        GEN O_2 = lideal_create(A, order, end_of_norm_2, gen_2);

        // we make sure the input ideal is a fixed point for the action of (R/2R)^*
        if (!lideal_equals(I_red, O_2)) {
            J = lideal_create(A, order, algmul(A,lideal_generator(I),end_of_norm_2), gmul(lideal_norm(I),gen_2));
        }

    }
    // else if (q == 3) {
    //     GEN automorphism = mkcol4(ghalf,ghalf,gen_0,gen_0);

    //     GEN I_red = lideal_create(A, order, lideal_generator(I), stoi(4));
    //     GEN O_4 = lideal_create(A, order, automorphism, stoi(4));


    //     output(algnorm(A,automorphism,0));
    //     output(lideal_norm(O_4));

    //     // we make sure the input ideal is a fixed point for the action of (R/2R)^*
    //     if (!lideal_equals(I_red, O_4)) {
    //         J = lideal_create(A, order, algmul(A,lideal_generator(I),automorphism), gmul(lideal_norm(I),stoi(4)));
    //     }

    //     output(lideal_norm(I));
    //     output(lideal_norm(J));

    //     GEN alpha = lideal_isom(I, J); // I1*alpha = I2
    //     if (!alpha) { printf("J is not isomorphic to input!\n"); }

    // }




    GEN klpt_sol = klpt_special_smooth_given_2e_J(J,fm);

    if (!klpt_sol) { avma = ltop; fprintf(stderr, "klpt_special_smooth_small_2e_input did not find a solution\n"); return NULL; }
    else return gerepilecopy(ltop, klpt_sol);

}


// assume basis of lideal_algebra(I) is of the form 1, i, j, ji, with i^2 = -q, j^2 = -p
GEN klpt_general_power_given_J(GEN I, GEN l, GEN J, GEN delta) {
    pari_sp ltop = avma;
    GEN A = lideal_algebra(I);
    GEN order = lideal_order(I);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);
    long q = itos_or_0(algnorm(A,mkcol4s(0,1,0,0),0));

    assert(q != 0);

    GEN N;
    N = lideal_norm(J);


    int l_is_square_N = Fp_issquare(l,N);
    int l_is_square_NI = Fp_issquare(l,lideal_norm(I));

    // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    long L1_exp = floor(dbllog2r(itor(p_div_N,10)) / dbllog2r(itor(l,10)) );
    GEN L1 = powis(l, L1_exp);

    GEN coeff, gamma, beta1 = NULL, beta1_norm;

    int beta_norm_is_square_N, parity = -1;
    unsigned int safety_ctr = 0;
    do {
        coeff = NULL;
        while((!coeff) && safety_ctr < 16) {
            L1 = gmul(L1, l);
            safety_ctr++;
            // if (q == 1) coeff = norm_equation_special(p, gmul(L1,N), 0, false);
            if (q == 1) coeff = norm_equation_special(p, gmul(L1,N), true , true);
            else coeff = norm_equation_special_q(p,q,gmul(L1,N));
        }
        if (!coeff) break;

        gamma = gtrans(coeff); // norm N*L1

        assert(gcmp(algnorm(A,gamma,0), gmul(L1,N)) == 0);

        // compute beta1 in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        beta1 = klpt_solve_beta(A, gamma, J, N);

        // make sure that there exists an e such that l^e*norm(beta) is a quadratic residue mod N
        if (beta1) {
            beta1_norm = algnorm(A,beta1,0);
            beta_norm_is_square_N = Fp_issquare(beta1_norm,N);

            if (beta_norm_is_square_N) {
                if (l_is_square_N) {}
                else if (!l_is_square_N) parity = 0;
            }
            else if (!beta_norm_is_square_N) {
                if (l_is_square_N) { beta1 = NULL; } // try again...
                else if (!l_is_square_N) parity = 1;
            }
        }
    }
    // sometimes a gamma is found such that gamma*j is in J, in which case beta1 is not found... a few repetitions should avoid that
    while ((!beta1) && (safety_ctr < 16)); // max 16 because we do not want L1 to get too big
    if (!beta1) {
        avma = ltop;
        // beta1 not found
        return NULL;
    }


    if (gcmp(ggcd(beta1_norm, N),gen_1) != 0) {
        // Problem with beta1
        avma = ltop;
        return NULL;
    }

    assert(alglatcontains(A, lideal_lattice(J), algmul(A,gamma,beta1), NULL));

    // compute beta2 in span(j,j*i) such that gamma*beta*delta/N - zeta in I with GCD(zeta,N(I)) = 1
    GEN beta2 = mkcol4s(1,0,0,0);
    if (gcmp(lideal_norm(I),gen_1) != 0)  beta2 = klpt_solve_beta_special(A, gamma, delta, N, I);

    if (!beta2) {
        avma = ltop;
        // beta2 not found
        return NULL;
    }

    GEN beta2_norm = algnorm(A,beta2,0);

    int beta_norm_is_square_NI = Fp_issquare(beta2_norm,lideal_norm(I));

    if (gcmp(ggcd(beta2_norm, N),gen_1) != 0) {
        // Problem with beta2
        avma = ltop;
        return NULL;
    }
    // recover beta as a CRT combination of beta1 and beta2
    GEN C0 = Z_chinese(gel(beta1,3), gel(beta2,3), N, lideal_norm(I));
    GEN D0 = Z_chinese(gel(beta1,4), gel(beta2,4), N, lideal_norm(I));

    GEN beta = mkcol4(gen_0,gen_0,C0,D0);

    GEN NNI = gmul(N,lideal_norm(I));
    GEN bound_L2 = gmul(gmul(gmul(p,NNI),gmul(NNI,NNI)),stoi(80)); // (p*N^3)*MARGIN

    long L2_exp = floor(dbllog2r(itor(bound_L2,10)) / dbllog2r(itor(l,10)));

    GEN beta_norm = algnorm(A,beta,0);

    int fail = 0;
    if (beta_norm_is_square_NI) {
        if (l_is_square_NI) {}
        else if (!l_is_square_NI) {
            if (parity == 1) { fail = 1; } // no way
            else { parity = 0; }
        }
    }
    else if (!beta_norm_is_square_NI) {
        if (l_is_square_NI) { fail = 1; } // no way
        else if (!l_is_square_NI) {
            if (parity == 0) { fail = 1; } // no way
            else { parity = 1; }
        }
    }

    if (fail) {
        avma = ltop;
        return NULL;
    }

    if ((parity != -1) && ((L2_exp % 2) != parity)) L2_exp++;

    GEN L2 = powis(l, L2_exp);

    for (int j = 0; j < 10; j++) {
        GEN lambda1 = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
        GEN lambda2 = gen_1;
        if (gcmp(lideal_norm(I),gen_1) != 0)
            lambda2 = Fp_sqrt(Fp_div(L2,beta_norm,lideal_norm(I)), lideal_norm(I));
        GEN lambda = Z_chinese(lambda1, lambda2, N, lideal_norm(I));
        GEN mu = klpt_strong_approximation(A, order, p, NNI, beta, gamma, L2, lambda);
        if (mu) {
            GEN alpha = algmul(A,gamma,mu);
            assert(alglatcontains(A, lideal_lattice(J), alpha, NULL));
            GEN generatorJ = lideal_generator_coprime(J,l);
            assert(alglatcontains(A, lideal_lattice(J), generatorJ, NULL));
            GEN generator = algmul(A,generatorJ,alg_conj(A,alpha));
            assert(alglatcontains(A, lideal_lattice(J), alg_conj(A,generator), NULL));
            GEN norm = gmul(L1,L2);

            assert(gcmp(gmul(N,L1), algnorm(A,gamma,0)) == 0);
            assert(gcmp(L2, algnorm(A,mu,0)) == 0);

            // N(alpha) = N*norm
            GEN newideal = lideal_create(A, order, generator, norm);
            return gerepilecopy(ltop, newideal);
        }

        L2 = gmul(L2,l);
        if (parity != -1) L2 = gmul(L2,l);
    }

    return NULL;
}

GEN klpt_general_power_given_J_fixed_norm(GEN I, GEN l, GEN J, GEN delta) {
    pari_sp ltop = avma;
    // clock_t t;

    GEN A = lideal_algebra(I);
    GEN order = lideal_order(I);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);

    GEN N;
    N = lideal_norm(J);

    // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    long L1_exp = floor(dbllog2r(itor(p_div_N,10)) / dbllog2r(itor(l,10)) )+13;
    GEN L1 = powis(l, L1_exp);


    GEN coeff, gamma, beta1, beta1_norm,beta2,beta2_norm,beta,beta_norm,mu;

    unsigned int safety_ctr = 0;

    GEN NNI = gmul(N,lideal_norm(I));

    long L2_exp;
    GEN L2;
    do {
      safety_ctr ++;
      // coeff = norm_equation_special(p,gmul(L1,N), 0, true);
      coeff = norm_equation_special_max_order(p,gmul(L1,N),false,true);
      if (!coeff) break;
      gamma = gtrans(coeff); // norm N*L1


      GEN X1;
      alglatcontains(A,order,gamma,&X1);
      GEN n1 =content(X1);

      //overestimate the size, so that after removing the constant factor, the lenght is exacltly signing_length
      L2_exp = signing_length - L1_exp + 2 + 2*Z_lval(ggcd(n1,L1),2);
      L2 = powis(l, L2_exp);


      //try a first time with gamma

      beta = mkcol4(gen_1,gen_0,gen_0,gen_0);
      // compute beta1 in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
      beta1 = klpt_solve_beta(A, gamma, J, N);
      if (!beta1) {
          avma = ltop;
          // beta1 not found
          return NULL;
      }
      beta1_norm = algnorm(A,beta1,0);

      if (gcmp(ggcd(beta1_norm, N),gen_1) != 0) {
          // Problem with beta1
          avma = ltop;
          return NULL;
      }
      // compute beta2 in span(j,j*i) such that gamma*beta*delta/N - zeta in I with GCD(zeta,N(I)) = 1

      if (Fp_issquare( Fp_mul(L2,beta1_norm,N), N )) {
        beta2 = mkcol4s(1,0,0,0);
        if (gcmp(lideal_norm(I),gen_1) != 0) beta2 = klpt_solve_beta_special(A, gamma, delta, N, I);



        if (!beta2) {
            avma = ltop;
            return NULL;
            // beta2 not found
        }

        beta2_norm = algnorm(A,beta2,0);


        if (gcmp(ggcd(beta2_norm, N),gen_1) != 0) {
            // Problem with beta2
            avma = ltop;
            return NULL;
        }
        if (!Fp_issquare( Fp_mul(L2,algnorm(A,beta2,0),lideal_norm(I)), lideal_norm(I))  ) {
          beta = NULL;
        }
      }
      else {
        beta = NULL;
      }
      // try with the conjugate of gamma
      if (!beta) {
        gamma = alg_conj(A,gamma);
        beta = mkcol4(gen_1,gen_0,gen_0,gen_0);
        // compute beta1 in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        beta1 = klpt_solve_beta(A, gamma, J, N);
        if (!beta1) {
            avma = ltop;
            // beta1 not found
            return NULL;
        }
        beta1_norm = algnorm(A,beta1,0);

        if (gcmp(ggcd(beta1_norm, N),gen_1) != 0) {
            // Problem with beta1
            avma = ltop;
            return NULL;
        }
        // compute beta2 in span(j,j*i) such that gamma*beta*delta/N - zeta in I with GCD(zeta,N(I)) = 1


        if (Fp_issquare( Fp_mul(L2,beta1_norm,N), N )) {
          beta2 = mkcol4s(1,0,0,0);
          if (gcmp(lideal_norm(I),gen_1) != 0) beta2 = klpt_solve_beta_special(A, gamma, delta, N, I);



          if (!beta2) {
              avma = ltop;
              return NULL;
              // beta2 not found
          }

          beta2_norm = algnorm(A,beta2,0);


          if (gcmp(ggcd(beta2_norm, N),gen_1) != 0) {
              // Problem with beta2
              avma = ltop;
              return NULL;
          }
          if (!Fp_issquare( Fp_mul(L2,algnorm(A,beta2,0),lideal_norm(I)), lideal_norm(I))  ) {
            beta = NULL;
          }
        }
        else {
          beta = NULL;
        }
      }


      if (beta) {
        // recover beta as a CRT combination of beta1 and beta2
        GEN C0 = Z_chinese(gel(beta1,3), gel(beta2,3), N, lideal_norm(I));
        GEN D0 = Z_chinese(gel(beta1,4), gel(beta2,4), N, lideal_norm(I));
        beta = mkcol4(gen_0,gen_0,C0,D0);
        beta_norm = algnorm(A,beta,0);
        GEN lambda1 = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
        GEN lambda2 = gen_1;
        if (gcmp(lideal_norm(I),gen_1) != 0)
            lambda2 = Fp_sqrt(Fp_div(L2,beta_norm,lideal_norm(I)), lideal_norm(I));
        GEN lambda = Z_chinese(lambda1, lambda2, N, lideal_norm(I));
        mu = klpt_strong_approximation(A, order, p, NNI, beta,gamma, L2, lambda);
      }
    } while ( !beta && !mu && (safety_ctr <100));
    if (!beta || !mu) {
        avma = ltop;
        return NULL;
        // beta2 not found
    }

    if (mu) {
        GEN alpha = algmul(A,gamma,mu);

        GEN generatorJ = lideal_generator_coprime(J,l);
        GEN generator = algmul(A,generatorJ,alg_conj(A,alpha));
        GEN norm = gmul(L1,L2);
        GEN newideal = lideal_create(A, order, generator, norm);
        return gerepilecopy(ltop, newideal);
    }

    return NULL;
}


GEN klpt_general_power_small_J(GEN I, GEN K, GEN l, GEN list_previous_NJ, GEN* NJ) {
    pari_sp ltop = avma;

    // find prime ideal J equivalent to K
    // delta in K and K*conj(delta)/N(K) = J
    GEN delta, J;

    J = lideal_equiv_prime_except(K,&delta,list_previous_NJ);

    GEN klpt_sol = klpt_general_power_given_J(I, l, J, delta);
    if (NJ && klpt_sol) {
        *NJ = lideal_norm(J);
        gerepileall(ltop, 2, &klpt_sol, NJ);
    }
    else if (NJ) {
        *NJ = gerepilecopy(ltop, lideal_norm(J));
    }
    else if (klpt_sol) {
        klpt_sol = gerepilecopy(ltop, klpt_sol);
    }
    return klpt_sol;
}


GEN klpt_general_power_random_J(GEN I, GEN K, GEN l, GEN bound_coeff_J) {
    pari_sp ltop = avma;
    //clock_t t;
    // find prime ideal equivalent to K
    // delta in K and K*conj(delta)/N(K) = J
    GEN delta, J;

    //t = tic();
    J = lideal_equiv_prime_random(K,&delta,bound_coeff_J);
    //TOC(t, "prime\t");
    if (!J) { avma = ltop; return NULL; }

    return klpt_general_power_given_J_fixed_norm(I, l, J, delta);
}



GEN eichler_norm_equation_special_smooth_given_J(GEN J0, GEN L0, GEN fm , GEN l, GEN delta) {
  pari_sp ltop = avma;

  GEN A = lideal_algebra(J0);
  GEN order = lideal_order(J0);

  GEN A_j = mkcol4s(0,0,1,0);
  GEN p = algnorm(A,A_j,0);

  GEN N=lideal_norm(J0);
  GEN J = J0;


  GEN alpha1 = lideal_generator(L0);

  alpha1 = algmul(A,delta,alpha1);
  alpha1 = algmul(A,alpha1,alg_conj(A,delta));


  GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,delta,0)),gen_0,gen_0,gen_0);
  alpha1 = algmul(A,alpha1,overnorm);

  GEN I1 = lideal_create(lideal_algebra(J),lideal_right_order(J),alpha1,l);
  GEN order1 = lideal_right_order(I1);


  // int size_N = dbllog2r(itor(N,10));
  // printf("eichler norm eq with size %d \n",size_N);


  bool needs_non_integer_sol = true;
  // try 1 times then abandon
  for (int n = 1; n <= 1; ++n) {
      GEN L1;
      L1= gen_1;


      GEN gamma = mkcol4(gen_1,gen_0,gen_0,gen_0);

      GEN fm_rem=fm;
      // compute beta in span(j,j*i) such that gamma*beta is in the eichler order Z+J
      GEN beta = klpt_solve_beta_special(A, gamma, gamma, gen_1, J);
      if (!beta) continue;
      GEN beta_norm = algnorm(A,beta,0);

      if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
          // printf("Problem with beta!\n");
          continue;
      }
      // generate L2 such that L2/beta_norm is a square modulo N
      int beta_norm_is_square = Fp_issquare(beta_norm,N), L2_is_square;
      GEN bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(100)), fm_2, L2; // (p*N^3)*MARGIN
      int size = dbllog2r(itor(bound_L2,10));
      _unused(size);



      int safety_ctr = 0;
      GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
      GEN betap = NULL;
      GEN beta_final = NULL;
      int countr = 0;
      int nope_countr = 0;
      do {
          countr ++;

          //first try
          do {
              fm_2 = famat_smallest(fm_rem, bound_L2);
              L2 = famat_prod(fm_2);

              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;

          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 1));
          //first few tries
          do {
              fm_2 = famat_random(fm_rem, bound_L2,2);
              L2 = famat_prod(fm_2);

              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;

          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 4));
          //final tries
          do {
              fm_2 = famat_random(fm_rem, bound_L2,0);
              L2 = famat_prod(fm_2);

              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;

          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 10));

          if (safety_ctr >= 10) { continue; } // not enough entropy in fm_rem to ensure a solution
          if ( gcmp(bound_L2,L2) == 1  ) {
            // printf("required amount of torsion: %d \n",size);
             continue; } // L2 is too big

          list_L2 = vec_append(list_L2,L2);



          // strong approximation
          GEN lambda = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
          if (needs_non_integer_sol) {
            GEN L2_mul_4 = gmul(L2,gpowgs(gen_2,2));
            GEN lambda_mul_2 = gmul(lambda,gen_2);
            betap = klpt_strong_approximation_divided_by_2(A, p, N, beta, L2_mul_4, lambda_mul_2);

          }
          else {
            betap = klpt_strong_approximation(A, order,p, N, beta, gamma, L2, lambda);
          }



          if (betap) {
            // checks the additional property on eichler orders
            // only works when 2 is coprime with NL
              GEN betap_plus_one = algadd(A,betap,mkcol4(gen_1,gen_0,gen_0,gen_0));
              int check = alglatcontains(A,order1,betap_plus_one,NULL);



              // GEN beta_temp = betap;
              // GEN alpha2 = algmul(A,alpha1,alg_conj(A,beta_temp));
              // GEN I2 = lideal_create(lideal_algebra(J),lideal_right_order(J),alpha2,l);
              // check = lideal_equals(I1,I2);
              //

                if ( ! check ) {
                  beta_final=betap;
                  return gerepilecopy(ltop, beta_final);
                }
                else {
                  nope_countr++;
                }

              }





      } while (!beta_final && (safety_ctr < 10) && (countr < 20) && (nope_countr <10) );

  }
  avma = ltop;
  return NULL;
}



GEN eichler_norm_equation_special_smooth_given_J_extremal(p_extremal_maximal_order extremal, GEN J0, GEN L0, GEN fm , GEN l, GEN delta) {
  pari_sp ltop = avma;
  // clock_t t;

  // float accumulated_time_ms = 0.;
  // t=tic();
  GEN A = lideal_algebra(J0);
  GEN order = extremal.order; // == lideal_order(J0);

  GEN p = algnorm(A,extremal.j,0);

  GEN N=lideal_norm(J0);
  GEN J = J0;
  // assert( lideal_equals() )
  GEN alpha1 = lideal_generator(L0);

  alpha1 = algmul(A,delta,alpha1);
  alpha1 = algmul(A,alpha1,alg_conj(A,delta));

  GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,delta,0)),gen_0,gen_0,gen_0);
  alpha1 = algmul(A,alpha1,overnorm);
  assert(alglatcontains(A,lideal_right_order(J),alpha1,NULL));

  GEN I1 = lideal_create(A,lideal_right_order(J),alpha1,l);


  // int size_N = dbllog2r(itor(N,10));

  bool needs_non_integer_sol = true;
  // try 2 times then abandon
  for (int n = 1; n <= 2; ++n) {
      GEN L1;
      L1= gen_1;
      // GEN coeff = NULL;
      // unsigned long ctr = 1;


      GEN gamma = mkcol4(gen_1,gen_0,gen_0,gen_0);

      GEN fm_rem=fm;
      // compute beta in span(j,j*i) such that gamma*beta is in the eichler order Z+J
      GEN beta_package = klpt_solve_beta_extremal(A, extremal, gamma, gamma, gen_1, J);
      GEN beta = gel(beta_package,1);
      GEN C0 = gel(beta_package,2);
      GEN D0 = gel(beta_package,3);
      if (!beta) continue;
      GEN beta_norm = algnorm(A,beta,0);
      if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
          // printf("Problem with beta!\n");
          continue;
      }
      // generate L2 such that L2/beta_norm is a square modulo N
      int beta_norm_is_square = Fp_issquare(beta_norm,N), L2_is_square;


      //this could be very dangerous
      GEN bound_L2, fm_2, L2;
      if (two_tors_height > 45 || n >1 ) {
        bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(100)); // (p*N^3)*MARGIN
      }
      else {
        bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(1000000000000)); // (p*N^4)*MARGIN
      }


      int size = dbllog2r(itor(bound_L2,10));
      _unused(size);
      // printf("size: %d\n",size);

      int safety_ctr = 0;
      GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
      GEN betap = NULL;
      GEN beta_final = NULL;
      int countr = 0;
      int nope_countr = 0;
      do {
          countr ++;


          do {
              fm_2 = famat_random(fm_rem, bound_L2,0);
              L2 = famat_prod(fm_2);

              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;

          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 10));

          if (safety_ctr >= 10) { continue; } // not enough entropy in fm_rem to ensure a solution
          if ( gcmp(bound_L2,L2) == 1  ) {
            // printf("required amount of torsion: %d \n",size);
             continue; } // L2 is too big

          list_L2 = vec_append(list_L2,L2);
          // strong approximation
          GEN lambda = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
          // if (!lambda) { continue; } // we already made sure this doesn't happen
          // t1 = tic();
          if (needs_non_integer_sol) {

            // todo: this only deals with the 2-part of the discriminant... sufficient to get full order when q = 1, but not q > 1

            GEN L2_mul_4 = gmul(L2,gpowgs(gen_2,2));
            GEN lambda_mul_2 = gmul(lambda,gen_2);
            // betap = klpt_strong_approximation_divided_by_2(A, p, N, beta, L2_mul_4, lambda_mul_2);
            betap = klpt_strong_approximation_divided_by_2_extremal(A, extremal, p, N, beta, C0, D0, L2_mul_4, lambda_mul_2);
            if (!betap) {
            }
            // if (betap && !alglatcontains(A,lideal_order(J0),betap,NULL)) {
            //   printf("betap  is not in O0 \n");
            // }
            // assert();
          }
          else {
            betap = klpt_strong_approximation_extremal(A, extremal, order,p, N, beta, C0, D0, gamma, L2, lambda);
          }
          // output(L2);
          // betap=klpt_strong_approximation_with_eichler_constraint(A, order, p, N, beta,gamma, L2, lambda,J,l ,L0);
          // float acc=toc(t1);
          // printf("SA t: %f ms \n",acc);
          // accumulated_time_ms += acc;
          if (betap) {
              // output(betaps);
              GEN beta_temp = algmul(A,gamma,betap);


              // checks the additional fact on eichler orders
              // only works when 2 is coprime with NL



              GEN alpha2 = algmul(A,alpha1,alg_conj(A,beta_temp));
              GEN I2 = lideal_create(lideal_algebra(J),lideal_right_order(J),alpha2,l);

              int check = lideal_equals(I1,I2);
              if (check) {
              }

                if ( ! check ) {
                    beta_final=beta_temp;
                    // printf("time: %f ms,   accumulated_time: %f ms  \n", toc(t),accumulated_time_ms);
                    // printf("success\n");
                    return gerepilecopy(ltop, beta_final);
                }
                else {
                  // printf("nope \n");
                  nope_countr++;
                }

              }





      } while (!beta_final && (safety_ctr < 10) && (countr < 20) && (nope_countr <10) );

  }
  avma = ltop;
  return NULL;
}

GEN eichler_norm_equation_special_smooth_given_fixed_J(GEN J0, GEN L0, GEN fm , GEN l, GEN delta) {
  pari_sp ltop = avma;
  // clock_t t,t1;
  printf("a \n");
  // float accumulated_time_ms = 0.;
  // t=tic();
  GEN A = lideal_algebra(J0);
  GEN order = lideal_order(J0);

  GEN A_j = mkcol4s(0,0,1,0);
  GEN p = algnorm(A,A_j,0);

  GEN N=lideal_norm(J0);

  GEN J = lideal_create(A, order, lideal_generator(J0), N);
  // GEN p_div_N = gadd(truedivii(p, N), gen_1);

  GEN alpha1 = lideal_generator_coprime(L0,gen_1);

  alpha1 = algmul(A,delta,alpha1);
  alpha1 = algmul(A,alpha1,alg_conj(A,delta));


  GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,delta,0)),gen_0,gen_0,gen_0);
  alpha1 = algmul(A,alpha1,overnorm);

  GEN I1 = lideal_create(lideal_algebra(J),lideal_right_order(J),alpha1,l);

  bool needs_non_integer_sol = true;

  // output(lideal_lattice(I_two));
  // try 10 times then abandon
  for (int n = 1; n <= 10; ++n) {
    printf("in loop 1 \n");
      GEN L1;
      L1= gen_1;
      // GEN coeff = NULL;
      // unsigned long ctr = 1;


      GEN gamma = mkcol4(gen_1,gen_0,gen_0,gen_0);

      GEN fm_rem=fm;
      // compute beta in span(j,j*i) such that gamma*beta is in the eichler order Z+J
      GEN beta = klpt_solve_beta_special(A, gamma, gamma, gen_1, J);
      if (!beta) continue;
      GEN beta_norm = algnorm(A,beta,0);

      if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
          // printf("Problem with beta!\n");
          continue;
      }
      // generate L2 such that L2/beta_norm is a square modulo N
      int beta_norm_is_square = Fp_issquare(beta_norm,N), L2_is_square;
      GEN bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(50)), fm_2, L2; // (p*N^3)*MARGIN

      int safety_ctr = 0;
      GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
      GEN betap = NULL;
      GEN beta_final = NULL;
      int countr = 0;
      int nope_countr = 0;
      do {
          countr ++;


          do {
              fm_2 = famat_random(fm_rem, bound_L2,0);
              L2 = famat_prod(fm_2);
              L2_is_square = Fp_issquare(L2,N);
              safety_ctr++;
          } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 30));

          if (safety_ctr >= 30) { continue; } // not enough entropy in fm_rem to ensure a solution
          if ( gcmp(bound_L2,L2) == 1  ) { continue; } // L2 is too big

          list_L2 = vec_append(list_L2,L2);

          // strong approximation
          GEN lambda = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
          // if (!lambda) { continue; } // we already made sure this doesn't happen
          // t1 = tic();
          if (needs_non_integer_sol) {
            GEN L2_mul_4 = gmul(L2,gpowgs(gen_2,2));
            GEN lambda_mul_2 = gmul(lambda,gen_2);
            betap = klpt_strong_approximation_divided_by_2(A, p, N, beta, L2_mul_4, lambda_mul_2);

            if (!alglatcontains(A,lideal_order(J0),betap,NULL)) {
              printf("betap  is not in O0 \n");
            }
            // assert();
          }
          else {
            betap = klpt_strong_approximation(A, order, p, N, beta,gamma,L2, lambda);
          }

          if (betap) {
              // output(betaps);
              // output(betap);
              GEN beta_temp = algmul(A,gamma,betap);




              GEN alpha2 = algmul(A,alpha1,alg_conj(A,beta_temp));
              GEN I2 = lideal_create(lideal_algebra(J),lideal_right_order(J),alpha2,l);

              int check = lideal_equals(I1,I2);

                if ( ! check ) {
                  beta_final=beta_temp;
                  // printf("time: %f ms,   accumulated_time: %f ms  \n", toc(t),accumulated_time_ms);
                  return gerepilecopy(ltop, beta_final);
                }
                else {
                  // printf("nope \n");
                  // output(lideal_generator_coprime(I1,gen_1));
                  // output(gmod(gel(beta_temp,1),gen_2));
                  // output(gmod(gel(beta_temp,2),gen_2));
                  // output(gmod(gel(beta_temp,3),gen_2));
                  // output(gmod(gel(beta_temp,4),gen_2));
                  nope_countr++;
                }

              }

              // GEN generator = algmul(A,lideal_generator_coprime(J0,norm),alg_conj(A,alpha));
              // GEN newideal = lideal_create(A, order, generator, norm);
              // printf("SA acc\t [%f ms]\n",  (accumulated_time_ms));
          // }



      } while (!beta_final && (safety_ctr < 40) && (countr < 40) && (nope_countr <20) );

  }
  avma = ltop;
  return NULL;
}


GEN eichler_norm_equation_special_smooth_small_J_fixed(eichler_package *e,GEN I, GEN L, GEN fm , GEN l, GEN list_previous_NJ) {
    pari_sp ltop = avma;
    //clock_t t;

    // find prime ideal J equivalent to I
    // delta in I and I*conj(delta)/N(I) = J
    GEN delta, J;

    //t = tic();
    J = lideal_equiv_prime_except(I,&delta,list_previous_NJ);
    if (!J) { avma = ltop; return NULL; }
    //TOC(t, "prime\t");

    // printf("log NJ = %ld\n", (long) dbllog2r(itor(NJ,10)));
    GEN sol = eichler_norm_equation_special_smooth_given_J(J, L, fm , l, delta);
    if (sol) {
      e->I = I;
      e-> L = L;
      e->beta = sol;
      e->J = J;
      e->delta = delta;
      e->gamma = mkcol4(gen_1,gen_0,gen_0,gen_0);
      // e->gamma = NULL;
      // return lideal_norm(NJ);
    }
    return sol;
    // *NJ = gerepilecopy(ltop, lideal_norm(J));

    // if (NJ && sol) {
    //     *NJ = lideal_norm(J);
    //     gerepileall(ltop, 2, &sol, NJ);
    // }
    // else if (NJ) {
    //     *NJ = gerepilecopy(ltop, lideal_norm(J));
    // }
    // else if (sol) {
    //     sol = gerepilecopy(ltop, sol);
    // }
    // return sol;
}


GEN eichler_norm_equation_special_smooth_random_J(eichler_package *e,GEN I, GEN L, GEN fm , GEN l, GEN bound_coeff_J ) {
    pari_sp ltop = avma;
    // printf("try new \n");
    //clock_t t;

    // find prime ideal J equivalent to I
    // delta in I and I*conj(delta)/N(I) = J
    GEN delta, J;

    //t = tic();
    J = lideal_equiv_prime_random(I,&delta,bound_coeff_J);
    //TOC(t, "prime\t");
    if (!J) { avma = ltop; return NULL; }

    GEN sol =  eichler_norm_equation_special_smooth_given_J(J, L, fm ,l, delta);

    //GEN sol = eichler_norm_equation_special_smooth_given_J(J, L, fm ,l, delta);
    if (sol) {
      e->I = I;
      e-> L = L;
      e->beta = sol;
      e->J = J;
      e->delta = delta;
      // e->gamma = NULL;
    }
    return sol;
}


GEN eichler_norm_equation_special_smooth_random_J_extremal(eichler_package *e, p_extremal_maximal_order extremal, GEN I, GEN L, GEN fm , GEN l, GEN bound_coeff_J ) {
    pari_sp ltop = avma;
    //clock_t t;

    // find prime ideal J equivalent to I
    // delta in I and I*conj(delta)/N(I) = J
    GEN delta, J;

    //t = tic();
    J = lideal_equiv_prime_random(I,&delta,bound_coeff_J);
    // assert(lideal_equals(J,lideal_mul(I,gdiv(alg_conj(lideal_algebra(J),delta),lideal_norm(I)) )));
    //TOC(t, "prime\t");
    if (!J) { avma = ltop; return NULL; }
    // p_extremal_maximal_order extremal = {lideal_order(J), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), 1};
    GEN sol =  eichler_norm_equation_special_smooth_given_J_extremal(extremal, J, L, fm ,l, delta);

    //GEN sol = eichler_norm_equation_special_smooth_given_J(J, L, fm ,l, delta);
    if (sol) {
      e->I = I;
      e-> L = L;
      e->beta = sol;
      e->J = J;
      e->delta = delta;
      e->gamma = NULL;
        gerepileall(ltop,6,&e->I,&e->L,&e->beta,&e->J,&e->delta,&sol);
    }
    return sol;
}


GEN eichler_norm_equation_special_smooth_small_J(eichler_package *e,GEN I, GEN L, GEN fm , GEN l) {
    pari_sp ltop = avma;
    //clock_t t;

    // find prime ideal J equivalent to I
    // delta in I and I*conj(delta)/N(I) = J
    GEN p = gneg(algdisc(lideal_algebra(I)));
    long issq = issquareall(p, &p); assert(issq);_unused(issq);

    GEN q = stoi(3);
    for (int i = 0; i < 40; ++i) {
        if (kronecker(gneg(q),p) == -1)  {
            // printf("we are in the alternate case \n");
            // output(q);
            p_extremal_maximal_order EMO = get_p_extremal_maximal_order(itos(q));

            // EMO = (p_extremal_maximal_order){global_setup.O0, global_setup.i, global_setup.j, 1};
            assert(EMO.order);
            // output(q);
            // assert(alglatcontains(lideal_algebra(I),lideal_right_order(I),lideal_generator(L),NULL));
            GEN I_new = connecting_ideal(lideal_algebra(I), EMO.order, lideal_right_order(I));
            // assert(alglatcontains(lideal_algebra(I),lideal_right_order(I_new),lideal_generator(L),NULL));

            // GEN smallest_norm = gel(lideal_short(I_new, NULL, gen_1),2);


            float bound_coeff_J= 5;
            int ctr =0;
            GEN sol = NULL;
            while (!sol && ctr < 20) {
                GEN sol = eichler_norm_equation_special_smooth_random_J_extremal(e, EMO, I_new, L, fm , l, stoi(bound_coeff_J) );
                if (sol) {
                    e->order = EMO;
                    gerepileall(ltop,9,&e->I,&e->L,&e->beta,&e->J,&e->delta,&sol,&e->order.order,&e->order.i,&e->order.j);
                    return sol;
                }
                ctr++;
                bound_coeff_J++;
            }
        }
        q = nextprime(gadd(q,gen_1));
    }

    avma = ltop;
    return NULL;
}




//the eichler_package contains two equivalent ideals I and J, a value delta such that
void eichler_norm_equation_special_smooth(eichler_package *e,GEN I, GEN L, GEN fm,GEN l) {
  // GEN NJ = gen_0;
  GEN sol = NULL;
  float bound_coeff_J= 5;
  int ctr =0;
  while (!sol && ctr < 20) {
    // start looking for the smallest equivalent primes, then try random ones
    if ((ctr < 1 && false) ) {
      sol = eichler_norm_equation_special_smooth_small_J(e,I, L, fm ,l);
    }
    else {
        sol = eichler_norm_equation_special_smooth_random_J(e,I, L, fm ,l, stoi(floor(bound_coeff_J)));
        bound_coeff_J += 1.;
    }
    // list_NJ = vec_append(list_NJ,NJ);
    ctr++;
  }
  //try with the extremal_order
  ctr =0;
  while (!sol && ctr < 1) {
      printf("have to use the extremal order \n");
      sol = eichler_norm_equation_special_smooth_small_J(e,I, L, fm ,l);
      // sol = eichler_norm_equation_special_smooth_small_J_old(e,I, L, fm ,l, NULL);
      ctr++;
  }

  if (!sol) {
    e->delta = NULL;
    printf("no solution found in eichler norm equation !!!! \n");
  }
  // else {
  //   // printf("found it after %d attempts \n",ctr);
  // }
}

//this is for the cases where we have a bad curve that would have possibly made the computation fail
void eichler_norm_equation_special_smooth_small_curve(eichler_package *e,GEN I, GEN L, GEN fm,GEN l) {
   // NJ = gen_0;
  GEN sol = NULL;
    int ctr =0;
    while (!sol && ctr < 3) {
        sol = eichler_norm_equation_special_smooth_small_J(e,I, L, fm ,l);
        // sol = eichler_norm_equation_special_smooth_small_J_old(e,I, L, fm ,l, NULL);
        ctr++;
    }

}

//this is for the first step where we know we have a very small prime equivalent ideal
void eichler_norm_equation_special_smooth_small_curve_fixed(eichler_package *e,GEN I, GEN L, GEN fm,GEN l) {
   // NJ = gen_0;
  GEN sol = NULL;
    int ctr =0;
    while (!sol && ctr < 3) {
        // sol = eichler_norm_equation_special_smooth_small_J(e,I, L, fm ,l);
        sol = eichler_norm_equation_special_smooth_small_J_fixed(e,I, L, fm ,l, NULL);
        ctr++;
    }
  if (!sol) {
    printf(" eichler norm equation with fixed equivalent ideal did not work ! \n");
  }
  else {
    // output(lideal_norm(e->J));
    //
  }
}

//possibly not correct ? (works only half the time)
// finds coeff C,D such that (C + D beta) sends ker L on ker K.
GEN find_coeff(GEN K, GEN L, GEN beta, GEN delta, GEN order_J) {
  pari_sp ltop = avma;

  GEN A = lideal_algebra(K);

  assert(alglatcontains(A,order_J,beta,NULL ));
  // assert(gequal(lideal_norm(K),lideal_norm(L)));
  // order_J = O_R(R)
  // let I = I(O_L(J), O_L(K)), connecting ideal, so
  // O_R(I) = O_L(K) = O_L(L)
  // O_L(I) = O_L(J)
  // we should have I*conj(delta)/N(I) = J

  // GEN N = gmul(gen_2,lideal_norm(K));
  GEN N = lideal_norm(K);

  GEN alphaL = lideal_generator(L);
  GEN alphaK = lideal_generator(K);

  GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,delta,0)),gen_0,gen_0,gen_0);

  alphaL = algmul(A,delta,alphaL);
  alphaL = algmul(A,alphaL,alg_conj(A,delta));
  alphaL = algmul(A,alphaL,overnorm);
  assert(alglatcontains(A,order_J,alphaL,NULL ));

  alphaK = algmul(A,delta,alphaK);
  alphaK = algmul(A,alphaK,alg_conj(A,delta));
  alphaK = algmul(A,alphaK,overnorm);
  assert(alglatcontains(A,order_J,alphaK,NULL ));

  GEN new_K = lideal_create_safe(A,order_J,alphaK,lideal_norm(K));
  // GEN new_L = lideal_create_safe(A,order_J,alphaL,lideal_norm(L));

  GEN beta_bar=alg_conj(A,beta);
  GEN K_basis = gmul(lideal_basis(new_K),lideal_scalar(new_K));
  GEN gamma = algmul(A,alphaL,beta_bar);


  GEN coord_alphaL,coord_gamma;
  GEN coord_K1,coord_K2,coord_K3,coord_K4;
  alglatcontains(A,order_J,alphaL,&coord_alphaL);
  alglatcontains(A,order_J,gamma,&coord_gamma);
  alglatcontains(A,order_J,gel(K_basis,1),&coord_K1);
  alglatcontains(A,order_J,gel(K_basis,2),&coord_K2);
  alglatcontains(A,order_J,gel(K_basis,3),&coord_K3);
  alglatcontains(A,order_J,gel(K_basis,4),&coord_K4);

  GEN matsys = mkmat2(coord_alphaL,coord_gamma);
  matsys = gconcat(matsys,mkmat4(coord_K1,coord_K2,coord_K3,coord_K4));

  // GEN matsys = gmul(mkmat2(alphaL,gamma),Q_denom(lideal_scalar(new_L)));
  // // output(matsys);
  // matsys= gconcat(matsys,gmul(K_basis,Q_remove_denom(K_scalar,NULL)));

  GEN ker;
  // if (mpodd(N)) { // prime case
  //     ker = matkermod(matsys, N, NULL);
  // }
  // // flag = 1 because integral entries
  // else {// power of 2 case
  //     ker = matkermod(matsys, gmul(N,Q_denom(K_scalar)), NULL);
  //   }

    ker = matkermod(matsys, N, NULL);

  unsigned long i = lg(ker)-1;
  GEN C0,D0;
  int check = 0;
  // int check;
  do {
      C0 = gel(gel(ker, i), 1);
      D0 = gel(gel(ker, i), 2);
      if (!((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0))) {
        check = 1;
        // GEN C0_alg = mkcol4(C0,gen_0,gen_0,gen_0);
        // GEN D0_alg = mkcol4(D0,gen_0,gen_0,gen_0);
        // GEN theta = algadd(A,algmul(A,D0_alg,beta),C0_alg);
        // GEN theta_bar = alg_conj(A,theta);
        // GEN K_test = lideal_create(A,order_J,algmul(A,alphaL,theta_bar),lideal_norm(K));
        // check = lideal_equals(new_K,K_test);
      }
      i--;
      if (i < 1) break;
  } while ((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0) || !( check) );
  if ((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0) || !( check) ) {
      avma = ltop;
      return NULL;
  }
  C0 = gmod(C0,lideal_norm(K));
  D0 = gmod(D0,lideal_norm(K));
  #ifndef NDEBUG
  GEN C0_alg = mkcol4(C0,gen_0,gen_0,gen_0);
  GEN D0_alg = mkcol4(D0,gen_0,gen_0,gen_0);
  GEN theta = algadd(A,algmul(A,D0_alg,beta),C0_alg);
  GEN theta_bar = alg_conj(A,theta);
  GEN K_test = lideal_create_safe(A,order_J,algmul(A,alphaL,theta_bar),lideal_norm(K));
  check = lideal_equals(new_K,K_test);
  assert(check);
  #endif

  GEN result = mkcol2(C0,D0);
  return gerepilecopy(ltop,result);
}


// I and K are both ideals is O0, of coprime norm
// Finds an l-path between the right order of I and the right order of I \cap K
// More precisely, finds an O0-ideal M of norm a power of l such that O_R(I \cap K) = O_R(I \cap M)
// Also assume that I has prime norm

GEN klpt_general_power(GEN I, GEN K, GEN l) {
    GEN list_NJ = const_vec(0, gen_0), NJ = gen_0;
    GEN klpt_sol = NULL;
    float bound_coeff_J = 5;
    int ctr = 0;

    while (!klpt_sol) {
        // start looking for the smallest equivalent primes, then try random ones
        //for signature we want random
        if ((ctr < 3) && false) klpt_sol = klpt_general_power_small_J(I, K, l, list_NJ, &NJ);
        else {
            klpt_sol = klpt_general_power_random_J(I, K, l, stoi(floor(bound_coeff_J)));
            bound_coeff_J += 1.;
        }
        list_NJ = vec_append(list_NJ,NJ);
        ctr++;
    }
    return klpt_sol;
}
