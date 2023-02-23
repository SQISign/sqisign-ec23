#define _unused(x) ((void)(x))
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "quaternion_tools.h"
#include "toolbox.h"
#include "precomputed.h"

GEN connecting_ideal(GEN A, GEN O1, GEN O2)  {
    pari_sp ltop = avma;
    GEN inter = alglatinter(A, O1, O2, NULL);
    GEN norm = alglatindex(A,inter,O1);
    GEN gens = gmul(alglat_get_primbasis(O2),gmul(alglat_get_scalar(O2),norm));



    GEN I1 = lideal_create(A, O1, gel(gens, 1), NULL);
    GEN I2 = lideal_create(A, O1, gel(gens, 2), NULL);
    GEN I3 = lideal_create(A, O1, gel(gens, 3), NULL);
    GEN I4 = lideal_create(A, O1, gel(gens, 4), NULL);
    GEN I5 = lideal_create(A, O1, mkcol4(norm,gen_0,gen_0,gen_0), NULL);


    GEN I = lideal_add(I1,lideal_add(I2,lideal_add(I3,lideal_add(I4,I5))));

    // check correctness


    // GEN inter2 = alglatinter(A, lideal_right_order(I), O2, NULL);



    GEN index;
    assert(alglatsubset(A,lideal_right_order(I), O2, &index));
    assert(gequal1(index));
    _unused(index);



    return gerepilecopy(ltop,I);
}


// returns the quadratic form of the left ideal equal to the norm form divided by
// the scalar alglat_get_scalar(lideal_lattice(lideal))
GEN alglat_gram(GEN A, GEN lat) {
    pari_sp ltop = avma;
    GEN basis = alglat_get_primbasis(lat);
    GEN gram = gmul(gtrans(basis), gmul(alg_gram(A),basis));
    return gerepilecopy(ltop,gram);
}


// returns an LLL-reduced basis of lat divided by
// the scalar alglat_get_scalar(lat)
GEN alglat_lll(GEN A, GEN lat) {
    pari_sp ltop = avma;
    GEN basis = alglat_get_primbasis(lat);
    GEN kerim = lllgramkerim(alglat_gram(A, lat));
    GEN red = gmul(basis,gel(kerim,2));
    return gerepilecopy(ltop,red);
}

// projects u along w
GEN alg_orth_proj(GEN A, GEN u, GEN w) {
    pari_sp ltop = avma;
    GEN normw = gmul(gtrans(w), gmul(alg_gram(A), w));
    GEN scalar = gmul(gtrans(u), gmul(alg_gram(A), w));
    scalar = gdiv(scalar,normw);
    GEN vec = gmul(w, scalar);
    GEN proj = gsub(u, vec);
    return gerepilecopy(ltop,proj);
}

GEN alg_list_orth_proj(GEN A, GEN input_list, GEN w) {
    pari_sp ltop = avma;

    GEN list = gcopy(input_list);
    long n = lg(list) - 1;

    for (int i = 1; i <= n; ++i) {
        gel(list,i) = alg_orth_proj(A,gel(list,i),w);
    }
    return gerepilecopy(ltop,list);
}

// returns the sublattice of lat orthogonal to the space spanned by W, plus N*W (because pari needs full rank)
GEN alglat_orthogonal_part_extremal(GEN A, GEN lat, GEN W_input, GEN N) {
    pari_sp ltop = avma;

    GEN W = gcopy(W_input);

    long n = lg(W) - 1;
    GEN list_proj = alglat_get_primbasis(lat);

    for (int i = 1; i <= n; ++i) {
        GEN w = gel(W,i);
        list_proj = alg_list_orth_proj(A,list_proj,w);
        for (int j = i+1; j <= n; ++j) {
            gel(W,j) = alg_orth_proj(A,gel(W,j),w);
        }
    }

    list_proj = gconcat(list_proj, gmul(W_input,N));
    GEN lat_proj = alglathnf(A,gmul(list_proj, alglat_get_scalar(lat)),0);

    // output(algnorm(A, gel(list_proj,1),0));
    // output(algnorm(A, gel(list_proj,2),0));
    // output(algnorm(A, gel(list_proj,3),0));
    // output(algnorm(A, gel(list_proj,4),0));
    // output(algnorm(A, gel(list_proj,5),0));
    // output(algnorm(A, gel(list_proj,6),0));

    GEN res = alglatinter(A, lat_proj, lat, NULL);
    return gerepilecopy(ltop,res);
}

GEN qfsolve_rand(GEN qf, GEN some_solution) {
    pari_sp ltop = avma;
    GEN x0 = some_solution;

    // compute a solution if none is provided
    if (!x0)
        x0 = qfsolve(qf);

    if (typ(x0) == t_INT) {
        // no solution
        return NULL;
    }

    long n = lg(qf) - 1;

    // random vector
    GEN r = cgetg(n+1, t_COL);
    for (int i = 1; i <= n; ++i) {
        gel(r,i) = stoi(random_Fl(2001));
        gel(r,i) = gsub(gel(r,i),stoi(1000));
    }

    // random solution = r*gram*r*x0-2*r*gram*x0*r

    GEN a = gmul(gmul(gtrans(r), gmul(qf, r)), x0);
    GEN b = gmul(gmul(gtrans(r), gmul(qf, x0)), r);
    GEN x = gsub(gsub(a,b),b);

    GEN gcd = content(x);

    if (gequal0(gcd))
        return qfsolve_rand(qf, x0);

    x = gdiv(x, gcd);
    return gerepilecopy(ltop,x);
}


p_extremal_maximal_order get_p_extremal_maximal_order_from_scratch(long long disc_abs) {
    pari_sp ltop = avma;

    GEN B = global_setup.B;
    GEN O0 = global_setup.O0;

    long q = 1; // global_setup.q;

    GEN qf = mkmat4(mkcol4(gmul(global_setup.p,stoi(q)),gen_0,gen_0,gen_0),
                    mkcol4(gen_0,global_setup.p,gen_0,gen_0),
                    mkcol4(gen_0,gen_0,stoi(q),gen_0),
                    mkcol4(gen_0,gen_0,gen_0,gneg(stoi(disc_abs))));

    GEN sol = qfsolve_rand(qf, NULL);

    // if (typ(sol) == t_INT) {
    if (!sol) {
        // no solution
        return (p_extremal_maximal_order){NULL,NULL,NULL,0};
    }

    // let theta = (a*i+b*j+c*k)/u, where sol = (a,b,c,u)
    // theta has discriminant -disc_abs
    // find an order containing theta
    GEN theta_u = mkcol4(gen_0, gel(sol,3), gel(sol,2), gel(sol,1));
    GEN u = gabs(gel(sol,4),0);

    assert(gcmp(algnorm(B, theta_u,0),gmul(gsqr(u),stoi(disc_abs))) == 0);
    assert(alglatcontains(B, O0, theta_u, NULL));

    GEN I = lideal_create(B, O0, theta_u, u); // I = O0*theta*u + O0*u

    // small correction if N(I) < u
    theta_u = gdiv(gmul(theta_u, lideal_norm(I)),u);
    u = lideal_norm(I);
    GEN theta = gdiv(theta_u,u);


    GEN OR = alglatrighttransporter(B, lideal_lattice(I), lideal_lattice(I));


    GEN W = mkmat2(mkcol4(gen_1,gen_0,gen_0,gen_0), theta);
    GEN lat_orth = alglat_orthogonal_part_extremal(B, OR, W, global_setup.p); // part of OR orthoginal to 1 and theta

    GEN lll = alglat_lll(B, lat_orth); // find shortest vector

    GEN frob_candidate = gmul(gel(lll,1), alglat_get_scalar(lat_orth));
    // output(algnorm(B, theta,0));
    // is it the Frobenius?
    if(gequal(algnorm(B, frob_candidate,0),global_setup.p)) {
        assert(alglatcontains(B, OR, mkcol4s(1,0,0,0), NULL));
        assert(alglatcontains(B, OR, theta, NULL));
        assert(alglatcontains(B, OR, frob_candidate, NULL));
        assert(alglatcontains(B, OR, algmul(B,theta,frob_candidate), NULL));

        assert(gequal(algnorm(B, theta, 0), stoi(disc_abs)));

        // orthogonal
        assert(gequal0(gmul(gtrans(theta), gmul(alg_gram(B), frob_candidate))));

        p_extremal_maximal_order extremal;
        extremal.order = OR;
        extremal.i = theta;
        extremal.j = frob_candidate;
        extremal.q = disc_abs;

        gerepileall(ltop, 3, &extremal.order, &extremal.i, &extremal.j);
        return extremal;
    }
    else {
        // no Frobenius, rerandomizing...
        // can only happen if #Cl(-d) ≠ 1
        return get_p_extremal_maximal_order_from_scratch(disc_abs);
    }

}


p_extremal_maximal_order get_p_extremal_maximal_order(long long disc_abs) {
    if (disc_abs < EXTREMAL_ORDERS_N) {
        if (global_setup.orders[disc_abs]) {
            p_extremal_maximal_order extremal;
            extremal.order = gcopy(global_setup.orders[disc_abs]->order);
            extremal.i = gcopy(global_setup.orders[disc_abs]->i);
            extremal.j = gcopy(global_setup.orders[disc_abs]->j);
            extremal.q = global_setup.orders[disc_abs]->q;
            return extremal;
        }
    }

    return get_p_extremal_maximal_order_from_scratch(disc_abs);
}
