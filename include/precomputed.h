#ifndef PRECOMPUTED_H
#define PRECOMPUTED_H

#include <pari/pari.h>
#include <stdbool.h>
#include "constants.h"
#include "mont.h"
#include "tedwards.h"
#include "quaternion_tools.h"

// Some useful field constants
extern fp minus_one_half;
extern fp minus_one_third;
// fp2_nonresidue()^((p-1)/2)
extern fp2 non_residue_p_minus_1_halves;

// each basis entry is a triple of the form P,Q,P+Q
extern proj torsion_basis[][3];
extern proj torsion_basis_sum[3];
extern point torsion_basis_ted_sum[3];
extern proj torsion_basis_twist[][3];
extern proj torsion_basis_twist_sum[3];
extern point torsion_basis_twist_ted_sum[3];
extern proj torsion_basis_two[3];

#define EXTREMAL_ORDERS_N 16

struct precomp_struct {
    // quaternion data

    GEN p; // the prime
    long q; // such that B = (-q,-p)
    GEN B; // the quaternion algebra
    GEN qf; // the quadratic form defined by the reduced norm with respect to the standard basis
    GEN O0; // the cannonical maximal order
    GEN one;
    GEN i;
    GEN j;
    GEN ji;
    // GEN torsion_fm; // factorisation matrix of the available torsion

    GEN O0_b1; // 1
    GEN O0_b2;
    GEN O0_b3;
    GEN O0_b4;

    GEN O0_to_standard;
    GEN standard_to_O0;

    GEN smooth_famat_for_klpt;

    // elliptic curve data

    proj E0;

    GEN *action_2, *action_3, *action_4;
    GEN *action_twist_2, *action_twist_3, *action_twist_4;
    GEN action_two_2, action_two_3, action_two_4;

    GEN gen_p_plus_fact, gen_p_minus_fact; // factorisation of p+1 and p-1
    GEN gen_p_plus_primary, gen_p_minus_primary; // primary decomposition (list of prime powers)
    GEN gen_odd_torsion;

    p_extremal_maximal_order* orders[EXTREMAL_ORDERS_N];

};


extern struct precomp_struct global_setup;

void init_precomputations();
void init_precomputations_generated();

#endif
