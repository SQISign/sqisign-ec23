
#ifndef QUATERNION_TOOLS_H
#define QUATERNION_TOOLS_H

#include <pari/pari.h>
#include <stdbool.h>

typedef struct p_extremal_maximal_order {
	GEN order;
	GEN i; // the element of small discriminant
	GEN j; // the element of norm p orthogonal to i
	long long q; // the discriminant of i
} p_extremal_maximal_order;

GEN connecting_ideal(GEN A, GEN O1, GEN O2);

GEN alglat_gram(GEN A, GEN lat);

GEN alglat_lll(GEN A, GEN lat);

GEN alg_orth_proj(GEN A, GEN u, GEN w);

GEN alg_list_orth_proj(GEN A, GEN input_list, GEN w);

GEN alglat_orthogonal_part_extremal(GEN A, GEN lat, GEN W_input, GEN N);

GEN qfsolve_rand(GEN qf, GEN some_solution);

p_extremal_maximal_order get_p_extremal_maximal_order_from_scratch(long long disc_abs);
p_extremal_maximal_order get_p_extremal_maximal_order(long long disc_abs);


#endif
