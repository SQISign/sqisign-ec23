
#ifndef KLPT_H
#define KLPT_H

#include <pari/pari.h>
#include "quaternion_tools.h"

//the eichler order is ZZ + J
// where I and L are as in the inputs of eichler_norm_equation_special
// beta is the endomorphism of smooth norm obtained through eichler norm norm_equation_special
// delta is in I and I*conj(delta)/N(I) = J
// gamma is an additional information
typedef struct eichler_package {
  p_extremal_maximal_order order;
  GEN I;
  GEN J;
  GEN L;
  GEN beta;
  GEN delta;
  GEN gamma;
} eichler_package;

// runs KLPT for the left ideal I in the special order of the quaternion algebra A
// the result is an ideal equivalent to I of norm dividing the integer whose factorisation matrix is fm
// Assumes the basis of A is 1, i, j, j*i, where i^2 = -1 and j^2 = -p
GEN klpt_special_smooth(GEN I, GEN fm);
// same as above, when I is of norm a small power of two (in which case one cannot find an equivalent prime ideal of small norm)
GEN klpt_special_smooth_small_2e_input(GEN I, GEN fm);

// runs norm equation solving in the eichler order corresonding to an ideal J equivalent to the input ideal I
// the result has norm dividing the integer whose factorisation matrix is fm
// we have an additional constraint depending on an other ideal L (with OL(L) isomorphic to OR(I)):
// we want the push forward of L by the result endomorphism to be independant of L (meaning that their kernel form a basis of the n(L) torsion)
// we assume that K has norm a power of l
// the output is contained in OR(J)
void eichler_norm_equation_special_smooth(eichler_package *e,GEN I, GEN L, GEN fm, GEN l);

//same as above but for a curve close to E0 in the bad situation where we need another extremal order
void eichler_norm_equation_special_smooth_small_curve(eichler_package *e,GEN I, GEN L, GEN fm,GEN l);
//same as above but in the case where the close curve is good
void eichler_norm_equation_special_smooth_small_curve_fixed(eichler_package *e,GEN I, GEN L, GEN fm,GEN l);




GEN klpt_general_power(GEN I, GEN K, GEN l);

//Compute a vector of two coefficients C,D such that the push forward of L by C +D beta is equal to K
GEN find_coeff(GEN K, GEN L, GEN beta, GEN delta, GEN order_J);

//ONLY WORKS when the norm of K and L is a power of 2
GEN find_coeff2(GEN K, GEN L, GEN beta, GEN delta, GEN order_J);


#endif
