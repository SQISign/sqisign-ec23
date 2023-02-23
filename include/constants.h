#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <assert.h>
#include "uintbig.h"

// The class of p mod 4, and its consequence for the order of curves
extern const long class_mod_4;
#define curve_order_is_p_plus_one (class_mod_4 == 3)

extern const long two_tors_height;
// The cofactor of 2^two_tors_height in p±1
extern const uintbig p_even_cofactor;
extern const long security_level;

//says if we need to use the two-torsion in the commitment
extern const bool need_even_commit;
//same as above but for the challenge
extern const bool need_even_chall;


//says if we can use the trick that do a bit of the translation from the challenge curve
extern const long len_tail;

// the signing isogeny has degree 2^signing_length
extern const long signing_length;
// we have the equality signin_length = two_tors_height * (signing_length_two_tors_height_step -1 ) + last_step_length
extern const long signing_length_two_tors_height_step;
extern const long last_step_length;

extern const char* p_str;
extern const char* all_the_torsion_str;

// The useful odd factors in p-1
extern const long p_minus_len;
extern const long p_minus_fact[];
extern const long p_minus_mult[];
// The cofactor of the useful odd torsion in p-1
extern const uintbig p_minus_odd_cofactor;

// The useful odd factors in p+1
extern const long p_plus_len;
extern const long p_plus_fact[];
extern const long p_plus_mult[];
// The cofactor of the useful odd torsion in p+1
extern const uintbig p_plus_odd_cofactor;

// Same as above, but along the curve/twist dichotomy
#define on_curve_len (curve_order_is_p_plus_one ? p_plus_len : p_minus_len)
#define on_curve_fact (curve_order_is_p_plus_one ? p_plus_fact : p_minus_fact)
#define on_curve_mult (curve_order_is_p_plus_one ? p_plus_mult : p_minus_mult)
#define on_curve_odd_cofactor (curve_order_is_p_plus_one ? p_plus_odd_cofactor : p_minus_odd_cofactor)
#define on_twist_len (!curve_order_is_p_plus_one ? p_plus_len : p_minus_len)
#define on_twist_fact (!curve_order_is_p_plus_one ? p_plus_fact : p_minus_fact)
#define on_twist_mult (!curve_order_is_p_plus_one ? p_plus_mult : p_minus_mult)
#define on_twist_odd_cofactor (!curve_order_is_p_plus_one ? p_plus_odd_cofactor : p_minus_odd_cofactor)

// the multiplicities to take to obtain log2(p) bits of torsion (for commitment)
extern const long p_minus_mult_com[];
extern const long p_plus_mult_com[];

// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
extern const long p_minus_mult_cha[];
extern const long p_plus_mult_cha[];

// inverse mapping of p_plus_fact and p_minus_fact
// Warning: unsafe if ell is not in the factors!
static inline long ell_to_index(long ell, bool *twist) {
  *twist = false;
  for (const long *f = on_curve_fact; *f <= ell; f++)
    if (*f == ell)
      return f - on_curve_fact;
  *twist = true;
  for (const long *f = on_twist_fact; *f <= ell; f++)
    if (*f == ell)
      return f - on_twist_fact;
  assert(0);
  return(0);
}
static inline long ell_to_e(long ell) {
  bool twist;
  int index = ell_to_index(ell, &twist);
  return (curve_order_is_p_plus_one != twist ? p_plus_mult : p_minus_mult)[index];
}

#endif
