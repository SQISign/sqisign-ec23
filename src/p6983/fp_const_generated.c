#include "precomputed.h"

// Field constants
fp minus_one_half =
  { 0x50ca4291ffffffffULL, 0xd1b004f94a5952c9ULL, 0xb25c76a437728f3bULL, 0x23091565b678a990ULL };
fp minus_one_third =
  { 0x1663036d55555554ULL, 0x6ceab14c6321c3b7ULL, 0x987b48daf498befaULL, 0x840c1c879df6376bULL };
// fp2_non_residue()^((p-1)/2)
fp2 non_residue_p_minus_1_halves = {
  { 0xbbda9829fdebd171ULL, 0xa9581a12615f83e2ULL, 0x9619f594cf016655ULL, 0x28733a4fcf36487aULL },
  { 0xc22d0587fce1ba29ULL, 0x66dc2998373bef38ULL, 0x3a552bb1523b611eULL, 0x8e31622a920dc180ULL }
};
