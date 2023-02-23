// p = 73743043621499797449074820543863456997944695372324032511999999999999999999999
//
// p-1 = 2 * 3^53 * 43 * 103^2 * 109 * 199 * 227 * 419 * 491 * 569 * 631 * 677 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 4283
// p+1 = 2^33 * 5^21 * 7^2 * 11 * 31 * 83 * 107 * 137 * 751 * 827 * 3691 * 4019 * 6983 * 517434778561 * 26602537156291

#include "constants.h"

const long class_mod_4 = 3;
const long two_tors_height = 33;
const bool need_even_commit = false;
const bool need_even_chall = false;
const long security_level = 128;

const long signing_length=1000 ;

const long signing_length_two_tors_height_step = 31;
const long last_step_length = 10;

//should be two_tors_height + delta
//but is a bit tricky right now
const long len_tail = 0;

const char* p_str =
  "73743043621499797449074820543863456997944695372324032511999999999999999999999";
const char* all_the_torsion_str =
  "5438036482562421961818827832407841404369503851311267868831048027176534732013273363338731215746134716355745303230004110609255351934976000000000000000000000";

const uintbig p_plus_odd_cofactor = { 0x68cd740600000000, 0x0016c5bcbd22f015, 0, 0 };
const uintbig p_minus_odd_cofactor = { 2, 0, 0, 0 };
const uintbig p_even_cofactor = { 0xa52ca964a8652149, 0x1bb9479de8d8027c,
				  0xdb3c54c8592e3b52, 0x51848ab2 };

#define M_LEN 19
const long p_minus_len = M_LEN;
const long p_minus_fact[M_LEN] =
  {  3, 43, 103, 109, 199, 227, 419, 491, 569, 631, 677, 857, 859, 883,
    1019, 1171, 1879, 2713, 4283 };
const long p_minus_mult[M_LEN] =
  { 53,  1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
       1,    1,    1,    1,    1 };




#define P_LEN 12
const long p_plus_len = P_LEN;
const long p_plus_fact[P_LEN] =
  {  5, 7, 11, 31, 83, 107, 137, 751, 827, 3691, 4019, 6983 };
const long p_plus_mult[P_LEN] =
  { 21, 2,  1,  1,  1,   1,   1,   1,   1,    1,    1,    1 };


// the multiplicities to take to obtain log2(p) bits of torsion (for commitment)
const long p_minus_mult_com[M_LEN] =
  { 0,  1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
       1,    1,    1,    1,    1 };
const long p_plus_mult_com[P_LEN] =
  { 0, 2,  1,  1,  1,   1,   1,   1,   1,    1,    1,    1 };


// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
const long p_minus_mult_cha[M_LEN] =
  { 53,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,    0,    0,    0,    0 };
const long p_plus_mult_cha[P_LEN] =
  { 21, 0,  0,  0,  0,   0,   0,   0,   0,    0,    0,    0 };
