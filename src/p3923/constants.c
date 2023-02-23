// p = 23759399264157352358673788613307970528646815114090876784643387662192449945599
//
// p-1 = 2 * 3^65 * 13 * 17 * 43 * 79 * 157 * 239 * 271 * 283 * 307 * 563 * 599 * 607 * 619 * 743 * 827 * 941 * 2357 * 10069
// p+1 = 2^65 * 5^2 * 7 * 11 * 19 * 29^2 * 37^2 * 47 * 197 * 263 * 281 * 461 * 521 * 3923 * 62731 * 96362257 * 3924006112952623

#include "constants.h"

const long class_mod_4 = 3;
const long two_tors_height = 65;

const long security_level = 127;
const long signing_length = 989;
const long signing_length_two_tors_height_step = 16;
const long last_step_length = 14;
const long len_tail = 14;
const bool need_even_commit = false;
const bool need_even_chall = true;

const char* p_str =
  "23759399264157352358673788613307970528646815114090876784643387662192449945599";
const char* all_the_torsion_str =
  "564509053393640936725750114155022155785978769656527351466904534115522598221337420839321402252432731293214917275172426130707490232370972580637058059468800";

const uintbig p_plus_odd_cofactor =
  {0x0000000000000000, 0x9ac144661f8dc6aa, 0x9949d317};
const uintbig p_minus_odd_cofactor =
  {0x4eaa};
const uintbig p_even_cofactor =
  {0xa9ac02fd374f2459, 0x5b988d31b19f81ed, 0x1a43abf56fae4a98};

#define M_LEN 18
const long p_minus_len = M_LEN;
const long p_minus_fact[M_LEN] =
  {3, 13, 17, 43, 79, 157, 239, 271, 283, 307, 563, 599, 607, 619, 743, 827, 941, 2357};
const long p_minus_mult[M_LEN] =
  {65, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

#define P_LEN 13
const long p_plus_len = P_LEN;
const long p_plus_fact[P_LEN] =
  {5, 7, 11, 19, 29, 37, 47, 197, 263, 281, 461, 521, 3923};
const long p_plus_mult[P_LEN] =
  {2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1};

// TODO: Adjust the parameters below!
//
// the multiplicities to take to obtain log2(p) bits of torsion when added with 2^two_tors_height (for commitment)
const long p_minus_mult_com[M_LEN] =
  {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const long p_plus_mult_com[P_LEN] =
  {2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 0};

// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
const long p_minus_mult_cha[M_LEN] =
  {65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const long p_plus_mult_cha[P_LEN] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
