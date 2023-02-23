#!/usr/bin/env sage

import itertools

def wordify(n, ws=64):
    ws = ws / 4
    h = ZZ(n).hex()
    return ", ".join('0x' + h[max(0,i-ws):i] for i in range(len(h), 0, -ws))

def fp_c_consts(p, smoothness=None, signing_length=None):
    p = ZZ(p)
    smoothness = smoothness or 2**12
    signing_length = signing_length or 1000
    pm = [(f,e) for f,e in factor(p-1) if 2 < f < smoothness]
    pp = [(f,e) for f,e in factor(p+1) if 2 < f < smoothness]
    twoval = (p - 2 + (p % 4)).valuation(2)
    twocof = (p - 2 + (p % 4)) >> twoval
    sign_steps = ceil(signing_length / twoval)
    return f"""
// p = { p }
//
// p-1 = { factor(p-1) }
// p+1 = { factor(p+1) }

#include "constants.h"

const long class_mod_4 = { p % 4 };
const long two_tors_height = { twoval };

const long security_level = { p.nbits() // 2 };
const long signing_length = { signing_length } ;
const long signing_length_two_tors_height_step = { sign_steps };
const long last_step_length = { signing_length - (sign_steps - 1) * twoval };

const char* p_str =
  "{p}";
const char* all_the_torsion_str =
  "{p**2 - 1}";

const uintbig p_plus_odd_cofactor =
  {{{ wordify((p+1) // prod(f**e for f,e in pp)) }}};
const uintbig p_minus_odd_cofactor =
  {{{ wordify((p-1) // prod(f**e for f,e in pm)) }}};
const uintbig p_even_cofactor =
  {{{ wordify(twocof) }}};

#define M_LEN { len(pm) }
const long p_minus_len = M_LEN;
const long p_minus_fact[M_LEN] =
  {{{ ", ".join(str(f) for f,_ in pm) }}};
const long p_minus_mult[M_LEN] =
  {{{ ", ".join(str(e) for _,e in pm) }}};

#define P_LEN { len(pp) }
const long p_plus_len = P_LEN;
const long p_plus_fact[P_LEN] =
  {{{ ", ".join(str(f) for f,_ in pp) }}};
const long p_plus_mult[P_LEN] =
  {{{ ", ".join(str(e) for _,e in pp) }}};

// TODO: Adjust the parameters below!
//
// the multiplicities to take to obtain log2(p) bits of torsion (for commitment)
const long p_minus_mult_com[M_LEN] =
  {{{ ", ".join(str(e) for _,e in pm) }}};
const long p_plus_mult_com[P_LEN] =
  {{{ ", ".join(str(e) for _,e in pp) }}};

// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
const long p_minus_mult_cha[M_LEN] =
  {{{ ", ".join(str(e) for _,e in pm) }}};
const long p_plus_mult_cha[P_LEN] =
  {{{ ", ".join(str(e) for _,e in pp) }}};
"""

def fp_asm_consts(p):
    return f"""
.global p
.global _p
p: _p:
    .quad { wordify(p) }

.inv_min_p_mod_r: /* -p^-1 mod 2^64 */
    .quad { wordify(Zmod(2**64)(-p)**-1) }

.global fp_1
.global _fp_1
fp_1: _fp_1: /* 2^256 mod p */
    .quad { wordify(2**256 % p) }

.r_squared_mod_p: /* (2^256)^2 mod p */
    .quad { wordify(2**512 % p) }

.p_minus_2:
    .quad { wordify(p - 2) }

.p_minus_1_halves:
    .quad { wordify((p - 1) // 2) }

/* Warning: this is specific to p = 3 mod 4 */
.p_plus_1_quarter:
    .quad { wordify((p + 1) // 4) }
"""

def fp2_consts(p):
    p = ZZ(p)
    if (p % 4 == 3):
        # Look for smallest non-quadratic residue in GF(p²)
        # (i+1 is a non-square iff p = 3 mod 8, so we don't test it)
        is_square = lambda a, b: GF(p)(a**2 + b**2).is_square()
        for a in itertools.count(2):
            for b in range(1, a):
                if not is_square(a, b):
                    re, im = ((a + b*GF(p)(-1).sqrt())**((p-1)//2)).polynomial()
                    return f"""
/* Arithmetic modulo X^2 + 1 */

// {a} + {b}i is a quadratic non-residue
fp2 fp2_non_residue() {{
  fp2 res;
  // TODO: fill implementation in
  return res;
}}
"""
    elif (p % 3 == 2):
        return f"""
/* Arithmetic modulo X^2 + 3 */

// √-3 is a quadratic non-residue
fp2 fp2_non_residue() {{
  return fp2_i;
}}
"""
    else:
        print ("// p = 1 mod 12 not supported")


####################

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate constants for prime')
    parser.add_argument('-s', '--smoothness', type=int, default=None,
                        help='smoothness bound')
    parser.add_argument('-l', '--sig-len', type=int, default=None,
                        help='length of isogeny signing chain')
    parser.add_argument('p', type=int, help='the prime')
    args = parser.parse_args()
    print('/********** Copy this to constants.c **********/')
    print(fp_c_consts(args.p, args.smoothness, args.sig_len))
    print('\n\n/********** Copy this at the top of fp.s **********/')
    print(fp_asm_consts(args.p))
    print('\n\n/********** Copy this at the top of fp2.c **********/')
    print(fp2_consts(args.p))
