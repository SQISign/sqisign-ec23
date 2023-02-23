#include "uintbig.h"

#define N_WORDS 4

void gentobig(uintbig *res, GEN a) {
    // assert(gsigne(a) >= 0);
    pari_sp ltop = avma;
    GEN b;
    res->c[0] = umodi2n(a,32);
    b = shifti(a, -32);
    res->c[0] += (umodi2n(b,32) << 32);
    b = shifti(b, -32);
    for (int i = 1; i < N_WORDS; i++) {
      res->c[i] = umodi2n(b,32);
      b = shifti(b, -32);
      res->c[i] += (umodi2n(b,32) << 32);
      b = shifti(b, -32);
    }
    //assert(isexactzero(b));
    avma = ltop;
}

GEN bigtogen(const uintbig *a) {
  pari_sp ltop = avma;
  GEN gena = gen_0;
  for (int i = N_WORDS - 1; i >= 0; i--) {
    GEN tmp = uu32toi(a->c[i] >> 32, a->c[i] % ((uint64_t)1 << 32));
    gena = gadd(mpshift(gena, 64), tmp);
  }
  return gerepilecopy(ltop, gena);
}
