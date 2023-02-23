#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include "mont.h"
#include "isogenies.h"

static inline void print_big(const uintbig *x) {
  int i = sizeof(uintbig) / 8 - 1;
  printf("0x%lx", x->c[i]);
  for (i--; i >= 0; i--)
    printf("%016lx", x->c[i]);
  printf("\n");
}


static inline void print_fp2(const fp2 *x) {
  printf("re ");
  print_big(&x->re.x);
  printf("im ");
  print_big(&x->im.x);
}

static inline void print_proj(const proj *P) {
  printf("X ");
  print_fp2(&P->x);
  printf("Z ");
  print_fp2(&P->z);
}

static inline unsigned long fp2_hash(fp2 x) {
  return (x.re.x.c[0]+3*x.re.x.c[1]+5*x.re.x.c[2]+7*x.re.x.c[3]
	  +11*x.im.x.c[0]+13*x.im.x.c[1]+17*x.im.x.c[2]+23*x.im.x.c[3]) % 100003;
}

static inline fp2 fp2_ratio(const fp2 *x, const fp2 *y) {
  fp2 tmp;
  tmp = *y;
  assert(!fp2_iszero(&tmp));
  fp2_inv(&tmp);
  fp2_mul2(&tmp, x);
  return tmp;
}

static inline void print_proj_hash(const proj *P){
  fp2 hash;
  fp2_add3(&hash,&P->x,&P->z);
  fp2_mul3(&hash,&hash,&P->x);
  printf("%ld", fp2_hash(hash));
}

static inline void proj2_print(proj2 x) {
  if (fp2_iszero(&x.z)) { printf("(infty)"); }
  else { printf("(%lu,%lu)", fp2_hash(fp2_ratio(&x.x,&x.z)), fp2_hash(fp2_ratio(&x.y,&x.z))); }
}

static inline void print_deg(isog_degree deg) {
  for (int i = 0; i < 30; i++)
    printf("%d,", degree_get(deg, i));
  printf("\n");
}
