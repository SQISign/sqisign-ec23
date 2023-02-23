#include <stdio.h>
#include "precomputed.h"

// Some useful field constants
fp minus_one_half;
fp minus_one_third;
// fp2_nonresidue()^((p-1)/2)
fp2 non_residue_p_minus_1_halves;

static void fp_code(const fp *x) {
  printf("{ 0x%lxULL, 0x%lxULL, 0x%lxULL, 0x%lxULL }",
	 x->x.c[0], x->x.c[1], x->x.c[2], x->x.c[3]);
}

int main(){
  uintbig pow;
  fp_add3(&minus_one_half, &fp_1, &fp_1);
  fp_neg1(&minus_one_half);
  fp_sub3(&minus_one_third, &minus_one_half, &fp_1);
  fp_inv(&minus_one_half);
  fp_inv(&minus_one_third);
  
  uintbig_div3_64(&pow, &p, 2); // (p-1)/2
  non_residue_p_minus_1_halves = fp2_non_residue();
  fp2_exp(&non_residue_p_minus_1_halves, &non_residue_p_minus_1_halves, &pow);
  assert(class_mod_4 == 3 || fp_iszero(&non_residue_p_minus_1_halves.im));
  
  printf("#include \"precomputed.h\"\n\n");
    
  printf("// Field constants\n");
  printf("fp minus_one_half =\n  ");
  fp_code(&minus_one_half);
  printf(";\nfp minus_one_third =\n  ");
  fp_code(&minus_one_third);
  printf(";\n// fp2_non_residue()^((p-1)/2)\n");
  printf("fp2 non_residue_p_minus_1_halves = {\n  ");
  fp_code(&non_residue_p_minus_1_halves.re);
  printf(",\n  ");
  fp_code(&non_residue_p_minus_1_halves.im);
  printf("\n};\n");
}
