#define _XOPEN_SOURCE
#define _unused(x) ((void)(x))

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "tedwards.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "sqisign.h"
#include "mont.h"
#include "quaternion_tools.h"



    // argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(800000000, 1<<18);
    init_precomputations();

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    GEN q = stoi(3);
    for (int i = 0; i < 100; ++i) {
        p_extremal_maximal_order EMO = get_p_extremal_maximal_order(itos(q));
        // output(q);
        // printf("kro %ld\n", kronecker(q,global_setup.p));
        if (EMO.order) {
            assert(kronecker(gneg(q),global_setup.p) == -1);
            GEN I = connecting_ideal(global_setup.B, global_setup.O0, EMO.order);
            assert(alglatcontains(global_setup.B,global_setup.O0,lideal_generator(I),NULL));
            assert(alglatcontains(global_setup.B,EMO.order,lideal_generator(I),NULL));
            _unused(I);
        }
        else {
            assert(kronecker(gneg(q),global_setup.p) == 1);
        }
        q = nextprime(gadd(q,gen_1));
    }

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);

    return 0;
}
