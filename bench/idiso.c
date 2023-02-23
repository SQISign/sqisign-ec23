#define _XOPEN_SOURCE

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <pari/pari.h>

#include "precomputed.h"
#include "sqisign.h"
#include "ideal.h"
#include "idiso.h"

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

int main(int argc, char **argv) {
  int keys = 10, samples = 10, seed = 1;

  int opt;
  while ((opt = getopt(argc, argv, "k:s:r:h")) != -1) {
    switch (opt) {
    case 'k':
      keys = atoi(optarg);
      break;
    case 's':
      samples = atoi(optarg);
      break;
    case 'r':
      seed = atoi(optarg);
      break;
    default:
      fprintf(stderr,
	      "Usage: %s [-k keys] [-s samples] [-r seed]\n",
	      argv[0]);
      exit(-1);
    }
  }

  pari_init(800000000, 1<<18);
  init_precomputations();

  setrand(mkintn(1, seed));
  srand48(seed);

  printf("### Sign\n");
  printf("# key\tcycles\t\tms\t\tlength\n");
  for (int kk = 0; kk < keys; kk++) {
    public_key pk;
    secret_key sk;
    keygen(&pk, &sk);
    GEN A = global_setup.B;

    // precomute last_I
    GEN beta = lideal_isom(sk.I_large_prime, sk.I_two); // I*beta = I_two;
    GEN last_I = lideal_create(A,lideal_right_order(sk.I_two),alg_conj(A,lideal_generator_coprime(sk.I_two,gen_1)),gpowgs(gen_2,two_tors_height));
    GEN alpha_last =lideal_generator_coprime(last_I,gen_1);
    alpha_last= algmul(A,alpha_last,alg_conj(A,beta));
    alpha_last= algmul(A,beta,alpha_last);
    alpha_last = algmul(A,alg_scalar(A,gdiv(gen_1,algnorm(A,beta,0))),alpha_last);

    last_I = lideal_create(A,lideal_right_order(sk.I_large_prime),alpha_last,gpowgs(gen_2,two_tors_height));
    for (int ind = 0; ind < samples; ind++) {
      // signature Sigma;
      two_walk_long sigma;
      init_trivial_two_walk_long(&sigma);
      uintbig zip[50];
      bool test;
      proj endpoint;
      GEN L_cut;
      GEN n;
      do{
          L_cut = lideal_random_2e(global_setup.B, global_setup.O0, 1000);
          assert(L_cut);
          alg_primitive(&n, A, lideal_order(sk.I_two), algmul(A, lideal_generator(L_cut), lideal_generator(sk.I_two)));
          // backtracking?
      } while(gcmp(n,gen_1) != 0);
      L_cut  =lideal_inter(L_cut,sk.I_large_prime);
      assert(alglatcontains(A,lideal_order(sk.I_large_prime),lideal_generator_coprime(L_cut,gen_1),NULL));
      clock_t t = -clock();
      uint64_t c = -rdtsc();
      ideal_to_isogeny_two_new(&sigma,sk.I_two,sk.I_large_prime,L_cut,last_I,&sk.phi_two,zip,&test,&endpoint);
      c += rdtsc();
      t += clock();
      free_two_walk_long(&sigma);

      printf("%d\t%" PRIu64 "\t%.3lf\t%ld\n", kk, c, 1000. * t / CLOCKS_PER_SEC, signing_length);
    }
  }

  return 0;
}
