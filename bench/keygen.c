#define _XOPEN_SOURCE

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h> 
#include <pari/pari.h>

#include "precomputed.h"
#include "sqisign.h"

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

int main(int argc, char **argv) {
  int samples = 100, seed = 1;
  uint64_t c_total = 0;
  float t_total = 0;
  
  int opt;
  while ((opt = getopt(argc, argv, "s:r:h")) != -1) {
    switch (opt) {
    case 's':
      samples = atoi(optarg);
      break;
    case 'r':
      seed = atoi(optarg);
      break;
    default:
      fprintf(stderr,
	      "Usage: %s [-s samples] [-r seed]\n",
	      argv[0]);
      exit(-1);
    }
  }
  
  pari_init(800000000, 1<<18);
  init_precomputations();
  
  setrand(mkintn(1, seed));
  srand48(seed);

  printf("### Keygen\n");
  printf("# cycles\tms\n");
  for (int i = 0; i < samples; i++) {
    public_key pk;
    secret_key sk;

    clock_t t = -clock();
    uint64_t c = -rdtsc();
    keygen(&pk, &sk);
    c += rdtsc();
    t += clock();

    printf("%" PRIu64 "\t%.3lf\n", c, 1000. * t / CLOCKS_PER_SEC);
    c_total += c;
    t_total += (float)(1000. * t / CLOCKS_PER_SEC);
  }
  printf("\nAvg\t\t%" PRIu64 "\t%.3lf\n", c_total/samples, t_total/samples);
    
  return 0;
}
