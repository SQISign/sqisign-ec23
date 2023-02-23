
#ifndef SQISIGN_H
#define SQISIGN_H

#include <pari/pari.h>
#include "idiso.h"

typedef struct public_key {
  proj E;
} public_key;

typedef struct secret_key {
  GEN I_large_prime;
  GEN I_two;
  GEN I_T;
  special_isogeny phi_T;
  two_walk_long phi_two;
} secret_key;

typedef struct signature {
  proj E_com;
  two_walk_long sigma;
} signature;

typedef struct compressed_signature {
  proj E_com;
  uintbig *zip;
} compressed_signature;


bool mont_power_dlp(uintbig *a,const proj *A, const proj *Q, const proj *P,const proj *PQ,long e);
void init_compressed_sig(compressed_signature *comp_sigma);
void free_compressed_sig(compressed_signature *comp_sigma);
void keygen(public_key *pk, secret_key *sk);
void commitment(GEN *coeff, GEN *I, GEN *I_even, odd_isogeny *phi_com,two_walk_long *phi_com_even);
void challenge(proj *E_cha, two_walk_long* phi_chall_even, const uintbig *m, const proj *E_com, const proj *basis_plus, const proj *basis_minus, GEN *dlog, proj *basis_two);
void response(two_walk_long *sigma, uintbig *zip, GEN coeff_ker_challenge_commitment, const secret_key *sk, const proj *basis_two, const proj *E_cha);
void sign(compressed_signature *comp_sigma, const secret_key *sk, const public_key *pk, const uintbig *m);
void decompress(two_walk *walk, proj *A, const uintbig *zip, long len,long last_step);
bool verif(compressed_signature *comp_sigma, const public_key *pk,const uintbig *m);
void response_new(two_walk_long *sigma, uintbig *zip, GEN coeff_ker_challenge_commitment, GEN I_com_even, const secret_key *sk, const proj *basis_two, const proj *E_cha, bool *test);
void sign_new(compressed_signature *comp_sigma, const secret_key *sk, const public_key *pk, const uintbig *m);
void decompress_new(two_walk *walk, proj *A, const uintbig *zip, long len,long last_step);
bool verif_new(compressed_signature *comp_sigma, const public_key *pk,const uintbig *m);

#endif
