#include "curve.h"

// Curve y² = x³ + x, with CM by √-1, j-invariant 1728
void init_curve(proj *E0) {
  E0->x = fp2_0;
  E0->z = fp2_1;
}

const int dist_trace = 0;
const int dist_norm = 1;
const int q_norm = 1;

// distorsion map on E0, twisted Edwards form
// (X:Y:Z:T) ↦ (iT : Z : Y : iX)
void ted0_dist(point *Q, const point *P) {
  point Pcopy = *P;
  fp2_mul3(&Q->x, &fp2_i, &Pcopy.t);
  Q->y = Pcopy.z;
  Q->z = Pcopy.y;
  fp2_mul3(&Q->t, &fp2_i, &Pcopy.x);
}

// distorsion map on E0, montgomery form
// (X:Y:Z) ↦ ( -X : iY : Z )
void mont0_dist(proj *Q, const proj *P) {
  fp2_neg2(&Q->x, &P->x);
  Q->z = P->z;
}

// distorsion map on E0, montgomery form
// just the same
void montxy0_dist(proj2 *Q, const proj2 *P) {
  fp2_neg2(&Q->x, &P->x);
  fp2_mul3(&Q->y, &fp2_i, &P->y);
  Q->z = P->z;
}

// frobenius map on E0, twisted Edwards form
// ( X : Y : Z : T ) ↦ ( X : Y : Z : T )^p
void ted0_frob(point *Q, const point *P) {
  fp2_frob2(&Q->x, &P->x);
  fp2_frob2(&Q->y, &P->y);
  fp2_frob2(&Q->z, &P->z);
  fp2_frob2(&Q->t, &P->t);
}

// frobenius map on the twist of E0, twisted Edwards form
// ( X : Y : Z : T ) ↦ ( τX^p : Y^p : Z^p : τT^p )
// where τ is the twisting factor
void ted0_frob_twist(point *Q, const point *P) {
  fp2_frob2(&Q->x, &P->x);
  fp2_frob2(&Q->y, &P->y);
  fp2_frob2(&Q->z, &P->z);
  fp2_frob2(&Q->t, &P->t);
  // fp2_mul2(&Q->x, &tau);
  // fp2_mul2(&Q->t, &tau);
  fp2_mul2(&Q->x, &non_residue_p_minus_1_halves);
  fp2_mul2(&Q->t, &non_residue_p_minus_1_halves);
}

// frobenius map on E0, montgomery form
void mont0_frob(proj *Q, const proj *P) {
  fp2_frob2(&Q->x, &P->x);
  fp2_frob2(&Q->z, &P->z);
}

// frobenius map on E0, montgomery form
void montxy0_frob(proj2 *Q, const proj2 *P) {
  fp2_frob2(&Q->x, &P->x);
  fp2_frob2(&Q->y, &P->y);
  fp2_frob2(&Q->z, &P->z);
}

// sqrt(-q) map in E0, Montgomery form
// Same as distortion map
void mont0_sqrt_minus_q(proj *Q, const proj *P) {
  mont0_dist(Q, P);
}

// sqrt(-q) map in E0, twisted Edwards form
void ted0_sqrt_minus_q(point *Q, const proj *E, const point *P) {
  (void)(E);
  ted0_dist(Q, P);
}
