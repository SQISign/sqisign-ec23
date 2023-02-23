#include "mont.h"
#include "tedwards.h"
#include "precomputed.h"

// well known supersingular curve (Montgomery form)
void init_curve(proj *E0);
// distorsion map on E0, twisted Edwards form
void ted0_dist(point *Q, const point *P);
// distorsion map on E0, Montgomery form
void mont0_dist(proj *Q, const proj *P);
// distorsion map on E0, Montgomery form
void montxy0_dist(proj2 *Q, const proj2 *P);
// frobenius map on E0, twisted Edwards form
void ted0_frob(point *Q, const point *P);
// frobenius map on the twist of E0, Montgomery form
void ted0_frob_twist(point *Q, const point *P);
// frobenius map on E0, Montgomery form
void mont0_frob(proj *Q, const proj *P);
// frobenius map on E0, Montgomery form
void montxy0_frob(proj2 *Q, const proj2 *P);
// sqrt(-q) map in E0, Montgomery form
void mont0_sqrt_minus_q(proj *Q, const proj *P);
// sqrt(-q) map in E0, twisted Edwards form
void ted0_sqrt_minus_q(point *Q, const proj *E, const point *P);

// Trace and norm of the distortion map
extern const int dist_trace;
extern const int dist_norm;
// q = 4n - t² or q = n - t²/4, where n = dist_norm and t =
// dist_trace, depending on whether t is odd or even
extern const int q_norm;
