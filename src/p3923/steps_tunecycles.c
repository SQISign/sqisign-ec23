int steps_guess(long long *bs,long long *gs,long long l)
{
  /* l=3: bs=0 gs=0 bench=5370 baseline=5264 */
  /* l=5: bs=0 gs=0 bench=7192 baseline=7292 */
  /* l=7: bs=0 gs=0 bench=10172 baseline=10394 */
  /* l=11: bs=0 gs=0 bench=16192 baseline=15832 */
  /* l=13: bs=0 gs=0 bench=17252 baseline=17106 */
  /* l=17: bs=0 gs=0 bench=21942 baseline=21636 */
  /* l=19: bs=0 gs=0 bench=26846 baseline=25372 */
  /* l=29: bs=0 gs=0 bench=40172 baseline=39962 */
  /* l=37: bs=0 gs=0 bench=49336 baseline=51482 */
  /* l=43: bs=0 gs=0 bench=54160 baseline=53914 */
  /* l=47: bs=0 gs=0 bench=62568 baseline=64248 */
  /* l=79: bs=0 gs=0 bench=97938 baseline=97492 */
  /* l=157: bs=0 gs=0 bench=191468 baseline=191866 */
  /* l=197: bs=0 gs=0 bench=265988 baseline=265970 */
  if (l <= 197) { *bs = 0; *gs = 0; return 1; }
  /* l=239: bs=8 gs=6 bench=327096 baseline=333366 */
  if (l <= 239) { *bs = 8; *gs = 6; return 1; }
  /* l=263: bs=2 gs=1 bench=341134 baseline=344964 */
  if (l <= 263) { *bs = 2; *gs = 1; return 1; }
  /* l=271: bs=10 gs=6 bench=336010 baseline=356494 */
  if (l <= 271) { *bs = 10; *gs = 6; return 1; }
  /* l=281: bs=0 gs=0 bench=359620 baseline=370664 */
  if (l <= 281) { *bs = 0; *gs = 0; return 1; }
  /* l=283: bs=10 gs=7 bench=346786 baseline=355810 */
  if (l <= 283) { *bs = 10; *gs = 7; return 1; }
  /* l=307: bs=12 gs=6 bench=371100 baseline=390810 */
  if (l <= 307) { *bs = 12; *gs = 6; return 1; }
  /* l=461: bs=14 gs=8 bench=511016 baseline=586776 */
  if (l <= 461) { *bs = 14; *gs = 8; return 1; }
  /* l=521: bs=16 gs=8 bench=567576 baseline=680830 */
  if (l <= 521) { *bs = 16; *gs = 8; return 1; }
  /* l=563: bs=14 gs=10 bench=571712 baseline=698498 */
  if (l <= 563) { *bs = 14; *gs = 10; return 1; }
  /* l=599: bs=16 gs=9 bench=610488 baseline=761310 */
  /* l=607: bs=16 gs=9 bench=632488 baseline=754220 */
  if (l <= 607) { *bs = 16; *gs = 9; return 1; }
  /* l=619: bs=14 gs=11 bench=632396 baseline=787534 */
  if (l <= 619) { *bs = 14; *gs = 11; return 1; }
  /* l=743: bs=16 gs=11 bench=750962 baseline=943926 */
  if (l <= 743) { *bs = 16; *gs = 11; return 1; }
  /* l=827: bs=16 gs=12 bench=830128 baseline=1062368 */
  if (l <= 827) { *bs = 16; *gs = 12; return 1; }
  /* l=941: bs=18 gs=13 bench=906868 baseline=1195968 */
  if (l <= 941) { *bs = 18; *gs = 13; return 1; }
  /* l=2357: bs=32 gs=18 bench=1860710 baseline=2993386 */
  if (l <= 2357) { *bs = 32; *gs = 18; return 1; }
  /* l=3923: bs=38 gs=25 bench=2854892 baseline=5077714 */
  if (l <= 3923) { *bs = 38; *gs = 25; return 1; }
  return 0;
}
