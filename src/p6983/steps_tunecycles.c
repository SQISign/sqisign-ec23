int steps_guess(long long *bs,long long *gs,long long l)
{
  /* l=3: bs=0 gs=0 bench=5336 baseline=5298 */
  /* l=5: bs=0 gs=0 bench=7412 baseline=7532 */
  /* l=7: bs=0 gs=0 bench=10772 baseline=10658 */
  /* l=11: bs=0 gs=0 bench=15734 baseline=16306 */
  /* l=31: bs=0 gs=0 bench=45792 baseline=44408 */
  /* l=43: bs=0 gs=0 bench=54128 baseline=55574 */
  /* l=83: bs=0 gs=0 bench=114044 baseline=116676 */
  /* l=103: bs=0 gs=0 bench=126682 baseline=126486 */
  /* l=107: bs=0 gs=0 bench=142382 baseline=146010 */
  if (l <= 107) { *bs = 0; *gs = 0; return 1; }
  /* l=109: bs=4 gs=2 bench=153568 baseline=145406 */
  if (l <= 109) { *bs = 4; *gs = 2; return 1; }
  /* l=137: bs=0 gs=0 bench=185730 baseline=180420 */
  /* l=199: bs=0 gs=0 bench=241696 baseline=242130 */
  /* l=227: bs=0 gs=0 bench=275396 baseline=275404 */
  if (l <= 227) { *bs = 0; *gs = 0; return 1; }
  /* l=419: bs=14 gs=7 bench=449808 baseline=505522 */
  if (l <= 419) { *bs = 14; *gs = 7; return 1; }
  /* l=491: bs=14 gs=8 bench=508118 baseline=592862 */
  if (l <= 491) { *bs = 14; *gs = 8; return 1; }
  /* l=569: bs=14 gs=10 bench=563110 baseline=691420 */
  if (l <= 569) { *bs = 14; *gs = 10; return 1; }
  /* l=631: bs=14 gs=11 bench=623960 baseline=761126 */
  if (l <= 631) { *bs = 14; *gs = 11; return 1; }
  /* l=677: bs=16 gs=10 bench=654780 baseline=815880 */
  if (l <= 677) { *bs = 16; *gs = 10; return 1; }
  /* l=751: bs=16 gs=11 bench=802244 baseline=1006642 */
  if (l <= 751) { *bs = 16; *gs = 11; return 1; }
  /* l=827: bs=20 gs=10 bench=858752 baseline=1116614 */
  if (l <= 827) { *bs = 20; *gs = 10; return 1; }
  /* l=857: bs=16 gs=13 bench=801850 baseline=1038378 */
  /* l=859: bs=16 gs=13 bench=804798 baseline=1046796 */
  if (l <= 859) { *bs = 16; *gs = 13; return 1; }
  /* l=883: bs=18 gs=12 bench=816334 baseline=1063794 */
  if (l <= 883) { *bs = 18; *gs = 12; return 1; }
  /* l=1019: bs=18 gs=14 bench=924364 baseline=1237068 */
  if (l <= 1019) { *bs = 18; *gs = 14; return 1; }
  /* l=1171: bs=22 gs=13 bench=1017432 baseline=1408884 */
  if (l <= 1171) { *bs = 22; *gs = 13; return 1; }
  /* l=1879: bs=30 gs=15 bench=1484564 baseline=2258266 */
  if (l <= 1879) { *bs = 30; *gs = 15; return 1; }
  /* l=2713: bs=32 gs=21 bench=1935288 baseline=3283824 */
  if (l <= 2713) { *bs = 32; *gs = 21; return 1; }
  /* l=3691: bs=38 gs=24 bench=2797162 baseline=4867644 */
  if (l <= 3691) { *bs = 38; *gs = 24; return 1; }
  /* l=4019: bs=40 gs=25 bench=2877540 baseline=5230412 */
  if (l <= 4019) { *bs = 40; *gs = 25; return 1; }
  /* l=4283: bs=38 gs=28 bench=3160004 baseline=5183308 */
  if (l <= 4283) { *bs = 38; *gs = 28; return 1; }
  /* l=6983: bs=62 gs=28 bench=4439882 baseline=8919712 */
  if (l <= 6983) { *bs = 62; *gs = 28; return 1; }
  return 0;
}
