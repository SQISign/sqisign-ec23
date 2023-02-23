#include "precomputed.h"

// Field constants
fp minus_one_half =
  { 0xfffffffffffffffdULL, 0xfa0811ef4bdada15ULL, 0x25934f2a29bd0b91ULL, 0x1d9607c09e15bf92ULL };
fp minus_one_third =
  { 0x5555555555555553ULL, 0x6d22b89dacc6a99fULL, 0x561d3d92923db3aaULL, 0x253bcd24092d5c1cULL };
// fp2_non_residue()^((p-1)/2)
fp2 non_residue_p_minus_1_halves = {
  { 0xa9048ce93ca6c160ULL, 0xe3d53238bf94dc76ULL, 0xe8fd778c83148dc6ULL, 0x1581f21244a1b9a4ULL },
  { 0x5345b0238bd071b8ULL, 0x9cca7ec6ef7a1394ULL, 0x233cd56fa3d9b138ULL, 0x1ae26e96d5ca280eULL }
};
