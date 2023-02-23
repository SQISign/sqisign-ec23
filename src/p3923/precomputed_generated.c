#include "precomputed.h"

// each basis entry is a triple of the form P,Q,P+Q
// this is initializing the point using the classical representation {0,...,p-1} for elements in GF(p).
// We don't use this representation for actual computation but rather the montgomery representation (the conversion is made in init_precomputations using the fp_enc function)
// hence the ***_uintbig[] defined below should not be used in any actual piece of code.

const uintbig torsion_basis_sum_uintbig[3][2][2] = 
{ { { { 0x213885624a6b7d44ULL, 0x4d6e1ef6903a099aULL, 0x299ed0105935f852ULL, 0x28331f2995680c56ULL },
  { 0xffc99c7e2ebb3649ULL, 0x72f6fff68ae29b4dULL, 0xe835b1af2f5739f3ULL, 0x2cb55b2fc8b5e159ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x94bb6ecff0316209ULL, 0xb834eb5afe4cb90eULL, 0x38e852f775acd0a9ULL, 0x13ab93e227fb0c9dULL },
  { 0xc27d475a3b45d76dULL, 0xd40b4a63075e4bbULL, 0xa1b785ab7f77bb04ULL, 0x7bcd661710ccf01ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xc1a26ad7fb4ed0d3ULL, 0xc4ca693d1a4246a3ULL, 0x91ecb5acf631c597ULL, 0x1d7075d6db4d712aULL },
  { 0x36eacbbf66a4380cULL, 0xb230708edeadcdf7ULL, 0xa14b18cdfb176465ULL, 0x7fa82bdce54e422ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };
const uintbig torsion_basis_twist_sum_uintbig[3][2][2] = 
{ { { { 0xfb8ed738ca8a9f8bULL, 0x46769fda2110def2ULL, 0x4cea29692e086f0ULL, 0x20d39eae298c2a88ULL },
  { 0xcc7504cb468eb45cULL, 0xedaa852c1ad6f766ULL, 0x233f01b069580665ULL, 0x1111c36121e04926ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x59dd381b2a6bc0d7ULL, 0x3b0ce6debadd47d4ULL, 0xb421ba7dfc3d868dULL, 0x23c96667ab2e7da0ULL },
  { 0x40df8da4e2aecadcULL, 0xcad76535a8b51d90ULL, 0x4e8e6c3025d6b098ULL, 0x32de005b981a8800ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6ff29f58a97b242dULL, 0xe58e29f51526caedULL, 0x4b624d7254b2e1ceULL, 0x15b280988ec3cbULL },
  { 0x83b66df857e42f33ULL, 0xff384f0a9195d22ULL, 0x484c65dce93c3709ULL, 0x299453f7583e0546ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };

const uintbig torsion_basis_ted_sum_uintbig[3][4][2] = 
{ { { { 0x534595cd82fbf653ULL, 0xd780d34a006c93deULL, 0xf6937526ea20bdf7ULL, 0x2186af25819d09e7ULL },
  { 0xdf17b43fadac4f12ULL, 0x25c6b4cf0f60c8a9ULL, 0x5f9d514061153d5fULL, 0x691a37bee6ba7faULL } },
  { { 0xe70ef4dfb906d0faULL, 0xb62330f402e998c6ULL, 0x111bd88ccd1602c1ULL, 0x17add22dc5b350eULL },
  { 0xc2b7e9798da4b9c2ULL, 0x882974e022a7f6eULL, 0xe2c9d6aa3a112048ULL, 0x31bdba7d584518c5ULL } },
  { { 0xcad0fa0d77ace963ULL, 0x17840abcd0a773bdULL, 0xa4979b33262641e6ULL, 0xa095f879a966e77ULL },
  { 0x8100109889790a13ULL, 0x3ea93c5d62930f36ULL, 0x13b0195ae4067566ULL, 0xa4317029b45cf0eULL } },
  { { 0xbe6cae2913ca8815ULL, 0x8de59cd0eae8162eULL, 0xf3c8f2f0bd4a625dULL, 0x24e6ddaaa3c2e7e2ULL },
  { 0xffedb7056eeaf8c6ULL, 0xb29b6e55bb108303ULL, 0x1438e9554fdbb500ULL, 0x13370bdd0c5fbef2ULL } } },
 { { { 0xc6fa78c7750785d8ULL, 0x79028194a47fc334ULL, 0x58e709dcc3ad806fULL, 0x268ba42b3d9c4af1ULL },
  { 0x48ba9fd5d356b15dULL, 0xef25eaea5a144369ULL, 0x53cd1a6da4b7f26ULL, 0x13039fed552d4fe1ULL } },
  { { 0x9fa6d2e305ce7cfULL, 0xfb326a9eddb4cc61ULL, 0x2dc8f62a20963883ULL, 0x2e58b71ccc4647d6ULL },
  { 0xd76404edd304ad5ULL, 0x5e3563b9ca3322f3ULL, 0x9915e8ed7f09929eULL, 0x1de3a15f4eadf14dULL } },
  { { 0xc478339460112ff8ULL, 0xeb160b67cecf87bbULL, 0xb97d6572e6702b07ULL, 0x32ecf871763d30f4ULL },
  { 0x803acd150a509fb5ULL, 0x6b2dfa5e5f013d72ULL, 0x557012dd82b62d98ULL, 0x2da36b1d28847849ULL } },
  { { 0xffd7214d4b2fc4acULL, 0x93095ccc1e3d3dcdULL, 0x26b6d98723a37efULL, 0x30f6e89c3bf0c463ULL },
  { 0x6f96e5d231a1c3b0ULL, 0x1ff733a9483f32c7ULL, 0xde47a0e4db03b5cULL, 0xb923520cb46866ULL } } },
 { { { 0xd4ff2e3a6c19c81fULL, 0x97ad1681c08fe6aULL, 0xea6d5d76fd1ca981ULL, 0x204bb34bedf1955bULL },
  { 0x85c0ede58ad370faULL, 0xa02d51df6568f371ULL, 0x785ddff0352b2fc6ULL, 0x2b3c39d43fad8e38ULL } },
  { { 0xa65eaea67e2d1a73ULL, 0x5b612668c5e71f88ULL, 0x502aeb473c953295ULL, 0xf7a685805a6e721ULL },
  { 0x1f3c1fad74259933ULL, 0x813dfdee9b907d21ULL, 0xbd6dad8558d2ddb7ULL, 0x19e34c5566e5e8bbULL } },
  { { 0xfcd18b5b67e401aeULL, 0x761b9c2929ad9232ULL, 0x656012145b806851ULL, 0x60967bb907046c6ULL },
  { 0x6b0307eb84afd57dULL, 0x6bce37a3803f481bULL, 0x44a8c04d0aa9d176ULL, 0xf8ad83c515ad243ULL } },
  { { 0x32ff1d0805833509ULL, 0x44ad7b1a46b45349ULL, 0xa7a6e09d60b620efULL, 0x2e483389137a6b9bULL },
  { 0x13b32c4666159c2bULL, 0xcff0affab26a56ffULL, 0xd7fd4d46ef4b5b2ULL, 0x13530c9e0b6eb3b2ULL } } } };
const uintbig torsion_basis_twist_ted_sum_uintbig[3][4][2] = 
{ { { { 0x2d330c494dcff4a6ULL, 0x718e854e6ab29bcdULL, 0x3c68e5701295560eULL, 0xc51cee79c309d59ULL },
  { 0x3dca338a1d9ba8d8ULL, 0x729c723e635796b0ULL, 0xec62bc140cb0aa5fULL, 0x27d2be3f8d48d805ULL } },
  { { 0xb0be24f0f7898bbaULL, 0x61b6973f93f9c15cULL, 0x2f0d2d8fcfef0620ULL, 0x1213ecd2224f4f3eULL },
  { 0x97c709c4cdc0a140ULL, 0xcf6a6a7369829275ULL, 0xad9245df1bfd665ULL, 0x2d6533f71cae94ebULL } },
  { { 0x117f3f2043f627d2ULL, 0xe19b16115fef2a82ULL, 0x4216ee4f4498d68bULL, 0x2d218211a466654ULL },
  { 0x27fd3b254003384cULL, 0xd02757439a3cef4fULL, 0x2ec2562aa1c89b81ULL, 0x32ee5a1e2f19d439ULL } },
  { { 0x4911459f55c4702aULL, 0x7f6b707b170e3228ULL, 0xda3b930779cdf047ULL, 0xfb874ae4c63112eULL },
  { 0x427435994573fdf6ULL, 0xae25d88f2fe1c1baULL, 0xfde53bc2a1ac8d91ULL, 0x948a6aedcbee6ceULL } } },
 { { { 0xac6b273a7a889686ULL, 0xb19544f148b91d5aULL, 0x679c3d6ff4c2e88aULL, 0x339a55d0e3029b1fULL },
  { 0xcfc2a5239d90f144ULL, 0x3dcb0b6ce808b4edULL, 0xf93031668398ed47ULL, 0x20d35f4715f51ffbULL } },
  { { 0xc77b7e493d887331ULL, 0x27e5d7c0978b17e9ULL, 0x583a03d7cfd714a6ULL, 0x322a3a3f86d5cccfULL },
  { 0x2f5fe30c5db7f17ULL, 0xafb8b2635d909bddULL, 0x7f9b68caf772f0dULL, 0x130834020dd7df87ULL } },
  { { 0x97dc405348a16fa9ULL, 0x42c991030cd00f44ULL, 0xea56e4d4d6b8501eULL, 0x2b94ce25e9919672ULL },
  { 0xacd02311b64fb381ULL, 0xf1d6fba716cbaaf1ULL, 0xf6eec5c766c02040ULL, 0x35f2645af96972fULL } },
  { { 0x39657e62799b0af6ULL, 0x68f69d55eb5cd9f5ULL, 0x7ad570e7b399fa81ULL, 0xc065dfb73b19b35ULL },
  { 0x1be42d849a55dc3ULL, 0x6bdc977120575d87ULL, 0x5ee1e9b1e3b59bc2ULL, 0x2ebe42828b0588c3ULL } } },
 { { { 0x2bfe2d4c8331fedaULL, 0xcdf8ec461a414625ULL, 0xfa70d89dd75b66e7ULL, 0x2faa138892ac3739ULL },
  { 0xfeee85eddb05321dULL, 0xbb1ad4f223cbd0beULL, 0x201c345159071eafULL, 0xcd598d0de6c7be6ULL } },
  { { 0x19c89e8871911b9eULL, 0xc27b587556d4cc8eULL, 0x7f87802d5cb844c0ULL, 0x16c27c1603ded39cULL },
  { 0xab4f338ad0dd3240ULL, 0xea1c0f0fcbff0aa2ULL, 0x4ed1148f67f33672ULL, 0x29e0c796464a6c4cULL } },
  { { 0xeb6cc69372e442e7ULL, 0xefe8e863549b14eaULL, 0x10373a0f52c2fecaULL, 0x94f073265a8b41ULL },
  { 0x90b5e99b55f118e8ULL, 0xe7187a3a5e837ff0ULL, 0x46e7634da9e88913ULL, 0x19a9b4b302197fd2ULL } },
  { { 0xe7608359ab1ce8c6ULL, 0x72c4ede08a03e06aULL, 0x8619e560bc25c0e7ULL, 0xc4149a4ade3c47ULL },
  { 0x3387ec62b8db48d4ULL, 0xf2246397eed10610ULL, 0xf0b17047e4d9dcfaULL, 0xad2f2e216339b8eULL } } } };

const uintbig torsion_basis_uintbig[13][3][2][2] = {
{ { { { 0xca79bcd5c5a7e63aULL, 0x313ab9262d9da56aULL, 0xc4ea66ce188eb8e0ULL, 0xcfd6cef9e574b35ULL },
  { 0x4a5043431fa1747ULL, 0xbe19cf075cb6d9a8ULL, 0x6b117cdd689eb9d8ULL, 0x3296d9eab225ad1eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xf1cf9f835f281c31ULL, 0x91db2668a96b17d2ULL, 0x26811ee0f703ae51ULL, 0xb0db510cf8d9fdbULL },
  { 0x32dfd0a965f3a615ULL, 0x2c9b52e82524afe0ULL, 0x7b1377341b6e5e62ULL, 0x2f100080762729a2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x11332d25495ba511ULL, 0xb6ab7bae2ac965d2ULL, 0x70c5c72f1c958f82ULL, 0x11076cdf4759281eULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xae96eac9f758e4feULL, 0xf7070b410c4840e4ULL, 0x2255699aa1b8c4cfULL, 0x2d4f53eaba693adbULL },
  { 0xa1f55f84fd00b02bULL, 0xe92639d027a0ac19ULL, 0x80bd7d89daa9eba3ULL, 0x2fc432be18268829ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x1b53499c06cd4bceULL, 0xe32b364225a92dfaULL, 0xef44d4ae2f317acULL, 0x2babb0a8c9429c8aULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5d3f867674ea8bf8ULL, 0x7fd58e55282461f1ULL, 0x36d3ce80ecf87087ULL, 0x224d38cd5b2de262ULL },
  { 0xcba36b72639843dcULL, 0x9e7ab46b5ade46deULL, 0x93f2633a3aa2ec0cULL, 0x3637d597a2b1075ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x9c27fe9307e69650ULL, 0x4c63762fad0a1031ULL, 0x97bec4728a222124ULL, 0x7282a5aabfba6c6ULL },
  { 0xcaf15c7a8ee3fc1bULL, 0x7007be22df5ba900ULL, 0xecb91cc67a9a8facULL, 0x278868d4eeb30aa8ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x83a7d6eaa93df142ULL, 0xb4dac9da4d5aece0ULL, 0xb27bf1560f52fafULL, 0x109ed66c1e3421a6ULL },
  { 0x69ec151294ae4d89ULL, 0x3b79cb7676f04b1cULL, 0xf6da1327df8eb335ULL, 0x2df3469c26521126ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xba51938efcaefc8bULL, 0x4df5af262f22757bULL, 0xd8424bb6af254ab3ULL, 0x2c083d2493efd087ULL },
  { 0xaff2de8d7475f7f9ULL, 0x41ecddfeca6e9f89ULL, 0x6fabe462419153baULL, 0x211004030b8bfa86ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xf39e9875bb64f725ULL, 0xee94ac63761b9818ULL, 0x8ab17d0cce189bc5ULL, 0x4a702f701afd28fULL },
  { 0x4caee109e0952f9dULL, 0xadd572b3265aff07ULL, 0x441e9f5ac712a3b6ULL, 0x19cbb80bce7a2e94ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xdc6b911866fc7af7ULL, 0x9365765cb0023ce8ULL, 0x2c2f37233bc57854ULL, 0x14f07d8b324826aaULL },
  { 0xdd3af4636832d325ULL, 0xba4fb795d8dfb700ULL, 0x324122295913bb39ULL, 0xc9fec8bbf8eb68dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9b09daa6c01fc174ULL, 0xb70e8a755de5d339ULL, 0x2fe06bec032b19edULL, 0x2b89f69340f45624ULL },
  { 0x9e38cd5c41c40a63ULL, 0x679b6f5ec3cdcf2aULL, 0xe5573817d5953d0dULL, 0x2031c87f7e568cffULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x55363246d2563c7bULL, 0xf250f04b327698f5ULL, 0xf7d2b27cf3c4fe1bULL, 0x248e707f549d0ee3ULL },
  { 0x4f752bad0c385447ULL, 0x7ac6d23c9f9824f1ULL, 0x2eb2c68f85f317fbULL, 0x8c5326ff5acad83ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x520363e398acee71ULL, 0x7570eeb7156eeb00ULL, 0x3e722494384cce44ULL, 0x1b82071778cb5769ULL },
  { 0xe2dd34b6a00e7453ULL, 0x2b0ea3f6c653174bULL, 0xd8ae2b7b7511e017ULL, 0x4eeed337998c16aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x155a26c5fd447708ULL, 0xcc4893c4dd29020fULL, 0xf4cc2a825e24242aULL, 0x2de8801419850bd6ULL },
  { 0xd216a9144a0c713eULL, 0xc655e566b6a25aa5ULL, 0xd35aa411241ad2abULL, 0x1a1ff06019032414ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x7d11c616abb1b8abULL, 0x59d0b0e5c1137ae9ULL, 0xbaf721cde4850946ULL, 0x15b43f6a39248929ULL },
  { 0xc984465a6c667304ULL, 0xa995b785c22f7090ULL, 0x8224e5dc647d44aULL, 0x306bebe8335d28c9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x2b3120313cfac646ULL, 0x3388c090e56e3ab0ULL, 0xcbf515a450ba86c7ULL, 0x2c4670ef7f13bd1ULL },
  { 0xee82b49cf28a8b3cULL, 0xda0bddc68d066175ULL, 0xc9618402e8f7c24aULL, 0xe698fbd56c2f937ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9f7e6cb8c98cd4f1ULL, 0xa115fe941eeb9ed4ULL, 0x78228e970c477260ULL, 0x2ec0928fda4458a3ULL },
  { 0x4f5026fa7cbe9bb1ULL, 0xaaf21d823c891921ULL, 0x806fa69fb50988b2ULL, 0x19620f9091a7c9b3ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x59193a45bc18b604ULL, 0xd8afde7d98d31f42ULL, 0x5b3d6aef3829f600ULL, 0xe6d1ecf9ef838baULL },
  { 0x9c69c1b388b24fe5ULL, 0xf8d3aa0823e8058bULL, 0x2dee6d22670e4bc6ULL, 0x121ce12ac0ce6b69ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x4ad5c62c6ec9574bULL, 0x6c189d2604465c66ULL, 0xcd82a0693e5078e2ULL, 0x309f01c02a314a9fULL },
  { 0x622cc254a9b05119ULL, 0x7faf5559b9f14706ULL, 0x3a1055ab938fed8ULL, 0x18795276b654e87ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0xa01ba98c99fae511ULL, 0x9d3e3677aa9941e3ULL, 0x227213cbf16b7601ULL, 0x31c2de9487fd7aabULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xe1f53768b047660cULL, 0x1fc6c063610766cfULL, 0xf154033ed0be4aa0ULL, 0x72ee1ce3a1aa06dULL },
  { 0x8476804dfa2ba689ULL, 0x1f04de71ea61e1adULL, 0x896f4d7d57109d7eULL, 0x2e109ed0eb0bd336ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb27ebc0209d367d4ULL, 0xb4b978bed2cf7e9ULL, 0x75e8326f62dccc80ULL, 0x108d5154ecc95e5bULL },
  { 0x630893082e942b00ULL, 0xae91edc0fc354183ULL, 0x12cdb31b981b1409ULL, 0x2275115676e8d5bcULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb514dfaf5cae8aeeULL, 0xeef464c4addd45bbULL, 0x9d87250690986b40ULL, 0x12a2d866b4fcf3f3ULL },
  { 0x92efc6bbed6565dULL, 0x5de4d5c8a61e0fcfULL, 0xa11123f6c70fce3aULL, 0x1e607b4a9e9619c4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x96879dd26ac60719ULL, 0x9a3b550080854b18ULL, 0xb7e5fcf0cd1dbcb4ULL, 0x91dd51dd320fb6bULL },
  { 0x3bfbe5687f88a3ddULL, 0x71bfb8f06038c174ULL, 0x493e9a95db0a63d5ULL, 0x3e505bbc9aaf8daULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x20845e936e22519eULL, 0xde91d3f2689c5ae0ULL, 0xaa8b7f6ccfc3c104ULL, 0x3002136fe6af5cc4ULL },
  { 0xf4cf2043d49aab5eULL, 0xc3e44e86135f58bULL, 0x67cbef18668ce22bULL, 0x21dd285b7738b2bcULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe1422de4a4b6707aULL, 0xef7e703025f2cd39ULL, 0x343e2722911d66a1ULL, 0x2415e4fb68fb9e60ULL },
  { 0x2de2b2e82b474f1fULL, 0x94de88bc3d24889eULL, 0xec273ddeccc939e0ULL, 0x169374fc9bc25fb0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x8eef168c01714995ULL, 0x8bb443e637a2a72eULL, 0x602fb5688579631cULL, 0x1ea156ee0435ca5fULL },
  { 0x53aa2817d23ad13cULL, 0xd9593c640dc52a7eULL, 0xf45128d37a33dabaULL, 0xcd4c4ade1dcbfdaULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x2359ec54fc4bc2e2ULL, 0x6f6c9f93e835547eULL, 0xb7981b906b2919c5ULL, 0x27994f19f2e2346fULL },
  { 0xd705056e4b81109dULL, 0xeebae00a1e4144adULL, 0xec49132a574e5660ULL, 0x1d9e076b64ece2b6ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xdd5fd172bc2c6e1dULL, 0xb3b4aadb8b617f83ULL, 0xf98c1bcd1f391fadULL, 0x40ffbcd4fc7eb7eULL },
  { 0x791d176eba21d864ULL, 0xc095b93b1c1a7175ULL, 0xa56190e084d3ec7eULL, 0x1b1e8dbc43045b50ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xb92761a80754184bULL, 0x1f22ca394560a178ULL, 0xc5dbe123c4b3ed45ULL, 0x219dc6ee073fe456ULL },
  { 0x7a8409a6f397dc52ULL, 0xa7a29fb8f335380fULL, 0x36c50bc5ecc25c44ULL, 0x1d78a8da440f199aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6a467d33a2395b84ULL, 0x307ce5a6fded6f89ULL, 0xaec74485740ccb01ULL, 0xd22cf311043ad05ULL },
  { 0x35eb1a2a7b966a2cULL, 0x2dbc566db3bed58dULL, 0x43d1e88906d27c9ULL, 0x1f9a6af1356598f0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xbb905eee2c36df66ULL, 0x4b7988047c900fe0ULL, 0x8f6460ee7c23a0a9ULL, 0x8dda1e3eef381bcULL },
  { 0xa44fcbd80a350a2bULL, 0xdf43fa56bc58ffebULL, 0x13030e53817b34f3ULL, 0x31e1a91b88449247ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xcadbb6b1f05a5efeULL, 0x3ed8026e89391bb8ULL, 0x26ed1ec271c31eccULL, 0x1e81e2704d1254f2ULL },
  { 0x351130f75be685e1ULL, 0xa8a241f25f6004ceULL, 0x84b1c6f9642f5649ULL, 0x17e8567f70d805cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5b4d5facb515b017ULL, 0x20cdb05a452b06d4ULL, 0x856502f2af2c0397ULL, 0x1ac08ef5f9322f2aULL },
  { 0x3be2de753264bca2ULL, 0x715abb832cb30248ULL, 0xb6a89491703078f1ULL, 0x241dd127dfa76b1bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x442b210f2eb73ddaULL, 0xf4854028ca719dd1ULL, 0x27c95f3eda5b8aa0ULL, 0x8737b7f01972088ULL },
  { 0xa2613faf633f2181ULL, 0xe89c5de20cc952adULL, 0x8f2918cc66f94cabULL, 0x2cccf4ae00fac777ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x4d0543fd3ea63fafULL, 0xff50033996610147ULL, 0x7b3ff9165f98025cULL, 0x28dcbead94c21cc1ULL },
  { 0x29e2b3a4c8d7aad2ULL, 0x769a2b0f87464751ULL, 0x477768831c4835c5ULL, 0x347dafcfa06b68eeULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x4babf341c17883acULL, 0xe7b61a2369018f0cULL, 0xfde49dfeb6af4b62ULL, 0x262aaca05b45967fULL },
  { 0x7d248b1201ed9916ULL, 0x49f1b6a6d49dca7ULL, 0xd819b2f02f807ed4ULL, 0x309f39c00e019cb9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xcaf0d30b32ec6532ULL, 0x7788a5848faf9123ULL, 0xabe3676db61bb964ULL, 0x1b009d5d93a9638eULL },
  { 0x2813b78a17849b0fULL, 0xc7e2ab3adfae80f0ULL, 0xf0fb7269448adc7cULL, 0x30d811c9ea0a2281ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
};

const uintbig torsion_basis_twist_uintbig[18][3][2][2] = {
{ { { { 0x6b4af976a18e4c91ULL, 0xf049b65364e0ad13ULL, 0x603011fd6eadfb42ULL, 0x178e4bd4945a4c63ULL },
  { 0x5a2e9c7429234fb7ULL, 0xc39f47cf7fcbd4c1ULL, 0xaae2526852621903ULL, 0x2fe22d2a70b67f49ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x23770becf5b7ccfcULL, 0x1296b356767f8e27ULL, 0x7e6375f209046790ULL, 0x2b65156c0ba18c88ULL },
  { 0x3fa89e69363ec2c1ULL, 0x76a4ba0236044008ULL, 0x9eef5d933955012ULL, 0xad3e58dfb423995ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x92e0fbba69339930ULL, 0x8abcec2b826b8ea7ULL, 0xf16e28b052230988ULL, 0x31f23c1c30ff97daULL },
  { 0x62b8581b963c11d0ULL, 0x989736279bc9dcf2ULL, 0x7d2aef6a0dbcecd0ULL, 0x19cbdcd1a26100aeULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x8c05944e43d6c593ULL, 0xe4a8c1de54a761b8ULL, 0xc0b1f39e0ad6ec9cULL, 0x2fbf916ddad90fd5ULL },
  { 0x14907cd925be78e3ULL, 0x2679d5aadd2fa163ULL, 0x5401a3b7ef3cd37bULL, 0x2d6dad001730949cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5d7132c3ab6f6d39ULL, 0x2da033d9747b0b00ULL, 0x3553cc4ad49cd4efULL, 0x11d135121840ee1aULL },
  { 0xe3ada297e2c9ffd1ULL, 0xc6abc1f2171f760cULL, 0x6e4c6cbfef998bfeULL, 0x1bdbb4587341b2baULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x52b94da958bd5f41ULL, 0xf53939d5ed6368c9ULL, 0xda211806f419830bULL, 0x15f5b8c9d422b5dbULL },
  { 0x14ff50b59a16c3ffULL, 0xcdd9e85d11dbd892ULL, 0xb4790fae2598238cULL, 0x1bedebfa53b3466aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x6ef7e4b8cd06ac36ULL, 0x27a4fb88f6373a3dULL, 0x6c0e69ed2297716eULL, 0x202d1746cb97f29dULL },
  { 0xc563b5eb9c856a87ULL, 0x6b0bfd237adfd63bULL, 0xf1c210195aa28819ULL, 0x15362aa6f1c46cd1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe4856285bd35f18eULL, 0xe0d736e6a3b1e059ULL, 0x1063e981f5d92b48ULL, 0x11c854fb3beabe98ULL },
  { 0x128be21dfed6fad4ULL, 0xa9640f7934441948ULL, 0x62931811503b280cULL, 0x59ce7f462107507ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xcc3ff172fb615499ULL, 0x9d9322079d8bed07ULL, 0x424723055043b464ULL, 0x51f363be80ab2bdULL },
  { 0x8e4760bbe54d20a7ULL, 0xb079155efcc1cd8ULL, 0x5dc0b22c3a239e0dULL, 0x201708b0f43a2a5bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x92c74ffc38ab0cb0ULL, 0x43f35a156c662584ULL, 0xf40ee6b014fe87daULL, 0x28b50c01fbde0806ULL },
  { 0x6829e924db55f4cbULL, 0xd466c8c2a9f6d5f3ULL, 0x3027ede51e20f263ULL, 0x12a35100d85bc9beULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xffb610f95cc4e06bULL, 0x8aff728422bbd3a6ULL, 0xf58d3b38e09a787aULL, 0x5de9a38f6b03c5eULL },
  { 0x607c0920d162347bULL, 0x7e8838664d62343eULL, 0xf18c7300be6ea4f0ULL, 0x107f0de2fa68a3f4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd5b2a4ad48fb3c5fULL, 0xc1b5f8ff511dbc62ULL, 0xb99bbeb61557f8faULL, 0x2e03e7498c34f5ceULL },
  { 0xa1020548366311dULL, 0x3f43b63d5cc802b6ULL, 0xb62c16c29e26ecaaULL, 0x263551878c7a4c21ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x450406ffa47e2d1cULL, 0xedf40466b2165288ULL, 0xbd0e1edcef0d3023ULL, 0x1ce89d9420b7518bULL },
  { 0x736ae2625a227056ULL, 0xad60a59e41a5197cULL, 0xe79a5bd0617e5f6bULL, 0x1ddb1dcbf7772508ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xdeb641fa771bd061ULL, 0x47f5ce17a3f3c596ULL, 0xf4bf74b5bd3c9142ULL, 0x27c9f25997305a29ULL },
  { 0xb27d59de7f542d0bULL, 0xde263b1553a8368fULL, 0x605a3d6978a64a47ULL, 0x2833006270cada1bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x3d7d56ee44f1cbcbULL, 0x70d938b11565ee51ULL, 0x87113b2798283165ULL, 0x10da4172adfd8d6fULL },
  { 0x600350d975b38d35ULL, 0xba5066ccae6a93c0ULL, 0x3e75fec1e1a831efULL, 0x241b901a7989b2b1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xc09f2c810d4b7b24ULL, 0xc025126af590d15cULL, 0x44e9e93cbf1f4ae9ULL, 0x5534a216af5440cULL },
  { 0x8b4753f6cd56a8b7ULL, 0xd41a573acfe20869ULL, 0x595cb0cb95e9b697ULL, 0x1f9bb0946f77ff6aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x8b3e0e986ffc0c6ULL, 0x31f9daa9bd04a712ULL, 0xdcc28686e4c996b7ULL, 0xd8db76999e2e978ULL },
  { 0x920daf7565db94e0ULL, 0xa4baceedef96c2beULL, 0x8295968ab2f39c4ULL, 0x2d18debd813862ddULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x1398f509379ce250ULL, 0xfaf544ac741ecf70ULL, 0xae8af97d4131b96ULL, 0xed3d79e431b12a7ULL },
  { 0x4aa2ce26dbc4dc02ULL, 0x65841eba2e6c7c1ULL, 0xae6233b3a0421745ULL, 0x17cbba3b9f02d310ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x194db6a89a5bbc9cULL, 0xe04d614df5da7ffULL, 0xf83f018bbb929656ULL, 0x236d50db2215d83aULL },
  { 0xcc3f786184e1f479ULL, 0x55378ba61351f4eaULL, 0xd39180703e132bc9ULL, 0x3470d5fdc23722dcULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x1672d9dbaa356104ULL, 0xc06353a0d1a915c7ULL, 0x63f77b92ed83eedULL, 0x1f6a5a2bac8cf928ULL },
  { 0xcbda71b8c4f8fdabULL, 0x939546e0c059874dULL, 0x8968e33c2732597cULL, 0x2f5a07a578318200ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5fb89aa85821e53bULL, 0x369d18dcefffb3f5ULL, 0xa892d360516f596bULL, 0x32ec1bbb6a3d6805ULL },
  { 0xcc39a85aa632b53dULL, 0xd09e484b96428cfcULL, 0x9dbaa428f08bf7fdULL, 0x13972de8a1cc4e8dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x18513e31826096a8ULL, 0xba0d9e48e2bc2f05ULL, 0xdac285b575ce9a02ULL, 0x2b4af7685843296cULL },
  { 0xde7ec00fc738be08ULL, 0x8ef8e923ecfd5d16ULL, 0x2cc50eccb6f048b5ULL, 0xd1862ffcc301b0aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xec3175b8712cb15ULL, 0x187c347be7c8cd61ULL, 0x57ab3416afdf78d8ULL, 0x1f416c700c072168ULL },
  { 0xc60f34836b18b476ULL, 0x8afce031c44e1249ULL, 0x803f61a44ca4b3a3ULL, 0x463d789a1569319ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x57f3ce09a5a18163ULL, 0x656c2adc2f257f6aULL, 0xe808fc4b460c02d8ULL, 0x2c2d9a6f112da095ULL },
  { 0x62cc68557b298647ULL, 0xc51dd0b0bd25faeeULL, 0x8a3bab2fd19d21faULL, 0x2b39ef6c1a20dc9bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xd91694a7138236feULL, 0x6f712e37d9d6941ULL, 0x310fcbecc4474c52ULL, 0x1a087c1ad8ecdc9aULL },
  { 0x74e9a14f9617b3dULL, 0x60123c6c920aecb5ULL, 0xca240d9235cf6bbULL, 0x6373f86f5b68ff9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x8a1d48cc51a1bbbcULL, 0x51e6ad945494bc0cULL, 0xe5c90bacac5b04b9ULL, 0x22353bb21f6378c8ULL },
  { 0x144bdad4ac4704cdULL, 0xfa82ded043b2fd49ULL, 0xa7d6ad4efd9fe193ULL, 0x9f9b976609eae95ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9989e032da0da5a8ULL, 0x2b1f91b3a36656d8ULL, 0xc1a590824afd02bfULL, 0x23c2e6e3b33debe0ULL },
  { 0x3df45389a10f19dULL, 0x95968c03125264a7ULL, 0x14676c8582a392d1ULL, 0x2d76c38da439ddb4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x9784cf65e853f066ULL, 0x63e455bd6fde5aaULL, 0xd0e1df84089ff756ULL, 0x286455b9272cc52fULL },
  { 0xe16dc88aeb4f1e01ULL, 0x442cab51083e3812ULL, 0xc319ced60e148f0fULL, 0x14a2c5ac7d7097e0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x726f803d3209209aULL, 0x3b95970748e0ff3dULL, 0xc5a16100cc2582a5ULL, 0x327264a7f6a698d2ULL },
  { 0x1c4d98a9d12dcaecULL, 0xf7e9f08d5f51375eULL, 0xf8a382964984a3fbULL, 0x1b7b3e8ab4b9aeb4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe1848507f3ccf195ULL, 0x6ad221e3be01a891ULL, 0x3f1e6eb81186847fULL, 0x104096a647f4111cULL },
  { 0x5f329a7f7b0ecd20ULL, 0xc2d04370e53edb5dULL, 0x3ba8b86fe6f7f518ULL, 0x2c7bf558d9b8366ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xe1d66ca5227fa153ULL, 0x95934b902fb8911bULL, 0x86f187e641ba2a76ULL, 0x1bf6dfbe2989d94ULL },
  { 0xb8cd2080dd494318ULL, 0x2429230de7a98f4cULL, 0xd053d61801850b3dULL, 0x1d7285c2c21edcddULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x453f48bb84f4d873ULL, 0xa2e59f2228ab3c09ULL, 0x7da0a247856dfb24ULL, 0x1669e2d41578dbdfULL },
  { 0xf5662d45e515e87aULL, 0x83cdf5135550fddcULL, 0x760b16c2e12ea49dULL, 0x1ecbcb5bed604b1aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7602edfc9ef6430dULL, 0x72974f07d2f566d3ULL, 0x39a79fd9742057c4ULL, 0x1ec225dac05ee8ceULL },
  { 0x5c1e86edd267d653ULL, 0xf4bbaf307b3f8178ULL, 0x9ecdc4d59311a04cULL, 0x2b30115fe56f9435ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xa044ca0a2759a73fULL, 0x6e2d44cb5bc08002ULL, 0xf659d08a185c1a30ULL, 0x284afd9034673afbULL },
  { 0xdd12264eb03652aeULL, 0xf9a50e4ec3bf468bULL, 0xa78acd97da2b28cdULL, 0x238213f945f72e3cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x89650cfabf036a9eULL, 0x70df3739af876595ULL, 0x44c67e5405109598ULL, 0x30e27d6ff2b6943ULL },
  { 0x7f9f09deedd5cc17ULL, 0x4556878f02a4f05aULL, 0x9867f402f21bde25ULL, 0x143d1eead62b0ca1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x4d074e482b9b833cULL, 0x92af4ff12a716154ULL, 0xc01b95602525f3caULL, 0x2c27cf8b546e90e6ULL },
  { 0x4cd087d7a2292de1ULL, 0xe886e8d860762ce6ULL, 0x616dc077c1d84268ULL, 0x11a8594425cb01e1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xb9226f138c688e00ULL, 0xc3d99fdf801d3f38ULL, 0x9025470fdc309d83ULL, 0x28372862c994c653ULL },
  { 0x79bf1237fd99dde3ULL, 0x4709a1bbd468e289ULL, 0xcd546b912f7298caULL, 0x30c6c4c1246bf49eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x80672efdd40afde7ULL, 0xa345618f1e0745d6ULL, 0x7e5bdf5ab9acb20ULL, 0x197ef607c69cb126ULL },
  { 0xda38f7f04011190dULL, 0x2d07efaf1616173bULL, 0xf2a480529abb1342ULL, 0x323a18675a1aa266ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xfd2a410c2de14fdaULL, 0xb595e0f5e405fe59ULL, 0xe06795f836f8c38bULL, 0x187867d219e27cc6ULL },
  { 0x855d0968bbacfd71ULL, 0xd47bfdbd4842a603ULL, 0xc646d674508660d9ULL, 0x2c316f3f1409ec14ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x6dbae15fc995cdd9ULL, 0xe659b3694c8f868aULL, 0xf29ee18d1487c00ULL, 0x6009a8aa2835a76ULL },
  { 0x4909a89a9b3085cdULL, 0x37fe7c731cd1e2aULL, 0xe6d6b331f617adbaULL, 0x1175022c9401c349ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa374e88b0cf57315ULL, 0xefe02cab4531f6e6ULL, 0x8b2029a0cf776cf1ULL, 0x1d2dabdf0a03700dULL },
  { 0xfd9fb15b0f8f9e2aULL, 0x2e79b441337e9482ULL, 0x20f3fc0dbe58ef05ULL, 0x1459caa0bcb4aaf2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd1ca533d4abfaea4ULL, 0xfb84ffaf7f19a03aULL, 0xc81c232f25471e29ULL, 0x2a65088f4a8329b6ULL },
  { 0x4250db03bc097092ULL, 0x57612ee21a31f768ULL, 0xa0c7e0f96fd1b7cbULL, 0x1c18fcb01632b637ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x7fcdf2747f2a0299ULL, 0x16133ae4df170ed7ULL, 0x7a84d70de4e51718ULL, 0xf8d289f50857c38ULL },
  { 0x7bab67267267deb9ULL, 0x590c33a8ab1b2fc7ULL, 0x7420651a5b71e8aeULL, 0x388e8cb0f9b0c41ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x64e6b8886f6dad06ULL, 0xa3f71f5c2499b01dULL, 0x49c1980d4f50fc51ULL, 0x162f287b48a88fc0ULL },
  { 0xffde665f2c342fd6ULL, 0x6e953febf8f56816ULL, 0xb0a388b5c5e5f3c0ULL, 0x242aad82c0362e0dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd37c677736970298ULL, 0xa4553cab58133b84ULL, 0xd297944d8e7a6c09ULL, 0x12e0c2cd77a9b919ULL },
  { 0x5cafb3fb5c0241b9ULL, 0x8fa9e4d029eb6fcaULL, 0xbdedff7247925ce3ULL, 0x18e0b04dfa0c3558ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xf2fa42443c411cf4ULL, 0x2ceff868701c94a9ULL, 0x65773229897cf301ULL, 0x180bb893b365a0faULL },
  { 0xbe9b8a3935d168daULL, 0xbb9daa0fbf25f2f9ULL, 0xfaafeaf387eacb4bULL, 0x134aa5ec5761f9caULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x3caaa3afec35452cULL, 0x81cb40c6460f38ecULL, 0x334266282554b98dULL, 0x2867f9bb1c1021ffULL },
  { 0x6997138a409abd1bULL, 0x6a25c46b367a1204ULL, 0xc0001e56bf8cabbaULL, 0x75a2199020245b5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x3d565217e24ed1c9ULL, 0xa18042a60a881fd7ULL, 0xeb5b10642931ed52ULL, 0x2d269c1a5d0a79bdULL },
  { 0x75eb97797848b0aULL, 0x72539ae69813f666ULL, 0x58412e6f96ffa192ULL, 0x22c0384be84beb41ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xdff1911c09aa9b18ULL, 0x4d4c6084733a3a54ULL, 0x57ce2503107baa5aULL, 0x2ba97f4cd6c6439aULL },
  { 0xceb1d7e9250e6b3dULL, 0xc95cb0738638447aULL, 0x7f6454fcb23e7031ULL, 0x19c45c3ff567dd19ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xc54cad2bda7d69feULL, 0x9d456cd9372d822aULL, 0xb5dd60b8aefa359fULL, 0x858e2a2866bef86ULL },
  { 0x4e135e67e14276aaULL, 0xd05dad50ab62c633ULL, 0x11ab109effb81af4ULL, 0x2b8be729653b7f3bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb8eb3c348b17ae69ULL, 0x6024d7495567b7ccULL, 0x50376891c7c6d133ULL, 0x1921a67b87bad5a9ULL },
  { 0x743024167053113cULL, 0xa92adcb5a2f0c6a9ULL, 0xa67ada2c22d3ecd5ULL, 0x2f2f66b2ee8d815aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x50ab0b306b991f5dULL, 0x30742a26540a8f3eULL, 0x8d131aec7e0f397cULL, 0x181569cff41296d3ULL },
  { 0x5567b228a31a1519ULL, 0x7285fa4339b40eb2ULL, 0xd3799ad1064c2abeULL, 0x36b953fc5ffd3f5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xcf886a27966da78fULL, 0xeddb7b9cd69a5869ULL, 0x94c799369d6c1158ULL, 0x63340ace3b1e87ULL },
  { 0x7cfc17c31636cff0ULL, 0x9b9549afb881530bULL, 0xc54fb5fbde134562ULL, 0x1bb584fcba09e985ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5ae4fef397f5a1acULL, 0x24959cfef91d2923ULL, 0xb1c70fc42828f9d3ULL, 0x2da57ca3fb964c59ULL },
  { 0x6c57d9700dd9cf8fULL, 0xe20ef12bcee90f1cULL, 0x3943f2cb03c75b96ULL, 0x1729522440285bc4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
};

const uintbig torsion_basis_two_uintbig[3][2][2] = 
{ { { { 0x8866ae14c4dfa062ULL, 0x2cba1eebb69fa792ULL, 0x2452d1852c43cd9ULL, 0xf939462ae560649ULL },
  { 0x37b25d4c1cde278cULL, 0x8a857aa744cd3febULL, 0x578dee317b487826ULL, 0xad3aac3b16bafc1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7938dad38f4426dcULL, 0x11b208f486a29ca7ULL, 0x6eb63309ceb9556aULL, 0x2bcf8c8507c3627fULL },
  { 0x96fa6bb5bc533c17ULL, 0x24d963c5018dbbcULL, 0xdd3706dbc23733d6ULL, 0x10477981f4d53065ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x674cabbf457d8385ULL, 0xfa47857bca10ed43ULL, 0x989df8ef746d18a6ULL, 0xdf789e99a0b79c8ULL },
  { 0x4fd3cbd004af08d9ULL, 0x4a016cf295beabd9ULL, 0xffe92b94a2f3d973ULL, 0x1b9d4af710c1a9d9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };


proj torsion_basis[13][3];
proj torsion_basis_sum[3];
point torsion_basis_ted_sum[3];
proj torsion_basis_twist[18][3];
proj torsion_basis_twist_sum[3];
point torsion_basis_twist_ted_sum[3];
proj torsion_basis_two[3];


void init_precomputations_generated() {
	global_setup.action_2 = malloc(13*sizeof(GEN));
	global_setup.action_3 = malloc(13*sizeof(GEN));
	global_setup.action_4 = malloc(13*sizeof(GEN));
	global_setup.action_twist_2 = malloc(18*sizeof(GEN));
	global_setup.action_twist_3 = malloc(18*sizeof(GEN));
	global_setup.action_twist_4 = malloc(18*sizeof(GEN));

	global_setup.action_2[0] = mkmat2(mkcol2(stoi(9ULL),stoi(18ULL)),mkcol2(stoi(1ULL),stoi(16ULL)));
	global_setup.action_3[0] = mkmat2(mkcol2(stoi(6ULL),stoi(10ULL)),mkcol2(stoi(12ULL),stoi(20ULL)));
	global_setup.action_4[0] = mkmat2(mkcol2(stoi(20ULL),stoi(0ULL)),mkcol2(stoi(23ULL),stoi(5ULL)));
	global_setup.action_2[1] = mkmat2(mkcol2(stoi(2ULL),stoi(1ULL)),mkcol2(stoi(2ULL),stoi(5ULL)));
	global_setup.action_3[1] = mkmat2(mkcol2(stoi(5ULL),stoi(1ULL)),mkcol2(stoi(1ULL),stoi(3ULL)));
	global_setup.action_4[1] = mkmat2(mkcol2(stoi(4ULL),stoi(5ULL)),mkcol2(stoi(1ULL),stoi(3ULL)));
	global_setup.action_2[2] = mkmat2(mkcol2(stoi(8ULL),stoi(3ULL)),mkcol2(stoi(4ULL),stoi(3ULL)));
	global_setup.action_3[2] = mkmat2(mkcol2(stoi(8ULL),stoi(2ULL)),mkcol2(stoi(5ULL),stoi(4ULL)));
	global_setup.action_4[2] = mkmat2(mkcol2(stoi(2ULL),stoi(6ULL)),mkcol2(stoi(3ULL),stoi(9ULL)));
	global_setup.action_2[3] = mkmat2(mkcol2(stoi(13ULL),stoi(4ULL)),mkcol2(stoi(5ULL),stoi(6ULL)));
	global_setup.action_3[3] = mkmat2(mkcol2(stoi(18ULL),stoi(15ULL)),mkcol2(stoi(10ULL),stoi(2ULL)));
	global_setup.action_4[3] = mkmat2(mkcol2(stoi(8ULL),stoi(13ULL)),mkcol2(stoi(17ULL),stoi(11ULL)));
	global_setup.action_2[4] = mkmat2(mkcol2(stoi(89ULL),stoi(198ULL)),mkcol2(stoi(784ULL),stoi(752ULL)));
	global_setup.action_3[4] = mkmat2(mkcol2(stoi(154ULL),stoi(103ULL)),mkcol2(stoi(155ULL),stoi(688ULL)));
	global_setup.action_4[4] = mkmat2(mkcol2(stoi(664ULL),stoi(739ULL)),mkcol2(stoi(134ULL),stoi(177ULL)));
	global_setup.action_2[5] = mkmat2(mkcol2(stoi(558ULL),stoi(1293ULL)),mkcol2(stoi(116ULL),stoi(811ULL)));
	global_setup.action_3[5] = mkmat2(mkcol2(stoi(963ULL),stoi(1184ULL)),mkcol2(stoi(1278ULL),stoi(407ULL)));
	global_setup.action_4[5] = mkmat2(mkcol2(stoi(777ULL),stoi(0ULL)),mkcol2(stoi(944ULL),stoi(592ULL)));
	global_setup.action_2[6] = mkmat2(mkcol2(stoi(44ULL),stoi(23ULL)),mkcol2(stoi(20ULL),stoi(3ULL)));
	global_setup.action_3[6] = mkmat2(mkcol2(stoi(29ULL),stoi(29ULL)),mkcol2(stoi(19ULL),stoi(19ULL)));
	global_setup.action_4[6] = mkmat2(mkcol2(stoi(21ULL),stoi(21ULL)),mkcol2(stoi(26ULL),stoi(26ULL)));
	global_setup.action_2[7] = mkmat2(mkcol2(stoi(104ULL),stoi(189ULL)),mkcol2(stoi(47ULL),stoi(93ULL)));
	global_setup.action_3[7] = mkmat2(mkcol2(stoi(81ULL),stoi(57ULL)),mkcol2(stoi(187ULL),stoi(117ULL)));
	global_setup.action_4[7] = mkmat2(mkcol2(stoi(33ULL),stoi(67ULL)),mkcol2(stoi(119ULL),stoi(164ULL)));
	global_setup.action_2[8] = mkmat2(mkcol2(stoi(181ULL),stoi(22ULL)),mkcol2(stoi(41ULL),stoi(82ULL)));
	global_setup.action_3[8] = mkmat2(mkcol2(stoi(178ULL),stoi(218ULL)),mkcol2(stoi(104ULL),stoi(86ULL)));
	global_setup.action_4[8] = mkmat2(mkcol2(stoi(53ULL),stoi(59ULL)),mkcol2(stoi(46ULL),stoi(210ULL)));
	global_setup.action_2[9] = mkmat2(mkcol2(stoi(47ULL),stoi(82ULL)),mkcol2(stoi(69ULL),stoi(234ULL)));
	global_setup.action_3[9] = mkmat2(mkcol2(stoi(98ULL),stoi(145ULL)),mkcol2(stoi(198ULL),stoi(184ULL)));
	global_setup.action_4[9] = mkmat2(mkcol2(stoi(48ULL),stoi(266ULL)),mkcol2(stoi(266ULL),stoi(233ULL)));
	global_setup.action_2[10] = mkmat2(mkcol2(stoi(164ULL),stoi(201ULL)),mkcol2(stoi(357ULL),stoi(297ULL)));
	global_setup.action_3[10] = mkmat2(mkcol2(stoi(347ULL),stoi(196ULL)),mkcol2(stoi(380ULL),stoi(115ULL)));
	global_setup.action_4[10] = mkmat2(mkcol2(stoi(59ULL),stoi(400ULL)),mkcol2(stoi(246ULL),stoi(402ULL)));
	global_setup.action_2[11] = mkmat2(mkcol2(stoi(65ULL),stoi(460ULL)),mkcol2(stoi(428ULL),stoi(456ULL)));
	global_setup.action_3[11] = mkmat2(mkcol2(stoi(232ULL),stoi(239ULL)),mkcol2(stoi(445ULL),stoi(290ULL)));
	global_setup.action_4[11] = mkmat2(mkcol2(stoi(439ULL),stoi(450ULL)),mkcol2(stoi(36ULL),stoi(82ULL)));
	global_setup.action_2[12] = mkmat2(mkcol2(stoi(1093ULL),stoi(1981ULL)),mkcol2(stoi(498ULL),stoi(2830ULL)));
	global_setup.action_3[12] = mkmat2(mkcol2(stoi(1025ULL),stoi(1829ULL)),mkcol2(stoi(2000ULL),stoi(2899ULL)));
	global_setup.action_4[12] = mkmat2(mkcol2(stoi(2040ULL),stoi(1937ULL)),mkcol2(stoi(3494ULL),stoi(1883ULL)));
	global_setup.action_twist_2[0] = mkmat2(mkcol2(strtoi("3111118678268948423337286873280"),strtoi("2824021227687250034949644573908")),mkcol2(strtoi("1067707206858390488616522664903"),strtoi("7189932782608589030636260394563")));
	global_setup.action_twist_3[0] = mkmat2(mkcol2(strtoi("7871770392801520228810155367325"),strtoi("3509556472256979032173534318859")),mkcol2(strtoi("267336979238764942205909212612"),strtoi("2429281068076017225163391900519")));
	global_setup.action_twist_4[0] = mkmat2(mkcol2(strtoi("3566491939805748734964185829497"),strtoi("9559049652056624285254006826070")),mkcol2(strtoi("6738156842295373753619190561561"),strtoi("6734559521071788719009361438346")));
	global_setup.action_twist_2[1] = mkmat2(mkcol2(stoi(1ULL),stoi(10ULL)),mkcol2(stoi(5ULL),stoi(12ULL)));
	global_setup.action_twist_3[1] = mkmat2(mkcol2(stoi(11ULL),stoi(0ULL)),mkcol2(stoi(7ULL),stoi(3ULL)));
	global_setup.action_twist_4[1] = mkmat2(mkcol2(stoi(3ULL),stoi(4ULL)),mkcol2(stoi(9ULL),stoi(10ULL)));
	global_setup.action_twist_2[2] = mkmat2(mkcol2(stoi(15ULL),stoi(11ULL)),mkcol2(stoi(15ULL),stoi(2ULL)));
	global_setup.action_twist_3[2] = mkmat2(mkcol2(stoi(7ULL),stoi(0ULL)),mkcol2(stoi(7ULL),stoi(11ULL)));
	global_setup.action_twist_4[2] = mkmat2(mkcol2(stoi(12ULL),stoi(2ULL)),mkcol2(stoi(0ULL),stoi(5ULL)));
	global_setup.action_twist_2[3] = mkmat2(mkcol2(stoi(21ULL),stoi(18ULL)),mkcol2(stoi(28ULL),stoi(22ULL)));
	global_setup.action_twist_3[3] = mkmat2(mkcol2(stoi(36ULL),stoi(5ULL)),mkcol2(stoi(36ULL),stoi(8ULL)));
	global_setup.action_twist_4[3] = mkmat2(mkcol2(stoi(28ULL),stoi(34ULL)),mkcol2(stoi(37ULL),stoi(15ULL)));
	global_setup.action_twist_2[4] = mkmat2(mkcol2(stoi(26ULL),stoi(59ULL)),mkcol2(stoi(22ULL),stoi(53ULL)));
	global_setup.action_twist_3[4] = mkmat2(mkcol2(stoi(75ULL),stoi(44ULL)),mkcol2(stoi(13ULL),stoi(5ULL)));
	global_setup.action_twist_4[4] = mkmat2(mkcol2(stoi(31ULL),stoi(17ULL)),mkcol2(stoi(48ULL),stoi(48ULL)));
	global_setup.action_twist_2[5] = mkmat2(mkcol2(stoi(92ULL),stoi(10ULL)),mkcol2(stoi(17ULL),stoi(65ULL)));
	global_setup.action_twist_3[5] = mkmat2(mkcol2(stoi(75ULL),stoi(11ULL)),mkcol2(stoi(102ULL),stoi(83ULL)));
	global_setup.action_twist_4[5] = mkmat2(mkcol2(stoi(70ULL),stoi(115ULL)),mkcol2(stoi(55ULL),stoi(87ULL)));
	global_setup.action_twist_2[6] = mkmat2(mkcol2(stoi(85ULL),stoi(15ULL)),mkcol2(stoi(60ULL),stoi(154ULL)));
	global_setup.action_twist_3[6] = mkmat2(mkcol2(stoi(72ULL),stoi(82ULL)),mkcol2(stoi(216ULL),stoi(168ULL)));
	global_setup.action_twist_4[6] = mkmat2(mkcol2(stoi(39ULL),stoi(169ULL)),mkcol2(stoi(61ULL),stoi(200ULL)));
	global_setup.action_twist_2[7] = mkmat2(mkcol2(stoi(214ULL),stoi(187ULL)),mkcol2(stoi(200ULL),stoi(57ULL)));
	global_setup.action_twist_3[7] = mkmat2(mkcol2(stoi(218ULL),stoi(217ULL)),mkcol2(stoi(186ULL),stoi(54ULL)));
	global_setup.action_twist_4[7] = mkmat2(mkcol2(stoi(134ULL),stoi(168ULL)),mkcol2(stoi(2ULL),stoi(137ULL)));
	global_setup.action_twist_2[8] = mkmat2(mkcol2(stoi(122ULL),stoi(68ULL)),mkcol2(stoi(10ULL),stoi(161ULL)));
	global_setup.action_twist_3[8] = mkmat2(mkcol2(stoi(167ULL),stoi(28ULL)),mkcol2(stoi(56ULL),stoi(117ULL)));
	global_setup.action_twist_4[8] = mkmat2(mkcol2(stoi(127ULL),stoi(52ULL)),mkcol2(stoi(215ULL),stoi(156ULL)));
	global_setup.action_twist_2[9] = mkmat2(mkcol2(stoi(199ULL),stoi(132ULL)),mkcol2(stoi(207ULL),stoi(108ULL)));
	global_setup.action_twist_3[9] = mkmat2(mkcol2(stoi(200ULL),stoi(143ULL)),mkcol2(stoi(165ULL),stoi(108ULL)));
	global_setup.action_twist_4[9] = mkmat2(mkcol2(stoi(180ULL),stoi(40ULL)),mkcol2(stoi(276ULL),stoi(127ULL)));
	global_setup.action_twist_2[10] = mkmat2(mkcol2(stoi(254ULL),stoi(409ULL)),mkcol2(stoi(174ULL),stoi(309ULL)));
	global_setup.action_twist_3[10] = mkmat2(mkcol2(stoi(158ULL),stoi(112ULL)),mkcol2(stoi(339ULL),stoi(406ULL)));
	global_setup.action_twist_4[10] = mkmat2(mkcol2(stoi(312ULL),stoi(267ULL)),mkcol2(stoi(501ULL),stoi(251ULL)));
	global_setup.action_twist_2[11] = mkmat2(mkcol2(stoi(105ULL),stoi(267ULL)),mkcol2(stoi(542ULL),stoi(494ULL)));
	global_setup.action_twist_3[11] = mkmat2(mkcol2(stoi(277ULL),stoi(462ULL)),mkcol2(stoi(534ULL),stoi(323ULL)));
	global_setup.action_twist_4[11] = mkmat2(mkcol2(stoi(349ULL),stoi(575ULL)),mkcol2(stoi(21ULL),stoi(250ULL)));
	global_setup.action_twist_2[12] = mkmat2(mkcol2(stoi(472ULL),stoi(382ULL)),mkcol2(stoi(572ULL),stoi(135ULL)));
	global_setup.action_twist_3[12] = mkmat2(mkcol2(stoi(366ULL),stoi(83ULL)),mkcol2(stoi(142ULL),stoi(242ULL)));
	global_setup.action_twist_4[12] = mkmat2(mkcol2(stoi(585ULL),stoi(508ULL)),mkcol2(stoi(290ULL),stoi(22ULL)));
	global_setup.action_twist_2[13] = mkmat2(mkcol2(stoi(423ULL),stoi(81ULL)),mkcol2(stoi(160ULL),stoi(196ULL)));
	global_setup.action_twist_3[13] = mkmat2(mkcol2(stoi(560ULL),stoi(27ULL)),mkcol2(stoi(293ULL),stoi(60ULL)));
	global_setup.action_twist_4[13] = mkmat2(mkcol2(stoi(14ULL),stoi(187ULL)),mkcol2(stoi(325ULL),stoi(605ULL)));
	global_setup.action_twist_2[14] = mkmat2(mkcol2(stoi(556ULL),stoi(226ULL)),mkcol2(stoi(128ULL),stoi(187ULL)));
	global_setup.action_twist_3[14] = mkmat2(mkcol2(stoi(77ULL),stoi(537ULL)),mkcol2(stoi(23ULL),stoi(667ULL)));
	global_setup.action_twist_4[14] = mkmat2(mkcol2(stoi(458ULL),stoi(542ULL)),mkcol2(stoi(40ULL),stoi(285ULL)));
	global_setup.action_twist_2[15] = mkmat2(mkcol2(stoi(743ULL),stoi(151ULL)),mkcol2(stoi(134ULL),stoi(84ULL)));
	global_setup.action_twist_3[15] = mkmat2(mkcol2(stoi(535ULL),stoi(416ULL)),mkcol2(stoi(346ULL),stoi(293ULL)));
	global_setup.action_twist_4[15] = mkmat2(mkcol2(stoi(690ULL),stoi(202ULL)),mkcol2(stoi(687ULL),stoi(137ULL)));
	global_setup.action_twist_2[16] = mkmat2(mkcol2(stoi(332ULL),stoi(137ULL)),mkcol2(stoi(473ULL),stoi(609ULL)));
	global_setup.action_twist_3[16] = mkmat2(mkcol2(stoi(297ULL),stoi(921ULL)),mkcol2(stoi(514ULL),stoi(645ULL)));
	global_setup.action_twist_4[16] = mkmat2(mkcol2(stoi(583ULL),stoi(799ULL)),mkcol2(stoi(886ULL),stoi(358ULL)));
	global_setup.action_twist_2[17] = mkmat2(mkcol2(stoi(2001ULL),stoi(1492ULL)),mkcol2(stoi(934ULL),stoi(356ULL)));
	global_setup.action_twist_3[17] = mkmat2(mkcol2(stoi(91ULL),stoi(2299ULL)),mkcol2(stoi(446ULL),stoi(2267ULL)));
	global_setup.action_twist_4[17] = mkmat2(mkcol2(stoi(1360ULL),stoi(1861ULL)),mkcol2(stoi(999ULL),stoi(997ULL)));
	global_setup.action_two_2 = mkmat2(mkcol2(strtoi("31479152875545029957"),stoi(480731893264431870ULL)),mkcol2(strtoi("26250023622549539277"),stoi(5414335271874073275ULL)));
	global_setup.action_two_3 = mkmat2(mkcol2(strtoi("31479152875545029957"),stoi(480731893264431870ULL)),mkcol2(strtoi("26250023622549539277"),stoi(5414335271874073275ULL)));
	global_setup.action_two_4 = mkmat2(mkcol2(strtoi("31479152875545029957"),stoi(480731893264431870ULL)),mkcol2(strtoi("26250023622549539277"),stoi(5414335271874073275ULL)));
	for (int i = 0; i < 3; ++i) {
		fp_enc( &(&(&torsion_basis_two[i])->x)->re, &(torsion_basis_two_uintbig[i][0][0]) );
		fp_enc( &(&(&torsion_basis_two[i])->x)->im, &(torsion_basis_two_uintbig[i][0][1]) );
		fp_enc( &(&(&torsion_basis_two[i])->z)->re, &(torsion_basis_two_uintbig[i][1][0]) );
		fp_enc( &(&(&torsion_basis_two[i])->z)->im, &(torsion_basis_two_uintbig[i][1][1]) );
	}
	for (int i = 0; i < 3; ++i) {
		fp_enc( &(&(&torsion_basis_sum[i])->x)->re, &(torsion_basis_sum_uintbig[i][0][0]) );
		fp_enc( &(&(&torsion_basis_sum[i])->x)->im, &(torsion_basis_sum_uintbig[i][0][1]) );
		fp_enc( &(&(&torsion_basis_sum[i])->z)->re, &(torsion_basis_sum_uintbig[i][1][0]) );
		fp_enc( &(&(&torsion_basis_sum[i])->z)->im, &(torsion_basis_sum_uintbig[i][1][1]) );
	}
	for (int i = 0; i < 3; ++i) {
		fp_enc( &(&(&torsion_basis_twist_sum[i])->x)->re, &(torsion_basis_twist_sum_uintbig[i][0][0]) );
		fp_enc( &(&(&torsion_basis_twist_sum[i])->x)->im, &(torsion_basis_twist_sum_uintbig[i][0][1]) );
		fp_enc( &(&(&torsion_basis_twist_sum[i])->z)->re, &(torsion_basis_twist_sum_uintbig[i][1][0]) );
		fp_enc( &(&(&torsion_basis_twist_sum[i])->z)->im, &(torsion_basis_twist_sum_uintbig[i][1][1]) );
	}
	for (int j=0;j<13;j++){
		for (int i = 0; i < 3; ++i) {
			fp_enc( &(&(&torsion_basis[j][i])->x)->re, &(torsion_basis_uintbig[j][i][0][0]) );
			fp_enc( &(&(&torsion_basis[j][i])->x)->im, &(torsion_basis_uintbig[j][i][0][1]) );
			fp_enc( &(&(&torsion_basis[j][i])->z)->re, &(torsion_basis_uintbig[j][i][1][0]) );
			fp_enc( &(&(&torsion_basis[j][i])->z)->im, &(torsion_basis_uintbig[j][i][1][1]) );
		}
	}
	for (int j=0;j<18;j++){
		for (int i = 0; i < 3; ++i) {
			fp_enc( &(&(&torsion_basis_twist[j][i])->x)->re, &(torsion_basis_twist_uintbig[j][i][0][0]) );
			fp_enc( &(&(&torsion_basis_twist[j][i])->x)->im, &(torsion_basis_twist_uintbig[j][i][0][1]) );
			fp_enc( &(&(&torsion_basis_twist[j][i])->z)->re, &(torsion_basis_twist_uintbig[j][i][1][0]) );
			fp_enc( &(&(&torsion_basis_twist[j][i])->z)->im, &(torsion_basis_twist_uintbig[j][i][1][1]) );
		}
	}
    	for (int i=0;i<3;i++){
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->x)->re, &(torsion_basis_ted_sum_uintbig[i][0][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->x)->im, &(torsion_basis_ted_sum_uintbig[i][0][1]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->y)->re, &(torsion_basis_ted_sum_uintbig[i][1][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->y)->im, &(torsion_basis_ted_sum_uintbig[i][1][1]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->z)->re, &(torsion_basis_ted_sum_uintbig[i][2][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->z)->im, &(torsion_basis_ted_sum_uintbig[i][2][1]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->t)->re, &(torsion_basis_ted_sum_uintbig[i][3][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->t)->im, &(torsion_basis_ted_sum_uintbig[i][3][1]) );
  	}
  	for (int i=0;i<3;i++){
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->re, &(torsion_basis_twist_ted_sum_uintbig[i][0][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->im, &(torsion_basis_twist_ted_sum_uintbig[i][0][1]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->re, &(torsion_basis_twist_ted_sum_uintbig[i][1][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->im, &(torsion_basis_twist_ted_sum_uintbig[i][1][1]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->re, &(torsion_basis_twist_ted_sum_uintbig[i][2][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->im, &(torsion_basis_twist_ted_sum_uintbig[i][2][1]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->re, &(torsion_basis_twist_ted_sum_uintbig[i][3][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->im, &(torsion_basis_twist_ted_sum_uintbig[i][3][1]) );
  	}

}

