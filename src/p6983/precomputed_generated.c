#include "precomputed.h"

// each basis entry is a triple of the form P,Q,P+Q
// this is initializing the point using the classical representation {0,...,p-1} for elements in GF(p).
// We don't use this representation for actual computation but rather the montgomery representation (the conversion is made in init_precomputations using the fp_enc function)
// hence the ***_uintbig[] defined below should not be used in any actual piece of code.

const uintbig torsion_basis_sum_uintbig[3][2][2] = 
{ { { { 0xb63776217b4e5dd1ULL, 0x3c5376f4c058feaaULL, 0xe2a1fb539d6e4cf8ULL, 0x4a04fcbc5704c188ULL },
  { 0x17e7398faaa45a5ULL, 0x93eb97fbafbf4902ULL, 0xb85797acfc9c2121ULL, 0x22cb2d480660b15aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd1ecb00768f1af92ULL, 0xe8b431417eb98117ULL, 0xb9f205ffd6b32bb3ULL, 0x31fbf1ff5fbd1353ULL },
  { 0x707aee78b0fa39a9ULL, 0xd167fabd26c0bd4eULL, 0xdeb30629657e14d9ULL, 0x598fdc8091da0b12ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xab10bb7b7cd11c5ULL, 0x35a7af95c502ecbdULL, 0x1ded8f87d07dd01eULL, 0x8da875d946d57e86ULL },
  { 0xdde788e5637db632ULL, 0x19144161949ac004ULL, 0x42702fab0830e6aaULL, 0x46c00f69c9bb7ebfULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };
const uintbig torsion_basis_twist_sum_uintbig[3][2][2] = 
{ { { { 0x6b5b454508d30648ULL, 0xd3a5099fcbe18cb8ULL, 0xbf6239f07ec3c681ULL, 0x9bdfd4e8c59f846ULL },
  { 0x9efc9353ae4b46edULL, 0xd737343671b8135aULL, 0xd11bd9c6ab85f3e0ULL, 0x704d5267925785d1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6535ce9a57a33501ULL, 0x1764614eac66189fULL, 0x7d3522fb9e0ca69aULL, 0x3733e7e86d088252ULL },
  { 0x9d17fffa3701f913ULL, 0x2489c9dd066cf90fULL, 0x9624c7ef84e85e4aULL, 0x48cf991577346510ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9959c82889c4ca40ULL, 0x9e95f59fd7909c9ULL, 0xfd64a7e72a240b1cULL, 0x4705c0687630a2ecULL },
  { 0x4bc7d3d1f40d6ab5ULL, 0xde00c53cf7c88785ULL, 0x6922d39e9c0e0133ULL, 0x1f6618f26040a0e4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };

const uintbig torsion_basis_ted_sum_uintbig[3][4][2] = 
{ { { { 0x6804b5759c368aedULL, 0xe303a762e5c63059ULL, 0x3747749cba6c7eb8ULL, 0x83a1c3263858f22bULL },
  { 0x972f8dfa3fb018e1ULL, 0x4f7a10add3cca860ULL, 0x2520741e501b4528ULL, 0x44e1bcd1b671ead7ULL } },
  { { 0x294e44430eb6ac6dULL, 0xf5943954a361ddf0ULL, 0x8272c9da9acfb736ULL, 0xc16610541a3f0ffULL },
  { 0xe8007786dd6ffa53ULL, 0x9f52d83d19582bd2ULL, 0xeb6ebc0bc4b324fdULL, 0x7f5580597135a0ccULL } },
  { { 0x2b7252e7c0f14258ULL, 0xd3a260971f5d6b5fULL, 0xc4098a8043cf8703ULL, 0x4c0b32a8b1417091ULL },
  { 0xcdb04b272659a527ULL, 0x56e2680567c72a30ULL, 0x8f9b7f09b476e70eULL, 0x6859168f6515a741ULL } },
  { { 0xec996133f8656755ULL, 0x719f0de5ec38c229ULL, 0xde291d6c854ae76bULL, 0x19c140b10b116e8bULL },
  { 0xd6d1626e88c24922ULL, 0x3bdc03d1fd9d9a64ULL, 0x9c62a1be57b14634ULL, 0x40ada2e32cdaa3feULL } } },
 { { { 0x63b111283cec28d5ULL, 0xf56e132c30d29ff5ULL, 0xc884b7d969e9ee22ULL, 0x855b69628bb924a2ULL },
  { 0xb782914c6280ebbcULL, 0x7a5d171c58937b8ULL, 0xd1721a75f2a8ac39ULL, 0x215f7725336872a7ULL } },
  { { 0x586f74c0fab5b7acULL, 0xd7516456e20bea2dULL, 0x70f67dccbf5966aeULL, 0x136b65e21593a9e1ULL },
  { 0xcc7108f10c5e7914ULL, 0x2b38391d61c7feb7ULL, 0x7755834ffeb890b1ULL, 0x56fb874cee0853e5ULL } },
  { { 0x7b3feec5ff977c5fULL, 0xb84decf3c92f4d5cULL, 0x278668ca7748ec5eULL, 0x444122199b83922eULL },
  { 0xfe35b35c88a7bb96ULL, 0x83caffd0026f3962ULL, 0x6147bc6ac93cb97dULL, 0x102bb6323de3ddf1ULL } },
  { { 0xb8d74e62eaedb4d8ULL, 0x92d0b11d8e86b73fULL, 0x1acc38a1c97acf8dULL, 0x3db922f554e07b90ULL },
  { 0x332b3896993a49dfULL, 0xdcd515b8347dd0b3ULL, 0x551e3a6f9d541ab8ULL, 0x8c31a1eed6452f54ULL } } },
 { { { 0xa79b3d1478c03566ULL, 0x17b24b9a6032309fULL, 0x8b4378aa58fbea42ULL, 0x4dbd519cc4d2d208ULL },
  { 0x1035f9a5a36e12c6ULL, 0x17704e60ad33e3f2ULL, 0xaf651c4f3c91b93eULL, 0x1beae496d61aa9bbULL } },
  { { 0xc8259e5766d90507ULL, 0x5c0c61bb85173057ULL, 0x57f85d8b3fe5759fULL, 0x85285bee2f8e68bbULL },
  { 0x81a9d77a04781a72ULL, 0xc9ae15ad570446f9ULL, 0x4d246a8a7c2d36d5ULL, 0x442e14f381a193f7ULL } },
  { { 0x9023a826623cf9f3ULL, 0x1c86ba281e27dd29ULL, 0xaefd62ec48d07185ULL, 0x4c5443fad53fafe3ULL },
  { 0x32a69cf2d628e9aeULL, 0xefb49d813bbbe8d1ULL, 0xd56142fe40033fd8ULL, 0x892dcc2c82443eccULL } },
  { { 0x4fc22902ecc3a1cfULL, 0xd12a2072b3546b74ULL, 0xe821037f2e25adccULL, 0x846bc7de4e102a6aULL },
  { 0xf79b071444efd2eaULL, 0x83c8addb8f5b1a00ULL, 0x84aedf83703148a0ULL, 0x6b122d769337e331ULL } } } };
const uintbig torsion_basis_twist_ted_sum_uintbig[3][4][2] = 
{ { { { 0x9ce4b08b107e14b7ULL, 0x29f962cf3b02e247ULL, 0x2b0422e3c4824391ULL, 0x94cecb4ff3dd3018ULL },
  { 0xb97e32b3ce7d5af1ULL, 0x6e4d76674cb8b554ULL, 0x34c1e50c3d76c9aULL, 0x9156368a0cc55a33ULL } },
  { { 0x16527173f66ba2d6ULL, 0xf67d3a16c524cd60ULL, 0xa9f4e311ea382f89ULL, 0x64c3f98973a3acf0ULL },
  { 0x444f7bc022011267ULL, 0xced7a217023e11caULL, 0x195eee068296f925ULL, 0xa13b792f78615f88ULL } },
  { { 0x2c5996fb6e398712ULL, 0x4a0d008b8df6136bULL, 0x7ff7a5a7a9010c66ULL, 0x4fc6b4b4d3e0f327ULL },
  { 0x5f23eb38395da2c4ULL, 0x2b45a524c4c45fe8ULL, 0x2bb05ecd7e3ea7e2ULL, 0x8db97566bba4756dULL } },
  { { 0x33efc447b5fefda9ULL, 0x48f921588cccf7a3ULL, 0x9d228b9e0c99b721ULL, 0x3e281b6e57e869f4ULL },
  { 0x8e051f56dfbc12dcULL, 0x5c84da01d7f494e5ULL, 0x2ffef855f2ebea4aULL, 0x810509176b47bf68ULL } } },
 { { { 0xbf3e2b64269a2031ULL, 0x26bcc59452e5c684ULL, 0x72f583fea1ddd13bULL, 0x5d3088d5fda634bbULL },
  { 0x5d342fc4fb369fdULL, 0x25860dc008662c5eULL, 0x5fbfb7277004973ULL, 0x9f9bd97c7c1b20c7ULL } },
  { { 0xd0231f1f65d6bb4bULL, 0x21e5bdaebbb3dec0ULL, 0x8d95aeb07f7a8333ULL, 0x44e0fb3e28c24d4bULL },
  { 0xc122145e757949e1ULL, 0x34cca5f4392f5099ULL, 0xe5b2fd2fc460803dULL, 0x3404ca0f227849efULL } },
  { { 0x5783e210d6964193ULL, 0x11103c16cac7d5b9ULL, 0x1ab6cf5210c1ab1ULL, 0x7b40a05af43d88f0ULL },
  { 0xadad4460464a0bd3ULL, 0xc065b2c36ed32513ULL, 0x830ba5709f8e2d1cULL, 0x7e833b6b919bcde9ULL } },
  { { 0xad00d9b1141daf92ULL, 0xa710b51507b12eeaULL, 0xbef195537b66c647ULL, 0x15e787e088626338ULL },
  { 0xfee87845eea15b8bULL, 0xaade9abdcf6ceba8ULL, 0xa02ae7933be7023eULL, 0x3224aaf374c8a87dULL } } },
 { { { 0x7ad8e03572e59189ULL, 0x72a7e3b9b00f6f53ULL, 0xe408f1844bbf535aULL, 0x9f3ef02399ca5da8ULL },
  { 0x2e340dceebb1ab91ULL, 0xa753e938bcb24dc1ULL, 0xbbc1b37fd3adbf89ULL, 0x67518760012526dbULL } },
  { { 0x41ad67a8714833e1ULL, 0x57f63169635fd931ULL, 0xa2c87bc305a65c1ULL, 0xca9f04148de18dfULL },
  { 0xb4e84c64a83c6d4cULL, 0x37638a584ccc39b2ULL, 0x187e6a37def39cd6ULL, 0x8ed73d76e15d2ceULL } },
  { { 0x2763e5329f037a1fULL, 0xb1e076b4ff39b5bfULL, 0xbc5022cb99948a7eULL, 0x9453e0f3646b238bULL },
  { 0xddd3b06d18ad968eULL, 0xf97b5c371e1cb456ULL, 0xd412991e99730d38ULL, 0x311f5bfd95dde764ULL } },
  { { 0xd3a5db1f7030a78dULL, 0x78fec307fc46a23bULL, 0x8a8a1dda73fec1d6ULL, 0x3b4a56868952f2afULL },
  { 0x2203e2927d65feefULL, 0x573dc407cb32bea5ULL, 0xdcac337a9b581b59ULL, 0x7ce7b7893b6600b7ULL } } } };

const uintbig torsion_basis_uintbig[12][3][2][2] = {
{ { { { 0xddd499f0af795a04ULL, 0x9d5debb8273280d4ULL, 0x49ba8bf942b0f3aULL, 0x235873646f1d25efULL },
  { 0x3107fa5063330763ULL, 0x9d12c09e79907fcULL, 0x7bab0c8bc5153baeULL, 0x2c511501d5a37825ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xaafe9d58a24ab637ULL, 0x90771182c2f5513cULL, 0x131ee6a606602208ULL, 0x8274b514a5a2efcdULL },
  { 0x2671262f2dd050b4ULL, 0xa5cf7dc5e3cbc51dULL, 0x39d89620299664c3ULL, 0x83a2543fb01a02d6ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd5c5adfb974ff63fULL, 0xef431ea19578eb28ULL, 0xdd45de8d082fde33ULL, 0x233b0024a94cf475ULL },
  { 0xaa44f02d4c9ca78bULL, 0xdeabd637ffd40bbaULL, 0x21e2106c76402b8eULL, 0x5625662aab384bd0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x4d703b3e9cbb983fULL, 0x4f439edce42a4e4ULL, 0x3c09638ae1d8022bULL, 0x64759e571a819eb1ULL },
  { 0xfe2740c1a4285e16ULL, 0xeaf5b434a1418c6eULL, 0x9ce53337dcbc9d11ULL, 0x6c348a0c58460c64ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5e477d17039543e6ULL, 0xc012885890bb2afdULL, 0xbb1899afce33aeb2ULL, 0x827a5d2e64b7c87cULL },
  { 0x64ad4c0854eb31fcULL, 0xcbb6121881d0e5aULL, 0xe4931573df2c34a2ULL, 0x76d080b74d3618bcULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd3e9e914a4989e1ULL, 0x2fc629aaeb32a88dULL, 0x86ba7661a1c5549eULL, 0xd1b3662144d9697ULL },
  { 0xc9960aa389a45410ULL, 0x9cb66db833fb87eULL, 0x32fe4607ca9e395cULL, 0x7b4cb5d499a28946ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xa7e679fb6c578a0ULL, 0xb7d0fec0a1da5705ULL, 0x8ec39d3e5e2f7111ULL, 0x6d63052f634297d6ULL },
  { 0xf2215ac381a6e6eaULL, 0x2cd8794c3df65353ULL, 0x38924fbe362c515bULL, 0xa2d6e473e068a0b3ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xc4174c167999e536ULL, 0x89881bec587622eaULL, 0x2d4e489e05d11ad6ULL, 0x30bc679a697206cbULL },
  { 0xce6f6686e24bfca2ULL, 0xf7ae4434daee8b64ULL, 0xf3963eb4ff445aadULL, 0x82428ad01211efe0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xbfdd0c825e3235ebULL, 0x68925d939fd0179dULL, 0x85e365add93000afULL, 0x719de821cb0258c9ULL },
  { 0x7f60d34b0e371245ULL, 0x1e94739de9d76535ULL, 0x55d804b35255d485ULL, 0x5fd3b55a753587c5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x2db6ae084ec4648dULL, 0x6a69aa2fdb0a6c0ULL, 0x8c9b11b7c9c5dc6cULL, 0x7dcf399353f22cbdULL },
  { 0x64a0a710fcbc77cULL, 0xb1d8f119a5309a8ULL, 0x302eb8c493a27f36ULL, 0x5bd174506f982b0dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa30fcc85c27b36a8ULL, 0x93b2b9b9f9ee718aULL, 0xbec5e4fb52eddf03ULL, 0x53a8d9f0af253d9bULL },
  { 0x8fb9261c8c57c447ULL, 0x601d2b884d596f34ULL, 0x5ceb0704c71fbb5ULL, 0x2ed37c99e83a450dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x62e1f2398af0974eULL, 0xdd7e9b3183c2239ULL, 0x62cec3bdba38fa19ULL, 0x23705837cb849942ULL },
  { 0xfe586453b676faa9ULL, 0x76adb69cd9e8c241ULL, 0xb36716e80752bc3cULL, 0xb8f664242343b18ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xcedd7ca96987a866ULL, 0xe21e2cb5ad4d4032ULL, 0x69e92c08d4c63839ULL, 0x762a9aa178ead5d3ULL },
  { 0x38aa9ce28a61e65eULL, 0x284072d758b4531ULL, 0xbe6c17bc16c9357eULL, 0x777cb73e1aa51029ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x42737cf8d2233628ULL, 0xd05f01079c8aea8aULL, 0xab25158ea0069b84ULL, 0x1ddff98c3d4f3fbcULL },
  { 0xd66ac267a48f7954ULL, 0x1c7d790f01fd0569ULL, 0x3bfe3fee89dc6c17ULL, 0x7a5329779777a0e8ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x16e7dbe079a07e60ULL, 0xd6a1c7b3c939452dULL, 0x1f850b968d0923e7ULL, 0x105a62904c7d9226ULL },
  { 0x38a7baee491a1a7eULL, 0xe446748483e93d04ULL, 0x1d7fd5c5e31e3bb9ULL, 0x4b1ace9ad9abd066ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xd8b4c27d634abcf3ULL, 0x28420186af46a617ULL, 0x96ab2b5b19675141ULL, 0x2e798f130b5b5b2ULL },
  { 0x3e20669bc780a412ULL, 0xf01d379763940856ULL, 0x3c84328055aa91b1ULL, 0x93bf647f68d71d72ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa5fb74f472b338d8ULL, 0x20d30614e3a64af7ULL, 0xb7a8b4c85b559dd3ULL, 0x15f4f515a700f88eULL },
  { 0x8dbcaa7774d97de2ULL, 0x3e8fce8812674333ULL, 0xdd84817c33e25523ULL, 0x418b38f97afcf76fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x293e78eec119e61eULL, 0x10e6db3296d41a64ULL, 0xd443f8e0e7f80a7eULL, 0x1a3a1a88f872f637ULL },
  { 0xe79b3d5774fae63eULL, 0x269397cb2ed4cd4bULL, 0x3046268727e5ea5fULL, 0x4e38e8bb0f63a305ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x530f936e53990dc2ULL, 0x6725d72cdc73a961ULL, 0x1dd7fac24e75afb6ULL, 0x710f16d942885a31ULL },
  { 0xea9f339d49299543ULL, 0xacf1ce826ef83648ULL, 0x145f9382f6bd6493ULL, 0x9eecbb9c305f1567ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x879ef67b89ad539ULL, 0x7370c14b651c981fULL, 0x956fdd37013439c3ULL, 0x1fd044d3b333456fULL },
  { 0xf44651149d9a0ae0ULL, 0xf3440f6da6d3496fULL, 0xa174df1bf8753d96ULL, 0x27b6cf04259a79d6ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd7e8e170f7b80c61ULL, 0xb1efd16cf7a59bd3ULL, 0x8fa46e0c8b0b0ULL, 0x86bcbb083f045d2cULL },
  { 0xdb3dee5c40f4537fULL, 0x5d101d553551698ULL, 0x92f2ead6000456b0ULL, 0x8191c9fa4a6589e1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xa90c0925db4d0a10ULL, 0x1843dcc6b4a67841ULL, 0x43bf4ae20f015da7ULL, 0x506e43fe19fb89b9ULL },
  { 0x3befe5f1c9d313c5ULL, 0x48d2e92fc8616598ULL, 0x11aeb250a537acc1ULL, 0x519fe3ffd12e0e70ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x563a62212174e873ULL, 0x54f71201f946687eULL, 0x6eb73c6a754c9898ULL, 0x9fabb8055e3fdedaULL },
  { 0x3220c9055c951f71ULL, 0xcead44d9f99019a2ULL, 0x99734520d1788008ULL, 0x116f205d93febe90ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x4fa169e9c0be761fULL, 0xe2ba05f95ffe3446ULL, 0x22c373588d86b37cULL, 0x863b8bf7f1ce8d60ULL },
  { 0x655d4210819e0fedULL, 0x9ee15f693920cfc4ULL, 0x49e710641e872c2aULL, 0x70dc287d5fafa49fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xee2d3485ed682c5bULL, 0x12fbaf9f0feefb91ULL, 0xf114a79720c73cbbULL, 0x1368359da4bebcULL },
  { 0x5dff46ca35704a5eULL, 0xbac78ae299ab14f3ULL, 0xf5b46d414ccaac75ULL, 0x2355239bc683800bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5de9d21faa0d8e6eULL, 0xac857217fc1b74d5ULL, 0x760552c5fdae0cccULL, 0x53f8168bf7a19578ULL },
  { 0x547d96877f28746cULL, 0x934e62b72a0dd73ULL, 0xaef81873fb0d2d43ULL, 0x3c721404abb9ebc3ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x77ec097c0d250aa7ULL, 0x1ea15090510c1b04ULL, 0x11273f558470859ULL, 0x59db2024034ba0d0ULL },
  { 0xf8074c874fa9ff4aULL, 0x3bbe056eea96cbb9ULL, 0x378d8d5ac4f17e62ULL, 0xeb03c162d8af414ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xcc6ff14a4948a491ULL, 0x8e04dc1dd750493fULL, 0x45a3715b43b33319ULL, 0x7998aa7f89799222ULL },
  { 0x4efc4ddb87288bafULL, 0x47d8613e480bc9dbULL, 0xc873a608520f6352ULL, 0x11c88a8313737ec8ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x18e8b3bd23342a9aULL, 0x207d52f6830d281cULL, 0x3658cbed341c115cULL, 0x5877cc12965e0923ULL },
  { 0x163074b2f7b4951aULL, 0x1db26566d7966ab4ULL, 0x276a0353ecead7bULL, 0x4e540870f4d77ecaULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x456adb008fd7e712ULL, 0x41fbb86cd6cdb747ULL, 0xfd3f09f721770416ULL, 0x842fc4f2a6e70753ULL },
  { 0x860f02125c96e411ULL, 0x76b7d04d5d0d5effULL, 0x22aea8db1b877855ULL, 0x832654ddefff8ff5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x720947c148a70e4aULL, 0x115a3844360f008eULL, 0x5d7d643f9dab6141ULL, 0x1a21db34009eb7b2ULL },
  { 0x36544c7a46c5ac68ULL, 0x3cf98373afce6e55ULL, 0x96c914ea8cee44f9ULL, 0x810dc0300424d6c6ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7008659d2e3e0c3aULL, 0x7d4fffa6c8177719ULL, 0xef46588c55d06aa5ULL, 0x8219a8fd30126644ULL },
  { 0xdb0ac7678e97d3b3ULL, 0x2755be7983df1718ULL, 0xfd993332f128b7daULL, 0x7cf93a69966cdb62ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7d9ef5f9a2753c15ULL, 0x24c5d55f54cbdee0ULL, 0x4d4fc6a7e91e0b9dULL, 0x293f57f715fc1d7bULL },
  { 0x5a929e7f489abb2eULL, 0x8052202c56621e0cULL, 0x98ae1cc03fc80782ULL, 0x809f9dc1b1f7ca3eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x6677bf56f75674baULL, 0x92a05768fa71171dULL, 0xa1565ae799b8f0eeULL, 0x15feb162b19f2492ULL },
  { 0x38b86ef6e2c8f313ULL, 0xb963e11e2042f114ULL, 0x270b16b7c46ed89cULL, 0xb8ffc6ce464d168ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xcfe8b3a295699794ULL, 0xca117967138f3e1eULL, 0x9202bb57668c7943ULL, 0x5af5455b3faab0d4ULL },
  { 0xf863d574b828ab05ULL, 0x28ab2beeb6141747ULL, 0x559ae6687da21468ULL, 0x509897fdcdce262eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xf10843ef92760c1dULL, 0xcb1881111858b443ULL, 0xca59d90b0c3faa6bULL, 0x4a4a96fc597958ccULL },
  { 0xb3ac1fc319b1e117ULL, 0x6bd0aeabe0f1ff92ULL, 0x517dfd5a3c277549ULL, 0x58cb34df8ea8f7bfULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
};

const uintbig torsion_basis_twist_uintbig[19][3][2][2] = {
{ { { { 0x9f29bc5b16e1ace1ULL, 0xbe1d189c6c9a86c8ULL, 0x975f24135e4ee5baULL, 0x847c6e79e0957225ULL },
  { 0x920df44435910134ULL, 0x14631eda71561da2ULL, 0x48296543ad481a31ULL, 0x2328b170c2f2fee5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x74fede480bc6178bULL, 0x7a4a1c79563a6e1ULL, 0x870ef8a81ca95d93ULL, 0x78c259b3f9a72f1dULL },
  { 0xe26047bf9cd04890ULL, 0xff97ddc873f20ab3ULL, 0x79b397d6af1d6319ULL, 0x315556913e4806b7ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9589aed9f616ad09ULL, 0xa1b80f128cdd185dULL, 0xe9420c4db6aff0c4ULL, 0x72b0f667a5725987ULL },
  { 0x33eea96e2d08b1d6ULL, 0x9e57271b99f6b6abULL, 0xac95983b6b2ddb9dULL, 0x36f2c9b7066e9600ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x38905e38c353aad3ULL, 0x12e01ecc3883e7ecULL, 0x209ca5db6fd0e8baULL, 0x77b512a2dd53962ULL },
  { 0x831a287195983c1bULL, 0xf8577e2951678a97ULL, 0xef7453521cb99f93ULL, 0x195f1e38d696cc93ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xf1457f098736f795ULL, 0x824eb77516bea085ULL, 0x8bc3c47e4810d374ULL, 0x7bb131823754b3c0ULL },
  { 0xa21340f6665287b1ULL, 0xc6b687ac721bd79cULL, 0xa6c0509d413debdfULL, 0x951fdcd710058d7aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x576c7aed61cdc383ULL, 0xda67128e6265573eULL, 0x80c400543896d0f4ULL, 0x378851162fa71951ULL },
  { 0x81b74793f38ad3f6ULL, 0xce5d820b40f28b51ULL, 0x8ca912ba2fbab842ULL, 0x67421ffe7b1b9a73ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x3c056447bf50760dULL, 0xfca85b24a784a2bULL, 0x25ef58692a648bceULL, 0x24cf689ef6444e48ULL },
  { 0x258559ee898e7dffULL, 0x5cb2ee326d715cc3ULL, 0x23827e08b8e98560ULL, 0x6d9f14c1df0c562aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb93bc827bab1a841ULL, 0x6c2b80ab675ba3bfULL, 0xdc34d9744f6741f7ULL, 0x92ba1c51b5415d36ULL },
  { 0x51317b40721b9908ULL, 0xfe292f495e96364fULL, 0xc3701d25b7f362bdULL, 0x691f2c25de413794ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xf15af76b8b2dd5bdULL, 0x2a2c7b87a2abe1b5ULL, 0xf77ebc2b9eeda91cULL, 0x10553e4140a5b1f8ULL },
  { 0x93199bf61ba4df66ULL, 0xe3a974c798f06779ULL, 0xea251979bb8357ccULL, 0x41c9e044e28216dfULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x1f81203aaa37c50fULL, 0x8b1d9a8c2583cb43ULL, 0x2c99f081b0a12452ULL, 0x964fd3abb3b955ccULL },
  { 0x1c551a18c000a1c5ULL, 0xd5e25e2fd9f99265ULL, 0xbdae29cee86e036aULL, 0x18da8235e9640f22ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xf209fa81eba20659ULL, 0xba8b0d2e57a0d4e1ULL, 0x1ec3e673968bc06dULL, 0xb71c7c1b4d96d11ULL },
  { 0xbfa03e6272d908eULL, 0xe38b34e41df0d421ULL, 0xb0688e20add95d76ULL, 0x4392dacbc4121d42ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x1e0eff5146ab38bcULL, 0x2e755d5260baf93dULL, 0x96e2bb08fcc6bd0aULL, 0x81255801e62e4937ULL },
  { 0x52b5fb8b705f4cbaULL, 0x1c405a25d8c775b5ULL, 0xdb39c088ff0a0760ULL, 0x9eb51b1a4a774069ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xc53b0c9be3f37b6cULL, 0x5221a961ad1cf746ULL, 0xaf838e015f9bd307ULL, 0x5576960b83f98e4fULL },
  { 0x1177fa9a6fab013dULL, 0x40ba91583afb12a0ULL, 0xbc5d8612d55f0bebULL, 0x3d109fb16004399aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x1a10356be875c4b5ULL, 0xc94fda2b6bbe6ea3ULL, 0xb9e09eeb85b38a9aULL, 0x1567ea2df8c5b9e4ULL },
  { 0x28f94db89320d302ULL, 0x5a288dad80c65d7dULL, 0x9127d8a4c6a0215cULL, 0x1356d5d54e4608cbULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7441b0e37037cd08ULL, 0x4d19a46d07fcfceULL, 0xa910baab6da00177ULL, 0x8cef7cf626fa4deaULL },
  { 0x6ed0d93bfc9d8074ULL, 0xe0627da277e5dd4cULL, 0x14f5500d07220d72ULL, 0x8148550b0ee23f37ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x33e310d453a31973ULL, 0x5b3ef60a01df752aULL, 0xa6fda7c0d6c4d996ULL, 0x32d965493d62b4bdULL },
  { 0xc8acdde383e4fb67ULL, 0x2abe1af6dbefc0fcULL, 0xb208a23309e571f3ULL, 0x8cf6a46d6d212c2dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb9857340ebb85179ULL, 0x79f383550a0a1cd6ULL, 0xdb7ce03f51f929d6ULL, 0xa1cc1e811ac362a6ULL },
  { 0x9bbaa4d8ac2a5a32ULL, 0xda593dba08dbe2f4ULL, 0xed5f20b8df20f4f8ULL, 0x33c0da35ad8eeb54ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x597574f26a03bbafULL, 0x85dda3ac061c7c0bULL, 0x9930627574c6033bULL, 0x5169f788866afe8cULL },
  { 0x80b54a60276e1f8cULL, 0x9933aec444f3777dULL, 0x140c87bbd061f6f8ULL, 0x9f2155e0f50212ffULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x1095ca15623e46ffULL, 0x37965fbc1caba63aULL, 0x7c3e2aec434b14a7ULL, 0x90de9eb1e0525db5ULL },
  { 0xdf31d5d0ed1a974fULL, 0x84857a5dc47aba24ULL, 0xe60ce08f35e84c8aULL, 0x41530db39d256847ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x1cecd6d5ca165a8fULL, 0xa18c1f28c82eda1ULL, 0xf2ab8c37acce9584ULL, 0x2234b5d71cdaca9dULL },
  { 0x1c6050fc3ff545adULL, 0x2391e1846eec1c41ULL, 0xf08358b0477456ccULL, 0x5d42c71bf191e51dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xd783385803a02902ULL, 0xa24333bfac164c9bULL, 0xe06f3d5eda2b20e7ULL, 0x3f029e9ff75e60d1ULL },
  { 0x1ef365b04bc2f502ULL, 0xf1487c1d43a5f311ULL, 0xc8da7330875d52deULL, 0x4c1cec34367501cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x2333248299ebafd2ULL, 0x9a115f75f79554b8ULL, 0xb165b1ca3b4142b7ULL, 0x66ae28be3607a511ULL },
  { 0xcbe321e5851f2ec7ULL, 0xf56f1aea68fb9ca4ULL, 0xce9162fe0edcc462ULL, 0x52140e6d52aad455ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x68e208ceb45a2613ULL, 0x8c992bc8b4c79c6cULL, 0x36501df567401011ULL, 0x7090aa7c4db81629ULL },
  { 0x83d1bec9d0e383daULL, 0xac4fc372d20fe5f4ULL, 0xf9f692f02a9b8d74ULL, 0x29edfdf28acd9665ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x77bdb9a73ab957d0ULL, 0x491c551bf27b0f2aULL, 0x5bbd6df24ab752a0ULL, 0x3971f4f82e7aa3c7ULL },
  { 0xac31ba65e055b565ULL, 0x97ddcbde0d484cd5ULL, 0x5b8b471862559289ULL, 0x4ab7d9ef40aae9c2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x8f8b656ad2647c55ULL, 0xed23c2417a09cb05ULL, 0xa7990764bda2fab4ULL, 0x48764beb2e7e6136ULL },
  { 0xcb84c39ac16b9a7aULL, 0x435683dc40904e7fULL, 0xb9088605df4d2a00ULL, 0x35f220f060c11649ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xfc8db01278ddb2dfULL, 0x48125226d3419360ULL, 0xe59fc324472fabe2ULL, 0x9c0ed815e3bc92a1ULL },
  { 0xdadabe961bb07b3eULL, 0xcaa137b86994e9c3ULL, 0xd2cb64f23b5ef78bULL, 0xa094886c1c7fc785ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5d7d38161fed43deULL, 0x28235131fb8c99c6ULL, 0x6f3eb5ac122cf985ULL, 0x311ea9c2f98745c4ULL },
  { 0x67456a55cf84d71eULL, 0x1ecfefc3adfce913ULL, 0xe21211556bc4f2c2ULL, 0x60eee2bcdedb881eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x284bb5164ab74ab0ULL, 0x1b137c1a518cf849ULL, 0xd61e6189ef084e03ULL, 0x75ded383e3331351ULL },
  { 0x31852098ec39d358ULL, 0x1f385bf6b2ae7f9bULL, 0x60c4ae9a3144c073ULL, 0x1517c11c7b354a51ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb659f0cbfaedf1dcULL, 0xdeb83c63edf2de25ULL, 0x1c8cc068194b3128ULL, 0x401c6dddcb218a60ULL },
  { 0xa120f7598ebc99dcULL, 0x44468d9957d389e4ULL, 0x40f14da222b86850ULL, 0x19a035c74173cbb3ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x8a6adb3d9ca0f215ULL, 0x6af4423037c47105ULL, 0x7c8e15c495fb4877ULL, 0x581bed27b605c9b9ULL },
  { 0xe3b3fb8853f537bbULL, 0x774bf5c3bccf0366ULL, 0x7a1f6f0ca69cb199ULL, 0x74c8b976c235e282ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x5a62869366df0945ULL, 0x75f7ff62a01c0476ULL, 0x24b2ec09c280a574ULL, 0x8e4a41f682a191dcULL },
  { 0xc833e27efd03d2c6ULL, 0x54c8769b9c781440ULL, 0xe340d8706a9ecfebULL, 0x9fe53f321a392d80ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x3fef3eabace988efULL, 0xd88d315af7fcdcb8ULL, 0x463f99f748e81aeeULL, 0x33c2ae0225b5e8d5ULL },
  { 0xc068145a36da51dbULL, 0xacdfd882e85047c8ULL, 0x1c89b5862cc7f407ULL, 0x4f99fb8e7b1fa337ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xc4a23edd7c1f9123ULL, 0xef685aa0f6e23476ULL, 0x2a48d41ad7ea794fULL, 0x66e8b16c27a2e83eULL },
  { 0x5841da96136fe006ULL, 0x1f85897b10511115ULL, 0xa455cc49268f4824ULL, 0x84acc669cc89a57eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x370e67e28a839773ULL, 0xe7cbad5b3e238697ULL, 0xc00dc48e7588ad38ULL, 0x253acb870435ea95ULL },
  { 0xe983596a8c82ebb3ULL, 0x9973dc9f578a3bceULL, 0x6862f9afb5c7aeecULL, 0x68c9338c7fd14ed5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa9c0793b5bf9d8f0ULL, 0xda161907eb83c92ULL, 0x123e92e56f29b6baULL, 0x84de0497fddbb95fULL },
  { 0xf03e9c232d194bc8ULL, 0xc61e4d73fa98dc4dULL, 0xa16fa227b46a98caULL, 0x9bec5debb81df0f7ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x2e0199dcf6e8a541ULL, 0x658ee32d7865aad6ULL, 0x659b6a416de7eb26ULL, 0x1ca45b4ebc83c580ULL },
  { 0x98d49899f0f5d96fULL, 0x47e88def9f19d03aULL, 0x1177c8b7ba06c52bULL, 0x8ed4f18f437ecfcdULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xbe1f9936c5e1e04eULL, 0x5921f882145cfb7cULL, 0x14891d5f7a91c4caULL, 0x344446c372f2dca8ULL },
  { 0x43f7dcaf57582a0bULL, 0x8838a28805ea70feULL, 0xc8428e546aa81daeULL, 0x6745862be080945aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb3eb3b9e5c78fc2bULL, 0xc88977660babdbfaULL, 0x9ac786da890765b3ULL, 0x596476c808b05452ULL },
  { 0x260bb90926ce3d2eULL, 0xe2e0c3518477287ULL, 0x4a697403bd7302f6ULL, 0x1f527716c8558919ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xda3fae4f32b383b0ULL, 0xa1ed338f44ced8b9ULL, 0x78c1933aa0c23757ULL, 0x44ae3ba3c27c2418ULL },
  { 0x3d85eac82ae35b8ULL, 0xe4347ced09f29a8fULL, 0x1df0fa95982b6e03ULL, 0x6900d08cccdc1158ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x9597333f66660b35ULL, 0x6c1dce9cb25e80a3ULL, 0x691a175c6ba845b6ULL, 0x3128b0dd6862e6a9ULL },
  { 0x9bc98d1a4867aafbULL, 0x70539d834f19e50ULL, 0xb88df08c8f9d0cffULL, 0x67bc0ee961711d89ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x97db5ced0ee9d4d6ULL, 0x473f5c07bb35847cULL, 0xa439021f16681152ULL, 0x9674f4ca59a242d7ULL },
  { 0xb3ec6e5d462fe84eULL, 0x11a73e00b834fb1cULL, 0xfa9171be621471a3ULL, 0x18421cd568039ff5ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6d3d5b6cad4a6f36ULL, 0x8368cb14165fcc87ULL, 0x1f87eed246a9f50cULL, 0x753deba048fcab1eULL },
  { 0xea4664d69eae1519ULL, 0xce242d4a3b3574f1ULL, 0xbda8479499ea17eaULL, 0x76e931a34333ab27ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xe400298cba30a0a1ULL, 0x1f056e6c9b7cdec5ULL, 0xa9eddf0ab400a130ULL, 0x9de25ff7a46db18dULL },
  { 0x5a278db2e83a3342ULL, 0xc1fea6c368d82f46ULL, 0x1b9f53142404c07aULL, 0xa1ea0119d8b5c938ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe42220559c91a0a6ULL, 0x8768166189a72ae4ULL, 0x4dcd843480ddf68fULL, 0x59512db5691e94dULL },
  { 0x144de4fe96d62ed1ULL, 0x7f69d31304a951f2ULL, 0xf4f99a5dee4c9fa3ULL, 0x4bf0f5f276556880ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa97b11be9ab9bc95ULL, 0x3468bce78caf7771ULL, 0x65cd24b7f16cf4fbULL, 0x13cd06df8c106730ULL },
  { 0x9da293dd506c5db7ULL, 0x6d0f8dcfa18555c7ULL, 0xc5d38e9bdcb0cd62ULL, 0x841b6fff4c1b8f69ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x6bbe8a4545d78376ULL, 0xe8a50aa79110fd75ULL, 0x29fe3f2ff71c28e3ULL, 0x52d36c79248669a8ULL },
  { 0xe2d0ccda2c246d34ULL, 0x7ec005b8bd9d9777ULL, 0x45ca4739be416a86ULL, 0x21dba43d5592a9ceULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x38d421aad05e4b39ULL, 0xb719017ee49c0a7dULL, 0xdd3508cce7690039ULL, 0x259025263172bc72ULL },
  { 0x873e305ac36a9ec4ULL, 0x6779158963c4f290ULL, 0x100c42fc04d04b66ULL, 0x66209baf984237d9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6ec42b5b7f924a13ULL, 0xbb514220242a0e90ULL, 0xd38e24a748d9cd45ULL, 0x27045bc1c567a673ULL },
  { 0x8cf9371f461688ccULL, 0xe838a2c9224421a0ULL, 0xdbef830db1c0dd8ULL, 0x6adb6f855c750962ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x21fdc21b1305ad54ULL, 0xdcdbfda7d91c0111ULL, 0x68f0ecc195de45ddULL, 0x8e24c0bbfc4f608dULL },
  { 0xddb1af3aeea8feb9ULL, 0x28cf2e753565871bULL, 0x456d036c3f570f01ULL, 0x6b222de2be7d4a70ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5a9a823d192effacULL, 0x36b7001c5fe5960aULL, 0x36c007cd7701ee53ULL, 0x37efe00d55d3c37aULL },
  { 0x19c3cd4f67994288ULL, 0x6b5a2337aff706cdULL, 0x3feeb30c436a2200ULL, 0x283dcc9a67a1e6a7ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x2395b6a3d93dbf89ULL, 0xf626785bc02927ccULL, 0x34a11eec727f189bULL, 0x7c5c33f65e35f4bbULL },
  { 0x34f1213eb7354d37ULL, 0x846a367c19771c94ULL, 0xd67bd63194415c8bULL, 0x5f5af5b6a7ea547dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x1f582c076563f3e2ULL, 0x6af97e7fd7132689ULL, 0x76884de0ea477044ULL, 0x2b18422a2083619ULL },
  { 0x6c996297e784f496ULL, 0x7fd771d69b2c815dULL, 0x322bae75a7c8dcbULL, 0x627bd2e2c6b25c39ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xdcbd721cbed9ae99ULL, 0x71addd6468f50f23ULL, 0xb2a501a75673edd1ULL, 0x1b33f1aa0820acaULL },
  { 0x6b7381e0e0b9808dULL, 0xbb06c2c2493f703dULL, 0x7c58aab18a0a2a3ULL, 0x5888c69221e2989dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x4163b37e1c089c49ULL, 0xb5a482b8de32b93eULL, 0xfd395654267a3e5fULL, 0x2913e823fa9cebc5ULL },
  { 0xb986eba34fc34cceULL, 0xd502342db7c2cc95ULL, 0x3cf77b8d15ca8ac3ULL, 0x3b7e668b128dd5c2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x28029c71bbaa8938ULL, 0x32c6aafc102f5720ULL, 0x7132273343074c96ULL, 0x30b2650cb448e0d7ULL },
  { 0xe38aa242b2a2bedbULL, 0x8aa2d0bd95d2568cULL, 0xe9aeb32066b16558ULL, 0x1564e4bec35e2719ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7547ce45f17762b6ULL, 0xdceb717f7c5e1dbdULL, 0x6f7239721baaa120ULL, 0x6320ff99b1932514ULL },
  { 0x6a6199eceb3ce7bdULL, 0x9f74c66e65c33e59ULL, 0xf3f1898e3b497ca2ULL, 0x338dca5f6b7e4d3aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x141240aab98e5bdfULL, 0x9c54ae3a355dc258ULL, 0x174e82d43144c4b6ULL, 0x9f6677fa16711e31ULL },
  { 0x1a28a97733f2bddaULL, 0x998565eb31f69ef8ULL, 0xf450cabd74663dbbULL, 0x6db3f8f5988bd6e8ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
};

const uintbig torsion_basis_two_uintbig[3][2][2] = 
{ { { { 0xf1fe96cc4eaf0775ULL, 0x14057ead15b9f0d8ULL, 0xa508db1478e635a8ULL, 0x9ee509352cbf2129ULL },
  { 0x5f78295f053ecdc2ULL, 0x659fdd3c830ea621ULL, 0x6141a8f51c191377ULL, 0x8e47b35c9916244fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x319c1e00de60157fULL, 0x841bddb428cf525fULL, 0xd7918cf202e2a0bULL, 0x6e4ab618eaddbf74ULL },
  { 0x4e08e27e7c3e665cULL, 0x37b564fd8333a5f4ULL, 0x228629f5691a695fULL, 0x9e390be952811cbdULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x2db272228241144dULL, 0xb8ab5e9985272d95ULL, 0x6f57fd4671aab293ULL, 0x92f3f8dd93584662ULL },
  { 0xf37819f49afefccdULL, 0xca0c5c3b735fcacfULL, 0xc33114ce1f0ccec6ULL, 0x4e815afd0fceb9a1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };


proj torsion_basis[12][3];
proj torsion_basis_sum[3];
point torsion_basis_ted_sum[3];
proj torsion_basis_twist[19][3];
proj torsion_basis_twist_sum[3];
point torsion_basis_twist_ted_sum[3];
proj torsion_basis_two[3];


void init_precomputations_generated() {
	global_setup.action_2 = malloc(12*sizeof(GEN));
	global_setup.action_3 = malloc(12*sizeof(GEN));
	global_setup.action_4 = malloc(12*sizeof(GEN));
	global_setup.action_twist_2 = malloc(19*sizeof(GEN));
	global_setup.action_twist_3 = malloc(19*sizeof(GEN));
	global_setup.action_twist_4 = malloc(19*sizeof(GEN));

	global_setup.action_2[0] = mkmat2(mkcol2(stoi(278938600389568ULL),stoi(476702643838875ULL)),mkcol2(stoi(315443836407367ULL),stoi(197898557813557ULL)));
	global_setup.action_3[0] = mkmat2(mkcol2(stoi(107003843221062ULL),stoi(153690750870483ULL)),mkcol2(stoi(43280178462046ULL),stoi(369833314982064ULL)));
	global_setup.action_4[0] = mkmat2(mkcol2(stoi(30297936706966ULL),stoi(87634623728094ULL)),mkcol2(stoi(385138224183876ULL),stoi(446539221496159ULL)));
	global_setup.action_2[1] = mkmat2(mkcol2(stoi(23ULL),stoi(32ULL)),mkcol2(stoi(11ULL),stoi(26ULL)));
	global_setup.action_3[1] = mkmat2(mkcol2(stoi(24ULL),stoi(38ULL)),mkcol2(stoi(19ULL),stoi(26ULL)));
	global_setup.action_4[1] = mkmat2(mkcol2(stoi(33ULL),stoi(40ULL)),mkcol2(stoi(23ULL),stoi(16ULL)));
	global_setup.action_2[2] = mkmat2(mkcol2(stoi(2ULL),stoi(2ULL)),mkcol2(stoi(3ULL),stoi(9ULL)));
	global_setup.action_3[2] = mkmat2(mkcol2(stoi(6ULL),stoi(8ULL)),mkcol2(stoi(10ULL),stoi(6ULL)));
	global_setup.action_4[2] = mkmat2(mkcol2(stoi(10ULL),stoi(6ULL)),mkcol2(stoi(9ULL),stoi(1ULL)));
	global_setup.action_2[3] = mkmat2(mkcol2(stoi(17ULL),stoi(19ULL)),mkcol2(stoi(19ULL),stoi(14ULL)));
	global_setup.action_3[3] = mkmat2(mkcol2(stoi(21ULL),stoi(2ULL)),mkcol2(stoi(7ULL),stoi(11ULL)));
	global_setup.action_4[3] = mkmat2(mkcol2(stoi(25ULL),stoi(26ULL)),mkcol2(stoi(1ULL),stoi(6ULL)));
	global_setup.action_2[4] = mkmat2(mkcol2(stoi(78ULL),stoi(2ULL)),mkcol2(stoi(70ULL),stoi(5ULL)));
	global_setup.action_3[4] = mkmat2(mkcol2(stoi(16ULL),stoi(36ULL)),mkcol2(stoi(21ULL),stoi(68ULL)));
	global_setup.action_4[4] = mkmat2(mkcol2(stoi(45ULL),stoi(39ULL)),mkcol2(stoi(63ULL),stoi(38ULL)));
	global_setup.action_2[5] = mkmat2(mkcol2(stoi(6ULL),stoi(47ULL)),mkcol2(stoi(88ULL),stoi(101ULL)));
	global_setup.action_3[5] = mkmat2(mkcol2(stoi(84ULL),stoi(24ULL)),mkcol2(stoi(84ULL),stoi(24ULL)));
	global_setup.action_4[5] = mkmat2(mkcol2(stoi(65ULL),stoi(95ULL)),mkcol2(stoi(40ULL),stoi(42ULL)));
	global_setup.action_2[6] = mkmat2(mkcol2(stoi(50ULL),stoi(135ULL)),mkcol2(stoi(86ULL),stoi(87ULL)));
	global_setup.action_3[6] = mkmat2(mkcol2(stoi(124ULL),stoi(8ULL)),mkcol2(stoi(80ULL),stoi(14ULL)));
	global_setup.action_4[6] = mkmat2(mkcol2(stoi(12ULL),stoi(98ULL)),mkcol2(stoi(88ULL),stoi(125ULL)));
	global_setup.action_2[7] = mkmat2(mkcol2(stoi(342ULL),stoi(547ULL)),mkcol2(stoi(150ULL),stoi(409ULL)));
	global_setup.action_3[7] = mkmat2(mkcol2(stoi(468ULL),stoi(433ULL)),mkcol2(stoi(163ULL),stoi(284ULL)));
	global_setup.action_4[7] = mkmat2(mkcol2(stoi(636ULL),stoi(30ULL)),mkcol2(stoi(185ULL),stoi(115ULL)));
	global_setup.action_2[8] = mkmat2(mkcol2(stoi(187ULL),stoi(164ULL)),mkcol2(stoi(523ULL),stoi(640ULL)));
	global_setup.action_3[8] = mkmat2(mkcol2(stoi(128ULL),stoi(500ULL)),mkcol2(stoi(510ULL),stoi(700ULL)));
	global_setup.action_4[8] = mkmat2(mkcol2(stoi(66ULL),stoi(723ULL)),mkcol2(stoi(519ULL),stoi(761ULL)));
	global_setup.action_2[9] = mkmat2(mkcol2(stoi(1215ULL),stoi(1233ULL)),mkcol2(stoi(988ULL),stoi(2476ULL)));
	global_setup.action_3[9] = mkmat2(mkcol2(stoi(687ULL),stoi(1741ULL)),mkcol2(stoi(2320ULL),stoi(3005ULL)));
	global_setup.action_4[9] = mkmat2(mkcol2(stoi(574ULL),stoi(3464ULL)),mkcol2(stoi(736ULL),stoi(3117ULL)));
	global_setup.action_2[10] = mkmat2(mkcol2(stoi(1354ULL),stoi(1402ULL)),mkcol2(stoi(54ULL),stoi(2665ULL)));
	global_setup.action_3[10] = mkmat2(mkcol2(stoi(3578ULL),stoi(3310ULL)),mkcol2(stoi(1505ULL),stoi(442ULL)));
	global_setup.action_4[10] = mkmat2(mkcol2(stoi(1752ULL),stoi(1313ULL)),mkcol2(stoi(163ULL),stoi(2267ULL)));
	global_setup.action_2[11] = mkmat2(mkcol2(stoi(4126ULL),stoi(6019ULL)),mkcol2(stoi(5874ULL),stoi(2857ULL)));
	global_setup.action_3[11] = mkmat2(mkcol2(stoi(5217ULL),stoi(311ULL)),mkcol2(stoi(160ULL),stoi(1767ULL)));
	global_setup.action_4[11] = mkmat2(mkcol2(stoi(3122ULL),stoi(5761ULL)),mkcol2(stoi(6479ULL),stoi(3861ULL)));
	global_setup.action_twist_2[0] = mkmat2(mkcol2(strtoi("9270895902870275720723784"),strtoi("17087778891795898109246687")),mkcol2(strtoi("16857494539085976434494894"),strtoi("10112349764809744176072939")));
	global_setup.action_twist_3[0] = mkmat2(mkcol2(strtoi("11028619819681929979714893"),strtoi("14765541424117899392974006")),mkcol2(strtoi("11167813334852078556352084"),strtoi("8354625847998089917081831")));
	global_setup.action_twist_4[0] = mkmat2(mkcol2(strtoi("9686131711198419447643943"),strtoi("15367862023993177766253752")),mkcol2(strtoi("12338580278705852834095995"),strtoi("9697113956481600449152780")));
	global_setup.action_twist_2[1] = mkmat2(mkcol2(stoi(34ULL),stoi(40ULL)),mkcol2(stoi(13ULL),stoi(9ULL)));
	global_setup.action_twist_3[1] = mkmat2(mkcol2(stoi(13ULL),stoi(34ULL)),mkcol2(stoi(15ULL),stoi(31ULL)));
	global_setup.action_twist_4[1] = mkmat2(mkcol2(stoi(10ULL),stoi(31ULL)),mkcol2(stoi(3ULL),stoi(33ULL)));
	global_setup.action_twist_2[2] = mkmat2(mkcol2(stoi(9516ULL),stoi(2858ULL)),mkcol2(stoi(8773ULL),stoi(1093ULL)));
	global_setup.action_twist_3[2] = mkmat2(mkcol2(stoi(9597ULL),stoi(10289ULL)),mkcol2(stoi(4911ULL),stoi(1013ULL)));
	global_setup.action_twist_4[2] = mkmat2(mkcol2(stoi(2711ULL),stoi(9169ULL)),mkcol2(stoi(1026ULL),stoi(7898ULL)));
	global_setup.action_twist_2[3] = mkmat2(mkcol2(stoi(25ULL),stoi(72ULL)),mkcol2(stoi(67ULL),stoi(84ULL)));
	global_setup.action_twist_3[3] = mkmat2(mkcol2(stoi(55ULL),stoi(23ULL)),mkcol2(stoi(77ULL),stoi(55ULL)));
	global_setup.action_twist_4[3] = mkmat2(mkcol2(stoi(52ULL),stoi(66ULL)),mkcol2(stoi(16ULL),stoi(57ULL)));
	global_setup.action_twist_2[4] = mkmat2(mkcol2(stoi(170ULL),stoi(73ULL)),mkcol2(stoi(13ULL),stoi(29ULL)));
	global_setup.action_twist_3[4] = mkmat2(mkcol2(stoi(126ULL),stoi(196ULL)),mkcol2(stoi(43ULL),stoi(74ULL)));
	global_setup.action_twist_4[4] = mkmat2(mkcol2(stoi(82ULL),stoi(116ULL)),mkcol2(stoi(99ULL),stoi(117ULL)));
	global_setup.action_twist_2[5] = mkmat2(mkcol2(stoi(122ULL),stoi(99ULL)),mkcol2(stoi(189ULL),stoi(105ULL)));
	global_setup.action_twist_3[5] = mkmat2(mkcol2(stoi(150ULL),stoi(190ULL)),mkcol2(stoi(147ULL),stoi(78ULL)));
	global_setup.action_twist_4[5] = mkmat2(mkcol2(stoi(165ULL),stoi(30ULL)),mkcol2(stoi(201ULL),stoi(62ULL)));
	global_setup.action_twist_2[6] = mkmat2(mkcol2(stoi(389ULL),stoi(179ULL)),mkcol2(stoi(215ULL),stoi(30ULL)));
	global_setup.action_twist_3[6] = mkmat2(mkcol2(stoi(104ULL),stoi(348ULL)),mkcol2(stoi(24ULL),stoi(316ULL)));
	global_setup.action_twist_4[6] = mkmat2(mkcol2(stoi(338ULL),stoi(34ULL)),mkcol2(stoi(35ULL),stoi(81ULL)));
	global_setup.action_twist_2[7] = mkmat2(mkcol2(stoi(366ULL),stoi(162ULL)),mkcol2(stoi(243ULL),stoi(125ULL)));
	global_setup.action_twist_3[7] = mkmat2(mkcol2(stoi(167ULL),stoi(34ULL)),mkcol2(stoi(15ULL),stoi(325ULL)));
	global_setup.action_twist_4[7] = mkmat2(mkcol2(stoi(213ULL),stoi(282ULL)),mkcol2(stoi(230ULL),stoi(278ULL)));
	global_setup.action_twist_2[8] = mkmat2(mkcol2(stoi(559ULL),stoi(349ULL)),mkcol2(stoi(391ULL),stoi(10ULL)));
	global_setup.action_twist_3[8] = mkmat2(mkcol2(stoi(453ULL),stoi(471ULL)),mkcol2(stoi(333ULL),stoi(117ULL)));
	global_setup.action_twist_4[8] = mkmat2(mkcol2(stoi(163ULL),stoi(276ULL)),mkcol2(stoi(80ULL),stoi(406ULL)));
	global_setup.action_twist_2[9] = mkmat2(mkcol2(stoi(401ULL),stoi(31ULL)),mkcol2(stoi(329ULL),stoi(230ULL)));
	global_setup.action_twist_3[9] = mkmat2(mkcol2(stoi(275ULL),stoi(605ULL)),mkcol2(stoi(95ULL),stoi(357ULL)));
	global_setup.action_twist_4[9] = mkmat2(mkcol2(stoi(271ULL),stoi(10ULL)),mkcol2(stoi(7ULL),stoi(360ULL)));
	global_setup.action_twist_2[10] = mkmat2(mkcol2(stoi(123ULL),stoi(513ULL)),mkcol2(stoi(84ULL),stoi(554ULL)));
	global_setup.action_twist_3[10] = mkmat2(mkcol2(stoi(398ULL),stoi(213ULL)),mkcol2(stoi(404ULL),stoi(280ULL)));
	global_setup.action_twist_4[10] = mkmat2(mkcol2(stoi(300ULL),stoi(589ULL)),mkcol2(stoi(665ULL),stoi(377ULL)));
	global_setup.action_twist_2[11] = mkmat2(mkcol2(stoi(505ULL),stoi(125ULL)),mkcol2(stoi(140ULL),stoi(352ULL)));
	global_setup.action_twist_3[11] = mkmat2(mkcol2(stoi(286ULL),stoi(579ULL)),mkcol2(stoi(267ULL),stoi(572ULL)));
	global_setup.action_twist_4[11] = mkmat2(mkcol2(stoi(406ULL),stoi(527ULL)),mkcol2(stoi(332ULL),stoi(451ULL)));
	global_setup.action_twist_2[12] = mkmat2(mkcol2(stoi(453ULL),stoi(743ULL)),mkcol2(stoi(836ULL),stoi(406ULL)));
	global_setup.action_twist_3[12] = mkmat2(mkcol2(stoi(12ULL),stoi(501ULL)),mkcol2(stoi(731ULL),stoi(848ULL)));
	global_setup.action_twist_4[12] = mkmat2(mkcol2(stoi(527ULL),stoi(594ULL)),mkcol2(stoi(155ULL),stoi(332ULL)));
	global_setup.action_twist_2[13] = mkmat2(mkcol2(stoi(846ULL),stoi(33ULL)),mkcol2(stoi(12ULL),stoi(37ULL)));
	global_setup.action_twist_3[13] = mkmat2(mkcol2(stoi(621ULL),stoi(844ULL)),mkcol2(stoi(148ULL),stoi(263ULL)));
	global_setup.action_twist_4[13] = mkmat2(mkcol2(stoi(450ULL),stoi(409ULL)),mkcol2(stoi(566ULL),stoi(433ULL)));
	global_setup.action_twist_2[14] = mkmat2(mkcol2(stoi(119ULL),stoi(868ULL)),mkcol2(stoi(249ULL),stoi(900ULL)));
	global_setup.action_twist_3[14] = mkmat2(mkcol2(stoi(337ULL),stoi(903ULL)),mkcol2(stoi(717ULL),stoi(683ULL)));
	global_setup.action_twist_4[14] = mkmat2(mkcol2(stoi(109ULL),stoi(248ULL)),mkcol2(stoi(628ULL),stoi(910ULL)));
	global_setup.action_twist_2[15] = mkmat2(mkcol2(stoi(140ULL),stoi(651ULL)),mkcol2(stoi(103ULL),stoi(1031ULL)));
	global_setup.action_twist_3[15] = mkmat2(mkcol2(stoi(837ULL),stoi(501ULL)),mkcol2(stoi(70ULL),stoi(335ULL)));
	global_setup.action_twist_4[15] = mkmat2(mkcol2(stoi(1152ULL),stoi(159ULL)),mkcol2(stoi(296ULL),stoi(19ULL)));
	global_setup.action_twist_2[16] = mkmat2(mkcol2(stoi(351ULL),stoi(875ULL)),mkcol2(stoi(1444ULL),stoi(1528ULL)));
	global_setup.action_twist_3[16] = mkmat2(mkcol2(stoi(638ULL),stoi(1456ULL)),mkcol2(stoi(1576ULL),stoi(1242ULL)));
	global_setup.action_twist_4[16] = mkmat2(mkcol2(stoi(151ULL),stoi(656ULL)),mkcol2(stoi(1691ULL),stoi(1728ULL)));
	global_setup.action_twist_2[17] = mkmat2(mkcol2(stoi(1280ULL),stoi(146ULL)),mkcol2(stoi(2213ULL),stoi(1433ULL)));
	global_setup.action_twist_3[17] = mkmat2(mkcol2(stoi(1333ULL),stoi(1665ULL)),mkcol2(stoi(1403ULL),stoi(1381ULL)));
	global_setup.action_twist_4[17] = mkmat2(mkcol2(stoi(1126ULL),stoi(2359ULL)),mkcol2(stoi(1064ULL),stoi(1587ULL)));
	global_setup.action_twist_2[18] = mkmat2(mkcol2(stoi(125ULL),stoi(2227ULL)),mkcol2(stoi(1295ULL),stoi(4158ULL)));
	global_setup.action_twist_3[18] = mkmat2(mkcol2(stoi(4179ULL),stoi(4239ULL)),mkcol2(stoi(2049ULL),stoi(105ULL)));
	global_setup.action_twist_4[18] = mkmat2(mkcol2(stoi(1577ULL),stoi(1336ULL)),mkcol2(stoi(3231ULL),stoi(2706ULL)));
	global_setup.action_two_2 = mkmat2(mkcol2(stoi(7058267855ULL),stoi(6646714874ULL)),mkcol2(stoi(938481723ULL),stoi(1531666737ULL)));
	global_setup.action_two_3 = mkmat2(mkcol2(stoi(7058267855ULL),stoi(6646714874ULL)),mkcol2(stoi(938481723ULL),stoi(1531666737ULL)));
	global_setup.action_two_4 = mkmat2(mkcol2(stoi(7058267855ULL),stoi(6646714874ULL)),mkcol2(stoi(938481723ULL),stoi(1531666737ULL)));
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
	for (int j=0;j<12;j++){
		for (int i = 0; i < 3; ++i) {
			fp_enc( &(&(&torsion_basis[j][i])->x)->re, &(torsion_basis_uintbig[j][i][0][0]) );
			fp_enc( &(&(&torsion_basis[j][i])->x)->im, &(torsion_basis_uintbig[j][i][0][1]) );
			fp_enc( &(&(&torsion_basis[j][i])->z)->re, &(torsion_basis_uintbig[j][i][1][0]) );
			fp_enc( &(&(&torsion_basis[j][i])->z)->im, &(torsion_basis_uintbig[j][i][1][1]) );
		}
	}
	for (int j=0;j<19;j++){
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

