/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#ifndef _GFPARIT_H
#define _GFPARIT_H

#include <stdint.h>
#include "intrin.h"

// -----------------------------------------------------------------------------
// radix-52 parameters and constants for High-Throughput (HT) implementations

#define HT_BRADIX 52                    // limb size 
#define HT_NWORDS 10                    // limb number
#define HT_BMASK  0xFFFFFFFFFFFFFULL    // 2^52 - 1
#define HT_MONTW  0x1301F632E294D       // the constant for Montgomery Multiplication

typedef __m512i htfe_t[HT_NWORDS];      // HT field element

// the prime p of the field
static const uint64_t ht_p[HT_NWORDS] = {
  0x1B90533C6C87B, 0xF457ACA8351B8, 0xF0B4F25C2721B, 0x5507516730CC1,
  0xDA7AAC6C567F3, 0xFBFCC69322C9C, 0x83AEDC88C425A, 0x5E3E4C4AB42D0,
  0xF89BFFC8AB0D1, 0x0065B48E8F740, };

// p * 2
static const uint64_t ht_pmul2[HT_NWORDS] = {
  0x3720A678D90F6, 0xE8AF59506A370, 0xE169E4B84E437, 0xAA0EA2CE61983,
  0xB4F558D8ACFE6, 0xF7F98D2645939, 0x075DB911884B5, 0xBC7C9895685A1,
  0xF137FF91561A2, 0x00CB691D1EE81, };
 
// R mod p = 2^520 mod p
static const uint64_t ht_montR[HT_NWORDS] = {
  0xA8EE9BFEFAA94, 0x5371A8DA66CDA, 0x78CE502D8F1AD, 0x199738693E81E,
  0x63663F7667FDE, 0x181C75DC7C56A, 0xBC1D37F29131E, 0xEB481412BEB74,
  0x97908B31B314E, 0x0025C95F2008E, };  

// R^2 mod p = (2^520)^2 mod p 
static const uint64_t ht_montR2[HT_NWORDS] = {
  0x70C9A15C8CEBF, 0x5F05D4936EAAF, 0xF7EA4ADD4639B, 0x746FDEA7066EB,
  0xE3AA0C3B19E9E, 0x3A0C5CE31A621, 0x114C5DF09803F, 0x702F11A883DAD,
  0x94ECD5F5EC3F3, 0x00034FA8BE69F, };

// -----------------------------------------------------------------------------
// radix-43 parameters and constants for Low-Latency (LL) implementations

#define LL_BRADIX 43                    // limb size
#define LL_NWORDS 12                    // limb number 
#define LL_VLIMBS 3                     // limb vector number (NWORDS/4)
#define LL_BMASK  0x7FFFFFFFFFFULL      // 2^43 - 1
#define LL_MONTW  0x1F632E294D          // the constant for Montgomery Multiplication

typedef __m512i llfe_t[LL_VLIMBS];      // LL field element

// the prime p of the field 
static const uint64_t ll_p[LL_NWORDS] = {
  0x10533C6C87B, 0x59506A37037, 0x709C86FD15E, 0x0660F85A792, 
  0x73550751673, 0x34F558D8ACF, 0x731A4C8B273, 0x6446212D7DF, 
  0x2B42D083AED, 0x61A2BC7C989, 0x03E26FFF22A, 0x032DA4747BA, };

// p * 2
static const uint64_t ll_pmul2[LL_NWORDS] = {
  0x20A678D90F6, 0x32A0D46E06E, 0x61390DFA2BD, 0x0CC1F0B4F25, 
  0x66AA0EA2CE6, 0x69EAB1B159F, 0x663499164E6, 0x488C425AFBF, 
  0x5685A1075DB, 0x434578F9312, 0x07C4DFFE455, 0x065B48E8F74, };

// R mod p = 2^516 mod p
static const uint64_t ll_montR[LL_NWORDS] = {
  0x72FE8F0ACC8, 0x0B6F6767762, 0x678AE874934, 0x00D931DD10C, 
  0x7AB6DB47E06, 0x39AA1E24F83, 0x03E40A41DF7, 0x550AD0E4504, 
  0x3D8F6B6CAD8, 0x3E928C8828A, 0x649E8022951, 0x00DE4DCCAEE, };

// R^2 mod p = (2^516)^2 mod p
static const uint64_t ll_montR2[LL_NWORDS] = {
  0x529F05814DE, 0x50C846940C7, 0x290B698A889, 0x0F91058C999, 
  0x50EDA6B9700, 0x30118B01CBE, 0x326643EEEC2, 0x4B7B05A78F3, 
  0x7377415838D, 0x5864018C56A, 0x633D0C02F46, 0x016D9B5D997, }; 

// -----------------------------------------------------------------------------
// radix-64 constants 

// prime p in radix-64
static const uint64_t u64_p[8] = {
  0x1b81b90533c6c87B, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
  0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf, };

// (p-2) in radix-64
static const uint64_t u64_psub2[8] = {
  0x1b81b90533c6c879, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
  0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf, };

// (p-1)/2 in radix-64
static const uint64_t u64_pdiv2[8] = {
  0x8dc0dc8299e3643d, 0xe1390dfa2bd6541a, 0xa8b398660f85a792, 0xd3d56362b3f9aa83,
  0x2d7dfe63499164e6, 0x5a16841d76e44621, 0xfe455868af1f2625, 0x32da4747ba07c4df, };

// -----------------------------------------------------------------------------
// (8x1)-way prime-field operations

void gfp_add_8x1w(htfe_t r, const htfe_t a, const htfe_t b);
void gfp_sub_8x1w(htfe_t r, const htfe_t a, const htfe_t b);
void gfp_mul_8x1w(htfe_t r, const htfe_t a, const htfe_t b);
void gfp_sqr_8x1w(htfe_t r, const htfe_t a);
void gfp_pow_8x1w(htfe_t r, const htfe_t a, const uint64_t *e);
void gfp_inv_8x1w(htfe_t r, const htfe_t a);
void gfp_rdcp_8x1w(htfe_t r, const htfe_t a);
void gfp_zero_8x1w(htfe_t r);
void gfp_num2mont_8x1w(htfe_t r, const htfe_t a);
void gfp_mont2num_8x1w(htfe_t r, const htfe_t a);
void gfp_copy_8x1w(htfe_t r, const htfe_t a);
void gfp_cmove_8x1w(htfe_t r, htfe_t a, const __m512i b);
void gfp_cswap_8x1w(htfe_t r, htfe_t a, const __m512i b);
__m512i gfp_issqr_8x1w(const htfe_t a);
__m512i gfp_iszero_8x1w(const htfe_t a);

// -----------------------------------------------------------------------------
// (2x4)-way prime-field operations

void gfp_addsubi_2x4w(llfe_t rs, const llfe_t ac, const llfe_t bd);
void gfp_subaddi_2x4w(llfe_t rs, const llfe_t ac, const llfe_t bd);
void gfp_addsubc_2x4w(llfe_t rs, const llfe_t ac, const llfe_t bd);
void gfp_subaddc_2x4w(llfe_t rs, const llfe_t ac, const llfe_t bd);
void gfp_mul_2x4w(llfe_t r, const llfe_t a, const llfe_t b);
void gfp_sqr_2x4w(llfe_t r, const llfe_t a);
void gfp_pow_2x4w(llfe_t r, const llfe_t a, const uint64_t *e);
void gfp_inv_2x4w(llfe_t r, const llfe_t a);
void gfp_rdcp_2x4w(llfe_t r, const llfe_t a);
void gfp_carryp_2x4w(llfe_t r);
void gfp_zero_2x4w(llfe_t r);
void gfp_num2mont_2x4w(llfe_t r, const llfe_t a);
void gfp_mont2num_2x4w(llfe_t r, const llfe_t a);
void gfp_copy_2x4w(llfe_t r, const llfe_t a);
void gfp_cmove_2x4w(llfe_t r, llfe_t a, const uint8_t b);
void gfp_cswap_2x4w(llfe_t r, llfe_t a, const uint8_t b);
uint8_t gfp_issqr_2x4w(const llfe_t a);
uint8_t gfp_iszero_2x4w(const llfe_t a);
void vec_permlh_2x4w(llfe_t r, const llfe_t a);
void vec_permzl_2x4w(llfe_t r, const llfe_t a);
void vec_permzh_2x4w(llfe_t r, const llfe_t a);
void vec_permll_2x4w(llfe_t r, const llfe_t a);
void vec_permhh_2x4w(llfe_t r, const llfe_t a);
void vec_blend_2x4w(llfe_t r, const llfe_t a, const llfe_t b, const uint8_t mask);

// -----------------------------------------------------------------------------
// other prime-field operations 

uint8_t u32_iseql(const uint32_t a, const uint32_t b);
void mpi43_carryp(uint64_t *a);

#endif
