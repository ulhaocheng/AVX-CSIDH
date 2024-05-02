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
// radix-29 parameters and constants for High-Throughput (HT) implementations

#define HT_BRADIX 29                    // limb size
#define HT_NWORDS 18                    // limb number 
#define HT_BMASK  0x1FFFFFFFUL          // 2^29 - 1
#define HT_MONTW 0x32E294D              // the constant for Montgomery Multiplication

typedef __m512i htfe_t[HT_NWORDS];      // HT field element

// the prime p of the field
static const uint64_t ht_p[HT_NWORDS] = {
  0x13C6C87B, 0x1C0DC829, 0x0B2A0D46, 0x0437E8AF, 0x14F25C27, 0x18660F85,
  0x141D459C, 0x18ACFE6A, 0x0DA7AAC6, 0x1499164E, 0x16BEFF31, 0x1B911884,
  0x02D083AE, 0x1F26255A, 0x0AC34578, 0x1137FF91, 0x0E8F740F, 0x00032DA4, };

// p * 2
static const uint32_t ht_pmul2[HT_NWORDS] = {
  0x078D90F6, 0x181B9053, 0x16541A8D, 0x086FD15E, 0x09E4B84E, 0x10CC1F0B,
  0x083A8B39, 0x1159FCD5, 0x1B4F558D, 0x09322C9C, 0x0D7DFE63, 0x17223109, 
  0x05A1075D, 0x1E4C4AB4, 0x15868AF1, 0x026FFF22, 0x1D1EE81F, 0x00065B48, };

// R mod p = 2^522 mod p
static const uint32_t ht_montR[HT_NWORDS] = {
  0x0BF7E1D5, 0x1944150E, 0x1DB05986, 0x0932B2DD, 0x044E5A15, 0x049DBF94, 
  0x055640F7, 0x1A92F0A2, 0x0B31E516, 0x06F67486, 0x07591D44, 0x00683014, 
  0x0B026CC6, 0x11020023, 0x0851A93B, 0x0B4C59FC, 0x0DF0AF96, 0x00018B87, };

// R^2 mod p = (2^522)^2 mod p 
static const uint32_t ht_montR2[HT_NWORDS] = {
  0x1C8CEBF0, 0x1B864D0A, 0x124DBAAB, 0x136BE0BA, 0x04ADD463, 0x1375FBF5,
  0x1BF7A9C1, 0x033D3CE8, 0x03AA0C3B, 0x118D310F, 0x1CE83173, 0x1BE13007,
  0x1AD114C5, 0x188D441E, 0x10FCDC0B, 0x1D9ABEBD, 0x0BE69F94, 0x0001A7D4, };

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

uint8_t u32_iseql(const uint32_t a, const uint32_t b);

// -----------------------------------------------------------------------------
// 1-way prime-field operations from [CCC+19] 

#define NWORDS 8

typedef uint64_t fp[NWORDS] __attribute__((aligned(64)));

extern const fp p;
extern const fp R_mod_p;
extern const fp R_squared_mod_p;	
extern const fp p_minus_1_halves;	

void fp_cswap(fp x, fp y, uint8_t c);
void fp_add(fp c, const fp a, const fp b);
void fp_sub(fp c, const fp a, const fp b);
void fp_mul(fp c, const fp a, const fp b);
void fp_sqr(fp b, const fp a);
void fp_inv(fp x);
uint8_t fp_issquare(fp const x);
void fp_random(fp x);

#define set_zero(x, NUM)	memset(x, 0, sizeof(uint64_t) * NUM);

#define set_one(x, NUM) {\
  int i;\
  x[0] = 1;   \
  for (i=1; i < NUM; i++)\
      x[i] = 0;\
}

#define copy(x, y, NUM)\
  memcpy(x, y, sizeof(uint64_t)*NUM);

int compare(uint64_t *x, uint64_t *y, int NUM);
int iszero(uint64_t *x, int NUM);
uint32_t isequal(uint32_t a, uint32_t b);
int32_t issmaller(int32_t x, int32_t y);

#endif
