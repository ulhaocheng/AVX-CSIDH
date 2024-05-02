/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#ifndef _TEDCURVE_H
#define _TEDCURVE_H

#include "gfparith.h"
#include "rng.h"
#include "utils.h"

// points for HT implementations 

// the projective point on twisted Edwards curve with y and z coordinates 
typedef struct projective_point {
  htfe_t y;  // projective y coordinate
  htfe_t z;  // projective z coordinate
} htpoint;

typedef htpoint htpoint_t[1]; 

// -----------------------------------------------------------------------------
// points for LL implementation 

// llpoint <y-coordinate | z-coordinate>
typedef llfe_t llpoint_t;

// -----------------------------------------------------------------------------
// constants for curve and isogeny arithmetic 

// the number of small primes l_i 
#define N 74  

#define HLMAX 294                       // (587+1)/2 = 294  

// small primes l_i
static int primeli[N] = { 
  349, 347, 337, 331, 317, 313, 311, 307, 293, 283, 281, 277, 271, 
  269, 263, 257, 251, 241, 239, 233, 229, 227, 223, 211, 199, 197, 
  193, 191, 181, 179, 173, 167, 163, 157, 151, 149, 139, 137, 131, 
  127, 113, 109, 107, 103, 101,  97,  89,  83,  79,  73,  71,  67,  
  61,  59,  53,  47,   43,  41,  37,  31,  29,  23,  19,  17,  13,
  11,   7,   5,   3,  587, 373, 367, 359, 353 };

// the bitlength of each l_i
static int bits_li[N] = { 
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 5, 
  5, 5, 5, 5, 4, 4, 3, 3, 2, 10, 9, 9, 9, 9 };  

// the shortest differential addition chains for each l_i (used in yMUL)
static uint32_t addc[N] = {
  0x231, 0x324,  0x10, 0x2D8, 0x140,  0x50,  0x14, 0x108, 0x101, 0x144, 0x148, 
  0x122, 0x141,  0x51, 0x12A, 0x109, 0x118, 0x1A2, 0x181,   0x0, 0x134, 0x194, 
  0x185, 0x191,   0x1, 0x198,  0x82,  0x88,  0x81,  0x8A,  0xC0,  0xA1,  0xD0,  
   0xC2,  0x98,  0xD4,  0xC5,  0x2B,  0x40,  0xE8,  0xE1,  0x4A,  0x60,   0xC,  
   0x68,  0x49,   0x0,  0x6C,   0x4,  0x14,  0x24,   0x9,  0x2C,  0x32,  0x36,
    0x1,  0x11,  0x18,   0x3,   0x8,   0xA,   0xD,   0x4,   0x5,   0x0,   0x1,
    0x1,   0x0,   0x0, 0x612, 0x268, 0x312,  0xD1, 0x352 };

// the bitlength of the addc[N] (used in yMUL)
static uint8_t addc_len[N] = {
  11, 11, 10, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
  10, 10, 10,  9, 10, 10, 10, 10,  9, 10,  9,  9,  9,  9,  9,  9,
   9,  9,  9,  9,  9,  9,  8,  9,  9,  8,  8,  8,  8,  8,  7,  8, 
   7,  7,  7,  7,  7,  7,  7,  6,  6,  6,  6,  5,  5,  5,  4,  4,
   3,  3,  2,  1,  0, 12, 11, 11, 11, 11 }; 

// -----------------------------------------------------------------------------
// (8x1)-way curve and isogeny operations 

void point_copy_8x1w(htpoint_t R, const htpoint_t P);
void point_cmove_8x1w(htpoint_t R, htpoint_t P, const __m512i b);
void point_cswap_8x1w(htpoint_t R, htpoint_t P, const __m512i b);
__m512i point_isinf_8x1w(const htpoint_t P);
void yDBL_8x1w(htpoint_t R, const htpoint_t P, const htpoint_t A);
void yADD_8x1w(htpoint_t R, const htpoint_t P, const htpoint_t Q, const htpoint_t PQ);
void yMUL_8x1w(htpoint_t R, const htpoint_t P, const htpoint_t A, const uint8_t k);
void elligator_8x1w(htpoint_t Tplus, htpoint_t Tminus, const htpoint_t A);
void yISOG_8x1w(htpoint R[], htpoint_t C, const htpoint_t P, const htpoint_t A, const uint8_t k);
void yEVAL_8x1w(htpoint_t R, const htpoint_t Q, const htpoint P[], const uint8_t k);

// -----------------------------------------------------------------------------
// (2x4)-way curve and isogeny operations 

void point_copy_2x4w(llpoint_t R, const llpoint_t P);
void point_cmove_2x4w(llpoint_t R, llpoint_t P, const uint8_t b);
void point_cswap_2x4w(llpoint_t R, llpoint_t P, const uint8_t b);
uint8_t point_isinf_2x4w(const llpoint_t P);
void yDBL_2x4w(llpoint_t R, const llpoint_t P, const llpoint_t A);
void yADD_2x4w(llpoint_t R, const llpoint_t P, const llpoint_t Q, const llpoint_t PQ);
void yMUL_2x4w(llpoint_t R, const llpoint_t P, const llpoint_t A, const uint8_t k);
void elligator_2x4w(llpoint_t Tplus, llpoint_t Tminus, const llpoint_t A);
void yISOG_2x4w(llpoint_t R[], llpoint_t C, const llpoint_t P, const llpoint_t A, const uint8_t k);
void yEVAL_2x4w(llpoint_t R, const llpoint_t Q, const llpoint_t P[], const uint8_t k);

#endif
