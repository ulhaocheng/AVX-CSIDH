/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "gfparith.h"

// (8x1)-way prime-field operations 

// field addition r = a + b mod 2p
// a, b in [0, 2p) -> r in [0, 2p)
void gfp_add_8x1w(htfe_t r, const htfe_t a, const htfe_t b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4],  b5  = b[5];
  __m512i b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15], b16 = b[16], b17 = b[17]; 
  __m512i r0, r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8;
  __m512i r9, r10, r11, r12, r13, r14, r15, r16, r17, smask;
  const __m512i vp0  = VSET1(ht_pmul2[0]),  vp1  = VSET1(ht_pmul2[1]);
  const __m512i vp2  = VSET1(ht_pmul2[2]),  vp3  = VSET1(ht_pmul2[3]);
  const __m512i vp4  = VSET1(ht_pmul2[4]),  vp5  = VSET1(ht_pmul2[5]);
  const __m512i vp6  = VSET1(ht_pmul2[6]),  vp7  = VSET1(ht_pmul2[7]);
  const __m512i vp8  = VSET1(ht_pmul2[8]),  vp9  = VSET1(ht_pmul2[9]);
  const __m512i vp10 = VSET1(ht_pmul2[10]), vp11 = VSET1(ht_pmul2[11]);
  const __m512i vp12 = VSET1(ht_pmul2[12]), vp13 = VSET1(ht_pmul2[13]);
  const __m512i vp14 = VSET1(ht_pmul2[14]), vp15 = VSET1(ht_pmul2[15]);
  const __m512i vp16 = VSET1(ht_pmul2[16]), vp17 = VSET1(ht_pmul2[17]);
  const __m512i vbmask = VSET1(HT_BMASK);

  // r = a + b
  r0  = VADD(a0,  b0);  r1  = VADD(a1,  b1);  r2  = VADD(a2,  b2);
  r3  = VADD(a3,  b3);  r4  = VADD(a4,  b4);  r5  = VADD(a5,  b5);
  r6  = VADD(a6,  b6);  r7  = VADD(a7,  b7);  r8  = VADD(a8,  b8);
  r9  = VADD(a9,  b9);  r10 = VADD(a10, b10); r11 = VADD(a11, b11);
  r12 = VADD(a12, b12); r13 = VADD(a13, b13); r14 = VADD(a14, b14);
  r15 = VADD(a15, b15); r16 = VADD(a16, b16); r17 = VADD(a17, b17);

  // r = a + b - 2p
  r0   = VSUB(r0,  vp0);  r1  = VSUB(r1,  vp1);  r2  = VSUB(r2,  vp2);
  r3   = VSUB(r3,  vp3);  r4  = VSUB(r4,  vp4);  r5  = VSUB(r5,  vp5);
  r6   = VSUB(r6,  vp6);  r7  = VSUB(r7,  vp7);  r8  = VSUB(r8,  vp8);
  r9   = VSUB(r9,  vp9);  r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);
  r12  = VSUB(r12, vp12); r13 = VSUB(r13, vp13); r14 = VSUB(r14, vp14);
  r15  = VSUB(r15, vp15); r16 = VSUB(r16, vp16); r17 = VSUB(r17, vp17);

  // check the current r is positive or negative 
  // carry propagation
  r1  = VADD(r1,  VSRA(r0,  HT_BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  HT_BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  HT_BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  HT_BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  HT_BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  HT_BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  HT_BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  HT_BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  HT_BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  HT_BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, HT_BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, HT_BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, HT_BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, HT_BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSRA(r14, HT_BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSRA(r15, HT_BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSRA(r16, HT_BRADIX)); r16 = VAND(r16, vbmask);

  // if r is positive, then the corresponding element in smask = 0; 
  // if r is negative, then the corresponding element in smask = all-1.
  smask = VSRA(r17, 63);
  // r = r + (2p & smask), add either 2p or 0 to the current r
  r0  = VADD(r0,  VAND(vp0,  smask)); r1  = VADD(r1,  VAND(vp1,  smask));
  r2  = VADD(r2,  VAND(vp2,  smask)); r3  = VADD(r3,  VAND(vp3,  smask));
  r4  = VADD(r4,  VAND(vp4,  smask)); r5  = VADD(r5,  VAND(vp5,  smask));
  r6  = VADD(r6,  VAND(vp6,  smask)); r7  = VADD(r7,  VAND(vp7,  smask));
  r8  = VADD(r8,  VAND(vp8,  smask)); r9  = VADD(r9,  VAND(vp9,  smask));
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask));
  r12 = VADD(r12, VAND(vp12, smask)); r13 = VADD(r13, VAND(vp13, smask));
  r14 = VADD(r14, VAND(vp14, smask)); r15 = VADD(r15, VAND(vp15, smask));
  r16 = VADD(r16, VAND(vp16, smask)); r17 = VADD(r17, VAND(vp17, smask));

  // carry propagation
  r1  = VADD(r1,  VSHR(r0,  HT_BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSHR(r1,  HT_BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSHR(r2,  HT_BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSHR(r3,  HT_BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSHR(r4,  HT_BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSHR(r5,  HT_BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSHR(r6,  HT_BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSHR(r7,  HT_BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSHR(r8,  HT_BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSHR(r9,  HT_BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSHR(r10, HT_BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSHR(r11, HT_BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSHR(r12, HT_BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSHR(r13, HT_BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSHR(r14, HT_BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSHR(r15, HT_BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSHR(r16, HT_BRADIX)); r16 = VAND(r16, vbmask);
  r17 = VAND(r17, vbmask);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// field subtraction r = a - b mod 2p 
// a, b in [0, 2p) -> r in [0, 2p)
void gfp_sub_8x1w(htfe_t r, const htfe_t a, const htfe_t b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4],  b5  = b[5];
  __m512i b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15], b16 = b[16], b17 = b[17]; 
  __m512i r0, r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8;
  __m512i r9, r10, r11, r12, r13, r14, r15, r16, r17, smask;
  const __m512i vp0  = VSET1(ht_pmul2[0]),  vp1  = VSET1(ht_pmul2[1]);
  const __m512i vp2  = VSET1(ht_pmul2[2]),  vp3  = VSET1(ht_pmul2[3]);
  const __m512i vp4  = VSET1(ht_pmul2[4]),  vp5  = VSET1(ht_pmul2[5]);
  const __m512i vp6  = VSET1(ht_pmul2[6]),  vp7  = VSET1(ht_pmul2[7]);
  const __m512i vp8  = VSET1(ht_pmul2[8]),  vp9  = VSET1(ht_pmul2[9]);
  const __m512i vp10 = VSET1(ht_pmul2[10]), vp11 = VSET1(ht_pmul2[11]);
  const __m512i vp12 = VSET1(ht_pmul2[12]), vp13 = VSET1(ht_pmul2[13]);
  const __m512i vp14 = VSET1(ht_pmul2[14]), vp15 = VSET1(ht_pmul2[15]);
  const __m512i vp16 = VSET1(ht_pmul2[16]), vp17 = VSET1(ht_pmul2[17]);
  const __m512i vbmask = VSET1(HT_BMASK);

  // r = a - b
  r0  = VSUB(a0,  b0);  r1  = VSUB(a1,  b1);  r2  = VSUB(a2,  b2);
  r3  = VSUB(a3,  b3);  r4  = VSUB(a4,  b4);  r5  = VSUB(a5,  b5);
  r6  = VSUB(a6,  b6);  r7  = VSUB(a7,  b7);  r8  = VSUB(a8,  b8);
  r9  = VSUB(a9,  b9);  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
  r12 = VSUB(a12, b12); r13 = VSUB(a13, b13); r14 = VSUB(a14, b14);
  r15 = VSUB(a15, b15); r16 = VSUB(a16, b16); r17 = VSUB(a17, b17);

  // check the current r is positive or negative 
  // carry propagation
  r1  = VADD(r1,  VSRA(r0,  HT_BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  HT_BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  HT_BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  HT_BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  HT_BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  HT_BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  HT_BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  HT_BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  HT_BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  HT_BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, HT_BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, HT_BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, HT_BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, HT_BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSRA(r14, HT_BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSRA(r15, HT_BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSRA(r16, HT_BRADIX)); r16 = VAND(r16, vbmask);

  // if r is positive, then the corresponding element in smask = 0; 
  // if r is negative, then the corresponding element in smask = all-1.
  smask = VSRA(r17, 63);
  // r = r + (2p & smask), add either 2p or 0 to the current r
  r0  = VADD(r0,  VAND(vp0,  smask)); r1  = VADD(r1,  VAND(vp1,  smask));
  r2  = VADD(r2,  VAND(vp2,  smask)); r3  = VADD(r3,  VAND(vp3,  smask));
  r4  = VADD(r4,  VAND(vp4,  smask)); r5  = VADD(r5,  VAND(vp5,  smask));
  r6  = VADD(r6,  VAND(vp6,  smask)); r7  = VADD(r7,  VAND(vp7,  smask));
  r8  = VADD(r8,  VAND(vp8,  smask)); r9  = VADD(r9,  VAND(vp9,  smask));
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask));
  r12 = VADD(r12, VAND(vp12, smask)); r13 = VADD(r13, VAND(vp13, smask));
  r14 = VADD(r14, VAND(vp14, smask)); r15 = VADD(r15, VAND(vp15, smask));
  r16 = VADD(r16, VAND(vp16, smask)); r17 = VADD(r17, VAND(vp17, smask));

  // carry propagation
  r1  = VADD(r1,  VSHR(r0,  HT_BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSHR(r1,  HT_BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSHR(r2,  HT_BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSHR(r3,  HT_BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSHR(r4,  HT_BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSHR(r5,  HT_BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSHR(r6,  HT_BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSHR(r7,  HT_BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSHR(r8,  HT_BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSHR(r9,  HT_BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSHR(r10, HT_BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSHR(r11, HT_BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSHR(r12, HT_BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSHR(r13, HT_BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSHR(r14, HT_BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSHR(r15, HT_BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSHR(r16, HT_BRADIX)); r16 = VAND(r16, vbmask);
  r17 = VAND(r17, vbmask);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// Montgomery multiplication r = a * b mod 2p
// multiplication (Karatsuba) separated with reduction (product-scanning)
// -> r in [0, 2p)
void gfp_mul_8x1w(htfe_t r, const htfe_t a, const htfe_t b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4],  b5  = b[5];
  __m512i b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15], b16 = b[16], b17 = b[17]; 
  __m512i z0,  z1,  z2,  z3,  z4,  z5,  z6,  z7,  z8,  z9,  z10, z11;
  __m512i z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23;
  __m512i z24, z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35;
  __m512i r0, r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8;
  __m512i r9, r10, r11, r12, r13, r14, r15, r16, r17, accu = VZERO;
  __m512i ta0, ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8;
  __m512i tb0, tb1, tb2, tb3, tb4, tb5, tb6, tb7, tb8;
  __m512i m0, m1,  m2,  m3,  m4,  m5,  m6,  m7,  m8;
  __m512i m9, m10, m11, m12, m13, m14, m15, m16, m17;
  const __m512i vp0  = VSET1(ht_p[0]),  vp1  = VSET1(ht_p[1]),  vp2  = VSET1(ht_p[2]);
  const __m512i vp3  = VSET1(ht_p[3]),  vp4  = VSET1(ht_p[4]),  vp5  = VSET1(ht_p[5]);
  const __m512i vp6  = VSET1(ht_p[6]),  vp7  = VSET1(ht_p[7]),  vp8  = VSET1(ht_p[8]);
  const __m512i vp9  = VSET1(ht_p[9]),  vp10 = VSET1(ht_p[10]), vp11 = VSET1(ht_p[11]);
  const __m512i vp12 = VSET1(ht_p[12]), vp13 = VSET1(ht_p[13]), vp14 = VSET1(ht_p[14]);
  const __m512i vp15 = VSET1(ht_p[15]), vp16 = VSET1(ht_p[16]), vp17 = VSET1(ht_p[17]);
  const __m512i vbmask = VSET1(HT_BMASK), vw = VSET1(HT_MONTW);

  // ---------------------------------------------------------------------------
  // compute zL(z0~z17) by aL(a0~a8) * bL(b0~b8)

  z0 = VMUL(a0, b0);

  z1 = VMUL(a0, b1); z1 = VADD(z1, VMUL(a1, b0));

  z2 = VMUL(a0, b2); z2 = VADD(z2, VMUL(a1, b1)); z2 = VADD(z2, VMUL(a2, b0));

  z3 = VMUL(a0, b3); z3 = VADD(z3, VMUL(a1, b2)); z3 = VADD(z3, VMUL(a2, b1)); 
  z3 = VADD(z3, VMUL(a3, b0));

  z4 = VMUL(a0, b4); z4 = VADD(z4, VMUL(a1, b3)); z4 = VADD(z4, VMUL(a2, b2)); 
  z4 = VADD(z4, VMUL(a3, b1)); z4 = VADD(z4, VMUL(a4, b0));

  z5 = VMUL(a0, b5); z5 = VADD(z5, VMUL(a1, b4)); z5 = VADD(z5, VMUL(a2, b3)); 
  z5 = VADD(z5, VMUL(a3, b2)); z5 = VADD(z5, VMUL(a4, b1)); z5 = VADD(z5, VMUL(a5, b0));

  z6 = VMUL(a0, b6); z6 = VADD(z6, VMUL(a1, b5)); z6 = VADD(z6, VMUL(a2, b4)); 
  z6 = VADD(z6, VMUL(a3, b3)); z6 = VADD(z6, VMUL(a4, b2)); z6 = VADD(z6, VMUL(a5, b1)); 
  z6 = VADD(z6, VMUL(a6, b0));

  z7 = VMUL(a0, b7); z7 = VADD(z7, VMUL(a1, b6)); z7 = VADD(z7, VMUL(a2, b5)); 
  z7 = VADD(z7, VMUL(a3, b4)); z7 = VADD(z7, VMUL(a4, b3)); z7 = VADD(z7, VMUL(a5, b2)); 
  z7 = VADD(z7, VMUL(a6, b1)); z7 = VADD(z7, VMUL(a7, b0));

  z8 = VMUL(a0, b8); z8 = VADD(z8, VMUL(a1, b7)); z8 = VADD(z8, VMUL(a2, b6)); 
  z8 = VADD(z8, VMUL(a3, b5)); z8 = VADD(z8, VMUL(a4, b4)); z8 = VADD(z8, VMUL(a5, b3)); 
  z8 = VADD(z8, VMUL(a6, b2)); z8 = VADD(z8, VMUL(a7, b1)); z8 = VADD(z8, VMUL(a8, b0));

  z9 = VMUL(a1, b8); z9 = VADD(z9, VMUL(a2, b7)); z9 = VADD(z9, VMUL(a3, b6)); 
  z9 = VADD(z9, VMUL(a4, b5)); z9 = VADD(z9, VMUL(a5, b4)); z9 = VADD(z9, VMUL(a6, b3)); 
  z9 = VADD(z9, VMUL(a7, b2)); z9 = VADD(z9, VMUL(a8, b1));

  z10 = VMUL(a2, b8); z10 = VADD(z10, VMUL(a3, b7)); z10 = VADD(z10, VMUL(a4, b6)); 
  z10 = VADD(z10, VMUL(a5, b5)); z10 = VADD(z10, VMUL(a6, b4)); z10 = VADD(z10, VMUL(a7, b3)); 
  z10 = VADD(z10, VMUL(a8, b2));

  z11 = VMUL(a3, b8); z11 = VADD(z11, VMUL(a4, b7)); z11 = VADD(z11, VMUL(a5, b6)); 
  z11 = VADD(z11, VMUL(a6, b5)); z11 = VADD(z11, VMUL(a7, b4)); z11 = VADD(z11, VMUL(a8, b3));

  z12 = VMUL(a4, b8); z12 = VADD(z12, VMUL(a5, b7)); z12 = VADD(z12, VMUL(a6, b6)); 
  z12 = VADD(z12, VMUL(a7, b5)); z12 = VADD(z12, VMUL(a8, b4));

  z13 = VMUL(a5, b8); z13 = VADD(z13, VMUL(a6, b7)); z13 = VADD(z13, VMUL(a7, b6)); 
  z13 = VADD(z13, VMUL(a8, b5));

  z14 = VMUL(a6, b8); z14 = VADD(z14, VMUL(a7, b7)); z14 = VADD(z14, VMUL(a8, b6));

  z15 = VMUL(a7, b8); z15 = VADD(z15, VMUL(a8, b7));

  z16 = VMUL(a8, b8);

  z17 = VZERO;

  // ---------------------------------------------------------------------------
  // compute zH(z18~z35) with aH(a9~a17)*bH(b9~b17)
  
  z18 = VMUL(a9, b9);

  z19 = VMUL(a9, b10); z19 = VADD(z19, VMUL(a10, b9));

  z20 = VMUL(a9, b11); z20 = VADD(z20, VMUL(a10, b10)); z20 = VADD(z20, VMUL(a11, b9));

  z21 = VMUL(a9, b12); z21 = VADD(z21, VMUL(a10, b11)); z21 = VADD(z21, VMUL(a11, b10)); 
  z21 = VADD(z21, VMUL(a12, b9));

  z22 = VMUL(a9, b13); z22 = VADD(z22, VMUL(a10, b12)); z22 = VADD(z22, VMUL(a11, b11)); 
  z22 = VADD(z22, VMUL(a12, b10)); z22 = VADD(z22, VMUL(a13, b9));

  z23 = VMUL(a9, b14); z23 = VADD(z23, VMUL(a10, b13)); z23 = VADD(z23, VMUL(a11, b12)); 
  z23 = VADD(z23, VMUL(a12, b11)); z23 = VADD(z23, VMUL(a13, b10)); z23 = VADD(z23, VMUL(a14, b9));

  z24 = VMUL(a9, b15); z24 = VADD(z24, VMUL(a10, b14)); z24 = VADD(z24, VMUL(a11, b13)); 
  z24 = VADD(z24, VMUL(a12, b12)); z24 = VADD(z24, VMUL(a13, b11)); z24 = VADD(z24, VMUL(a14, b10)); 
  z24 = VADD(z24, VMUL(a15, b9));

  z25 = VMUL(a9, b16); z25 = VADD(z25, VMUL(a10, b15)); z25 = VADD(z25, VMUL(a11, b14)); 
  z25 = VADD(z25, VMUL(a12, b13)); z25 = VADD(z25, VMUL(a13, b12)); z25 = VADD(z25, VMUL(a14, b11)); 
  z25 = VADD(z25, VMUL(a15, b10)); z25 = VADD(z25, VMUL(a16, b9));

  z26 = VMUL(a9, b17); z26 = VADD(z26, VMUL(a10, b16)); z26 = VADD(z26, VMUL(a11, b15)); 
  z26 = VADD(z26, VMUL(a12, b14)); z26 = VADD(z26, VMUL(a13, b13)); z26 = VADD(z26, VMUL(a14, b12)); 
  z26 = VADD(z26, VMUL(a15, b11)); z26 = VADD(z26, VMUL(a16, b10)); z26 = VADD(z26, VMUL(a17, b9));

  z27 = VMUL(a10, b17); z27 = VADD(z27, VMUL(a11, b16)); z27 = VADD(z27, VMUL(a12, b15)); 
  z27 = VADD(z27, VMUL(a13, b14)); z27 = VADD(z27, VMUL(a14, b13)); z27 = VADD(z27, VMUL(a15, b12)); 
  z27 = VADD(z27, VMUL(a16, b11)); z27 = VADD(z27, VMUL(a17, b10));

  z28 = VMUL(a11, b17); z28 = VADD(z28, VMUL(a12, b16)); z28 = VADD(z28, VMUL(a13, b15)); 
  z28 = VADD(z28, VMUL(a14, b14)); z28 = VADD(z28, VMUL(a15, b13)); z28 = VADD(z28, VMUL(a16, b12)); 
  z28 = VADD(z28, VMUL(a17, b11));

  z29 = VMUL(a12, b17); z29 = VADD(z29, VMUL(a13, b16)); z29 = VADD(z29, VMUL(a14, b15)); 
  z29 = VADD(z29, VMUL(a15, b14)); z29 = VADD(z29, VMUL(a16, b13)); z29 = VADD(z29, VMUL(a17, b12));

  z30 = VMUL(a13, b17); z30 = VADD(z30, VMUL(a14, b16)); z30 = VADD(z30, VMUL(a15, b15)); 
  z30 = VADD(z30, VMUL(a16, b14)); z30 = VADD(z30, VMUL(a17, b13));

  z31 = VMUL(a14, b17); z31 = VADD(z31, VMUL(a15, b16)); z31 = VADD(z31, VMUL(a16, b15)); 
  z31 = VADD(z31, VMUL(a17, b14));

  z32 = VMUL(a15, b17); z32 = VADD(z32, VMUL(a16, b16)); z32 = VADD(z32, VMUL(a17, b15));

  z33 = VMUL(a16, b17); z33 = VADD(z33, VMUL(a17, b16));

  z34 = VMUL(a17, b17);

  z35 = VZERO;

  // ---------------------------------------------------------------------------
  // ta (ta0~ta8) = aL(a0~a8) + aH(a9~a17)

  ta0 = VADD(a0, a9);  ta1 = VADD(a1, a10); ta2 = VADD(a2, a11);
  ta3 = VADD(a3, a12); ta4 = VADD(a4, a13); ta5 = VADD(a5, a14);
  ta6 = VADD(a6, a15); ta7 = VADD(a7, a16); ta8 = VADD(a8, a17);

  // tb (tb0~tb8) = bL(b0~b8) + bH(b9~b17)

  tb0 = VADD(b0, b9);  tb1 = VADD(b1, b10); tb2 = VADD(b2, b11);
  tb3 = VADD(b3, b12); tb4 = VADD(b4, b13); tb5 = VADD(b5, b14);
  tb6 = VADD(b6, b15); tb7 = VADD(b7, b16); tb8 = VADD(b8, b17);

  // ---------------------------------------------------------------------------
  // zM = ta * tb - zL - zH

  m0 = VMUL(ta0, tb0);

  m1 = VMUL(ta0, tb1); m1 = VADD(m1, VMUL(ta1, tb0));

  m2 = VMUL(ta0, tb2); m2 = VADD(m2, VMUL(ta1, tb1)); m2 = VADD(m2, VMUL(ta2, tb0));

  m3 = VMUL(ta0, tb3); m3 = VADD(m3, VMUL(ta1, tb2)); m3 = VADD(m3, VMUL(ta2, tb1)); 
  m3 = VADD(m3, VMUL(ta3, tb0));

  m4 = VMUL(ta0, tb4); m4 = VADD(m4, VMUL(ta1, tb3)); m4 = VADD(m4, VMUL(ta2, tb2)); 
  m4 = VADD(m4, VMUL(ta3, tb1)); m4 = VADD(m4, VMUL(ta4, tb0));

  m5 = VMUL(ta0, tb5); m5 = VADD(m5, VMUL(ta1, tb4)); m5 = VADD(m5, VMUL(ta2, tb3)); 
  m5 = VADD(m5, VMUL(ta3, tb2)); m5 = VADD(m5, VMUL(ta4, tb1)); m5 = VADD(m5, VMUL(ta5, tb0));

  m6 = VMUL(ta0, tb6); m6 = VADD(m6, VMUL(ta1, tb5)); m6 = VADD(m6, VMUL(ta2, tb4)); 
  m6 = VADD(m6, VMUL(ta3, tb3)); m6 = VADD(m6, VMUL(ta4, tb2)); m6 = VADD(m6, VMUL(ta5, tb1)); 
  m6 = VADD(m6, VMUL(ta6, tb0));

  m7 = VMUL(ta0, tb7); m7 = VADD(m7, VMUL(ta1, tb6)); m7 = VADD(m7, VMUL(ta2, tb5)); 
  m7 = VADD(m7, VMUL(ta3, tb4)); m7 = VADD(m7, VMUL(ta4, tb3)); m7 = VADD(m7, VMUL(ta5, tb2)); 
  m7 = VADD(m7, VMUL(ta6, tb1)); m7 = VADD(m7, VMUL(ta7, tb0));

  m8 = VMUL(ta0, tb8); m8 = VADD(m8, VMUL(ta1, tb7)); m8 = VADD(m8, VMUL(ta2, tb6)); 
  m8 = VADD(m8, VMUL(ta3, tb5)); m8 = VADD(m8, VMUL(ta4, tb4)); m8 = VADD(m8, VMUL(ta5, tb3)); 
  m8 = VADD(m8, VMUL(ta6, tb2)); m8 = VADD(m8, VMUL(ta7, tb1)); m8 = VADD(m8, VMUL(ta8, tb0));

  m9 = VMUL(ta1, tb8); m9 = VADD(m9, VMUL(ta2, tb7)); m9 = VADD(m9, VMUL(ta3, tb6)); 
  m9 = VADD(m9, VMUL(ta4, tb5)); m9 = VADD(m9, VMUL(ta5, tb4)); m9 = VADD(m9, VMUL(ta6, tb3)); 
  m9 = VADD(m9, VMUL(ta7, tb2)); m9 = VADD(m9, VMUL(ta8, tb1));

  m10 = VMUL(ta2, tb8); m10 = VADD(m10, VMUL(ta3, tb7)); m10 = VADD(m10, VMUL(ta4, tb6)); 
  m10 = VADD(m10, VMUL(ta5, tb5)); m10 = VADD(m10, VMUL(ta6, tb4)); m10 = VADD(m10, VMUL(ta7, tb3)); 
  m10 = VADD(m10, VMUL(ta8, tb2));

  m11 = VMUL(ta3, tb8); m11 = VADD(m11, VMUL(ta4, tb7)); m11 = VADD(m11, VMUL(ta5, tb6)); 
  m11 = VADD(m11, VMUL(ta6, tb5)); m11 = VADD(m11, VMUL(ta7, tb4)); m11 = VADD(m11, VMUL(ta8, tb3));

  m12 = VMUL(ta4, tb8); m12 = VADD(m12, VMUL(ta5, tb7)); m12 = VADD(m12, VMUL(ta6, tb6)); 
  m12 = VADD(m12, VMUL(ta7, tb5)); m12 = VADD(m12, VMUL(ta8, tb4));

  m13 = VMUL(ta5, tb8); m13 = VADD(m13, VMUL(ta6, tb7)); m13 = VADD(m13, VMUL(ta7, tb6)); 
  m13 = VADD(m13, VMUL(ta8, tb5));

  m14 = VMUL(ta6, tb8); m14 = VADD(m14, VMUL(ta7, tb7)); m14 = VADD(m14, VMUL(ta8, tb6));

  m15 = VMUL(ta7, tb8); m15 = VADD(m15, VMUL(ta8, tb7));

  m16 = VMUL(ta8, tb8);

  m17 = VZERO;

  m0  = VSUB(m0,  VADD(z0,  z18)); m1  = VSUB(m1,  VADD(z1,  z19)); 
  m2  = VSUB(m2,  VADD(z2,  z20)); m3  = VSUB(m3,  VADD(z3,  z21)); 
  m4  = VSUB(m4,  VADD(z4,  z22)); m5  = VSUB(m5,  VADD(z5,  z23)); 
  m6  = VSUB(m6,  VADD(z6,  z24)); m7  = VSUB(m7,  VADD(z7,  z25)); 
  m8  = VSUB(m8,  VADD(z8,  z26)); m9  = VSUB(m9,  VADD(z9,  z27)); 
  m10 = VSUB(m10, VADD(z10, z28)); m11 = VSUB(m11, VADD(z11, z29)); 
  m12 = VSUB(m12, VADD(z12, z30)); m13 = VSUB(m13, VADD(z13, z31)); 
  m14 = VSUB(m14, VADD(z14, z32)); m15 = VSUB(m15, VADD(z15, z33)); 
  m16 = VSUB(m16, VADD(z16, z34)); m17 = VSUB(m17, VADD(z17, z35));

  // ---------------------------------------------------------------------------
  // z = z + zM

  z9  = VADD(z9,  m0);  z10 = VADD(z10, m1);  z11 = VADD(z11, m2); 
  z12 = VADD(z12, m3);  z13 = VADD(z13, m4);  z14 = VADD(z14, m5); 
  z15 = VADD(z15, m6);  z16 = VADD(z16, m7);  z17 = VADD(z17, m8); 
  z18 = VADD(z18, m9);  z19 = VADD(z19, m10); z20 = VADD(z20, m11); 
  z21 = VADD(z21, m12); z22 = VADD(z22, m13); z23 = VADD(z23, m14); 
  z24 = VADD(z24, m15); z25 = VADD(z25, m16); z26 = VADD(z26, m17);

  // ---------------------------------------------------------------------------
  // 1st loop of Montgomery reduction

  accu = VADD(accu, VAND(z0, vbmask)); z1 = VADD(z1, VSHR(z0, HT_BRADIX)); r0 = VAND(accu, vbmask); 
  r0 = VAND(vbmask, VMUL(r0, vw)); accu = VADD(accu, VMUL(r0, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp1)); accu = VADD(accu, VAND(z1, vbmask)); z2 = VADD(z2, VSHR(z1, HT_BRADIX)); 
  r1 = VAND(accu, vbmask); r1 = VAND(vbmask, VMUL(r1, vw)); accu = VADD(accu, VMUL(r1, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp2)); accu = VADD(accu, VMUL(r1, vp1)); accu = VADD(accu, VAND(z2, vbmask)); 
  z3 = VADD(z3, VSHR(z2, HT_BRADIX)); r2 = VAND(accu, vbmask); r2 = VAND(vbmask, VMUL(r2, vw)); 
  accu = VADD(accu, VMUL(r2, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp3)); accu = VADD(accu, VMUL(r1, vp2)); accu = VADD(accu, VMUL(r2, vp1)); 
  accu = VADD(accu, VAND(z3, vbmask)); z4 = VADD(z4, VSHR(z3, HT_BRADIX)); r3 = VAND(accu, vbmask); 
  r3 = VAND(vbmask, VMUL(r3, vw)); accu = VADD(accu, VMUL(r3, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp4)); accu = VADD(accu, VMUL(r1, vp3)); accu = VADD(accu, VMUL(r2, vp2)); 
  accu = VADD(accu, VMUL(r3, vp1)); accu = VADD(accu, VAND(z4, vbmask)); z5 = VADD(z5, VSHR(z4, HT_BRADIX)); 
  r4 = VAND(accu, vbmask); r4 = VAND(vbmask, VMUL(r4, vw)); accu = VADD(accu, VMUL(r4, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp5)); accu = VADD(accu, VMUL(r1, vp4)); accu = VADD(accu, VMUL(r2, vp3)); 
  accu = VADD(accu, VMUL(r3, vp2)); accu = VADD(accu, VMUL(r4, vp1)); accu = VADD(accu, VAND(z5, vbmask)); 
  z6 = VADD(z6, VSHR(z5, HT_BRADIX)); r5 = VAND(accu, vbmask); r5 = VAND(vbmask, VMUL(r5, vw)); 
  accu = VADD(accu, VMUL(r5, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp6)); accu = VADD(accu, VMUL(r1, vp5)); accu = VADD(accu, VMUL(r2, vp4)); 
  accu = VADD(accu, VMUL(r3, vp3)); accu = VADD(accu, VMUL(r4, vp2)); accu = VADD(accu, VMUL(r5, vp1)); 
  accu = VADD(accu, VAND(z6, vbmask)); z7 = VADD(z7, VSHR(z6, HT_BRADIX)); r6 = VAND(accu, vbmask); 
  r6 = VAND(vbmask, VMUL(r6, vw)); accu = VADD(accu, VMUL(r6, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp7)); accu = VADD(accu, VMUL(r1, vp6)); accu = VADD(accu, VMUL(r2, vp5)); 
  accu = VADD(accu, VMUL(r3, vp4)); accu = VADD(accu, VMUL(r4, vp3)); accu = VADD(accu, VMUL(r5, vp2)); 
  accu = VADD(accu, VMUL(r6, vp1)); accu = VADD(accu, VAND(z7, vbmask)); z8 = VADD(z8, VSHR(z7, HT_BRADIX)); 
  r7 = VAND(accu, vbmask); r7 = VAND(vbmask, VMUL(r7, vw)); accu = VADD(accu, VMUL(r7, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp8)); accu = VADD(accu, VMUL(r1, vp7)); accu = VADD(accu, VMUL(r2, vp6)); 
  accu = VADD(accu, VMUL(r3, vp5)); accu = VADD(accu, VMUL(r4, vp4)); accu = VADD(accu, VMUL(r5, vp3)); 
  accu = VADD(accu, VMUL(r6, vp2)); accu = VADD(accu, VMUL(r7, vp1)); accu = VADD(accu, VAND(z8, vbmask)); 
  z9 = VADD(z9, VSHR(z8, HT_BRADIX)); r8 = VAND(accu, vbmask); r8 = VAND(vbmask, VMUL(r8, vw)); 
  accu = VADD(accu, VMUL(r8, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp9)); accu = VADD(accu, VMUL(r1, vp8)); accu = VADD(accu, VMUL(r2, vp7)); 
  accu = VADD(accu, VMUL(r3, vp6)); accu = VADD(accu, VMUL(r4, vp5)); accu = VADD(accu, VMUL(r5, vp4)); 
  accu = VADD(accu, VMUL(r6, vp3)); accu = VADD(accu, VMUL(r7, vp2)); accu = VADD(accu, VMUL(r8, vp1)); 
  accu = VADD(accu, VAND(z9, vbmask)); z10 = VADD(z10, VSHR(z9, HT_BRADIX)); r9 = VAND(accu, vbmask); 
  r9 = VAND(vbmask, VMUL(r9, vw)); accu = VADD(accu, VMUL(r9, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp10)); accu = VADD(accu, VMUL(r1, vp9)); accu = VADD(accu, VMUL(r2, vp8)); 
  accu = VADD(accu, VMUL(r3, vp7)); accu = VADD(accu, VMUL(r4, vp6)); accu = VADD(accu, VMUL(r5, vp5)); 
  accu = VADD(accu, VMUL(r6, vp4)); accu = VADD(accu, VMUL(r7, vp3)); accu = VADD(accu, VMUL(r8, vp2)); 
  accu = VADD(accu, VMUL(r9, vp1)); accu = VADD(accu, VAND(z10, vbmask)); z11 = VADD(z11, VSHR(z10, HT_BRADIX)); 
  r10 = VAND(accu, vbmask); r10 = VAND(vbmask, VMUL(r10, vw)); accu = VADD(accu, VMUL(r10, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp11)); accu = VADD(accu, VMUL(r1, vp10)); accu = VADD(accu, VMUL(r2, vp9)); 
  accu = VADD(accu, VMUL(r3, vp8)); accu = VADD(accu, VMUL(r4, vp7)); accu = VADD(accu, VMUL(r5, vp6)); 
  accu = VADD(accu, VMUL(r6, vp5)); accu = VADD(accu, VMUL(r7, vp4)); accu = VADD(accu, VMUL(r8, vp3)); 
  accu = VADD(accu, VMUL(r9, vp2)); accu = VADD(accu, VMUL(r10, vp1)); accu = VADD(accu, VAND(z11, vbmask)); 
  z12 = VADD(z12, VSHR(z11, HT_BRADIX)); r11 = VAND(accu, vbmask); r11 = VAND(vbmask, VMUL(r11, vw)); 
  accu = VADD(accu, VMUL(r11, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp12)); accu = VADD(accu, VMUL(r1, vp11)); accu = VADD(accu, VMUL(r2, vp10)); 
  accu = VADD(accu, VMUL(r3, vp9)); accu = VADD(accu, VMUL(r4, vp8)); accu = VADD(accu, VMUL(r5, vp7)); 
  accu = VADD(accu, VMUL(r6, vp6)); accu = VADD(accu, VMUL(r7, vp5)); accu = VADD(accu, VMUL(r8, vp4)); 
  accu = VADD(accu, VMUL(r9, vp3)); accu = VADD(accu, VMUL(r10, vp2)); accu = VADD(accu, VMUL(r11, vp1)); 
  accu = VADD(accu, VAND(z12, vbmask)); z13 = VADD(z13, VSHR(z12, HT_BRADIX)); r12 = VAND(accu, vbmask); 
  r12 = VAND(vbmask, VMUL(r12, vw)); accu = VADD(accu, VMUL(r12, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp13)); accu = VADD(accu, VMUL(r1, vp12)); accu = VADD(accu, VMUL(r2, vp11)); 
  accu = VADD(accu, VMUL(r3, vp10)); accu = VADD(accu, VMUL(r4, vp9)); accu = VADD(accu, VMUL(r5, vp8)); 
  accu = VADD(accu, VMUL(r6, vp7)); accu = VADD(accu, VMUL(r7, vp6)); accu = VADD(accu, VMUL(r8, vp5)); 
  accu = VADD(accu, VMUL(r9, vp4)); accu = VADD(accu, VMUL(r10, vp3)); accu = VADD(accu, VMUL(r11, vp2)); 
  accu = VADD(accu, VMUL(r12, vp1)); accu = VADD(accu, VAND(z13, vbmask)); z14 = VADD(z14, VSHR(z13, HT_BRADIX)); 
  r13 = VAND(accu, vbmask); r13 = VAND(vbmask, VMUL(r13, vw)); accu = VADD(accu, VMUL(r13, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp14)); accu = VADD(accu, VMUL(r1, vp13)); accu = VADD(accu, VMUL(r2, vp12)); 
  accu = VADD(accu, VMUL(r3, vp11)); accu = VADD(accu, VMUL(r4, vp10)); accu = VADD(accu, VMUL(r5, vp9)); 
  accu = VADD(accu, VMUL(r6, vp8)); accu = VADD(accu, VMUL(r7, vp7)); accu = VADD(accu, VMUL(r8, vp6)); 
  accu = VADD(accu, VMUL(r9, vp5)); accu = VADD(accu, VMUL(r10, vp4)); accu = VADD(accu, VMUL(r11, vp3)); 
  accu = VADD(accu, VMUL(r12, vp2)); accu = VADD(accu, VMUL(r13, vp1)); accu = VADD(accu, VAND(z14, vbmask)); 
  z15 = VADD(z15, VSHR(z14, HT_BRADIX)); r14 = VAND(accu, vbmask); r14 = VAND(vbmask, VMUL(r14, vw)); 
  accu = VADD(accu, VMUL(r14, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp15)); accu = VADD(accu, VMUL(r1, vp14)); accu = VADD(accu, VMUL(r2, vp13)); 
  accu = VADD(accu, VMUL(r3, vp12)); accu = VADD(accu, VMUL(r4, vp11)); accu = VADD(accu, VMUL(r5, vp10)); 
  accu = VADD(accu, VMUL(r6, vp9)); accu = VADD(accu, VMUL(r7, vp8)); accu = VADD(accu, VMUL(r8, vp7)); 
  accu = VADD(accu, VMUL(r9, vp6)); accu = VADD(accu, VMUL(r10, vp5)); accu = VADD(accu, VMUL(r11, vp4)); 
  accu = VADD(accu, VMUL(r12, vp3)); accu = VADD(accu, VMUL(r13, vp2)); accu = VADD(accu, VMUL(r14, vp1)); 
  accu = VADD(accu, VAND(z15, vbmask)); z16 = VADD(z16, VSHR(z15, HT_BRADIX)); r15 = VAND(accu, vbmask); 
  r15 = VAND(vbmask, VMUL(r15, vw)); accu = VADD(accu, VMUL(r15, vp0)); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp16)); accu = VADD(accu, VMUL(r1, vp15)); accu = VADD(accu, VMUL(r2, vp14)); 
  accu = VADD(accu, VMUL(r3, vp13)); accu = VADD(accu, VMUL(r4, vp12)); accu = VADD(accu, VMUL(r5, vp11)); 
  accu = VADD(accu, VMUL(r6, vp10)); accu = VADD(accu, VMUL(r7, vp9)); accu = VADD(accu, VMUL(r8, vp8)); 
  accu = VADD(accu, VMUL(r9, vp7)); accu = VADD(accu, VMUL(r10, vp6)); accu = VADD(accu, VMUL(r11, vp5)); 
  accu = VADD(accu, VMUL(r12, vp4)); accu = VADD(accu, VMUL(r13, vp3)); accu = VADD(accu, VMUL(r14, vp2)); 
  accu = VADD(accu, VMUL(r15, vp1)); accu = VADD(accu, VAND(z16, vbmask)); z17 = VADD(z17, VSHR(z16, HT_BRADIX)); 
  r16 = VAND(accu, vbmask); r16 = VAND(vbmask, VMUL(r16, vw)); accu = VADD(accu, VMUL(r16, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r0, vp17)); accu = VADD(accu, VMUL(r1, vp16)); accu = VADD(accu, VMUL(r2, vp15)); 
  accu = VADD(accu, VMUL(r3, vp14)); accu = VADD(accu, VMUL(r4, vp13)); accu = VADD(accu, VMUL(r5, vp12)); 
  accu = VADD(accu, VMUL(r6, vp11)); accu = VADD(accu, VMUL(r7, vp10)); accu = VADD(accu, VMUL(r8, vp9)); 
  accu = VADD(accu, VMUL(r9, vp8)); accu = VADD(accu, VMUL(r10, vp7)); accu = VADD(accu, VMUL(r11, vp6)); 
  accu = VADD(accu, VMUL(r12, vp5)); accu = VADD(accu, VMUL(r13, vp4)); accu = VADD(accu, VMUL(r14, vp3)); 
  accu = VADD(accu, VMUL(r15, vp2)); accu = VADD(accu, VMUL(r16, vp1)); accu = VADD(accu, VAND(z17, vbmask)); 
  z18 = VADD(z18, VSHR(z17, HT_BRADIX)); r17 = VAND(accu, vbmask); r17 = VAND(vbmask, VMUL(r17, vw)); 
  accu = VADD(accu, VMUL(r17, vp0)); accu = VSHR(accu, HT_BRADIX);

  // ---------------------------------------------------------------------------
  // 2nd loop of Montgomery reduction

  accu = VADD(accu, VMUL(r1, vp17)); accu = VADD(accu, VMUL(r2, vp16)); accu = VADD(accu, VMUL(r3, vp15)); 
  accu = VADD(accu, VMUL(r4, vp14)); accu = VADD(accu, VMUL(r5, vp13)); accu = VADD(accu, VMUL(r6, vp12)); 
  accu = VADD(accu, VMUL(r7, vp11)); accu = VADD(accu, VMUL(r8, vp10)); accu = VADD(accu, VMUL(r9, vp9)); 
  accu = VADD(accu, VMUL(r10, vp8)); accu = VADD(accu, VMUL(r11, vp7)); accu = VADD(accu, VMUL(r12, vp6)); 
  accu = VADD(accu, VMUL(r13, vp5)); accu = VADD(accu, VMUL(r14, vp4)); accu = VADD(accu, VMUL(r15, vp3)); 
  accu = VADD(accu, VMUL(r16, vp2)); accu = VADD(accu, VMUL(r17, vp1)); accu = VADD(accu, VAND(z18, vbmask)); 
  z19 = VADD(z19, VSHR(z18, HT_BRADIX)); 
  r0 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r2, vp17)); accu = VADD(accu, VMUL(r3, vp16)); accu = VADD(accu, VMUL(r4, vp15)); 
  accu = VADD(accu, VMUL(r5, vp14)); accu = VADD(accu, VMUL(r6, vp13)); accu = VADD(accu, VMUL(r7, vp12)); 
  accu = VADD(accu, VMUL(r8, vp11)); accu = VADD(accu, VMUL(r9, vp10)); accu = VADD(accu, VMUL(r10, vp9)); 
  accu = VADD(accu, VMUL(r11, vp8)); accu = VADD(accu, VMUL(r12, vp7)); accu = VADD(accu, VMUL(r13, vp6)); 
  accu = VADD(accu, VMUL(r14, vp5)); accu = VADD(accu, VMUL(r15, vp4)); accu = VADD(accu, VMUL(r16, vp3)); 
  accu = VADD(accu, VMUL(r17, vp2)); accu = VADD(accu, VAND(z19, vbmask)); z20 = VADD(z20, VSHR(z19, HT_BRADIX)); 
  r1 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r3, vp17)); accu = VADD(accu, VMUL(r4, vp16)); accu = VADD(accu, VMUL(r5, vp15)); 
  accu = VADD(accu, VMUL(r6, vp14)); accu = VADD(accu, VMUL(r7, vp13)); accu = VADD(accu, VMUL(r8, vp12)); 
  accu = VADD(accu, VMUL(r9, vp11)); accu = VADD(accu, VMUL(r10, vp10)); accu = VADD(accu, VMUL(r11, vp9)); 
  accu = VADD(accu, VMUL(r12, vp8)); accu = VADD(accu, VMUL(r13, vp7)); accu = VADD(accu, VMUL(r14, vp6)); 
  accu = VADD(accu, VMUL(r15, vp5)); accu = VADD(accu, VMUL(r16, vp4)); accu = VADD(accu, VMUL(r17, vp3)); 
  accu = VADD(accu, VAND(z20, vbmask)); z21 = VADD(z21, VSHR(z20, HT_BRADIX)); 
  r2 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r4, vp17)); accu = VADD(accu, VMUL(r5, vp16)); accu = VADD(accu, VMUL(r6, vp15)); 
  accu = VADD(accu, VMUL(r7, vp14)); accu = VADD(accu, VMUL(r8, vp13)); accu = VADD(accu, VMUL(r9, vp12)); 
  accu = VADD(accu, VMUL(r10, vp11)); accu = VADD(accu, VMUL(r11, vp10)); accu = VADD(accu, VMUL(r12, vp9)); 
  accu = VADD(accu, VMUL(r13, vp8)); accu = VADD(accu, VMUL(r14, vp7)); accu = VADD(accu, VMUL(r15, vp6)); 
  accu = VADD(accu, VMUL(r16, vp5)); accu = VADD(accu, VMUL(r17, vp4)); accu = VADD(accu, VAND(z21, vbmask)); 
  z22 = VADD(z22, VSHR(z21, HT_BRADIX)); 
  r3 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r5, vp17)); accu = VADD(accu, VMUL(r6, vp16)); accu = VADD(accu, VMUL(r7, vp15)); 
  accu = VADD(accu, VMUL(r8, vp14)); accu = VADD(accu, VMUL(r9, vp13)); accu = VADD(accu, VMUL(r10, vp12)); 
  accu = VADD(accu, VMUL(r11, vp11)); accu = VADD(accu, VMUL(r12, vp10)); accu = VADD(accu, VMUL(r13, vp9)); 
  accu = VADD(accu, VMUL(r14, vp8)); accu = VADD(accu, VMUL(r15, vp7)); accu = VADD(accu, VMUL(r16, vp6)); 
  accu = VADD(accu, VMUL(r17, vp5)); accu = VADD(accu, VAND(z22, vbmask)); z23 = VADD(z23, VSHR(z22, HT_BRADIX)); 
  r4 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r6, vp17)); accu = VADD(accu, VMUL(r7, vp16)); accu = VADD(accu, VMUL(r8, vp15)); 
  accu = VADD(accu, VMUL(r9, vp14)); accu = VADD(accu, VMUL(r10, vp13)); accu = VADD(accu, VMUL(r11, vp12)); 
  accu = VADD(accu, VMUL(r12, vp11)); accu = VADD(accu, VMUL(r13, vp10)); accu = VADD(accu, VMUL(r14, vp9)); 
  accu = VADD(accu, VMUL(r15, vp8)); accu = VADD(accu, VMUL(r16, vp7)); accu = VADD(accu, VMUL(r17, vp6)); 
  accu = VADD(accu, VAND(z23, vbmask)); z24 = VADD(z24, VSHR(z23, HT_BRADIX)); 
  r5 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r7, vp17)); accu = VADD(accu, VMUL(r8, vp16)); accu = VADD(accu, VMUL(r9, vp15)); 
  accu = VADD(accu, VMUL(r10, vp14)); accu = VADD(accu, VMUL(r11, vp13)); accu = VADD(accu, VMUL(r12, vp12)); 
  accu = VADD(accu, VMUL(r13, vp11)); accu = VADD(accu, VMUL(r14, vp10)); accu = VADD(accu, VMUL(r15, vp9)); 
  accu = VADD(accu, VMUL(r16, vp8)); accu = VADD(accu, VMUL(r17, vp7)); accu = VADD(accu, VAND(z24, vbmask)); 
  z25 = VADD(z25, VSHR(z24, HT_BRADIX)); 
  r6 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r8, vp17)); accu = VADD(accu, VMUL(r9, vp16)); accu = VADD(accu, VMUL(r10, vp15)); 
  accu = VADD(accu, VMUL(r11, vp14)); accu = VADD(accu, VMUL(r12, vp13)); accu = VADD(accu, VMUL(r13, vp12)); 
  accu = VADD(accu, VMUL(r14, vp11)); accu = VADD(accu, VMUL(r15, vp10)); accu = VADD(accu, VMUL(r16, vp9)); 
  accu = VADD(accu, VMUL(r17, vp8)); accu = VADD(accu, VAND(z25, vbmask)); z26 = VADD(z26, VSHR(z25, HT_BRADIX)); 
  r7 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r9, vp17)); accu = VADD(accu, VMUL(r10, vp16)); accu = VADD(accu, VMUL(r11, vp15)); 
  accu = VADD(accu, VMUL(r12, vp14)); accu = VADD(accu, VMUL(r13, vp13)); accu = VADD(accu, VMUL(r14, vp12)); 
  accu = VADD(accu, VMUL(r15, vp11)); accu = VADD(accu, VMUL(r16, vp10)); accu = VADD(accu, VMUL(r17, vp9)); 
  accu = VADD(accu, VAND(z26, vbmask)); z27 = VADD(z27, VSHR(z26, HT_BRADIX)); 
  r8 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r10, vp17)); accu = VADD(accu, VMUL(r11, vp16)); accu = VADD(accu, VMUL(r12, vp15)); 
  accu = VADD(accu, VMUL(r13, vp14)); accu = VADD(accu, VMUL(r14, vp13)); accu = VADD(accu, VMUL(r15, vp12)); 
  accu = VADD(accu, VMUL(r16, vp11)); accu = VADD(accu, VMUL(r17, vp10)); accu = VADD(accu, VAND(z27, vbmask)); 
  z28 = VADD(z28, VSHR(z27, HT_BRADIX)); 
  r9 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r11, vp17)); accu = VADD(accu, VMUL(r12, vp16)); accu = VADD(accu, VMUL(r13, vp15)); 
  accu = VADD(accu, VMUL(r14, vp14)); accu = VADD(accu, VMUL(r15, vp13)); accu = VADD(accu, VMUL(r16, vp12)); 
  accu = VADD(accu, VMUL(r17, vp11)); accu = VADD(accu, VAND(z28, vbmask)); z29 = VADD(z29, VSHR(z28, HT_BRADIX)); 
  r10 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r12, vp17)); accu = VADD(accu, VMUL(r13, vp16)); accu = VADD(accu, VMUL(r14, vp15)); 
  accu = VADD(accu, VMUL(r15, vp14)); accu = VADD(accu, VMUL(r16, vp13)); accu = VADD(accu, VMUL(r17, vp12)); 
  accu = VADD(accu, VAND(z29, vbmask)); z30 = VADD(z30, VSHR(z29, HT_BRADIX)); 
  r11 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r13, vp17)); accu = VADD(accu, VMUL(r14, vp16)); accu = VADD(accu, VMUL(r15, vp15)); 
  accu = VADD(accu, VMUL(r16, vp14)); accu = VADD(accu, VMUL(r17, vp13)); accu = VADD(accu, VAND(z30, vbmask)); 
  z31 = VADD(z31, VSHR(z30, HT_BRADIX)); 
  r12 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r14, vp17)); accu = VADD(accu, VMUL(r15, vp16)); accu = VADD(accu, VMUL(r16, vp15)); 
  accu = VADD(accu, VMUL(r17, vp14)); accu = VADD(accu, VAND(z31, vbmask)); z32 = VADD(z32, VSHR(z31, HT_BRADIX)); 
  r13 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r15, vp17)); accu = VADD(accu, VMUL(r16, vp16)); accu = VADD(accu, VMUL(r17, vp15)); 
  accu = VADD(accu, VAND(z32, vbmask)); z33 = VADD(z33, VSHR(z32, HT_BRADIX)); 
  r14 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r16, vp17)); accu = VADD(accu, VMUL(r17, vp16)); accu = VADD(accu, VAND(z33, vbmask)); 
  z34 = VADD(z34, VSHR(z33, HT_BRADIX)); 
  r15 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r17, vp17)); accu = VADD(accu, VAND(z34, vbmask)); z35 = VADD(z35, VSHR(z34, HT_BRADIX)); 
  r16 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  r17 = VADD(accu, z35);

  // ---------------------------------------------------------------------------

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// Montgomery squaring r = a^2 mod 2p
// squaring (product-scanning) interleaved with reduction (product-scanning)
// -> r in [0, 2p)
void gfp_sqr_8x1w(htfe_t r, const htfe_t a)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i z0,  z1,  z2,  z3,  z4,  z5,  z6,  z7,  z8,  z9,  z10, z11;
  __m512i z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23;
  __m512i z24, z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35;
  __m512i r0, r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8;
  __m512i r9, r10, r11, r12, r13, r14, r15, r16, r17, accu = VZERO;
  const __m512i vp0  = VSET1(ht_p[0]),  vp1  = VSET1(ht_p[1]),  vp2  = VSET1(ht_p[2]);
  const __m512i vp3  = VSET1(ht_p[3]),  vp4  = VSET1(ht_p[4]),  vp5  = VSET1(ht_p[5]);
  const __m512i vp6  = VSET1(ht_p[6]),  vp7  = VSET1(ht_p[7]),  vp8  = VSET1(ht_p[8]);
  const __m512i vp9  = VSET1(ht_p[9]),  vp10 = VSET1(ht_p[10]), vp11 = VSET1(ht_p[11]);
  const __m512i vp12 = VSET1(ht_p[12]), vp13 = VSET1(ht_p[13]), vp14 = VSET1(ht_p[14]);
  const __m512i vp15 = VSET1(ht_p[15]), vp16 = VSET1(ht_p[16]), vp17 = VSET1(ht_p[17]);
  const __m512i vbmask = VSET1(HT_BMASK), vw = VSET1(HT_MONTW);

  // ---------------------------------------------------------------------------
  // 1st loop of integer squaring

  z0 = VMUL(a0, a0); 

  z1 = VMUL(a0, a1); 
  z1 = VADD(z1, z1);

  z2 = VMUL(a0, a2); 
  z2 = VADD(z2, z2);
  z2 = VADD(z2, VMUL(a1, a1)); 

  z3 = VMUL(a0, a3); z3 = VADD(z3, VMUL(a1, a2)); 
  z3 = VADD(z3, z3);

  accu = VADD(accu, VAND(z0, vbmask)); z1 = VADD(z1, VSHR(z0, HT_BRADIX)); r0 = VAND(accu, vbmask); 
  r0 = VAND(vbmask, VMUL(r0, vw)); accu = VADD(accu, VMUL(r0, vp0)); accu = VSHR(accu, HT_BRADIX);

  z4 = VMUL(a0, a4); z4 = VADD(z4, VMUL(a1, a3)); 
  z4 = VADD(z4, z4);
  z4 = VADD(z4, VMUL(a2, a2)); 

  accu = VADD(accu, VMUL(r0, vp1)); accu = VADD(accu, VAND(z1, vbmask)); z2 = VADD(z2, VSHR(z1, HT_BRADIX)); 
  r1 = VAND(accu, vbmask); r1 = VAND(vbmask, VMUL(r1, vw)); accu = VADD(accu, VMUL(r1, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  z5 = VMUL(a0, a5); z5 = VADD(z5, VMUL(a1, a4)); z5 = VADD(z5, VMUL(a2, a3));
  z5 = VADD(z5, z5);

  accu = VADD(accu, VMUL(r0, vp2)); accu = VADD(accu, VMUL(r1, vp1)); accu = VADD(accu, VAND(z2, vbmask)); 
  z3 = VADD(z3, VSHR(z2, HT_BRADIX)); r2 = VAND(accu, vbmask); r2 = VAND(vbmask, VMUL(r2, vw)); 
  accu = VADD(accu, VMUL(r2, vp0)); accu = VSHR(accu, HT_BRADIX);

  z6 = VMUL(a0, a6); z6 = VADD(z6, VMUL(a1, a5)); z6 = VADD(z6, VMUL(a2, a4));
  z6 = VADD(z6, z6);
  z6 = VADD(z6, VMUL(a3, a3)); 

  accu = VADD(accu, VMUL(r0, vp3)); accu = VADD(accu, VMUL(r1, vp2)); accu = VADD(accu, VMUL(r2, vp1)); 
  accu = VADD(accu, VAND(z3, vbmask)); z4 = VADD(z4, VSHR(z3, HT_BRADIX)); r3 = VAND(accu, vbmask); 
  r3 = VAND(vbmask, VMUL(r3, vw)); accu = VADD(accu, VMUL(r3, vp0)); accu = VSHR(accu, HT_BRADIX);

  z7 = VMUL(a0, a7); z7 = VADD(z7, VMUL(a1, a6)); z7 = VADD(z7, VMUL(a2, a5));
  z7 = VADD(z7, VMUL(a3, a4));
  z7 = VADD(z7, z7);

  accu = VADD(accu, VMUL(r0, vp4)); accu = VADD(accu, VMUL(r1, vp3)); accu = VADD(accu, VMUL(r2, vp2)); 
  accu = VADD(accu, VMUL(r3, vp1)); accu = VADD(accu, VAND(z4, vbmask)); z5 = VADD(z5, VSHR(z4, HT_BRADIX)); 
  r4 = VAND(accu, vbmask); r4 = VAND(vbmask, VMUL(r4, vw)); accu = VADD(accu, VMUL(r4, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  z8 = VMUL(a0, a8); z8 = VADD(z8, VMUL(a1, a7)); z8 = VADD(z8, VMUL(a2, a6));
  z8 = VADD(z8, VMUL(a3, a5)); 
  z8 = VADD(z8, z8);
  z8 = VADD(z8, VMUL(a4, a4)); 

  accu = VADD(accu, VMUL(r0, vp5)); accu = VADD(accu, VMUL(r1, vp4)); accu = VADD(accu, VMUL(r2, vp3)); 
  accu = VADD(accu, VMUL(r3, vp2)); accu = VADD(accu, VMUL(r4, vp1)); accu = VADD(accu, VAND(z5, vbmask)); 
  z6 = VADD(z6, VSHR(z5, HT_BRADIX)); r5 = VAND(accu, vbmask); r5 = VAND(vbmask, VMUL(r5, vw)); 
  accu = VADD(accu, VMUL(r5, vp0)); accu = VSHR(accu, HT_BRADIX);

  z9 = VMUL(a0, a9); z9 = VADD(z9, VMUL(a1, a8)); z9 = VADD(z9, VMUL(a2, a7));
  z9 = VADD(z9, VMUL(a3, a6)); z9 = VADD(z9, VMUL(a4, a5)); 
  z9 = VADD(z9, z9);

  accu = VADD(accu, VMUL(r0, vp6)); accu = VADD(accu, VMUL(r1, vp5)); accu = VADD(accu, VMUL(r2, vp4)); 
  accu = VADD(accu, VMUL(r3, vp3)); accu = VADD(accu, VMUL(r4, vp2)); accu = VADD(accu, VMUL(r5, vp1)); 
  accu = VADD(accu, VAND(z6, vbmask)); z7 = VADD(z7, VSHR(z6, HT_BRADIX)); r6 = VAND(accu, vbmask); 
  r6 = VAND(vbmask, VMUL(r6, vw)); accu = VADD(accu, VMUL(r6, vp0)); accu = VSHR(accu, HT_BRADIX);

  z10 = VMUL(a0, a10); z10 = VADD(z10, VMUL(a1, a9)); z10 = VADD(z10, VMUL(a2, a8));
  z10 = VADD(z10, VMUL(a3, a7)); z10 = VADD(z10, VMUL(a4, a6)); 
  z10 = VADD(z10, z10);
  z10 = VADD(z10, VMUL(a5, a5));

  accu = VADD(accu, VMUL(r0, vp7)); accu = VADD(accu, VMUL(r1, vp6)); accu = VADD(accu, VMUL(r2, vp5)); 
  accu = VADD(accu, VMUL(r3, vp4)); accu = VADD(accu, VMUL(r4, vp3)); accu = VADD(accu, VMUL(r5, vp2)); 
  accu = VADD(accu, VMUL(r6, vp1)); accu = VADD(accu, VAND(z7, vbmask)); z8 = VADD(z8, VSHR(z7, HT_BRADIX)); 
  r7 = VAND(accu, vbmask); r7 = VAND(vbmask, VMUL(r7, vw)); accu = VADD(accu, VMUL(r7, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  z11 = VMUL(a0, a11); z11 = VADD(z11, VMUL(a1, a10)); z11 = VADD(z11, VMUL(a2, a9));
  z11 = VADD(z11, VMUL(a3, a8)); z11 = VADD(z11, VMUL(a4, a7)); z11 = VADD(z11, VMUL(a5, a6));
  z11 = VADD(z11, z11);

  accu = VADD(accu, VMUL(r0, vp8)); accu = VADD(accu, VMUL(r1, vp7)); accu = VADD(accu, VMUL(r2, vp6)); 
  accu = VADD(accu, VMUL(r3, vp5)); accu = VADD(accu, VMUL(r4, vp4)); accu = VADD(accu, VMUL(r5, vp3)); 
  accu = VADD(accu, VMUL(r6, vp2)); accu = VADD(accu, VMUL(r7, vp1)); accu = VADD(accu, VAND(z8, vbmask)); 
  z9 = VADD(z9, VSHR(z8, HT_BRADIX)); r8 = VAND(accu, vbmask); r8 = VAND(vbmask, VMUL(r8, vw)); 
  accu = VADD(accu, VMUL(r8, vp0)); accu = VSHR(accu, HT_BRADIX);

  z12 = VMUL(a0, a12); z12 = VADD(z12, VMUL(a1, a11)); z12 = VADD(z12, VMUL(a2, a10));
  z12 = VADD(z12, VMUL(a3, a9)); z12 = VADD(z12, VMUL(a4, a8)); z12 = VADD(z12, VMUL(a5, a7));
  z12 = VADD(z12, z12);
  z12 = VADD(z12, VMUL(a6, a6));

  accu = VADD(accu, VMUL(r0, vp9)); accu = VADD(accu, VMUL(r1, vp8)); accu = VADD(accu, VMUL(r2, vp7)); 
  accu = VADD(accu, VMUL(r3, vp6)); accu = VADD(accu, VMUL(r4, vp5)); accu = VADD(accu, VMUL(r5, vp4)); 
  accu = VADD(accu, VMUL(r6, vp3)); accu = VADD(accu, VMUL(r7, vp2)); accu = VADD(accu, VMUL(r8, vp1)); 
  accu = VADD(accu, VAND(z9, vbmask)); z10 = VADD(z10, VSHR(z9, HT_BRADIX)); r9 = VAND(accu, vbmask); 
  r9 = VAND(vbmask, VMUL(r9, vw)); accu = VADD(accu, VMUL(r9, vp0)); accu = VSHR(accu, HT_BRADIX);

  z13 = VMUL(a0, a13); z13 = VADD(z13, VMUL(a1, a12)); z13 = VADD(z13, VMUL(a2, a11));
  z13 = VADD(z13, VMUL(a3, a10)); z13 = VADD(z13, VMUL(a4, a9)); z13 = VADD(z13, VMUL(a5, a8));
  z13 = VADD(z13, VMUL(a6, a7));
  z13 = VADD(z13, z13); 

  accu = VADD(accu, VMUL(r0, vp10)); accu = VADD(accu, VMUL(r1, vp9)); accu = VADD(accu, VMUL(r2, vp8)); 
  accu = VADD(accu, VMUL(r3, vp7)); accu = VADD(accu, VMUL(r4, vp6)); accu = VADD(accu, VMUL(r5, vp5)); 
  accu = VADD(accu, VMUL(r6, vp4)); accu = VADD(accu, VMUL(r7, vp3)); accu = VADD(accu, VMUL(r8, vp2)); 
  accu = VADD(accu, VMUL(r9, vp1)); accu = VADD(accu, VAND(z10, vbmask)); z11 = VADD(z11, VSHR(z10, HT_BRADIX)); 
  r10 = VAND(accu, vbmask); r10 = VAND(vbmask, VMUL(r10, vw)); accu = VADD(accu, VMUL(r10, vp0)); 
  accu = VSHR(accu, HT_BRADIX);

  z14 = VMUL(a0, a14); z14 = VADD(z14, VMUL(a1, a13)); z14 = VADD(z14, VMUL(a2, a12));
  z14 = VADD(z14, VMUL(a3, a11)); z14 = VADD(z14, VMUL(a4, a10)); z14 = VADD(z14, VMUL(a5, a9));
  z14 = VADD(z14, VMUL(a6, a8)); 
  z14 = VADD(z14, z14);
  z14 = VADD(z14, VMUL(a7, a7)); 

  accu = VADD(accu, VMUL(r0, vp11)); accu = VADD(accu, VMUL(r1, vp10)); accu = VADD(accu, VMUL(r2, vp9)); 
  accu = VADD(accu, VMUL(r3, vp8)); accu = VADD(accu, VMUL(r4, vp7)); accu = VADD(accu, VMUL(r5, vp6)); 
  accu = VADD(accu, VMUL(r6, vp5)); accu = VADD(accu, VMUL(r7, vp4)); accu = VADD(accu, VMUL(r8, vp3)); 
  accu = VADD(accu, VMUL(r9, vp2)); accu = VADD(accu, VMUL(r10, vp1)); accu = VADD(accu, VAND(z11, vbmask)); 
  z12 = VADD(z12, VSHR(z11, HT_BRADIX)); r11 = VAND(accu, vbmask); r11 = VAND(vbmask, VMUL(r11, vw)); 
  accu = VADD(accu, VMUL(r11, vp0)); accu = VSHR(accu, HT_BRADIX);

  z15 = VMUL(a0, a15); z15 = VADD(z15, VMUL(a1, a14)); z15 = VADD(z15, VMUL(a2, a13)); 
  z15 = VADD(z15, VMUL(a3, a12)); z15 = VADD(z15, VMUL(a4, a11)); z15 = VADD(z15, VMUL(a5, a10)); 
  z15 = VADD(z15, VMUL(a6, a9)); z15 = VADD(z15, VMUL(a7, a8)); 
  z15 = VADD(z15, z15); 

  accu = VADD(accu, VMUL(r0, vp12)); accu = VADD(accu, VMUL(r1, vp11)); accu = VADD(accu, VMUL(r2, vp10)); 
  accu = VADD(accu, VMUL(r3, vp9)); accu = VADD(accu, VMUL(r4, vp8)); accu = VADD(accu, VMUL(r5, vp7)); 
  accu = VADD(accu, VMUL(r6, vp6)); accu = VADD(accu, VMUL(r7, vp5)); accu = VADD(accu, VMUL(r8, vp4)); 
  accu = VADD(accu, VMUL(r9, vp3)); accu = VADD(accu, VMUL(r10, vp2)); accu = VADD(accu, VMUL(r11, vp1)); 
  accu = VADD(accu, VAND(z12, vbmask)); z13 = VADD(z13, VSHR(z12, HT_BRADIX)); r12 = VAND(accu, vbmask); 
  r12 = VAND(vbmask, VMUL(r12, vw)); accu = VADD(accu, VMUL(r12, vp0)); accu = VSHR(accu, HT_BRADIX);

  z16 = VMUL(a0, a16); z16 = VADD(z16, VMUL(a1, a15)); z16 = VADD(z16, VMUL(a2, a14)); 
  z16 = VADD(z16, VMUL(a3, a13)); z16 = VADD(z16, VMUL(a4, a12)); z16 = VADD(z16, VMUL(a5, a11)); 
  z16 = VADD(z16, VMUL(a6, a10)); z16 = VADD(z16, VMUL(a7, a9)); 
  z16 = VADD(z16, z16);
  z16 = VADD(z16, VMUL(a8, a8)); 

  accu = VADD(accu, VMUL(r0, vp13)); accu = VADD(accu, VMUL(r1, vp12)); accu = VADD(accu, VMUL(r2, vp11)); 
  accu = VADD(accu, VMUL(r3, vp10)); accu = VADD(accu, VMUL(r4, vp9)); accu = VADD(accu, VMUL(r5, vp8)); 
  accu = VADD(accu, VMUL(r6, vp7)); accu = VADD(accu, VMUL(r7, vp6)); accu = VADD(accu, VMUL(r8, vp5)); 
  accu = VADD(accu, VMUL(r9, vp4)); accu = VADD(accu, VMUL(r10, vp3)); accu = VADD(accu, VMUL(r11, vp2)); 
  accu = VADD(accu, VMUL(r12, vp1)); accu = VADD(accu, VAND(z13, vbmask)); z14 = VADD(z14, VSHR(z13, HT_BRADIX)); 
  r13 = VAND(accu, vbmask); r13 = VAND(vbmask, VMUL(r13, vw)); accu = VADD(accu, VMUL(r13, vp0)); 
  accu = VSHR(accu, HT_BRADIX);
  
  z17 = VMUL(a0, a17); z17 = VADD(z17, VMUL(a1, a16)); z17 = VADD(z17, VMUL(a2, a15)); 
  z17 = VADD(z17, VMUL(a3, a14)); z17 = VADD(z17, VMUL(a4, a13)); z17 = VADD(z17, VMUL(a5, a12)); 
  z17 = VADD(z17, VMUL(a6, a11)); z17 = VADD(z17, VMUL(a7, a10)); z17 = VADD(z17, VMUL(a8, a9)); 
  z17 = VADD(z17, z17); 

  accu = VADD(accu, VMUL(r0, vp14)); accu = VADD(accu, VMUL(r1, vp13)); accu = VADD(accu, VMUL(r2, vp12)); 
  accu = VADD(accu, VMUL(r3, vp11)); accu = VADD(accu, VMUL(r4, vp10)); accu = VADD(accu, VMUL(r5, vp9)); 
  accu = VADD(accu, VMUL(r6, vp8)); accu = VADD(accu, VMUL(r7, vp7)); accu = VADD(accu, VMUL(r8, vp6)); 
  accu = VADD(accu, VMUL(r9, vp5)); accu = VADD(accu, VMUL(r10, vp4)); accu = VADD(accu, VMUL(r11, vp3)); 
  accu = VADD(accu, VMUL(r12, vp2)); accu = VADD(accu, VMUL(r13, vp1)); accu = VADD(accu, VAND(z14, vbmask)); 
  z15 = VADD(z15, VSHR(z14, HT_BRADIX)); r14 = VAND(accu, vbmask); r14 = VAND(vbmask, VMUL(r14, vw)); 
  accu = VADD(accu, VMUL(r14, vp0)); accu = VSHR(accu, HT_BRADIX);

  // ---------------------------------------------------------------------------
  // 2nd loop of integer squaring + Montgomery reduction

  z18 = VMUL(a1, a17); z18 = VADD(z18, VMUL(a2, a16)); z18 = VADD(z18, VMUL(a3, a15)); 
  z18 = VADD(z18, VMUL(a4, a14)); z18 = VADD(z18, VMUL(a5, a13)); z18 = VADD(z18, VMUL(a6, a12)); 
  z18 = VADD(z18, VMUL(a7, a11)); z18 = VADD(z18, VMUL(a8, a10)); 
  z18 = VADD(z18, z18);
  z18 = VADD(z18, VMUL(a9, a9)); 

  accu = VADD(accu, VMUL(r0, vp15)); accu = VADD(accu, VMUL(r1, vp14)); accu = VADD(accu, VMUL(r2, vp13)); 
  accu = VADD(accu, VMUL(r3, vp12)); accu = VADD(accu, VMUL(r4, vp11)); accu = VADD(accu, VMUL(r5, vp10)); 
  accu = VADD(accu, VMUL(r6, vp9)); accu = VADD(accu, VMUL(r7, vp8)); accu = VADD(accu, VMUL(r8, vp7)); 
  accu = VADD(accu, VMUL(r9, vp6)); accu = VADD(accu, VMUL(r10, vp5)); accu = VADD(accu, VMUL(r11, vp4)); 
  accu = VADD(accu, VMUL(r12, vp3)); accu = VADD(accu, VMUL(r13, vp2)); accu = VADD(accu, VMUL(r14, vp1)); 
  accu = VADD(accu, VAND(z15, vbmask)); z16 = VADD(z16, VSHR(z15, HT_BRADIX)); r15 = VAND(accu, vbmask); 
  r15 = VAND(vbmask, VMUL(r15, vw)); accu = VADD(accu, VMUL(r15, vp0)); accu = VSHR(accu, HT_BRADIX);

  z19 = VMUL(a2, a17); z19 = VADD(z19, VMUL(a3, a16)); z19 = VADD(z19, VMUL(a4, a15)); 
  z19 = VADD(z19, VMUL(a5, a14)); z19 = VADD(z19, VMUL(a6, a13)); z19 = VADD(z19, VMUL(a7, a12)); 
  z19 = VADD(z19, VMUL(a8, a11)); z19 = VADD(z19, VMUL(a9, a10)); 
  z19 = VADD(z19, z19);

  accu = VADD(accu, VMUL(r0, vp16)); accu = VADD(accu, VMUL(r1, vp15)); accu = VADD(accu, VMUL(r2, vp14)); 
  accu = VADD(accu, VMUL(r3, vp13)); accu = VADD(accu, VMUL(r4, vp12)); accu = VADD(accu, VMUL(r5, vp11)); 
  accu = VADD(accu, VMUL(r6, vp10)); accu = VADD(accu, VMUL(r7, vp9)); accu = VADD(accu, VMUL(r8, vp8)); 
  accu = VADD(accu, VMUL(r9, vp7)); accu = VADD(accu, VMUL(r10, vp6)); accu = VADD(accu, VMUL(r11, vp5)); 
  accu = VADD(accu, VMUL(r12, vp4)); accu = VADD(accu, VMUL(r13, vp3)); accu = VADD(accu, VMUL(r14, vp2)); 
  accu = VADD(accu, VMUL(r15, vp1)); accu = VADD(accu, VAND(z16, vbmask)); z17 = VADD(z17, VSHR(z16, HT_BRADIX)); 
  r16 = VAND(accu, vbmask); r16 = VAND(vbmask, VMUL(r16, vw)); accu = VADD(accu, VMUL(r16, vp0)); 
  accu = VSHR(accu, HT_BRADIX);


  z20 = VMUL(a3, a17); z20 = VADD(z20, VMUL(a4, a16)); z20 = VADD(z20, VMUL(a5, a15)); 
  z20 = VADD(z20, VMUL(a6, a14)); z20 = VADD(z20, VMUL(a7, a13)); z20 = VADD(z20, VMUL(a8, a12)); 
  z20 = VADD(z20, VMUL(a9, a11)); 
  z20 = VADD(z20, z20);
  z20 = VADD(z20, VMUL(a10, a10)); 

  accu = VADD(accu, VMUL(r0, vp17)); accu = VADD(accu, VMUL(r1, vp16)); accu = VADD(accu, VMUL(r2, vp15)); 
  accu = VADD(accu, VMUL(r3, vp14)); accu = VADD(accu, VMUL(r4, vp13)); accu = VADD(accu, VMUL(r5, vp12)); 
  accu = VADD(accu, VMUL(r6, vp11)); accu = VADD(accu, VMUL(r7, vp10)); accu = VADD(accu, VMUL(r8, vp9)); 
  accu = VADD(accu, VMUL(r9, vp8)); accu = VADD(accu, VMUL(r10, vp7)); accu = VADD(accu, VMUL(r11, vp6)); 
  accu = VADD(accu, VMUL(r12, vp5)); accu = VADD(accu, VMUL(r13, vp4)); accu = VADD(accu, VMUL(r14, vp3)); 
  accu = VADD(accu, VMUL(r15, vp2)); accu = VADD(accu, VMUL(r16, vp1)); accu = VADD(accu, VAND(z17, vbmask)); 
  z18 = VADD(z18, VSHR(z17, HT_BRADIX)); r17 = VAND(accu, vbmask); r17 = VAND(vbmask, VMUL(r17, vw)); 
  accu = VADD(accu, VMUL(r17, vp0)); accu = VSHR(accu, HT_BRADIX);

  z21 = VMUL(a4, a17); z21 = VADD(z21, VMUL(a5, a16)); z21 = VADD(z21, VMUL(a6, a15)); 
  z21 = VADD(z21, VMUL(a7, a14)); z21 = VADD(z21, VMUL(a8, a13)); z21 = VADD(z21, VMUL(a9, a12)); 
  z21 = VADD(z21, VMUL(a10, a11)); 
  z21 = VADD(z21, z21);

  accu = VADD(accu, VMUL(r1, vp17)); accu = VADD(accu, VMUL(r2, vp16)); accu = VADD(accu, VMUL(r3, vp15)); 
  accu = VADD(accu, VMUL(r4, vp14)); accu = VADD(accu, VMUL(r5, vp13)); accu = VADD(accu, VMUL(r6, vp12)); 
  accu = VADD(accu, VMUL(r7, vp11)); accu = VADD(accu, VMUL(r8, vp10)); accu = VADD(accu, VMUL(r9, vp9)); 
  accu = VADD(accu, VMUL(r10, vp8)); accu = VADD(accu, VMUL(r11, vp7)); accu = VADD(accu, VMUL(r12, vp6)); 
  accu = VADD(accu, VMUL(r13, vp5)); accu = VADD(accu, VMUL(r14, vp4)); accu = VADD(accu, VMUL(r15, vp3)); 
  accu = VADD(accu, VMUL(r16, vp2)); accu = VADD(accu, VMUL(r17, vp1)); accu = VADD(accu, VAND(z18, vbmask)); 
  z19 = VADD(z19, VSHR(z18, HT_BRADIX)); 
  r0 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z22 = VMUL(a5, a17); z22 = VADD(z22, VMUL(a6, a16)); z22 = VADD(z22, VMUL(a7, a15)); 
  z22 = VADD(z22, VMUL(a8, a14)); z22 = VADD(z22, VMUL(a9, a13)); z22 = VADD(z22, VMUL(a10, a12)); 
  z22 = VADD(z22, z22);
  z22 = VADD(z22, VMUL(a11, a11)); 

  accu = VADD(accu, VMUL(r2, vp17)); accu = VADD(accu, VMUL(r3, vp16)); accu = VADD(accu, VMUL(r4, vp15)); 
  accu = VADD(accu, VMUL(r5, vp14)); accu = VADD(accu, VMUL(r6, vp13)); accu = VADD(accu, VMUL(r7, vp12)); 
  accu = VADD(accu, VMUL(r8, vp11)); accu = VADD(accu, VMUL(r9, vp10)); accu = VADD(accu, VMUL(r10, vp9)); 
  accu = VADD(accu, VMUL(r11, vp8)); accu = VADD(accu, VMUL(r12, vp7)); accu = VADD(accu, VMUL(r13, vp6)); 
  accu = VADD(accu, VMUL(r14, vp5)); accu = VADD(accu, VMUL(r15, vp4)); accu = VADD(accu, VMUL(r16, vp3)); 
  accu = VADD(accu, VMUL(r17, vp2)); accu = VADD(accu, VAND(z19, vbmask)); z20 = VADD(z20, VSHR(z19, HT_BRADIX)); 
  r1 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z23 = VMUL(a6, a17); z23 = VADD(z23, VMUL(a7, a16)); z23 = VADD(z23, VMUL(a8, a15)); 
  z23 = VADD(z23, VMUL(a9, a14)); z23 = VADD(z23, VMUL(a10, a13)); z23 = VADD(z23, VMUL(a11, a12)); 
  z23 = VADD(z23, z23);

  accu = VADD(accu, VMUL(r3, vp17)); accu = VADD(accu, VMUL(r4, vp16)); accu = VADD(accu, VMUL(r5, vp15)); 
  accu = VADD(accu, VMUL(r6, vp14)); accu = VADD(accu, VMUL(r7, vp13)); accu = VADD(accu, VMUL(r8, vp12)); 
  accu = VADD(accu, VMUL(r9, vp11)); accu = VADD(accu, VMUL(r10, vp10)); accu = VADD(accu, VMUL(r11, vp9)); 
  accu = VADD(accu, VMUL(r12, vp8)); accu = VADD(accu, VMUL(r13, vp7)); accu = VADD(accu, VMUL(r14, vp6)); 
  accu = VADD(accu, VMUL(r15, vp5)); accu = VADD(accu, VMUL(r16, vp4)); accu = VADD(accu, VMUL(r17, vp3)); 
  accu = VADD(accu, VAND(z20, vbmask)); z21 = VADD(z21, VSHR(z20, HT_BRADIX)); 
  r2 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z24 = VMUL(a7, a17); z24 = VADD(z24, VMUL(a8, a16)); z24 = VADD(z24, VMUL(a9, a15)); 
  z24 = VADD(z24, VMUL(a10, a14)); z24 = VADD(z24, VMUL(a11, a13)); 
  z24 = VADD(z24, z24);
  z24 = VADD(z24, VMUL(a12, a12)); 


  accu = VADD(accu, VMUL(r4, vp17)); accu = VADD(accu, VMUL(r5, vp16)); accu = VADD(accu, VMUL(r6, vp15)); 
  accu = VADD(accu, VMUL(r7, vp14)); accu = VADD(accu, VMUL(r8, vp13)); accu = VADD(accu, VMUL(r9, vp12)); 
  accu = VADD(accu, VMUL(r10, vp11)); accu = VADD(accu, VMUL(r11, vp10)); accu = VADD(accu, VMUL(r12, vp9)); 
  accu = VADD(accu, VMUL(r13, vp8)); accu = VADD(accu, VMUL(r14, vp7)); accu = VADD(accu, VMUL(r15, vp6)); 
  accu = VADD(accu, VMUL(r16, vp5)); accu = VADD(accu, VMUL(r17, vp4)); accu = VADD(accu, VAND(z21, vbmask)); 
  z22 = VADD(z22, VSHR(z21, HT_BRADIX)); 
  r3 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z25 = VMUL(a8, a17); z25 = VADD(z25, VMUL(a9, a16)); z25 = VADD(z25, VMUL(a10, a15)); 
  z25 = VADD(z25, VMUL(a11, a14)); z25 = VADD(z25, VMUL(a12, a13)); 
  z25 = VADD(z25, z25);

  accu = VADD(accu, VMUL(r5, vp17)); accu = VADD(accu, VMUL(r6, vp16)); accu = VADD(accu, VMUL(r7, vp15)); 
  accu = VADD(accu, VMUL(r8, vp14)); accu = VADD(accu, VMUL(r9, vp13)); accu = VADD(accu, VMUL(r10, vp12)); 
  accu = VADD(accu, VMUL(r11, vp11)); accu = VADD(accu, VMUL(r12, vp10)); accu = VADD(accu, VMUL(r13, vp9)); 
  accu = VADD(accu, VMUL(r14, vp8)); accu = VADD(accu, VMUL(r15, vp7)); accu = VADD(accu, VMUL(r16, vp6)); 
  accu = VADD(accu, VMUL(r17, vp5)); accu = VADD(accu, VAND(z22, vbmask)); z23 = VADD(z23, VSHR(z22, HT_BRADIX)); 
  r4 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z26 = VMUL(a9, a17); z26 = VADD(z26, VMUL(a10, a16)); z26 = VADD(z26, VMUL(a11, a15)); 
  z26 = VADD(z26, VMUL(a12, a14)); 
  z26 = VADD(z26, z26);
  z26 = VADD(z26, VMUL(a13, a13)); 

  accu = VADD(accu, VMUL(r6, vp17)); accu = VADD(accu, VMUL(r7, vp16)); accu = VADD(accu, VMUL(r8, vp15)); 
  accu = VADD(accu, VMUL(r9, vp14)); accu = VADD(accu, VMUL(r10, vp13)); accu = VADD(accu, VMUL(r11, vp12)); 
  accu = VADD(accu, VMUL(r12, vp11)); accu = VADD(accu, VMUL(r13, vp10)); accu = VADD(accu, VMUL(r14, vp9)); 
  accu = VADD(accu, VMUL(r15, vp8)); accu = VADD(accu, VMUL(r16, vp7)); accu = VADD(accu, VMUL(r17, vp6)); 
  accu = VADD(accu, VAND(z23, vbmask)); z24 = VADD(z24, VSHR(z23, HT_BRADIX)); 
  r5 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z27 = VMUL(a10, a17); z27 = VADD(z27, VMUL(a11, a16)); z27 = VADD(z27, VMUL(a12, a15)); 
  z27 = VADD(z27, VMUL(a13, a14)); 
  z27 = VADD(z27, z27);

  accu = VADD(accu, VMUL(r7, vp17)); accu = VADD(accu, VMUL(r8, vp16)); accu = VADD(accu, VMUL(r9, vp15)); 
  accu = VADD(accu, VMUL(r10, vp14)); accu = VADD(accu, VMUL(r11, vp13)); accu = VADD(accu, VMUL(r12, vp12)); 
  accu = VADD(accu, VMUL(r13, vp11)); accu = VADD(accu, VMUL(r14, vp10)); accu = VADD(accu, VMUL(r15, vp9)); 
  accu = VADD(accu, VMUL(r16, vp8)); accu = VADD(accu, VMUL(r17, vp7)); accu = VADD(accu, VAND(z24, vbmask)); 
  z25 = VADD(z25, VSHR(z24, HT_BRADIX)); 
  r6 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z28 = VMUL(a11, a17); z28 = VADD(z28, VMUL(a12, a16)); z28 = VADD(z28, VMUL(a13, a15)); 
  z28 = VADD(z28, z28);
  z28 = VADD(z28, VMUL(a14, a14)); 

  accu = VADD(accu, VMUL(r8, vp17)); accu = VADD(accu, VMUL(r9, vp16)); accu = VADD(accu, VMUL(r10, vp15)); 
  accu = VADD(accu, VMUL(r11, vp14)); accu = VADD(accu, VMUL(r12, vp13)); accu = VADD(accu, VMUL(r13, vp12)); 
  accu = VADD(accu, VMUL(r14, vp11)); accu = VADD(accu, VMUL(r15, vp10)); accu = VADD(accu, VMUL(r16, vp9)); 
  accu = VADD(accu, VMUL(r17, vp8)); accu = VADD(accu, VAND(z25, vbmask)); z26 = VADD(z26, VSHR(z25, HT_BRADIX)); 
  r7 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z29 = VMUL(a12, a17); z29 = VADD(z29, VMUL(a13, a16)); z29 = VADD(z29, VMUL(a14, a15)); 
  z29 = VADD(z29, z29);

  accu = VADD(accu, VMUL(r9, vp17)); accu = VADD(accu, VMUL(r10, vp16)); accu = VADD(accu, VMUL(r11, vp15)); 
  accu = VADD(accu, VMUL(r12, vp14)); accu = VADD(accu, VMUL(r13, vp13)); accu = VADD(accu, VMUL(r14, vp12)); 
  accu = VADD(accu, VMUL(r15, vp11)); accu = VADD(accu, VMUL(r16, vp10)); accu = VADD(accu, VMUL(r17, vp9)); 
  accu = VADD(accu, VAND(z26, vbmask)); z27 = VADD(z27, VSHR(z26, HT_BRADIX)); 
  r8 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z30 = VMUL(a13, a17); z30 = VADD(z30, VMUL(a14, a16)); 
  z30 = VADD(z30, z30);
  z30 = VADD(z30, VMUL(a15, a15)); 

  accu = VADD(accu, VMUL(r10, vp17)); accu = VADD(accu, VMUL(r11, vp16)); accu = VADD(accu, VMUL(r12, vp15)); 
  accu = VADD(accu, VMUL(r13, vp14)); accu = VADD(accu, VMUL(r14, vp13)); accu = VADD(accu, VMUL(r15, vp12)); 
  accu = VADD(accu, VMUL(r16, vp11)); accu = VADD(accu, VMUL(r17, vp10)); accu = VADD(accu, VAND(z27, vbmask)); 
  z28 = VADD(z28, VSHR(z27, HT_BRADIX)); 
  r9 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z31 = VMUL(a14, a17); z31 = VADD(z31, VMUL(a15, a16)); 
  z31 = VADD(z31, z31);

  accu = VADD(accu, VMUL(r11, vp17)); accu = VADD(accu, VMUL(r12, vp16)); accu = VADD(accu, VMUL(r13, vp15)); 
  accu = VADD(accu, VMUL(r14, vp14)); accu = VADD(accu, VMUL(r15, vp13)); accu = VADD(accu, VMUL(r16, vp12)); 
  accu = VADD(accu, VMUL(r17, vp11)); accu = VADD(accu, VAND(z28, vbmask)); z29 = VADD(z29, VSHR(z28, HT_BRADIX)); 
  r10 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z32 = VMUL(a15, a17); 
  z32 = VADD(z32, z32);
  z32 = VADD(z32, VMUL(a16, a16)); 

  accu = VADD(accu, VMUL(r12, vp17)); accu = VADD(accu, VMUL(r13, vp16)); accu = VADD(accu, VMUL(r14, vp15)); 
  accu = VADD(accu, VMUL(r15, vp14)); accu = VADD(accu, VMUL(r16, vp13)); accu = VADD(accu, VMUL(r17, vp12)); 
  accu = VADD(accu, VAND(z29, vbmask)); z30 = VADD(z30, VSHR(z29, HT_BRADIX)); 
  r11 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z33 = VMUL(a16, a17); 
  z33 = VADD(z33, z33);

  accu = VADD(accu, VMUL(r13, vp17)); accu = VADD(accu, VMUL(r14, vp16)); accu = VADD(accu, VMUL(r15, vp15)); 
  accu = VADD(accu, VMUL(r16, vp14)); accu = VADD(accu, VMUL(r17, vp13)); accu = VADD(accu, VAND(z30, vbmask)); 
  z31 = VADD(z31, VSHR(z30, HT_BRADIX)); 
  r12 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  z34 = VMUL(a17, a17);

  z35 = VZERO;

  accu = VADD(accu, VMUL(r14, vp17)); accu = VADD(accu, VMUL(r15, vp16)); accu = VADD(accu, VMUL(r16, vp15)); 
  accu = VADD(accu, VMUL(r17, vp14)); accu = VADD(accu, VAND(z31, vbmask)); z32 = VADD(z32, VSHR(z31, HT_BRADIX)); 
  r13 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r15, vp17)); accu = VADD(accu, VMUL(r16, vp16)); accu = VADD(accu, VMUL(r17, vp15)); 
  accu = VADD(accu, VAND(z32, vbmask)); z33 = VADD(z33, VSHR(z32, HT_BRADIX)); 
  r14 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r16, vp17)); accu = VADD(accu, VMUL(r17, vp16)); accu = VADD(accu, VAND(z33, vbmask)); 
  z34 = VADD(z34, VSHR(z33, HT_BRADIX)); 
  r15 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  accu = VADD(accu, VMUL(r17, vp17)); accu = VADD(accu, VAND(z34, vbmask)); z35 = VADD(z35, VSHR(z34, HT_BRADIX)); 
  r16 = VAND(accu, vbmask); accu = VSHR(accu, HT_BRADIX);

  r17 = VADD(accu, z35);

  // ---------------------------------------------------------------------------

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// field exponentiation r = a^e mod 2p
// the exponent e is a *public* parameter and is the *same* for all 8 instances
// -> r in [0, 2p)
void gfp_pow_8x1w(htfe_t r, const htfe_t a, const uint64_t *e)
{
  htfe_t b;
  uint64_t t;
  int i, j;
  
  for (i = 0; i < HT_NWORDS; i++) {
    r[i] = VSET1(ht_montR[i]);
    b[i] = a[i];    
  }

  for (i = 0; i < 8; i++) {
    t = e[i];
    for (j = 0; j < 64; j++, t>>=1) {
      if (t&1) gfp_mul_8x1w(r, r, b);
      gfp_sqr_8x1w(b, b);
    }
  }
}

// field multiplicative inversion r = a^(-1) = a^(p-2) mod 2p
// -> r in [0, 2p)
void gfp_inv_8x1w(htfe_t r, const htfe_t a)
{
  gfp_pow_8x1w(r, a, u64_psub2);
}

// reduce the field element from [0, 2p) to [0, p)
// a in [0, 2p) -> r in [0, p)
void gfp_rdcp_8x1w(htfe_t r, const htfe_t a)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i r0, r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8;
  __m512i r9, r10, r11, r12, r13, r14, r15, r16, r17, smask;
  const __m512i vp0  = VSET1(ht_p[0]),  vp1  = VSET1(ht_p[1]);
  const __m512i vp2  = VSET1(ht_p[2]),  vp3  = VSET1(ht_p[3]);
  const __m512i vp4  = VSET1(ht_p[4]),  vp5  = VSET1(ht_p[5]);
  const __m512i vp6  = VSET1(ht_p[6]),  vp7  = VSET1(ht_p[7]);
  const __m512i vp8  = VSET1(ht_p[8]),  vp9  = VSET1(ht_p[9]);
  const __m512i vp10 = VSET1(ht_p[10]), vp11 = VSET1(ht_p[11]);
  const __m512i vp12 = VSET1(ht_p[12]), vp13 = VSET1(ht_p[13]);
  const __m512i vp14 = VSET1(ht_p[14]), vp15 = VSET1(ht_p[15]);
  const __m512i vp16 = VSET1(ht_p[16]), vp17 = VSET1(ht_p[17]);
  const __m512i vbmask = VSET1(HT_BMASK);

  // r = a - p
  r0  = VSUB(a0,  vp0);  r1  = VSUB(a1,  vp1);  r2  = VSUB(a2,  vp2);
  r3  = VSUB(a3,  vp3);  r4  = VSUB(a4,  vp4);  r5  = VSUB(a5,  vp5);
  r6  = VSUB(a6,  vp6);  r7  = VSUB(a7,  vp7);  r8  = VSUB(a8,  vp8);
  r9  = VSUB(a9,  vp9);  r10 = VSUB(a10, vp10); r11 = VSUB(a11, vp11);
  r12 = VSUB(a12, vp12); r13 = VSUB(a13, vp13); r14 = VSUB(a14, vp14);
  r15 = VSUB(a15, vp15); r16 = VSUB(a16, vp16); r17 = VSUB(a17, vp17);

  // check the current r is positive or negative 
  // carry propagation
  r1  = VADD(r1,  VSRA(r0,  HT_BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  HT_BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  HT_BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  HT_BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  HT_BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  HT_BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  HT_BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  HT_BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  HT_BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  HT_BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, HT_BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, HT_BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, HT_BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, HT_BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSRA(r14, HT_BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSRA(r15, HT_BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSRA(r16, HT_BRADIX)); r16 = VAND(r16, vbmask);

  // if r is positive, then the corresponding element in smask = 0; 
  // if r is negative, then the corresponding element in smask = all-1.
  smask = VSRA(r17, 63);
  // r = r + (2p & smask), add either 2p or 0 to the current r
  r0  = VADD(r0,  VAND(vp0,  smask)); r1  = VADD(r1,  VAND(vp1,  smask));
  r2  = VADD(r2,  VAND(vp2,  smask)); r3  = VADD(r3,  VAND(vp3,  smask));
  r4  = VADD(r4,  VAND(vp4,  smask)); r5  = VADD(r5,  VAND(vp5,  smask));
  r6  = VADD(r6,  VAND(vp6,  smask)); r7  = VADD(r7,  VAND(vp7,  smask));
  r8  = VADD(r8,  VAND(vp8,  smask)); r9  = VADD(r9,  VAND(vp9,  smask));
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask));
  r12 = VADD(r12, VAND(vp12, smask)); r13 = VADD(r13, VAND(vp13, smask));
  r14 = VADD(r14, VAND(vp14, smask)); r15 = VADD(r15, VAND(vp15, smask));
  r16 = VADD(r16, VAND(vp16, smask)); r17 = VADD(r17, VAND(vp17, smask));

  // carry propagation
  r1  = VADD(r1,  VSHR(r0,  HT_BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSHR(r1,  HT_BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSHR(r2,  HT_BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSHR(r3,  HT_BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSHR(r4,  HT_BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSHR(r5,  HT_BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSHR(r6,  HT_BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSHR(r7,  HT_BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSHR(r8,  HT_BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSHR(r9,  HT_BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSHR(r10, HT_BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSHR(r11, HT_BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSHR(r12, HT_BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSHR(r13, HT_BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSHR(r14, HT_BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSHR(r15, HT_BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSHR(r16, HT_BRADIX)); r16 = VAND(r16, vbmask);
  r17 = VAND(r17, vbmask);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// set the field element r to be 0
// -> r = 0
void gfp_zero_8x1w(htfe_t r)
{
  r[0]  = VZERO; r[1]  = VZERO; r[2]  = VZERO;
  r[3]  = VZERO; r[4]  = VZERO; r[5]  = VZERO;
  r[6]  = VZERO; r[7]  = VZERO; r[8]  = VZERO;
  r[9]  = VZERO; r[10] = VZERO; r[11] = VZERO;
  r[12] = VZERO; r[13] = VZERO; r[14] = VZERO;
  r[15] = VZERO; r[16] = VZERO; r[17] = VZERO;
}

// convert an integer from number domain to Montgomery domain r = a * R^2 mod 2p
// -> r in [0, 2p)
void gfp_num2mont_8x1w(htfe_t r, const htfe_t a)
{
  htfe_t vR2;
  int i;

  for (i = 0; i < HT_NWORDS; i++) vR2[i] = VSET1(ht_montR2[i]);
  gfp_mul_8x1w(r, a, vR2);
}

// convert an integer from Montgomery domain to number domain r = a * 1 mod 2p
// -> r in [0, p)
void gfp_mont2num_8x1w(htfe_t r, const htfe_t a)
{
  htfe_t vone;

  gfp_zero_8x1w(vone);
  vone[0] = VSET1(1);
  gfp_mul_8x1w(r, a, vone);
  gfp_rdcp_8x1w(r, r);
}

// copy field element a to r
void gfp_copy_8x1w(htfe_t r, const htfe_t a)
{
  r[0]  = a[0];  r[1]  = a[1];  r[2]  = a[2];
  r[3]  = a[3];  r[4]  = a[4];  r[5]  = a[5];
  r[6]  = a[6];  r[7]  = a[7];  r[8]  = a[8];
  r[9]  = a[9];  r[10] = a[10]; r[11] = a[11];
  r[12] = a[12]; r[13] = a[13]; r[14] = a[14];
  r[15] = a[15]; r[16] = a[16]; r[17] = a[17];
}

// conditional move
// move field element a to r if b == 1; do not move if b == 0
void gfp_cmove_8x1w(htfe_t r, htfe_t a, const __m512i b)
{
  __m512i r0  = r[0],  r1  = r[1],  r2  = r[2],  r3  = r[3],  r4  = r[4],  r5  = r[5];
  __m512i r6  = r[6],  r7  = r[7],  r8  = r[8],  r9  = r[9],  r10 = r[10], r11 = r[11];
  __m512i r12 = r[12], r13 = r[13], r14 = r[14], r15 = r[15], r16 = r[16], r17 = r[17];
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i x0, x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8;
  __m512i x9, x10, x11, x12, x13, x14, x15, x16, x17;
  const __m512i xmask = VSUB(VZERO, VAND(b, VSET1(1)));

  x0  = VXOR(r0,  a0);  x1  = VXOR(r1,  a1);  x2  = VXOR(r2,  a2);
  x3  = VXOR(r3,  a3);  x4  = VXOR(r4,  a4);  x5  = VXOR(r5,  a5);
  x6  = VXOR(r6,  a6);  x7  = VXOR(r7,  a7);  x8  = VXOR(r8,  a8);
  x9  = VXOR(r9,  a9);  x10 = VXOR(r10, a10); x11 = VXOR(r11, a11);
  x12 = VXOR(r12, a12); x13 = VXOR(r13, a13); x14 = VXOR(r14, a14);
  x15 = VXOR(r15, a15); x16 = VXOR(r16, a16); x17 = VXOR(r17, a17);

  x0  = VAND(x0,  xmask); x1  = VAND(x1,  xmask); x2  = VAND(x2,  xmask);
  x3  = VAND(x3,  xmask); x4  = VAND(x4,  xmask); x5  = VAND(x5,  xmask);
  x6  = VAND(x6,  xmask); x7  = VAND(x7,  xmask); x8  = VAND(x8,  xmask);
  x9  = VAND(x9,  xmask); x10 = VAND(x10, xmask); x11 = VAND(x11, xmask);
  x12 = VAND(x12, xmask); x13 = VAND(x13, xmask); x14 = VAND(x14, xmask);
  x15 = VAND(x15, xmask); x16 = VAND(x16, xmask); x17 = VAND(x17, xmask);

  r0  = VXOR(r0,  x0);   r1  = VXOR(r1,  x1);   r2  = VXOR(r2,  x2);
  r3  = VXOR(r3,  x3);   r4  = VXOR(r4,  x4);   r5  = VXOR(r5,  x5);
  r6  = VXOR(r6,  x6);   r7  = VXOR(r7,  x7);   r8  = VXOR(r8,  x8);
  r9  = VXOR(r9,  x9);   r10 = VXOR(r10, x10);  r11 = VXOR(r11, x11);
  r12 = VXOR(r12, x12);  r13 = VXOR(r13, x13);  r14 = VXOR(r14, x14);
  r15 = VXOR(r15, x15);  r16 = VXOR(r16, x16);  r17 = VXOR(r17, x17);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// conditional swap
// swap field elements a and r if b == 1; do not swap if b == 0
void gfp_cswap_8x1w(htfe_t r, htfe_t a, const __m512i b)
{
  __m512i r0  = r[0],  r1  = r[1],  r2  = r[2],  r3  = r[3],  r4  = r[4],  r5  = r[5];
  __m512i r6  = r[6],  r7  = r[7],  r8  = r[8],  r9  = r[9],  r10 = r[10], r11 = r[11];
  __m512i r12 = r[12], r13 = r[13], r14 = r[14], r15 = r[15], r16 = r[16], r17 = r[17];
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i x0, x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8;
  __m512i x9, x10, x11, x12, x13, x14, x15, x16, x17;
  const __m512i xmask = VSUB(VZERO, b);

  x0  = VXOR(r0,  a0);  x1  = VXOR(r1,  a1);  x2  = VXOR(r2,  a2);
  x3  = VXOR(r3,  a3);  x4  = VXOR(r4,  a4);  x5  = VXOR(r5,  a5);
  x6  = VXOR(r6,  a6);  x7  = VXOR(r7,  a7);  x8  = VXOR(r8,  a8);
  x9  = VXOR(r9,  a9);  x10 = VXOR(r10, a10); x11 = VXOR(r11, a11);
  x12 = VXOR(r12, a12); x13 = VXOR(r13, a13); x14 = VXOR(r14, a14);
  x15 = VXOR(r15, a15); x16 = VXOR(r16, a16); x17 = VXOR(r17, a17);

  x0  = VAND(x0,  xmask); x1  = VAND(x1,  xmask); x2  = VAND(x2,  xmask);
  x3  = VAND(x3,  xmask); x4  = VAND(x4,  xmask); x5  = VAND(x5,  xmask);
  x6  = VAND(x6,  xmask); x7  = VAND(x7,  xmask); x8  = VAND(x8,  xmask);
  x9  = VAND(x9,  xmask); x10 = VAND(x10, xmask); x11 = VAND(x11, xmask);
  x12 = VAND(x12, xmask); x13 = VAND(x13, xmask); x14 = VAND(x14, xmask);
  x15 = VAND(x15, xmask); x16 = VAND(x16, xmask); x17 = VAND(x17, xmask);

  r0  = VXOR(r0,  x0);   r1  = VXOR(r1,  x1);   r2  = VXOR(r2,  x2);
  r3  = VXOR(r3,  x3);   r4  = VXOR(r4,  x4);   r5  = VXOR(r5,  x5);
  r6  = VXOR(r6,  x6);   r7  = VXOR(r7,  x7);   r8  = VXOR(r8,  x8);
  r9  = VXOR(r9,  x9);   r10 = VXOR(r10, x10);  r11 = VXOR(r11, x11);
  r12 = VXOR(r12, x12);  r13 = VXOR(r13, x13);  r14 = VXOR(r14, x14);
  r15 = VXOR(r15, x15);  r16 = VXOR(r16, x16);  r17 = VXOR(r17, x17);

  a0  = VXOR(a0,  x0);   a1  = VXOR(a1,  x1);   a2  = VXOR(a2,  x2);
  a3  = VXOR(a3,  x3);   a4  = VXOR(a4,  x4);   a5  = VXOR(a5,  x5);
  a6  = VXOR(a6,  x6);   a7  = VXOR(a7,  x7);   a8  = VXOR(a8,  x8);
  a9  = VXOR(a9,  x9);   a10 = VXOR(a10, x10);  a11 = VXOR(a11, x11);
  a12 = VXOR(a12, x12);  a13 = VXOR(a13, x13);  a14 = VXOR(a14, x14);
  a15 = VXOR(a15, x15);  a16 = VXOR(a16, x16);  a17 = VXOR(a17, x17);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;

  a[0]  = a0;  a[1]  = a1;  a[2]  = a2;  a[3]  = a3;  a[4]  = a4;  a[5]  = a5;
  a[6]  = a6;  a[7]  = a7;  a[8]  = a8;  a[9]  = a9;  a[10] = a10; a[11] = a11;
  a[12] = a12; a[13] = a13; a[14] = a14; a[15] = a15; a[16] = a16; a[17] = a17;
}

// check whether the field element is a square in the field  
// if a is     a square in the field, return 1; 
// if a is not a square in the field, return 0.
__m512i gfp_issqr_8x1w(const htfe_t a)
{  
  htfe_t t0, t1;
  __m512i r = VZERO;
  int i;

  // compute r = a^((p-1)/2) - 1 
  // r is 0 if a is a square; r is non-0 if a is not a square
  gfp_copy_8x1w(t0, a);                 // t0 = a
  gfp_pow_8x1w(t1, t0, u64_pdiv2);      // t1 = a^((p-1)/2)
  gfp_rdcp_8x1w(t1, t1);                // t1 in [0, p) and strictly radix-29 now 
  for(i = 0; i < HT_NWORDS; i++)        // r = t1 - 1
    r = VOR(r, VSUB(t1[i], VSET1(ht_montR[i])));

  r = VAND(r, VSET1(HT_BMASK));         // r is 29-bit now 
  r = VSUB(VZERO, r);                   // r = 0 - r
  r = VSHR(r, 63);                      // r = r >> 63

  // now r == 0 means "a" is a square; r == 1 means "a" is not a square 
  // so we need to XOR r
  r = VXOR(r, VSET1(1));

  return r;
}

// check whether the field element is zero 
// if a is     zero, return 1; 
// if a is not zero, return 0.
__m512i gfp_iszero_8x1w(const htfe_t a)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i t = VZERO;

  // t = OR all the limbs of a
  t = VOR(t, a0);   t = VOR(t, a1);   t = VOR(t, a2);  
  t = VOR(t, a3);   t = VOR(t, a4);   t = VOR(t, a5);  
  t = VOR(t, a6);   t = VOR(t, a7);   t = VOR(t, a8);  
  t = VOR(t, a9);   t = VOR(t, a10);  t = VOR(t, a11);  
  t = VOR(t, a12);  t = VOR(t, a13);  t = VOR(t, a14);  
  t = VOR(t, a15);  t = VOR(t, a16);  t = VOR(t, a17);

  // if a == 0 then t is all-0; otherwise t is non-0
  t = VAND(t, VSET1(HT_BMASK));
  t = VSUB(VZERO, t); 
  t = VSHR(t, 63);
  t = VXOR(t, VSET1(1));

  return t;
}

// -----------------------------------------------------------------------------
// other prime-field operations 

// check whether two 32-bit unsinged integers are equal
// if two integers are     equal, return 1;
// if two integers are not equal, return 0.
uint8_t u32_iseql(const uint32_t a, const uint32_t b)
{
  uint32_t r = 0;
  unsigned char *ta = (unsigned char *)&a;
  unsigned char *tb = (unsigned char *)&b;
  r = (ta[0] ^ tb[0]) | (ta[1] ^ tb[1]) | (ta[2] ^ tb[2]) |  (ta[3] ^ tb[3]);
  r = (-r);
  r = r >> 31;
  return (uint8_t)(1-r);

}

//------------------------------------------------------------------------------
// some functions used in [CCC+19] implementation

int compare(uint64_t *x, uint64_t *y, int NUM)
{
	int i;

	for (i = NUM-1; i >= 0; i--) if (x[i] != y[i]) return x[i] > y[i] ? 1 : -1; 
	return 0;
}

int iszero(uint64_t *x, int NUM)
{
  int i;
  uint32_t r = 0;

  for (i = 0; i < NUM; i++) r |= x[i];
  return (!r);
}

uint32_t isequal(uint32_t a, uint32_t b)
{
  uint32_t r = 0;
  unsigned char *ta = (unsigned char *)&a;
  unsigned char *tb = (unsigned char *)&b;
  r = (ta[0] ^ tb[0]) | (ta[1] ^ tb[1]) | (ta[2] ^ tb[2]) |  (ta[3] ^ tb[3]);
  r = (-r);
  r = r >> 31;
  return (uint32_t)(1-r);
}

int32_t issmaller(int32_t x, int32_t y)
{
	int32_t xy = x ^ y;
	int32_t c = x - y;
	c ^= xy & (c ^ x);
	return (c >> 31);
}

