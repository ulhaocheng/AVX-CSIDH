/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "tedcurve.h"
#include <string.h>

// In this file, we use capital "X" and "Z" to denote the projective X and Z 
// coordinates on Montgomery curve while use lowercase "y" and "z" to denote 
// the projective y and z coordinates on twisted Edwards curve.
// We use capital "A24plus", "A24minus", "C24", "A", "C" to denote the 
// coefficients of Montgomery curve while use lowercase "a", "d", "a-d" to 
// denote the coeffcients of twisted Edwards curve. 

// As explained in [CCC+19, Sect 4.1.2],
// A24plus  = A+2C = a 
// A24minus = A-2C = d 
// C24      = 4C   = a-d  
// therefore, 4A = 2*(A24plus+A24minus) = 2*(a+d) and A' = A/C = 4A/4C = 2(a+d)/(a-d).

// -----------------------------------------------------------------------------

// (8x1)-way curve and isogeny operations 

// copy point P to R
void point_copy_8x1w(htpoint_t R, const htpoint_t P)
{
  gfp_copy_8x1w(R->y, P->y);
  gfp_copy_8x1w(R->z, P->z);
}

// conditional moving points 
// move points P to R if b == 1; do not swap if b == 0
void point_cmove_8x1w(htpoint_t R, htpoint_t P, const __m512i b)
{
  gfp_cmove_8x1w(R->y, P->y, b);
  gfp_cmove_8x1w(R->z, P->z, b);
}

// conditional swapping points 
// swap points R and P if b == 1; do not swap if b == 0
void point_cswap_8x1w(htpoint_t R, htpoint_t P, const __m512i b)
{
  gfp_cswap_8x1w(R->y, P->y, b);
  gfp_cswap_8x1w(R->z, P->z, b);
}

// check if a is exactly equal to p (a must in strictly radix-29)
// only if a == p, r will be 0; otherwise r will be non-0
static void gfp_checkp_8x1w(htfe_t r, const htfe_t a)
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

  // r = a - p
  r0  = VSUB(a0,  vp0);  r1  = VSUB(a1,  vp1);  r2  = VSUB(a2,  vp2);
  r3  = VSUB(a3,  vp3);  r4  = VSUB(a4,  vp4);  r5  = VSUB(a5,  vp5);
  r6  = VSUB(a6,  vp6);  r7  = VSUB(a7,  vp7);  r8  = VSUB(a8,  vp8);
  r9  = VSUB(a9,  vp9);  r10 = VSUB(a10, vp10); r11 = VSUB(a11, vp11);
  r12 = VSUB(a12, vp12); r13 = VSUB(a13, vp13); r14 = VSUB(a14, vp14);
  r15 = VSUB(a15, vp15); r16 = VSUB(a16, vp16); r17 = VSUB(a17, vp17);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; r[16] = r16; r[17] = r17;
}

// check whether the point is a point at infinity 
// if P is     a point at infinity, return 1; 
// if P is not a point at infinity, return 0.
// P's coordinates in [0, 2p) ->
__m512i point_isinf_8x1w(const htpoint_t P)
{
  htfe_t t0, t1;
  __m512i r;

  // Checking a point whether is the infinity on Montgomery curve is to check 
  // whether the projective Z-coordinate is 0. Therefore we here just need to 
  // check if y-z == 0. 
  gfp_sub_8x1w(t0, P->y, P->z);    // t0 = Z = y-z in [0, 2p)
  // Now there is a special case that t0 might be exactly equal to p.

  // So we need to check if t0 is exactly equal to p. 
  // Here we do not use "rdcp" function, but develop a specialized "checkp" 
  // function to save the cost of carry propagation.
  gfp_checkp_8x1w(t1, t0);         // if t0 == p then t1 = 0; otherwise, t1 = non-0

  r = VOR(gfp_iszero_8x1w(t0), gfp_iszero_8x1w(t1));

  return r;      
}

// (8x1)-way y-coordinate doubling R = [2]P on twisted Edwards curve,
// which is very similar to x-coordinate doubling on Montgomery curve. 
// NOTE: A->y = A24plus = a, A->z = C24 = a-d
void yDBL_8x1w(htpoint_t R, const htpoint_t P, const htpoint_t A)
{
  htfe_t t0, t1;

  // Here we use Montgomery curve projective X and Z coordinates in the comments 
  // for easy understanding.

  gfp_sqr_8x1w(t0, P->y);         // t0 = (X-Z)^2
  gfp_sqr_8x1w(t1, P->z);         // t1 = (X+Z)^2
  gfp_mul_8x1w(R->z, A->z, t0);   // zR = C24*(X-Z)^2
  gfp_mul_8x1w(R->y, R->z, t1);   // yR = C24*(X-Z)^2*(X+Z)^2
  gfp_sub_8x1w(t1, t1, t0);       // t1 = (X+Z)^2-(X-Z)^2 = 4XZ
  gfp_mul_8x1w(t0, A->y, t1);     // t0 = A24plus*4XZ
  gfp_add_8x1w(R->z, R->z, t0);   // zR = C24*(X-Z)^2+A24plus*4XZ
  gfp_mul_8x1w(t0, R->z, t1);     // t0 = [C24*(X-Z)^2+A24plus*4XZ]*4XZ 

  // Now yR is actually X2 of xDBL; t0 is actually Z2 of xDBL
  // We convert them to twsited Edwards curve projective y and z coordinates. 
  gfp_add_8x1w(R->z, R->y, t0);   // zR = X2 + Z2
  gfp_sub_8x1w(R->y, R->y, t0);   // yR = X2 - Z2
}

// (8x1)-way y-coordinate addition R = P + Q on twisted Edwards curve,
// which is very similar to x-coordinate addition on Montgomery curve. 
// NOTE: A->y = A24plus = a, A->z = C24 = a-d.
void yADD_8x1w(htpoint_t R, const htpoint_t P, const htpoint_t Q, const htpoint_t PQ)
{
  htfe_t t0, t1, t2, t3;

  // Here we use Montgomery curve projective X and Z coordinates in the comments 
  // for easy understanding.

  gfp_add_8x1w(t2, PQ->z, PQ->y); // t2 = 2X1 
  gfp_sub_8x1w(t3, PQ->z, PQ->y); // t3 = 2Z1
  gfp_mul_8x1w(t0, P->z, Q->y);   // t0 = (X2+Z2)*(X3-Z3)
  gfp_mul_8x1w(t1, P->y, Q->z);   // t1 = (X2-Z2)*(X3+Z3)
  gfp_sub_8x1w(R->z, t0, t1);     // zR = (X2+Z2)*(X3-Z3)-(X2-Z2)*(X3+Z3)
  gfp_add_8x1w(R->y, t0, t1);     // yR = (X2+Z2)*(X3-Z3)+(X2-Z2)*(X3+Z3)
  gfp_sqr_8x1w(R->z, R->z);       // zR = [(X2+Z2)*(X3-Z3)-(X2-Z2)*(X3+Z3)]^2
  gfp_sqr_8x1w(R->y, R->y);       // yR = [(X2+Z2)*(X3-Z3)+(X2-Z2)*(X3+Z3)]^2
  gfp_mul_8x1w(t0, R->y, t3);     // t0 = 2Z1*[(X2+Z2)*(X3-Z3)+(X2-Z2)*(X3+Z3)]^2
  gfp_mul_8x1w(t1, R->z, t2);     // t1 = 2X1*[(X2+Z2)*(X3-Z3)-(X2-Z2)*(X3+Z3)]^2

  // Now t0 is actually X5 of xADD; t1 is actually Z5 of xADD
  // We convert them to twsited Edwards curve projective y and z coordinates. 
  gfp_sub_8x1w(R->y, t0, t1);     // yR = X5 - Z5
  gfp_add_8x1w(R->z, t0, t1);     // zR = X5 + Z5
}

// (8x1)-way y-coordinate scalar multiplication R = [k]P on twsited Edwards curve.
// The scalar k is a *public* parameter and the same for all 8 instances.
// NOTE: A->y = A24plus = a, A->z = C24 = a-d.
// Please see details about yMUL in [CCC+19, Sect 4.2].
void yMUL_8x1w(htpoint_t R, const htpoint_t P, const htpoint_t A, const uint8_t k)
{
  htpoint T[3], Q, TQ;
  uint32_t m = addc[k];
  __m512i inf;
  int i, f_inf;

  point_copy_8x1w(&T[0], P);            // T0 = P   
  yDBL_8x1w(&T[1], P, A);               // T1 = [2]P
  yADD_8x1w(&T[2], &T[1], &T[0], P);    // T2 = [3]P

  for (i = 0; i < addc_len[k]; i++) {

    inf = point_isinf_8x1w(&T[m&1]);    // constant-time check if T_{m&1} is the point at infinity 
    f_inf = VORRDC(inf);                // f_inf == 1 if any one of 8 instances has the point at infinity
    
    // The following if-else branch only depends on the randomness. 
    // (this if-else branch is also kept in the original [CCC+19] x64 implementation)
    // A fully constant-time scalar multiplication can be implemented with a classic constant-time Montgomery ladder.

    // f_inf == 0 means there is no infinity point in 8 instances,  
    // then we just perform yADD for all instances.
    if (!f_inf) yADD_8x1w(&Q, &T[2], &T[(m&1)^1], &T[m&1]);
    // f_inf != 0 means at least one of eight instances has the point at infinity,
    // then we perform both yDBL and yADD, and later use CMOVE to get the correct values.
    else {
      yDBL_8x1w(&TQ, &T[2], A);         
      yADD_8x1w(&Q, &T[2], &T[(m&1)^1], &T[m&1]); 
      // for the instances having the point at infinity, we move the TQ to Q
      point_cmove_8x1w(&Q, &TQ, inf);  
    }
    
    point_copy_8x1w(&T[0], &T[(m&1)^1]);
    point_copy_8x1w(&T[1], &T[2]);
    point_copy_8x1w(&T[2], &Q);
    m >>= 1;
  } 
  point_copy_8x1w(R, &T[2]);
} 

// compare two mpi64 integers 
static int mpi64_compare(uint64_t *x, const uint64_t *y)
{
  int i;

  for (i = 7; i >= 0; i--) 
    if (x[i] != y[i]) return x[i] > y[i] ? 1 : -1; 
  
  return 0;
}

// (8x1)-way elligator.
// NOTE: A->y = A24plus = a, A->z = C24 = a-d.
// Please see details about this elligator in [CCC+19, Sect 3].
void elligator_8x1w(htpoint_t Tplus, htpoint_t Tminus, const htpoint_t A)
{
  htfe_t tc, t0, t1, c0, c1, alpha, vu, vR;;
  uint64_t u64[8][8];
  uint32_t u32[8][16], u[8][HT_NWORDS];
  __m512i m;
  int i;
 
  gfp_zero_8x1w(Tplus->y);                        // yT+ = 0
  gfp_zero_8x1w(Tminus->y);                       // yT- = 0
  gfp_zero_8x1w(alpha);                           // alpha = 0

  for (i = 0; i < 8; i++) {
    // generate eight random values u0...u7 
    do {
      fp_random(u64[i]);
    } while (mpi64_compare(u64[i], u64_pdiv2) > 0);// repeat if u > (p-1)/2 
    mpi_conv_64to32(u32[i], u64[i]);               // convert to radix-32
    mpi_conv_32to29(u[i], u32[i], HT_NWORDS, 16);  // convert to radix-29      
  }

  for (i = 0; i < HT_NWORDS; i++) { 
    vu[i] = set_vector(u[7][i], u[6][i], u[5][i], u[4][i], u[3][i], u[2][i], u[1][i], u[0][i]);
    vR[i] = VSET1(ht_montR[i]);
  }

  // Here we use Montgomery curve ceoffcients A and C in the comments 
  // for easy understanding.

  gfp_num2mont_8x1w(vu, vu);                      // convert u to Montgomery domain 
  gfp_sqr_8x1w(Tplus->z, vu);                     // zT+ = u^2
  gfp_add_8x1w(c0, Tplus->z, vR);                 // c0 = u^2+1
  gfp_sub_8x1w(tc, Tplus->z, vR);                 // tc = u^2-1
  gfp_mul_8x1w(c1, A->z, tc);                     //!c1 = 4*C*(u^2-1) = C*(u^2-1)
  gfp_sub_8x1w(Tminus->z, A->y, A->z);            // zT- = d = A-2C 
  gfp_add_8x1w(Tminus->z, Tminus->z, A->y);       // zT- = a+d = 2A
  gfp_add_8x1w(Tminus->z, Tminus->z, Tminus->z);  //!zT- = 2*(a+d) = 4A = A
  gfp_mul_8x1w(t0, Tminus->z, c1);                // t0 = A*C*(u^2-1)
  gfp_sqr_8x1w(t1, Tminus->z);                    // t1 = A^2
  gfp_mul_8x1w(t1, t1, Tplus->z);                 // t1 = (Au)^2
  gfp_sqr_8x1w(tc, c1);                           // tc = (C*(u^2-1))^2
  gfp_add_8x1w(t1, t1, tc);                       // t1 = (Au)^2+(C*(u^2-1))^2
  gfp_mul_8x1w(tc, t0, t1);                       // tc = A*C*(u^2-1)*[(Au)^2+(C*(u^2-1))^2]  
  gfp_rdcp_8x1w(tc, tc);                          // reduce tc to [0, p)
  gfp_cmove_8x1w(alpha, vu, gfp_iszero_8x1w(tc)); // alpha = u if A == 0; alpha = 0 otherwise
  gfp_mul_8x1w(c0, alpha, c0);                    // c0 = alpha*(u^2+1) 
  gfp_mul_8x1w(alpha, alpha, c1);                 // alpha = alpha*C*(u^2-1)
  gfp_copy_8x1w(Tplus->y, Tminus->z);             // yT+ = A
  gfp_sub_8x1w(Tminus->y, Tminus->y, Tminus->z);  // yT- = -A
  gfp_mul_8x1w(Tminus->y, Tminus->y, Tplus->z);   // yT- = -A*u^2
  gfp_add_8x1w(Tplus->y, Tplus->y, alpha);        // yT+ = A+alpha*C*(u^2-1)
  gfp_sub_8x1w(Tminus->y, Tminus->y, alpha);      // yT- = -A*u^2-alpha*C*(u^2-1)
  gfp_add_8x1w(tc, tc, c0);                       // tc = A*C*(u^2-1)*[(A*u)^2+(C*(u^2-1))^2] + alpha*(u^2+1)
  
  m = gfp_issqr_8x1w(tc);                         // legendre_symbol
  gfp_cswap_8x1w(Tplus->y, Tminus->y, VXOR(m, VSET1(1)));

  // We convert them to twsited Edwards curve projective y and z coordinates.
  gfp_add_8x1w(Tplus->z, Tplus->y, c1);           // zT+ = A+alpha*C*(u^2-1) + C*(u^2-1)
  gfp_sub_8x1w(Tplus->y, Tplus->y, c1);           // yT+ = A+alpha*C*(u^2-1) - C*(u^2-1) 
  gfp_add_8x1w(Tminus->z, Tminus->y, c1);         // zT- = -A*u^2-alpha*C*(u^2-1) + C*(u^2-1)
  gfp_sub_8x1w(Tminus->y, Tminus->y, c1);         // yT- = -A*u^2-alpha*C*(u^2-1) - C*(u^2-1) 
}

// (8x1)-way y-coordinate isogeny computation on twisted Edwards curve.
// NOTE: A->y = A24plus = a, A->z = C24 = a-d.
void yISOG_8x1w(htpoint R[], htpoint_t C, const htpoint_t P, const htpoint_t A, const uint8_t k)
{
  uint8_t mask;
  int lbits = bits_li[k], i, l = primeli[k], s = l>>1; 
  htfe_t t0, t1, td, By0, By1, Bz0, Bz1;

  gfp_copy_8x1w(t0, A->y);              // t0 = a
  gfp_sub_8x1w(td, A->y, A->z);         // td = a - (a-d) = d
  gfp_copy_8x1w(t1, td);                // t1 = d

  gfp_copy_8x1w(By0, P->y);             // By0 = yP
  gfp_copy_8x1w(Bz0, P->z);             // Bz0 = zP

  point_copy_8x1w(&R[0], P);            // P
  yDBL_8x1w(&R[1], P, A);               // [2]P

  for (i = 2; i < s; i++) {
    gfp_mul_8x1w(By0, By0, R[i-1].y);   // By0 = By0 * yR_{i-1}
    gfp_mul_8x1w(Bz0, Bz0, R[i-1].z);   // Bz0 = Bz0 * zR_{i-1}
    yADD_8x1w(&R[i], &R[i-1], P, &R[i-2]);
  }

  gfp_mul_8x1w(By1, By0, R[s-1].y);     // By1 = By0 * yR_{s-1}
  gfp_mul_8x1w(Bz1, Bz0, R[s-1].z);     // Bz1 = Bz0 * zR_{s-1}
  mask = u32_iseql(l, 3)^1;
  gfp_cswap_8x1w(By0, By1, VSET1(mask));           
  gfp_cswap_8x1w(Bz0, Bz1, VSET1(mask));

  // left-to-right computing a^l and d^l
  lbits -= 1;
  for (i = 1; i <= lbits; i++) {
    gfp_sqr_8x1w(t0, t0);                   
    gfp_sqr_8x1w(t1, t1);
    if ((l>>(lbits-i)) & 1) {
      gfp_mul_8x1w(t0, t0, A->y);
      gfp_mul_8x1w(t1, t1, td);
    }
  }

  gfp_sqr_8x1w(By0, By0); gfp_sqr_8x1w(By0, By0); gfp_sqr_8x1w(By0, By0);
  gfp_sqr_8x1w(Bz0, Bz0); gfp_sqr_8x1w(Bz0, Bz0); gfp_sqr_8x1w(Bz0, Bz0);
  
  gfp_mul_8x1w(C->y, t0, Bz0);
  gfp_mul_8x1w(C->z, t1, By0);
  gfp_sub_8x1w(C->z, C->y, C->z);       // z coordinate stores a-d
}

// (8x1)-way y-coordinate isogeny evaluation on twisted Edwards curve.
void yEVAL_8x1w(htpoint_t R, const htpoint_t Q, const htpoint P[], const uint8_t k)
{
  htfe_t t0, t1, s0, s1;
  htpoint T;
  int i, s;

  point_copy_8x1w(&T, Q);      

  gfp_mul_8x1w(s0, T.y, P[0].z);
  gfp_mul_8x1w(s1, T.z, P[0].y);
  gfp_add_8x1w(R->y, s0, s1);
  gfp_sub_8x1w(R->z, s0, s1);

  s = primeli[k]>>1;
  for (i = 1; i < s; i++) {
    gfp_mul_8x1w(s0, T.y, P[i].z);
    gfp_mul_8x1w(s1, T.z, P[i].y);
    gfp_add_8x1w(t0, s0, s1);
    gfp_sub_8x1w(t1, s0, s1);
    gfp_mul_8x1w(R->y, R->y, t0);
    gfp_mul_8x1w(R->z, R->z, t1);
  }

  gfp_sqr_8x1w(R->y, R->y);
  gfp_sqr_8x1w(R->z, R->z);
  gfp_add_8x1w(t0, T.z, T.y);
  gfp_sub_8x1w(t1, T.z, T.y);
  gfp_mul_8x1w(t0, R->y, t0);
  gfp_mul_8x1w(t1, R->z, t1);
  gfp_sub_8x1w(R->y, t0 ,t1);
  gfp_add_8x1w(R->z, t0, t1);
}

// -----------------------------------------------------------------------------
// 1-way curve and isogeny operations from [CCC+19] code with no modification

uint8_t point_isinf(const proj P)
{
  fp tmp;

  fp_sub(tmp, P[0], P[1]);				
  return iszero(tmp, NWORDS);	
};

void point_copy(proj Q, const proj P)
{
  copy(Q[0], P[0], NWORDS);
  copy(Q[1], P[1], NWORDS);
};

void point_cswap(proj R, proj P, const uint8_t b)
{
  fp_cswap(R[0], P[0], b&1);
  fp_cswap(R[1], P[1], b&1);
}

void yDBL(proj Q, const proj P, const proj A)
{
  fp tmp0, tmp1;

  // yDBL is performed on the isomorphic Montgomery curve
  fp_sqr(tmp0, P[0]);
  fp_sqr(tmp1, P[1]);

  fp_mul(Q[1], A[1], tmp0);
  fp_mul(Q[0], Q[1], tmp1);
  fp_sub(tmp1, tmp1, tmp0);
  fp_mul(tmp0, A[0], tmp1);
  fp_add(Q[1], Q[1], tmp0);
  fp_mul(tmp0, Q[1], tmp1);

  // the result should be mapped into the twisted Edwards curve
  fp_add(Q[1], Q[0], tmp0);
  fp_sub(Q[0], Q[0], tmp0);
};

void yADD(proj R, const proj P, const proj Q, const proj PQ)
{
  fp tmp0, tmp1, xD, zD;

  // the difference should be into the isomorphic Montgomery curve
  fp_add(xD, PQ[1], PQ[0]);
  fp_sub(zD, PQ[1], PQ[0]);

  fp_mul(tmp0, P[1], Q[0]);
  fp_mul(tmp1, P[0], Q[1]);

  fp_sub(R[1], tmp0, tmp1);
  fp_add(R[0], tmp0, tmp1);

  fp_sqr(R[1], R[1]);
  fp_sqr(R[0], R[0]);

  fp_mul(tmp0, R[0], zD);
  fp_mul(tmp1, R[1], xD);

  // the result should be mapped into the twisted Edwards curve
  fp_sub(R[0], tmp0, tmp1);
  fp_add(R[1], tmp0, tmp1);
};

void yMUL(proj Q, const proj P, const proj A, const uint8_t i)
{
	proj R[3], T;
  uint32_t tmp = addc[i];
  int j;

  // R0 = P; R1 = 2P; R2 = 3P
	point_copy(R[0], P);
	yDBL(R[1], P, A);
	yADD(R[2], R[1], R[0], P);

  // main loop
	for (j = 0; j < addc_len[i]; j++)
	{
    if (point_isinf(R[tmp&1]))
      yDBL(T, R[2], A);
    else
      yADD(T, R[2], R[(tmp & 0x1) ^ 0x1], R[tmp & 0x1]);

		point_copy(R[0], R[(tmp&1)^1]);
		point_copy(R[1], R[2]);
		point_copy(R[2], T);
    tmp >>= 1;
	};

	point_copy(Q, R[2]);	
};

void elligator(proj T_plus, proj T_minus, const proj A)
{
  fp u, tmp, u2_plus_1, Cu2_minus_1, tmp_0, tmp_1, alpha, beta;

  // T+ <- 0 and T- <- 0
  set_zero(T_plus[0], NWORDS);			
  set_zero(T_minus[0], NWORDS);			

	fp_random(u);
	while ( compare(u, (uint64_t *)p_minus_1_halves, NWORDS) > 0) fp_random(u);

	fp_mul(u, u, R_squared_mod_p);	// mapping u into the Montgomery domain

	set_zero(alpha, NWORDS);			// 0
	fp_add(beta, alpha, u);						// u
	
	fp_sqr(T_plus[1], u);						// u^2
	fp_add(u2_plus_1, T_plus[1], R_mod_p);		// u^2 + 1
	fp_sub(tmp, T_plus[1], R_mod_p);			// u^2 - 1
	fp_mul(Cu2_minus_1, A[1], tmp);					// C' * (u^2 - 1)

	// The goal is to evaluate in the projective Montgomery curve isomorphic to A
	fp_sub(T_minus[1], A[0], A[1]);					// A' := 2 * (a + d) and C' := (a - d) are
	fp_add(T_minus[1], T_minus[1], A[0]);			// projective constants of the isomorphic
	fp_add(T_minus[1], T_minus[1], T_minus[1]);		// Montgomery curve

	fp_mul(tmp_0, T_minus[1], Cu2_minus_1);			// A' * C' * (u^2 - 1)

	fp_sqr(tmp_1, T_minus[1]);				// (A')^2
	fp_mul(tmp_1, tmp_1, T_plus[1]);		// (A' * u)^2
	fp_sqr(tmp, Cu2_minus_1);				// [C' * (u^2 - 1)]^2
	fp_add(tmp_1, tmp_1, tmp);				// (A' * u)^2 + [C' * (u^2 - 1)]^2
	
	fp_mul(tmp, tmp_0, tmp_1);				// {A' * C' * (u^2 - 1)} * {(A' * u)^2 + [C' * (u^2 - 1)]^2} =? { [C' * (u^2 - 1)]^2 * w}^2
	
	//
	fp_cswap(alpha, beta ,iszero(tmp, NWORDS));	// alpha = 0 if A' = 0; alpha = u otherwise
  // fp_cswap(alpha, beta, fp_iszero(tmp));
	fp_mul(u2_plus_1, alpha, u2_plus_1);					// u2_plus_1 = 0 if A' != 0; u2_plus_1 = u^3 + u
	fp_mul(alpha, alpha, Cu2_minus_1);						// alpha * C' * (u^2 - 1)

	// the projective y-coordinate of T_{+} or T_{-}
	fp_add(T_plus[0], T_plus[0], T_minus[1]);			// [C' * (u^2 - 1) * v] = A'
	// the projective y-coordinate of T_{-} or T_{+}
	fp_sub(T_minus[0], T_minus[0], T_minus[1]);			// -[C' * (u^2 - 1) * v] = -A'
	fp_mul(T_minus[0], T_minus[0], T_plus[1]);			// -[ (C' * v) + A'] * (u^2) = -A' * (u^2)

	fp_add(T_plus[0], T_plus[0], alpha);				//  A' + alpha*C'*(u^2-1)
	fp_sub(T_minus[0], T_minus[0], alpha);				// -A' * (u^2) - alpha*C'*(u^2-1)

	fp_add(tmp, tmp, u2_plus_1);						// Now, if A'=0 then tmp = u^3 + u
	// Only one legendre symbol computation is required
	uint8_t legendre_symbol = fp_issquare(tmp) &  0x1;
	fp_cswap(T_plus[0], T_minus[0], legendre_symbol ^ 1);		// constant-time swap for determining T_{+} or T_{-}.
	
	// Finally, we mapping the points into the Edward's curves
	// T_{+}
	fp_add(T_plus[1], T_plus[0], Cu2_minus_1);
	fp_sub(T_plus[0], T_plus[0], Cu2_minus_1);
	// T_{-}
	fp_add(T_minus[1], T_minus[0], Cu2_minus_1);
	fp_sub(T_minus[0], T_minus[0], Cu2_minus_1);
};

void yISOG(proj Pk[], proj C, const proj P, const proj A, const uint8_t i)
{
  uint8_t mask;
  int64_t bits_l = bits_li[i];
  uint64_t j, l = primeli[i],	s = l>>1;

	// ---
	fp By[2], Bz[2], tmp_0, tmp_1, tmp_d;

	copy(tmp_0, A[0], NWORDS);		// a
	fp_sub(tmp_d, A[0], A[1]);			// d
	copy(tmp_1, tmp_d, NWORDS);

	copy(By[0], P[0], NWORDS); copy(By[1], P[0], NWORDS);
	copy(Bz[0], P[1], NWORDS); copy(Bz[1], P[1], NWORDS);

	point_copy(Pk[0], P);				// P
	yDBL(Pk[1], P, A);				// [2]P
		
	for(j = 2; j < s; j++)
	{
		fp_mul(By[0], By[0], Pk[j - 1][0]);
		fp_mul(Bz[0], Bz[0], Pk[j - 1][1]);
		yADD(Pk[j], Pk[j - 1], P, Pk[j - 2]);	// [j + 1]P
	};

	mask = isequal(l, 3) ^ 1;		// If l = 3 then we keep with the current values of By[0] and Bz[0]. This ask is done in constant-time
	fp_mul(By[1], By[0], Pk[s - 1][0]);	// This an extra cost for a degree-3 construction
	fp_mul(Bz[1], Bz[0], Pk[s - 1][1]);	// This an extra cost for a degree-3 construction
	fp_cswap(By[0], By[1], mask);		// constant-time swap: dummy or not dummy, that is the question.
	fp_cswap(Bz[0], Bz[1], mask);		// constant-time swap: dummy or not dummy, that is the question.

	// left-to-right method for computing a^l and d^l
	bits_l -= 1;
	for(j = 1; j <= bits_l; j++)
	{
		fp_sqr(tmp_0, tmp_0);
		fp_sqr(tmp_1, tmp_1);
		if( ( (l >> (bits_l - j)) & 1 ) != 0)
		{
			fp_mul(tmp_0, tmp_0, A[0]);
			fp_mul(tmp_1, tmp_1, tmp_d);
		};
	};

	for(j = 0; j < 3; j++)
	{
		fp_sqr(By[0], By[0]);
		fp_sqr(Bz[0], Bz[0]);
	};

	fp_mul(C[0], tmp_0, Bz[0]);
	fp_mul(C[1], tmp_1, By[0]);
	fp_sub(C[1], C[0], C[1]);
};

void yEVAL(proj R, const proj Q, const proj Pk[], const uint8_t i)
{
	int j;
	fp tmp_0, tmp_1, s_0, s_1;

	proj tmp_Q;
	point_copy(tmp_Q, Q);	// This is for allowing Q <- image of Q

	// Evaluating Q
	fp_mul(s_0, tmp_Q[0], Pk[0][1]);
	fp_mul(s_1, tmp_Q[1], Pk[0][0]);
	// Mapping R into the isomorphic Montgomery curve
	fp_add(R[0], s_0, s_1);
	fp_sub(R[1], s_0, s_1);

	uint64_t s = (primeli[i] >> 1);
	for(j = 1; j < s; j++)
	{
		// Evaluating Q
		fp_mul(s_0, tmp_Q[0], Pk[j][1]);
		fp_mul(s_1, tmp_Q[1], Pk[j][0]);
		fp_add(tmp_0, s_0, s_1);
		fp_sub(tmp_1, s_0, s_1);
		fp_mul(R[0], R[0], tmp_0);
		fp_mul(R[1], R[1], tmp_1);
	};

	fp_sqr(R[0], R[0]);
	fp_sqr(R[1], R[1]);
	// Mapping Q into the isomorphic Montgomery curve
	fp_add(tmp_0, tmp_Q[1], tmp_Q[0]);
	fp_sub(tmp_1, tmp_Q[1], tmp_Q[0]);
	fp_mul(tmp_0, R[0], tmp_0);
	fp_mul(tmp_1, R[1], tmp_1);
	// Mapping R into the Edwards curve
	fp_sub(R[0], tmp_0, tmp_1);
	fp_add(R[1], tmp_0, tmp_1);
};
