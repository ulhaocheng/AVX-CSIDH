/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "tedcurve.h"

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

// check if a is exactly equal to p (a must in strictly radix-52)
// only if a == p, r will be 0; otherwise r will be non-0
static void gfp_checkp_8x1w(htfe_t r, const htfe_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8], a9 = a[9]; 
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
  const __m512i vp0 = VSET1(ht_p[0]), vp1 = VSET1(ht_p[1]);
  const __m512i vp2 = VSET1(ht_p[2]), vp3 = VSET1(ht_p[3]);
  const __m512i vp4 = VSET1(ht_p[4]), vp5 = VSET1(ht_p[5]);
  const __m512i vp6 = VSET1(ht_p[6]), vp7 = VSET1(ht_p[7]);
  const __m512i vp8 = VSET1(ht_p[8]), vp9 = VSET1(ht_p[9]);

  // r = a - p
  r0 = VSUB(a0, vp0); r1 = VSUB(a1, vp1); r2 = VSUB(a2, vp2); r3 = VSUB(a3, vp3);
  r4 = VSUB(a4, vp4); r5 = VSUB(a5, vp5); r6 = VSUB(a6, vp6); r7 = VSUB(a7, vp7);
  r8 = VSUB(a8, vp8); r9 = VSUB(a9, vp9);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4;  
  r[5] = r5; r[6] = r6; r[7] = r7; r[8] = r8; r[9] = r9; 
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
  uint64_t u64[8][8], u[8][HT_NWORDS];
  __m512i m;
  int i;
 
  gfp_zero_8x1w(Tplus->y);                        // yT+ = 0
  gfp_zero_8x1w(Tminus->y);                       // yT- = 0
  gfp_zero_8x1w(alpha);                           // alpha = 0

  for (i = 0; i < 8; i++) {
    // generate eight random values u0...u7 
    do {
      mpi64_random(u64[i]);
    } while (mpi64_compare(u64[i], u64_pdiv2) > 0);// repeat if u > (p-1)/2 
    mpi_conv_64to52(u[i], u64[i], HT_NWORDS, 8);  // convert to radix-52      
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

// (2x4)-way curve and isogeny operations 

// copy point P to R
void point_copy_2x4w(llpoint_t R, const llpoint_t P)
{
  gfp_copy_2x4w(R, P);
}

// conditional moving points 
// move points P to R if b == 1; do not swap if b == 0
void point_cmove_2x4w(llpoint_t R, llpoint_t P, const uint8_t b)
{
  gfp_cmove_2x4w(R, P, b);
}

// conditional swapping points 
// swap points R and P if b == 1; do not swap if b == 0
void point_cswap_2x4w(llpoint_t R, llpoint_t P, const uint8_t b)
{
  gfp_cswap_2x4w(R, P, b);
}

// check if a is exactly equal to p (a must in strictly radix-43)
// only if a == p, r will be 0; otherwise r will be non-0
static void gfp_checkp_2x4w(llfe_t r, const llfe_t a)
{
  const __m512i p0 = VSET(0, 0, 0, 0, ll_p[9] , ll_p[6], ll_p[3], ll_p[0]);
  const __m512i p1 = VSET(0, 0, 0, 0, ll_p[10], ll_p[7], ll_p[4], ll_p[1]);
  const __m512i p2 = VSET(0, 0, 0, 0, ll_p[11], ll_p[8], ll_p[5], ll_p[2]);
  const __m512i a0 = a[0], a1 = a[1], a2 = a[2];
  __m512i r0, r1, r2;

  // r = a - p
  r0 = VSUB(a0, p0); r1 = VSUB(a1, p1); r2 = VSUB(a2, p2);

  r[0] = r0; r[1] = r1; r[2] = r2;
}

// check whether the point is a point at infinity 
// if P is     the point at infinity, return 1; 
// if P is not the point at infinity, return 0.
// P's coordinates in [0, 2p) ->
uint8_t point_isinf_2x4w(const llpoint_t P)
{
  llfe_t y, z, t0, t1;
  uint8_t r;
  
  vec_permzh_2x4w(y, P);           // 0 | y
  vec_permzl_2x4w(z, P);           // 0 | z
  gfp_addsubc_2x4w(t0, y, z);      // t0 = 0 | Z=y-z in [0, 2p) and strictly radix-43
  // Now there is a special case that t0 might be exactly equal to p.

  // So we need to check if t0 is exactly equal to p. 
  // Here we do not use "rdcp" function, but develop a specialized "checkp" 
  // function to save the cost of carry propagation.
  gfp_checkp_2x4w(t1, t0);         // if t0 == p then t1 = 0; otherwise, t1 = non-0
  r = gfp_iszero_2x4w(t0) | gfp_iszero_2x4w(t1);

  return r;                      
}

// (2x4)-way y-coordinate doubling R = [2]P on twisted Edwards curve,
// which is very similar to x-coordinate doubling on Montgomery curve. 
// NOTE: A = <A24plus | C24> = <a | a-d>
void yDBL_2x4w(llpoint_t R, const llpoint_t P, const llpoint_t A)
{
  llfe_t t0, t1, t2, t3;

  // Here we use Montgomery curve projective X and Z coordinates in the comments 
  // for easy understanding.

  gfp_sqr_2x4w(t0, P);                  // t0 = (X-Z)^2     | (X+Z)^2
  
  vec_permzh_2x4w(t1, t0);              // t1 = 0           | (X-Z)^2
  gfp_addsubi_2x4w(t2, t0, t1);         // t2 = (X-Z)^2     | (X+Z)^2-(X-Z)^2=4XZ
  
  vec_permlh_2x4w(t2, t2);              // t2 = 4XZ         | (X-Z)^2
  gfp_mul_2x4w(t1, A, t2);              // t1 = A24plus*4XZ | C24*(X-Z)^2
  
  vec_permzh_2x4w(t3, t1);              // t3 = 0           | A24plus*4XZ
  gfp_subaddi_2x4w(t3, t1, t3);         // t3 = A24plus*4XZ | C24*(X-Z)^2+A24plus*4XZ

  vec_blend_2x4w(t2, t2, t1, 0x0F);     // t2 = 4XZ                           | C24*(X-Z)^2           
  vec_permlh_2x4w(t3, t3);              // t3 = C24*(X-Z)^2+A24plus*4XZ       | A24plus*4XZ
  vec_blend_2x4w(t0, t3, t0, 0x0F);     // t0 = C24*(X-Z)^2+A24plus*4XZ       | (X+Z)^2
  gfp_mul_2x4w(t0, t2, t0);             // t0 = [C24*(X-Z)^2+A24plus*4XZ]*4XZ | C24*(X-Z)^2*(X+Z)^2
  
  // now t0 contains Z2 | X2 of xDBL 
  vec_permlh_2x4w(t1, t0);              // t1 = X2    | Z2
  gfp_subaddc_2x4w(R, t1, t0);          // R  = X2-Z2 | X2+Z2
}

// (2x4)-way y-coordinate addition R = P + Q on twisted Edwards curve,
// which is very similar to x-coordinate addition on Montgomery curve. 
// NOTE: A = <A24plus | C24> = <a | a-d>
void yADD_2x4w(llpoint_t R, const llpoint_t P, const llpoint_t Q, const llpoint_t PQ)
{
  llfe_t t0, t1, t2;

  // Here we use Montgomery curve projective X and Z coordinates in the comments 
  // for easy understanding.

  vec_permlh_2x4w(t0, PQ);              // t0 = X1+Z1           | X1-Z1
  gfp_addsubi_2x4w(t0, PQ, t0);         // t0 = 2X1             | 2Z1

  vec_permlh_2x4w(t1, Q);               // t1 = X3+Z3           | X3-Z3
  gfp_mul_2x4w(t1, t1, P);              // t1 = (X2-Z2)*(X3+Z3) | (X2+Z2)*(X3-Z3)

  vec_permlh_2x4w(t2, t1);              // t2 = (X2+Z2)*(X3-Z3) | (X2-Z2)*(X3+Z3)
  gfp_subaddi_2x4w(t2, t1, t2);         // t2 = (X2-Z2)*(X3+Z3)-(X2+Z2)*(X3-Z3) | (X2+Z2)*(X3-Z3)+(X2-Z2)*(X3+Z3)

  gfp_sqr_2x4w(t2, t2);                 // t2 = [(X2-Z2)*(X3+Z3)-(X2+Z2)*(X3-Z3)]^2 | [(X2+Z2)*(X3-Z3)+(X2-Z2)*(X3+Z3)]^2
  gfp_mul_2x4w(t0, t0, t2);             // t0 = 2X1*[(X2+Z2)*(X3-Z3)-(X2-Z2)*(X3+Z3)]^2 | 2Z1*[(X2-Z2)*(X3+Z3)+(X2+Z2)*(X3-Z3)]^2 

  // now t0 contains Z5 | X5 of xADD
  vec_permlh_2x4w(t1, t0);              // t1 = X5    | Z5
  gfp_subaddc_2x4w(R, t1, t0);          // R  = X5-Z5 | X5+Z5
}

// (2x4)-way y-coordinate scalar multiplication R = [k]P on twsited Edwards curve.
// The scalar k is a *public* parameter.
// NOTE: A = <A24plus | C24> = <a | a-d>
// Please see details about yMUL in [CCC+19, Sect 4.2].
void yMUL_2x4w(llpoint_t R, const llpoint_t P, const llpoint_t A, const uint8_t k)
{
  llpoint_t T[3], Q;
  uint32_t m = addc[k];
  int i;

  point_copy_2x4w(T[0], P);             // T0 = P
  yDBL_2x4w(T[1], P, A);                // T1 = [2]P
  yADD_2x4w(T[2], T[1], T[0], P);       // T2 = [3]P
  
  for (i = 0; i < addc_len[k]; i++) {

    // The following if-else branch only depends on the randomness. 
    // (this if-else branch is also kept in the original [CCC+19] x64 implementation)
    // A fully constant-time scalar multiplication can be implemented with a classic constant-time Montgomery ladder.

    if (point_isinf_2x4w(T[m&1])) yDBL_2x4w(Q, T[2], A);
    else yADD_2x4w(Q, T[2], T[(m&1)^1], T[m&1]);

    point_copy_2x4w(T[0], T[(m&1)^1]);
    point_copy_2x4w(T[1], T[2]);
    point_copy_2x4w(T[2], Q);
    m >>= 1;
  }
  point_copy_2x4w(R, T[2]);
}

// (2x4)-way elligator.
// NOTE: A = <A24plus | C24> = <a | a-d>
// Please see details about this elligator in [CCC+19, Sect 3].
void elligator_2x4w(llpoint_t Tplus, llpoint_t Tminus, const llpoint_t A)
{
  llfe_t t0, t1, t2, t3, t4, t5, t6, vu, vone, vR2, alpha;
  uint64_t u64[8], u[LL_NWORDS];
  uint8_t m;
  int i;

  gfp_zero_2x4w(alpha);                 // alpha = 0 | 0

  // to get the projective coefficient A of the isomorphic Montgomery curve 
  vec_permlh_2x4w(t3, A);               //  t3 = a-d     | a
  gfp_addsubc_2x4w(t1, t3, A);          //  t1 = -d      | d
  gfp_subaddc_2x4w(t0, t1, t3);         //  t0 = a-2d    | a+d
  gfp_subaddc_2x4w(t6, t0, t0);         //  t6 = 2(a-2d) | 2(a+d)=A in Montgomery curve

  vec_permll_2x4w(t6, t6);

  do {
    mpi64_random(u64);
  } while (mpi64_compare(u64, u64_pdiv2) > 0);// repeat if u > (p-1)/2 
  mpi_conv_64to43(u, u64, LL_NWORDS, 8);  // convert to radix-43      

  // vu = u' | u'
  vu[0] = set_vector(u[9] , u[6], u[3], u[0], u[9] , u[6], u[3], u[0]);
  vu[1] = set_vector(u[10], u[7], u[4], u[1], u[10], u[7], u[4], u[1]);
  vu[2] = set_vector(u[11], u[8], u[5], u[2], u[11], u[8], u[5], u[2]); 
  // vone = 1 | 1
  vone[0] = VSET(ll_montR[9] , ll_montR[6], ll_montR[3], ll_montR[0], ll_montR[9] , ll_montR[6], ll_montR[3], ll_montR[0]);
  vone[1] = VSET(ll_montR[10], ll_montR[7], ll_montR[4], ll_montR[1], ll_montR[10], ll_montR[7], ll_montR[4], ll_montR[1]);
  vone[2] = VSET(ll_montR[11], ll_montR[8], ll_montR[5], ll_montR[2], ll_montR[11], ll_montR[8], ll_montR[5], ll_montR[2]);
  // vR2 = mont_R2 | mont_R2
  vR2[0] = VSET(ll_montR2[9] , ll_montR2[6], ll_montR2[3], ll_montR2[0], ll_montR2[9] , ll_montR2[6], ll_montR2[3], ll_montR2[0]);
  vR2[1] = VSET(ll_montR2[10], ll_montR2[7], ll_montR2[4], ll_montR2[1], ll_montR2[10], ll_montR2[7], ll_montR2[4], ll_montR2[1]);
  vR2[2] = VSET(ll_montR2[11], ll_montR2[8], ll_montR2[5], ll_montR2[2], ll_montR2[11], ll_montR2[8], ll_montR2[5], ll_montR2[2]);

  vec_blend_2x4w(t1, vu, t6, 0x0F);     //  t1 = u'  | A
  vec_blend_2x4w(t2, vR2, A, 0x0F);     //  t2 = R2  | C
  gfp_mul_2x4w(t1, t1, t2);             //  t1 = u   | AC
  vec_permhh_2x4w(vu, t1);              //  vu = u   | u

  vec_blend_2x4w(t0, vu, t6, 0x0F);     //  t0 = u   | A
  gfp_sqr_2x4w(t0, t0);                 //  t0 = u^2 | A^2

  vec_permhh_2x4w(t4, t0);              //  t4 = u^2 | u^2
  gfp_addsubc_2x4w(t2, t4, vone);       // !t2 = u^2+1 | u^2-1

  vec_blend_2x4w(t3, t3, t1, 0x0F);     //  t3 = C        | AC
  vec_permll_2x4w(t1, t2);              //  t1 = u^2-1    | u^2-1
  gfp_mul_2x4w(t1, t3, t1);             //  t1 = C(u^2-1) | AC(u^2-1)

  vec_blend_2x4w(t3, t1, t0, 0x0F);     //  t3 = C(u^2-1)     | A^2
  vec_blend_2x4w(t0, t1, t4, 0x0F);     //  t0 = C(u^2-1)     | u^2
  gfp_mul_2x4w(t0, t3, t0);             //  t0 = [C(u^2-1)]^2 | (Au)^2

  vec_blend_2x4w(t3, t0, alpha, 0x0F);  //  t3 = [C(u^2-1)]^2           | 0
  vec_permll_2x4w(t0, t0);              //  t0 = (Au)^2                 | (Au)^2
  vec_blend_2x4w(t0, t0, t6, 0x0F);     //  t0 = (Au)^2                 | A
  gfp_addsubc_2x4w(t0, t3, t0);         //  t0 = [C(u^2-1)]^2 + (Au)^2  | -A

  vec_permlh_2x4w(t5, t1);              //  t5 = AC(u^2-1)             | C(u^2-1)
  vec_blend_2x4w(t4, t5, t4, 0x0F);     // !t4 = AC(u^2-1)             | u^2
  gfp_mul_2x4w(t5, t0, t4);             //  t5 = t                     | -Au^2

  vec_permhh_2x4w(t0, t5);              //  t0 = t | t
  gfp_rdcp_2x4w(t0, t0);                // !t0 in [0, p)

  gfp_cmove_2x4w(alpha, vu, gfp_iszero_2x4w(t0)); // alpha = u | u or 0 | 0

  vec_permhh_2x4w(t1, t1);              // !t1 = C(u^2-1) | C(u^2-1)
  vec_blend_2x4w(t2, t2, t1, 0x0F);     //  t2 = u^2+1    | C(u^2-1)
  gfp_mul_2x4w(t2, t2, alpha);          // !t2 = alpha*(u^2+1) | alpha*C(u^2-1)

  vec_blend_2x4w(t5, t6, t5, 0x0F);     // !t5 = A                | -Au^2
  vec_permll_2x4w(t3, t2);              //  t3 = alpha*C(u^2-1)   | alpha*C(u^2-1)
  gfp_addsubc_2x4w(t3, t5, t3);         // !t3 = A+alpha*C(u^2-1) | -Au^2-alpha*C(u^2-1)

  gfp_addsubc_2x4w(t2, t0, t2);         //  t2 = alpha*(u^2+1) + t| 
  vec_permzh_2x4w(t2, t2);              // !t2 = 0 | alpha*(u^2+1) + t

  m = gfp_issqr_2x4w(t2);               // legendre_symbol

  vec_permhh_2x4w(Tplus, t3);           //   T+ = A+alpha*C(u^2-1)     | A+alpha*C(u^2-1)
  vec_permll_2x4w(Tminus, t3);          //   T- = -Au^2-alpha*C(u^2-1) | -Au^2-alpha*C(u^2-1)

  gfp_subaddc_2x4w(Tplus, Tplus, t1);   //   T+ = A+alpha*C(u^2-1)-C(u^2-1)     | A+alpha*C(u^2-1)+C(u^2-1) 
  gfp_subaddc_2x4w(Tminus, Tminus, t1); //   T- = -Au^2-alpha*C(u^2-1)-C(u^2-1) | -Au^2-alpha*C(u^2-1)+C(u^2-1)

  gfp_cswap_2x4w(Tplus, Tminus, m^1);
}

// (2x4)-way y-coordinate isogeny computation on twisted Edwards curve.
// NOTE: A = <A24plus | C24> = <a | a-d>
void yISOG_2x4w(llpoint_t R[], llpoint_t C, const llpoint_t P, const llpoint_t A, const uint8_t k)
{
  uint8_t mask;
  int lbits = bits_li[k];
  uint32_t i, l = primeli[k], s = l>>1; 
  llpoint_t T0, T1, AD;
  llfe_t td, t0;

  vec_permlh_2x4w(t0, A);               // t0 = a-d  | a
  gfp_addsubi_2x4w(td, t0, A);          // td = 2a-d | d

  point_copy_2x4w(T0, P);               // By0 = yP, Bz0 = zP

  point_copy_2x4w(R[0], P);             // P 
  yDBL_2x4w(R[1], P, A);                // [2]P

  for (i = 2; i < s; i++) {
    gfp_mul_2x4w(T0, T0, R[i-1]);       // By0 * yR_{i-1} | Bz0 * zR_{i-1}
    yADD_2x4w(R[i], R[i-1], P, R[i-2]);
  }

  gfp_mul_2x4w(T1, T0, R[s-1]);         // By0 * yR_{s-1} | Bz0 * zR_{s-1}
  mask = u32_iseql(l, 3)^1;
  gfp_cswap_2x4w(T0, T1, mask);       

  vec_blend_2x4w(AD, A, td, 0x0F);      // AD = a | d
  point_copy_2x4w(t0, AD);              // t0 = a | d

  // left-to-right computing a^l and d^l
  lbits -= 1;
  for (i = 1; i <= lbits; i++) {
    gfp_sqr_2x4w(t0, t0);
    if ((l>>(lbits-i)) & 1) gfp_mul_2x4w(t0, t0, AD);
  }
  gfp_sqr_2x4w(T0, T0); gfp_sqr_2x4w(T0, T0); gfp_sqr_2x4w(T0, T0);

  vec_permlh_2x4w(T0, T0);              // T0 = Bz0 | By0
  gfp_mul_2x4w(T1, t0, T0);             // T1 = a'  | d'

  vec_permlh_2x4w(t0, T1);              // t0 = d'    | a'
  gfp_addsubc_2x4w(t0, t0, T1);         // t0 = d'+a' | a'-d'

  vec_blend_2x4w(C, T1, t0, 0x0F);      // C = a' | a'-d'
}

// (2x4)-way y-coordinate isogeny evaluation on twisted Edwards curve.
void yEVAL_2x4w(llpoint_t R, const llpoint_t Q, const llpoint_t P[], const uint8_t k)
{
  llpoint_t T, U, V, t0;
  int i, s;

  point_copy_2x4w(T, Q);

  vec_permlh_2x4w(t0, P[0]);            // t0 = zP0 | yP0
  gfp_mul_2x4w(U, T, t0);               // U  = yQ*zP0 | zQ*yP0
  vec_permlh_2x4w(t0, U);               // t0 = zQ*yP0 |yQ*zP0
  gfp_addsubi_2x4w(R, t0, U);           // R  = yQ*zP0+zT*yP0 | yT*zP0-zT*yP0

  s = primeli[k]>>1;
  for (i = 1; i < s; i++) {
    vec_permlh_2x4w(t0, P[i]);          // t0 = zPi | yPi
    gfp_mul_2x4w(U, T, t0);             // U  = yQ*zPi | zQ*yPi

    vec_permlh_2x4w(t0, U);             // t0 = zQ*yPi | yQ*zPi
    gfp_addsubi_2x4w(V, t0, U);         // V  = yQ*zPi+zQ*yPi | yQ*zPi-yQ*zPi

    gfp_mul_2x4w(R, R, V);
  }

  gfp_sqr_2x4w(R, R);

  vec_permlh_2x4w(t0, T);               // t0 = zQ | yQ
  gfp_addsubc_2x4w(V, T, t0);           // V  = yQ+zQ | zQ-yQ

  gfp_mul_2x4w(V, R, V);      

  vec_permlh_2x4w(t0, V);               // t0 = zV | yV
  gfp_subaddc_2x4w(R, V, t0);
}

