/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "action.h"
#include <string.h>
#include <stdio.h>

// cmove for small 7-bit integer 
static void u8_cmove(uint8_t *r, const uint8_t a, const uint8_t b)
{
  uint8_t x;

  x = *r ^ a;
  x = x & (-b);
  *r = *r ^ x;
}

// look up the secret exponent in constant time 
static __m512i querye_8x1w(size_t pos, const __m512i e[])
{
  __m512i r = e[0], x;
  uint8_t b;
  size_t i;

  for (i = 0; i < N; i++) {
    b = u32_iseql(i, pos);
    x = VXOR(r, e[i]);
    x = VAND(x, VSET1(-b));
    r = VXOR(r, x);
  }
  return r;
}

// look up the secret exponent in constant time (from [CCC+19] code)
static uint32_t lookup(size_t pos, uint8_t const priv[])
{
  int b;
  uint8_t r = priv[0];
  for(size_t i = 1; i < N; i++)
  {
    b = isequal(i, pos);
    u8_cmove(&r, priv[i], b);
  }
  return r;
}

// check whether the small integer is zero
static __m512i u8_iszero_8x1w(const __m512i a)
{
  __m512i r;

  r = VSUB(VZERO, a);
  r = VSHR(r, 63);
  return VXOR(VAND(r, VSET1(1)), VSET1(1));
}

// The unbatched component which is a low-latency CSIDH class group action 
// will be executed sequentially 8 times. (from [CCC+19] code)
static void action_1w(proj C, const uint8_t *sk, const proj A, const uint8_t* visocnt)
{
  uint8_t ba[N], sizeba = 0, compba[N], sicoba = 0, lastiso;
  uint8_t e[N], bc, ec = 0, fnsh[N] = { 0 }, isocnt[N] = { 0 };
  proj A0, T0, T1, G0, G1, K[HLMAX];
  int total = 0, si, i, j, sum = 0;

  // Initialize variables for computing CSIDH class group action.
  point_copy(A0, A);                          // initialize curve
  memcpy(isocnt, visocnt, N);                 // initialize the isogeny counter
  for (i = 0; i < N; i++) sum += isocnt[i];   // calculate how many isogenies need to be computed 
  for (i = 0; i < N; i++) e[i] = sk[i];       // initialize the secret exponents

  // Since there are only few isogeny computations in the unbatched component, 
  // we will not use SIMBA here         
  // initialize BA and COMPBA                         
  for (i = 0; i < N; i++) {                   
    if (isocnt[i] == 0) {
      compba[sicoba] = i;
      sicoba += 1;
      fnsh[i] = 1;
    }
    else {
      lastiso = i;
      ba[sizeba] = i;
      sizeba += 1;
    }
  }

  // the main loop 
  while (total < sum) {
    elligator(T1, T0, A0);
    yDBL(T0, T0, A0);
    yDBL(T0, T0, A0);
    yDBL(T1, T1, A0);
    yDBL(T1, T1, A0);

    for (i = 0; i < sicoba; i++) {
      yMUL(T0, T0, A0, compba[i]);    
      yMUL(T1, T1, A0, compba[i]);
    }

    for (i = 0; i < sizeba; i++) {
      if (fnsh[ba[i]]) continue;
      else {
        ec = lookup(ba[i], e);
        point_cswap(T0, T1, ec&1);
        point_copy(G0, T0);
        point_copy(G1, T1);

        for (j = i+1; j < sizeba; j++)
          if (!fnsh[ba[j]]) yMUL(G0, G0, A0, ba[j]);

        if ((!point_isinf(G0)) && (!point_isinf(G1))) {
          bc = (uint8_t) (isequal(ec>>1, 0))&1;
          yISOG(K, A0, G0, A0, ba[i]);

          if (ba[i] != lastiso) {
            yEVAL(T0, T0, K, ba[i]);
            yEVAL(T1, T1, K, ba[i]);

            yMUL(T1, T1, A0, ba[i]);         
          }

          e[ba[i]] = ((((ec>>1) - (bc^1))^bc) << 1) ^ ((ec&1)^bc);
          isocnt[ba[i]] -= 1;
          total += 1;
        }
        else {
          yMUL(T1, T1, A0, ba[i]);
        }

        point_cswap(T0, T1, ec&1);

        if (isocnt[ba[i]] == 0 ) {
          fnsh[ba[i]] = 1;
          compba[sicoba] = ba[i];
          sicoba += 1;
        }
      } 
    }
  }
  point_copy(C, A0);
}

// The batched component. 
static void action_8x1w(htpoint_t C, __m512i* visocnt, __m512i* e, const __m512i *sk, const htpoint_t A)
{
  uint8_t ba[NUMBA][MAXSIZEBA], sizeba[NUMBA], compba[NUMBA][N], sicoba[NUMBA], lastiso[NUMBA];
  uint8_t fnsh[N] = { 0 }, isocnt[N] = { 0 };
  __m512i ec = VZERO, inf, bc, vone = VSET1(1), t, x;  
  htpoint A0, A1, T0, T1, G0, G1, K[HLMAX];
  int count = 0, total = 0, i, j, m = 0, numba = NUMBA, n_inf;

  // Initialize SIMBA variables. 
  memcpy(ba, BATCHES, N);
  memcpy(sizeba, SIZEBA, NUMBA);
  memcpy(compba, COMPBA, NUMBA*N);
  memcpy(sicoba, SICOBA, NUMBA);
  memcpy(lastiso, LASTISO, NUMBA);

  // Initialize variables for computing CSIDH class group action.
  point_copy_8x1w(&A0, A);
  memcpy(isocnt, B, N);

  while (total < NUMISO) {
    m = (m+1) % numba;                  // public parameter SIMBA index

    if (count == ROUND*NUMBA) {         // merge the BAs when required rounds finished
      m = 0;                    
      sicoba[0] = 0;          
      sizeba[0] = 0;
      numba = 1;

      for (i = 0; i < N; i++) {
        if (!isocnt[i]) {
          compba[0][sicoba[0]] = i;
          sicoba[0] += 1;
        }
        else {
          lastiso[0] = i;
          ba[0][sizeba[0]] = i;
          sizeba[0] += 1;
        }
      }
    }

    elligator_8x1w(&T1, &T0, &A0);
    yDBL_8x1w(&T0, &T0, &A0); 
    yDBL_8x1w(&T0, &T0, &A0);
    yDBL_8x1w(&T1, &T1, &A0);
    yDBL_8x1w(&T1, &T1, &A0); 

    for (i = 0; i < sicoba[m]; i++) {
      yMUL_8x1w(&T0, &T0, &A0, compba[m][i]); 
      yMUL_8x1w(&T1, &T1, &A0, compba[m][i]);
    }

    for (i = 0; i < sizeba[m]; i++) {
      if (fnsh[ba[m][i]]) continue;
      else {
        ec = querye_8x1w(ba[m][i], e);
        point_cswap_8x1w(&T0, &T1,  VAND(ec, VSET1(1)));
        point_copy_8x1w(&G0, &T0);
        point_copy_8x1w(&G1, &T0);

        for (j = i+1; j < sizeba[m]; j++) 
          if (!fnsh[ba[m][j]]) yMUL_8x1w(&G0, &G0, &A0, ba[m][j]);

        // combined
        inf = point_isinf_8x1w(&G0);
        inf = VOR(inf, point_isinf_8x1w(&G1));
        n_inf = VADDRDC(inf);

        if (n_inf <= 5) {
          bc = u8_iszero_8x1w(VSHR(ec, 1));
          //////////////////////////////////////////////////////////
          // extra-dummy
          htpoint A1, A2, T2, T3, T4, T5;

          point_copy_8x1w(&A1, &A0);
          point_copy_8x1w(&A2, &A0);
          point_copy_8x1w(&T2, &T0);
          point_copy_8x1w(&T3, &T1);
          point_copy_8x1w(&T4, &T0);
          point_copy_8x1w(&T5, &T1);

          yISOG_8x1w(K, &A1, &G0, &A0, ba[m][i]);

          if (ba[m][i] != lastiso[m]) {

            yEVAL_8x1w(&T2, &T2, K, ba[m][i]);
            yEVAL_8x1w(&T3, &T3, K, ba[m][i]);
            
            yMUL_8x1w(&T3, &T3, &A1, ba[m][i]);
          }
          t = VSHR(ec, 1);
          t = VSUB(t, VXOR(bc, vone));
          t = VXOR(t, bc);
          t = VSHL(t, 1);

          t = VXOR(t, VXOR(VAND(ec, vone), bc));
          x = VSUB(inf, VSET1(1));        // not infinity x = 1; infinity x = 0
          t = VXOR(t, e[ba[m][i]]);       // t = t ^ e
          e[ba[m][i]] = VXOR(e[ba[m][i]], VAND(x, t));  // e = e ^ ((t^e)&x)

          isocnt[ba[m][i]] -= 1;
          total += 1;        
          
          point_cswap_8x1w(&T0, &T2, VXOR(inf, VSET1(1)));
          point_cswap_8x1w(&T1, &T3, VXOR(inf, VSET1(1)));
          point_cswap_8x1w(&A0, &A1, VXOR(inf, VSET1(1)));

          yMUL_8x1w(&T4, &T4, &A2, ba[m][i]);
          yMUL_8x1w(&T5, &T5, &A2, ba[m][i]);
          point_cswap_8x1w(&T0, &T4, inf);
          point_cswap_8x1w(&T1, &T5, inf);

          // combined
          visocnt[ba[m][i]] = VSUB(visocnt[ba[m][i]], VXOR(inf, VSET1(1)));
        }
        else {
          yMUL_8x1w(&T0, &T0, &A0, ba[m][i]);
          yMUL_8x1w(&T1, &T1, &A0, ba[m][i]);
        }

        point_cswap_8x1w(&T0, &T1, VAND(ec, vone));
        
        if (!isocnt[ba[m][i]]) {
          fnsh[ba[m][i]] = 1;
          compba[m][sicoba[m]] = ba[m][i];
          sicoba[m] += 1;
        }
        /////////////////////////////////
      }
    }
    count += 1;
  }
  point_copy_8x1w(C, &A0);
}

// The complete CSIDH group action using the combined method. 
void action(htpoint_t C, const __m512i *sk, const htpoint_t A)
{
  __m512i visocnt[N], e[N];
  htpoint A0;
  int i, j;

  // Initialize the vectorized isogeny counter.
  for (i = 0; i < N; i++) visocnt[i] = VSET1(B[i]);

  // Initialize the exponent vector.
  for (i = 0; i < N; i++) e[i] = sk[i];

  // Execute the batched component. 
  action_8x1w(&A0, visocnt, e, sk, A);

  // ---------------------------------------------------------------------------
  // the batched component has completed at this moment
  // we now extract the variables of each instance for the low-latency component

  uint8_t llisocnt[8][N], lle[8][N];
  proj llC[8], llA[8];
  uint32_t a29[8][HT_NWORDS] = {0}, a32[8][16] = {0};
  uint32_t ad29[8][HT_NWORDS] = {0}, ad32[8][16] = {0};
  uint64_t a64[8][8] = {0}, ad64[8][8] = {0}, one[8] = {1};

  // Extract the isogeny counter of each instance.
  for (i = 0; i < N; i++) {
    llisocnt[0][i] = ((uint64_t *)&visocnt[i])[0];
    llisocnt[1][i] = ((uint64_t *)&visocnt[i])[1];
    llisocnt[2][i] = ((uint64_t *)&visocnt[i])[2];
    llisocnt[3][i] = ((uint64_t *)&visocnt[i])[3];
    llisocnt[4][i] = ((uint64_t *)&visocnt[i])[4];
    llisocnt[5][i] = ((uint64_t *)&visocnt[i])[5];
    llisocnt[6][i] = ((uint64_t *)&visocnt[i])[6];
    llisocnt[7][i] = ((uint64_t *)&visocnt[i])[7];
  }

  // convert coefficients from Montgomery domain to number domain,
  // because the Montgomery domain (R = 2^512) of low-latency component is different 
  // from the Montgomery domain (R' = 2^522) of high-throughput component 
  gfp_mont2num_8x1w(A0.y, A0.y);
  gfp_mont2num_8x1w(A0.z, A0.z);

  // extract the curve coefficient for each instance
  get_channel_8x1w(a29[0], A0.y, 0); get_channel_8x1w(a29[1], A0.y, 1);
  get_channel_8x1w(a29[2], A0.y, 2); get_channel_8x1w(a29[3], A0.y, 3);
  get_channel_8x1w(a29[4], A0.y, 4); get_channel_8x1w(a29[5], A0.y, 5);
  get_channel_8x1w(a29[6], A0.y, 6); get_channel_8x1w(a29[7], A0.y, 7);

  get_channel_8x1w(ad29[0], A0.z, 0); get_channel_8x1w(ad29[1], A0.z, 1);
  get_channel_8x1w(ad29[2], A0.z, 2); get_channel_8x1w(ad29[3], A0.z, 3);
  get_channel_8x1w(ad29[4], A0.z, 4); get_channel_8x1w(ad29[5], A0.z, 5);
  get_channel_8x1w(ad29[6], A0.z, 6); get_channel_8x1w(ad29[7], A0.z, 7);

  // convert from radix-29 to radix-64
  for (i = 0; i < 8; i++) {
    mpi_conv_29to32(a32[i], a29[i], 16, HT_NWORDS);
    mpi_conv_29to32(ad32[i], ad29[i], 16, HT_NWORDS);
  }
  for (i = 0; i < 8; i++) {
    mpi_conv_32to64(a64[i], a32[i]);
    mpi_conv_32to64(ad64[i], ad32[i]);
  }

  // form the coefficient for the unbatched component
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      llA[i][0][j]= a64[i][j];
      llA[i][1][j]= ad64[i][j];
    }
    fp_mul(llA[i][0], llA[i][0], R_squared_mod_p);  // convert to Montgomery domain
    fp_mul(llA[i][1], llA[i][1], R_squared_mod_p);  // convert to Montgomery domain
  }

  // extract the secret exponent for each instance
  for (i = 0; i < N; i++) {
    lle[0][i] = ((uint64_t *)&e[i])[0];
    lle[1][i] = ((uint64_t *)&e[i])[1];
    lle[2][i] = ((uint64_t *)&e[i])[2];
    lle[3][i] = ((uint64_t *)&e[i])[3];
    lle[4][i] = ((uint64_t *)&e[i])[4];
    lle[5][i] = ((uint64_t *)&e[i])[5];
    lle[6][i] = ((uint64_t *)&e[i])[6];
    lle[7][i] = ((uint64_t *)&e[i])[7];
  }

  // ---------------------------------------------------------------------------
  // perform sequentially low-latency implementation for each instance

  for (i = 0; i < 8; i++) action_1w(llC[i], lle[i], llA[i], llisocnt[i]);

  // ---------------------------------------------------------------------------
  // the entire CSIDH class group action has completed at this moment

  for (i = 0; i < 8; i++) {
    fp_mul(llC[i][0], llC[i][0], one);
    fp_mul(llC[i][1], llC[i][1], one);
  }

  for(i = 0; i < 8; i++) {
    mpi_conv_64to32(a32[i], llC[i][0]);
    mpi_conv_64to32(ad32[i], llC[i][1]);
    mpi_conv_32to29(a29[i], a32[i], HT_NWORDS, 16);
    mpi_conv_32to29(ad29[i], ad32[i], HT_NWORDS, 16);
  }

  // form the final (8x1)-way result
  for (i = 0; i < HT_NWORDS; i++) {
    C->y[i] = set_vector(a29[7][i], a29[6][i], a29[5][i], a29[4][i], a29[3][i], a29[2][i], a29[1][i], a29[0][i]);
    C->z[i] = set_vector(ad29[7][i], ad29[6][i], ad29[5][i], ad29[4][i], ad29[3][i], ad29[2][i], ad29[1][i], ad29[0][i]);
  }
  gfp_num2mont_8x1w(C->y, C->y);
  gfp_num2mont_8x1w(C->z, C->z);
}

// The functions below are from [CCC+19] code for generating the secret key.

static void cmov(int8_t *r, const int8_t a, uint32_t b)
{
  uint32_t t;
  b = -b; /* Now b is either 0 or 0xffffffff */
  t = (*r ^ a) & b;
  *r ^= t;
}

void random_sk(uint8_t *sk)
{
  uint8_t i, tmp, r;
  int8_t exp, sgn;
  for(i = 0; i < N; i++)
  {

    r = B[i] & 0x1;						// B_i mod 2

    // exp is randomly selected from |[ 0, B ]|
    randombytes(&tmp, 1);
    while ( issmaller((int32_t)B[i], (int32_t)tmp) == -1 )	// constant-time comparison
      randombytes(&tmp, 1);

    exp = (int8_t)tmp;

    // Mapping integers from |[ 0, B |] into
    //                                      |[ -B/2, B/2]| if B is even, or
    //                                      |[ -(B+1)/2, (B-1)/2 ]| if B is odd.
    exp = ( (exp << 1) - (B[i] + r) ) >> 1;

    // Mapping into the set |[-B, B]|.
    exp = (exp << 1) + r;
    sgn = exp >> 7;	// sign of exp

    // Next, to write  key[i] = e || ((1 + sgn)/2)
    cmov(&exp, -exp, sgn == -1);
    sk[i] = (exp << 1) ^ (1 & (1 + sgn));
  }
}
