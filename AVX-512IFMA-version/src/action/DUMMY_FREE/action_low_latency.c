/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "action_low_latency.h"
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
static uint8_t querye_2x4w(size_t pos, uint8_t const e[])
{
  uint8_t r = e[0], b;
  size_t i;

  for(i = 1; i < N; i++)
  {
    b = u32_iseql(i, pos);
    u8_cmove(&r, e[i], b);
  }
  return r;
}

// check whether the small 7-bit integer is zero 
static uint8_t u8_iszero(const uint8_t a)
{
  uint8_t r;

  r = 0-a;
  r = r >> 7;
  return (r&1)^1;
}

// low-latency CSIDH class group action based on [CCC+19] 
void action_2x4w(llpoint_t C, const uint8_t *sk, const llpoint_t A)
{
  uint8_t ba[NUMBA][MAXSIZEBA], sizeba[NUMBA], compba[NUMBA][N], sicoba[NUMBA], lastiso[NUMBA];
  uint8_t e[N], bc, ec = 0, fnsh[N] = { 0 }, isocnt[N] = { 0 };
  llpoint_t A0, T0, T1, T2, T3, G0, G1, K[HLMAX], Z;
  int count = 0, total = 0, si, i, j, m = 0, numba = NUMBA;

  // Initialize SIMBA variables. 
  memcpy(ba, BATCHES, N);
  memcpy(sizeba, SIZEBA, NUMBA);
  memcpy(compba, COMPBA, NUMBA*N);
  memcpy(sicoba, SICOBA, NUMBA);
  memcpy(lastiso, LASTISO, NUMBA);

  // Initialize variables for computing CSIDH class group action.
  point_copy_2x4w(A0, A);
  memcpy(isocnt, B, N);
  for (i = 0; i < N; i++) e[i] = sk[i]; // initialize the secret exponent

  // the main loop
  while (total < NUMISO) {
    m = (m+1) % numba;                  // public parameter SIMBA index

    if (count == NUMBA*ROUND) {         // merge the BAs when required rounds finished
      m = 0;                    
      sicoba[0] = 0;          
      sizeba[0] = 0;
      numba = 1;

      for (i = 0; i < N; i++) {
        if (isocnt[i] == 0) {
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

    elligator_2x4w(T1, T0, A0);
    yDBL_2x4w(T0, T0, A0);
    yDBL_2x4w(T0, T0, A0);
    yDBL_2x4w(T1, T1, A0);
    yDBL_2x4w(T1, T1, A0);

    for (i = 0; i < sicoba[m]; i++) {
      yMUL_2x4w(T0, T0, A0, compba[m][i]); 
      yMUL_2x4w(T1, T1, A0, compba[m][i]);
    }

    for (i = 0; i < sizeba[m]; i++) {
      if (fnsh[ba[m][i]]) continue;
      else {
        ec = querye_2x4w(ba[m][i], e);
        point_cswap_2x4w(T0, T1, ec&1);
        point_copy_2x4w(G0, T0);
        point_copy_2x4w(G1, T1);

        for (j = i+1; j < sizeba[m]; j++)
          if (!fnsh[ba[m][j]]) yMUL_2x4w(G0, G0, A0, ba[m][j]);

        if ((!point_isinf_2x4w(G0)) && (!point_isinf_2x4w(G1))) {
          bc = u8_iszero(ec>>1);
          yISOG_2x4w(K, A0, G0, A0, ba[m][i]);

          if (ba[m][i] != lastiso[m]) {
            yEVAL_2x4w(T0, T0, K, ba[m][i]);
            yEVAL_2x4w(T1, T1, K, ba[m][i]);

            yMUL_2x4w(T1, T1, A0, ba[m][i]);  
          }
          
          e[ba[m][i]] = ((((ec>>1) - (bc^1))^bc) << 1) ^ ((ec&1)^bc);
          isocnt[ba[m][i]] -= 1;
          total += 1;
        }
        else {
          yMUL_2x4w(T1, T1, A0, ba[m][i]);
        }

        point_cswap_2x4w(T0, T1, ec&1);

        if (isocnt[ba[m][i]] == 0 ) {
          fnsh[ba[m][i]] = 1;
          compba[m][sicoba[m]] = ba[m][i];
          sicoba[m] += 1;
        }
      } 
    }
    count += 1;
  }
  point_copy_2x4w(C, A0);
}

// The functions below are from [CCC+19] code for generating the secret key.

static void cmov(int8_t *r, const int8_t a, uint32_t b)
{
  uint32_t t;
  b = -b; /* Now b is either 0 or 0xffffffff */
  t = (*r ^ a) & b;
  *r ^= t;
}

// constant-time comparison: -1 if x < y, 0 otherwise.
static int32_t issmaller(int32_t x, int32_t y)
{
  int32_t xy = x ^ y;
  int32_t c = x - y;
  c ^= xy & (c ^ x);
  return (c >> 31);
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
