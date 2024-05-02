/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "utils.h"

void mpi_print(const char *c, const uint32_t *a, int len)
{
  int i;

  printf("%s", c);
  for (i = len-1; i > 0; i--) printf("%08X", a[i]);
  printf("%08X\n", a[0]);
}

void mpi_conv_29to32(uint32_t *r, const uint32_t *a, int rlen, int alen)
{
  int i, j, bits_in_word, bits_to_shift;
  uint32_t word;

  i = j = 0;
  bits_in_word = bits_to_shift = 0;
  word = 0;
  while ((i < rlen) && (j < alen)) {
    word |= (a[j] << bits_in_word);
    bits_to_shift = (32 - bits_in_word);
    bits_in_word += 29;
    if (bits_in_word >= 32) {
      r[i++] = word;
      word = ((bits_to_shift > 0) ? (a[j] >> bits_to_shift) : 0);
      bits_in_word = ((bits_to_shift > 0) ? (29 - bits_to_shift) : 0);
    }
    j++;
  }
  if (i < rlen) r[i++] = word;
  for (; i < rlen; i++) r[i] = 0;
}

void mpi_conv_32to29(uint32_t *r, const uint32_t *a, int rlen, int alen)
{
  int i, j, shr_pos, shl_pos;
  uint32_t word, temp;

  i = j = 0;
  shr_pos=32; shl_pos=0;
  temp = 0;
  while ((i < rlen) && (j < alen)) {
    word = ((temp >> shr_pos) | (a[j] << shl_pos));
    r[i] = (word & HT_BMASK);
    shr_pos-=3, shl_pos+=3;
    if ((shr_pos > 0) && (shl_pos < 32)) temp = a[j++];
    if (shr_pos <= 0) shr_pos += 32;
    if (shl_pos >= 32) shl_pos -= 32;
    // Any shift past 31 is undefined!
    if (shr_pos == 32) temp = 0;
    i++;
  }
  if (i < rlen) r[i++] = ((temp >> shr_pos) & HT_BMASK);
  for (; i < rlen; i++) r[i] = 0;
}

void mpi_conv_32to64(uint64_t *r, const uint32_t *a)
{
  r[0] = (uint64_t) a[0] | ((uint64_t) a[1]<<32); 
  r[1] = (uint64_t) a[2] | ((uint64_t) a[3]<<32); 
  r[2] = (uint64_t) a[4] | ((uint64_t) a[5]<<32); 
  r[3] = (uint64_t) a[6] | ((uint64_t) a[7]<<32); 
  r[4] = (uint64_t) a[8] | ((uint64_t) a[9]<<32); 
  r[5] = (uint64_t) a[10] | ((uint64_t) a[11]<<32); 
  r[6] = (uint64_t) a[12] | ((uint64_t) a[13]<<32); 
  r[7] = (uint64_t) a[14] | ((uint64_t) a[15]<<32); 
}

void mpi_conv_64to32(uint32_t *r, const uint64_t *a)
{
  r[0] = (uint32_t) a[0];
  r[1] = (uint32_t)(a[0]>>32);
  r[2] = (uint32_t) a[1];
  r[3] = (uint32_t)(a[1]>>32);
  r[4] = (uint32_t) a[2];
  r[5] = (uint32_t)(a[2]>>32);
  r[6] = (uint32_t) a[3];
  r[7] = (uint32_t)(a[3]>>32);
  r[8] = (uint32_t) a[4];
  r[9] = (uint32_t)(a[4]>>32);
  r[10] = (uint32_t) a[5];
  r[11] = (uint32_t)(a[5]>>32);
  r[12] = (uint32_t) a[6];
  r[13] = (uint32_t)(a[6]>>32);
  r[14] = (uint32_t) a[7];
  r[15] = (uint32_t)(a[7]>>32);
}

__m512i set_vector(const uint32_t a7, const uint32_t a6, const uint32_t a5, const uint32_t a4, 
                   const uint32_t a3, const uint32_t a2, const uint32_t a1, const uint32_t a0)
{
  __m512i r;

  ((uint64_t *)&r)[0] = a0; ((uint64_t *)&r)[1] = a1;
  ((uint64_t *)&r)[2] = a2; ((uint64_t *)&r)[3] = a3;
  ((uint64_t *)&r)[4] = a4; ((uint64_t *)&r)[5] = a5;
  ((uint64_t *)&r)[6] = a6; ((uint64_t *)&r)[7] = a7;

  return r;
}

void get_channel_8x1w(uint32_t *r, const htfe_t a, const int ch) 
{
  int i;

  for(i = 0; i < HT_NWORDS; i++){
    r[i] = (uint32_t) ((uint64_t *)&a[i])[ch];
  }
}
