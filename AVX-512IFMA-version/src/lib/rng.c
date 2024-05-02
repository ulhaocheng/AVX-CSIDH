/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "rng.h"

void randombytes(void *r, size_t len)
{
  static int fd = -1;
  ssize_t n;

  if ((fd < 0) && (0 > (fd = open("/dev/urandom", O_RDONLY)))) exit(1);
  for (size_t i = 0; i < len; i += n)
    if (0 >= (n = read(fd, (char *)r + i, len-i))) exit(2);
}

// generate a random radix-64 field element
void mpi64_random(uint64_t *x)
{
  int i;

  while (1) {
    randombytes(x, sizeof(uint64_t)*8);
    uint64_t m = ((uint64_t) 1 << 63) - 1;  // prime p is 511-bit 
    x[7] &= m;

    for (i = 7; i >= 0; i--)
      if (x[i] < u64_p[i]) return;
      else if (x[i] > u64_p[i]) break;
  }
}
