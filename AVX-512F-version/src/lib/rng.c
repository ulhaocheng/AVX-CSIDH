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
