/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#ifndef _RNG_H
#define _RNG_H

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include "gfparith.h"

void randombytes(void *r, size_t len);
void mpi64_random(uint64_t *r);

#endif
