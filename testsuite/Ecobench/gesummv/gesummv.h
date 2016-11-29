#ifndef CONVOLUTION2D_H
#define CONVOLUTION2D_H

#include "common.h"

void gesummv(DATA_TYPE A[SIZE_AB], DATA_TYPE B[SIZE_AB], DATA_TYPE x[SIZE_XY], DATA_TYPE y[SIZE_XY], DATA_TYPE tmp[SIZE_XY]);

#endif // End of CONVOLUTION2D_H
