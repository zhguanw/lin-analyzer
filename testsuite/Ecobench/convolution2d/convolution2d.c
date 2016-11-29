#include "convolution2d.h"

void convolution2d(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE]) {
	int i, j;
	DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;
#ifdef FLOAT_TYPE
	c11 = +0.2f;  c21 = +0.5f;  c31 = -0.8f;
	c12 = -0.3f;  c22 = +0.6f;  c32 = -0.9f;
	c13 = +0.4f;  c23 = +0.7f;  c33 = +0.10f;
#else // Double type
	c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
	c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
	c13 = +0.4;  c23 = +0.7;  c33 = +0.10;
#endif


	for (i = 1; i < NI - 1; ++i) {
		for (j = 1; j < NJ - 1; ++j) {
			B[i*NJ + j] = c11 * A[(i - 1)*NJ + (j - 1)]  +  c12 * A[(i + 0)*NJ + (j - 1)]  +  c13 * A[(i + 1)*NJ + (j - 1)]
				+ c21 * A[(i - 1)*NJ + (j + 0)]  +  c22 * A[(i + 0)*NJ + (j + 0)]  +  c23 * A[(i + 1)*NJ + (j + 0)]
				+ c31 * A[(i - 1)*NJ + (j + 1)]  +  c32 * A[(i + 0)*NJ + (j + 1)]  +  c33 * A[(i + 1)*NJ + (j + 1)];
		}
	}
}
