#include "convolution3d.h"

void convolution3d(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE]) {
	int i, j, k;
	DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;
#ifdef FLOAT_TYPE
	c11 = +2.0f;  c21 = +5.0f;  c31 = -8.0f;
	c12 = -3.0f;  c22 = +6.0f;  c32 = -9.0f;
	c13 = +4.0f;  c23 = +7.0f;  c33 = +10.0f;
#else // Double type
	c11 = +2.0;  c21 = +5.0;  c31 = -8.0;
	c12 = -3.0;  c22 = +6.0;  c32 = -9.0;
	c13 = +4.0;  c23 = +7.0;  c33 = +10.0;
#endif

	for (i = 1; i < NI - 1; ++i) {
		for (j = 1; j < NJ - 1; ++j) {
			for (k = 1; k < NK -1; ++k) {
				//printf("i:%d\nj:%d\nk:%d\n", i, j, k);
				B[i*(NK * NJ) + j*NK + k] = c11 * A[(i - 1)*(NK * NJ) + (j - 1)*NK + (k - 1)]  +  c13 * A[(i + 1)*(NK * NJ) + (j - 1)*NK + (k - 1)]
					     +   c21 * A[(i - 1)*(NK * NJ) + (j - 1)*NK + (k - 1)]  +  c23 * A[(i + 1)*(NK * NJ) + (j - 1)*NK + (k - 1)]
					     +   c31 * A[(i - 1)*(NK * NJ) + (j - 1)*NK + (k - 1)]  +  c33 * A[(i + 1)*(NK * NJ) + (j - 1)*NK + (k - 1)]
					     +   c12 * A[(i + 0)*(NK * NJ) + (j - 1)*NK + (k + 0)]  +  c22 * A[(i + 0)*(NK * NJ) + (j + 0)*NK + (k + 0)]
					     +   c32 * A[(i + 0)*(NK * NJ) + (j + 1)*NK + (k + 0)]  +  c11 * A[(i - 1)*(NK * NJ) + (j - 1)*NK + (k + 1)]
					     +   c13 * A[(i + 1)*(NK * NJ) + (j - 1)*NK + (k + 1)]  +  c21 * A[(i - 1)*(NK * NJ) + (j + 0)*NK + (k + 1)]
					     +   c23 * A[(i + 1)*(NK * NJ) + (j + 0)*NK + (k + 1)]  +  c31 * A[(i - 1)*(NK * NJ) + (j + 1)*NK + (k + 1)]
					     +   c33 * A[(i + 1)*(NK * NJ) + (j + 1)*NK + (k + 1)];
			}
		}
	}
}
