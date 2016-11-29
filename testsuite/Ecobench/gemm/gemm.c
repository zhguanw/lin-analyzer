#include "gemm.h"

void gemm(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B], DATA_TYPE C[SIZE_C])
{
	int i,j,k;

	for (i = 0; i < NI; i++) {
		for (j = 0; j < NJ; j++) {
			C[i*NJ + j] *= BETA;

			for (k = 0; k < NK; ++k) {
				C[i*NJ + j] += ALPHA * A[i*NK + k] * B[k*NJ + j];
			}
		}
	}

}
