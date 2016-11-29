#include "syrk.h"

void syrk(DATA_TYPE A[SIZE], DATA_TYPE C[SIZE])
{
	int i, j, k;

	/*  C := alpha*A*A' + beta*C */
//	for (i = 0; i < N; i++) {
//		for (j = 0; j < N; j++) {
//			C[i*M + j] *= BETA;
//		}
//	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i*M + j] *= BETA;
			for (k = 0; k < M; k++) {
				C[i*N + j] += ALPHA * A[i*M + k] * A[j*M + k];
			}
		}
	}

}
