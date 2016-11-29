#include "syr2k.h"

void syr2k(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE], DATA_TYPE C[SIZE])
{
	int i, j, k;

//  	for (i = 0; i < N; i++) {
//   		for (j = 0; j < N; j++) {
//			C[i*N + j] *= BETA;
//		}
//	}

  	for (i = 0; i < N; i++) {
   		for (j = 0; j < N; j++) {
			C[i*N + j] *= BETA;
			for (k = 0; k < M; k++) {
	  			C[i*N + j] += ALPHA * A[i*M + k] * B[j*M + k];
	 		 	C[i*N + j] += ALPHA * B[i*M + k] * A[j*M + k];
			}
		}
	}

}
