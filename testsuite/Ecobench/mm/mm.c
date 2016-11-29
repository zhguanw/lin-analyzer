#include "mm.h"

void mm(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B], DATA_TYPE C[SIZE_C]) {

	int i, j, k;
	
  	for (i = 0; i < NI; i++) {
		for (j = 0; j < NJ; j++) {
			DATA_TYPE tmpSum1 = 0.0f;
			for (k = 0; k < NK; ++k) {
				tmpSum1 += A[i*NK + k] * B[k*NJ + j];
			}
			C[i*NJ + j] = tmpSum1;
		}
	}

}
