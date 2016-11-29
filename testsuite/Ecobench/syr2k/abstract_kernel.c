/**********************************************************************
Author: Guanwen Zhong
Associated Filename: abstract_kernel.c
Purpose: The abstract kernel wrapper for 2D Convolution kernel
Origin: 2D Convolution benchmarch from Polybench GPU v1.0

This benchmark is modified from cpu-based 2D Convolution application in
Polybench GPU v1.0. For more information, please check:
http://web.cse.ohio-state.edu/~pouchet/software/polybench/GPU/index.html

**********************************************************************
Reserved.
© Copyright

**********************************************************************/

#include "abstract_kernel.h"

DATA_TYPE A[SIZE], B[SIZE], HW_C[SIZE];

void init(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE], DATA_TYPE C[SIZE]) {

	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i*N + j] = ((DATA_TYPE) i*j + 2) / N;
		}

		for (j = 0; j < M; j++) {
			A[i*N + j] = ((DATA_TYPE) i*j) / N;
			B[i*N + j] = ((DATA_TYPE) i*j + 1) / N;
		}
	}

}

void abstract_kernel (void * arguments) {

	init(A, B, HW_C);

	syr2k( A, B, HW_C);

}

void syr2k_sw(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE], DATA_TYPE C[SIZE])
{
	int i, j, k;

  	for (i = 0; i < N; i++) {
   		for (j = 0; j < N; j++) {
			C[i*N + j] *= BETA;
		}
	}

  	for (i = 0; i < N; i++) {
   		for (j = 0; j < N; j++) {
			for (k = 0; k < M; k++) {
	  			C[i*N + j] += ALPHA * A[i*M + k] * B[j*M + k];
	 		 	C[i*N + j] += ALPHA * B[i*M + k] * A[j*M + k];
			}
		}
	}
}

int compare_result(char * input_fileName, char * hw_result) {

	DATA_TYPE* C = (DATA_TYPE*) malloc(SIZE * sizeof(DATA_TYPE));
	init(A, B, C);

	// Run software version
	syr2k_sw( A, B, C);

	int i;
	for (i=0; i<SIZE; i++) {
		if(fabs(HW_C[i] - C[i]) > 0.00001){
			free(C);
			return -1; // Failure
		}
	}

	free(C);
	return 0; // Success
}
