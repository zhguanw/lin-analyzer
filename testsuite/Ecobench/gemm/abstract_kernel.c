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

DATA_TYPE A[SIZE_A], B[SIZE_B], HW_C[SIZE_C];

void init(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B], DATA_TYPE C[SIZE_C]) {

	int i, j, k;

  	for (i = 0; i < NI; i++) {
    	for (j = 0; j < NK; j++) {
      		A[i*NK + j] = ((DATA_TYPE) i*j) / NI;
		}
	}

  	for (i = 0; i < NK; i++) {
    	for (j = 0; j < NJ; j++) {
      		B[i*NJ + j] = ((DATA_TYPE) i*j + 1) / NJ;
		}
	}

  	for (i = 0; i < NI; i++) {
    	for (j = 0; j < NJ; j++) {
      		C[i*NJ + j] = ((DATA_TYPE) i*j + 2) / NJ;
		}
	}
}

void abstract_kernel (void * arguments) {

	init(A, B, HW_C);

	gemm( A, B, HW_C);

}

void gemm_sw(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B], DATA_TYPE C[SIZE_C])
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

int compare_result(char * input_fileName, char * hw_result) {

	DATA_TYPE* C = (DATA_TYPE*) malloc(SIZE_C * sizeof(DATA_TYPE));
	init(A, B, C);

	// Run software version
	gemm_sw( A, B, C);

	int i;
	for (i=0; i<SIZE_C; i++) {
		if(fabs(HW_C[i] - C[i]) > 0.00001){
			free(C);
			return -1; // Failure
		}
	}

	free(C);
	return 0; // Success
}
