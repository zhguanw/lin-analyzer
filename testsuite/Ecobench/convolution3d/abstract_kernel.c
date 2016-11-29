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

DATA_TYPE A[SIZE], HW_B[SIZE];

void init(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE]) {

	int i, j, k;

	for (i = 0; i < NI; ++i) {
		for (j = 0; j < NJ; ++j) {
			for (k = 0; k < NK; ++k) {
				A[i*(NK * NJ) + j*NK + k] = i % 12 + 2 * (j % 7) + 3 * (k % 13);
				B[i*(NK * NJ) + j*NK + k] = (DATA_TYPE) 0;
			}
		}
	}
}

void abstract_kernel (void * arguments) {

	init(A, HW_B);

	convolution3d( A, HW_B);

}

void conv3D_sw(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE])
{
	int i, j, k;
	DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;

	c11 = +2;  c21 = +5;  c31 = -8;
	c12 = -3;  c22 = +6;  c32 = -9;
	c13 = +4;  c23 = +7;  c33 = +10;

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

int compare_result(char * input_fileName, char * hw_result) {

	DATA_TYPE* B = (DATA_TYPE*) malloc(SIZE * sizeof(DATA_TYPE));
	init(A, B);
	// Run software version
	conv3D_sw( A, B );

	int i;
	for (i=0; i<SIZE; i++) {
		if(fabs(HW_B[i] - B[i]) > 0.00001){
			free(B);
			return -1; // Failure
		}
	}

	free(B);
	return 0; // Success
}
