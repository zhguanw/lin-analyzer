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
	int i, j;

	for (i = 0; i < NI; ++i) {
		for (j = 0; j < NJ; ++j) {
			A[i*NJ + j] = (float)rand()/RAND_MAX;
			B[i*NJ + j] = 0.0f;
       	}
    }
}

void abstract_kernel (void * arguments) {

	init(A, HW_B);

	convolution2d( A, HW_B);

}

void conv2D_sw(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE])
{
	int i, j;
	DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;

	c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
	c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
	c13 = +0.4;  c23 = +0.7;  c33 = +0.10;


	for (i = 1; i < NI - 1; ++i) // 0
	{
		for (j = 1; j < NJ - 1; ++j) // 1
		{
			B[i*NJ + j] = c11 * A[(i - 1)*NJ + (j - 1)]  +  c12 * A[(i + 0)*NJ + (j - 1)]  +  c13 * A[(i + 1)*NJ + (j - 1)]
				+ c21 * A[(i - 1)*NJ + (j + 0)]  +  c22 * A[(i + 0)*NJ + (j + 0)]  +  c23 * A[(i + 1)*NJ + (j + 0)]
				+ c31 * A[(i - 1)*NJ + (j + 1)]  +  c32 * A[(i + 0)*NJ + (j + 1)]  +  c33 * A[(i + 1)*NJ + (j + 1)];
		}
	}
}

int compare_result(char * input_fileName, char * hw_result) {
	int i, j;
	DATA_TYPE* B = (DATA_TYPE*) malloc(SIZE * sizeof(DATA_TYPE));

	for (i = 0; i < NI; ++i) {
		for (j = 0; j < NJ; ++j) {
			B[i*NJ + j] = (DATA_TYPE) 0;
       	}
    }
	
	// Run software version
	conv2D_sw( A, B );

	for (i=0; i<SIZE; i++) {
		if(fabs(HW_B[i] - B[i]) > 0.00001){
			free(B);
			return -1; // Failure
		}
	}

	free(B);
	return 0; // Success
}
