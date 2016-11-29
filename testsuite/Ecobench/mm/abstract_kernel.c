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

DATA_TYPE hw_A[SIZE_A], hw_B[SIZE_B], hw_C[SIZE_C];
DATA_TYPE sw_A[SIZE_A], sw_B[SIZE_B], sw_C[SIZE_C];

void init(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B]) {

	int i, j;

	for (i = 0; i < NI; i++) {
		for (j = 0; j < NK; j++) {
			A[i*NI + j] = ((DATA_TYPE) i*j) / NI;
		}
	}

	for (i = 0; i < NK; i++) {
		for (j = 0; j < NJ; j++) {
			B[i*NK + j] = ((DATA_TYPE) i*(j+1)) / NJ;
		}
	}
	
}

void abstract_kernel (void * arguments) {

	init(hw_A, hw_B);
	mm( hw_A, hw_B, hw_C );

}

void mm2_sw(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B], DATA_TYPE C[SIZE_C])
{

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

int compare_result(char * input_fileName, char * hw_result) {

	init(sw_A, sw_B);

	// Run software version
	mm2_sw( sw_A, sw_B, sw_C);

	int i;
	for (i=0; i<SIZE_C; i++) {
		if(fabs(hw_C[i] - sw_C[i]) > 0.00001){
			return -1; // Failure
		}
	}

	return 0; // Success
}
