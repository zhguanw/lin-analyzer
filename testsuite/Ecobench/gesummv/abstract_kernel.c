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

DATA_TYPE hw_a[SIZE_AB], hw_b[SIZE_AB], hw_x[SIZE_XY], hw_y[SIZE_XY], hw_temp[SIZE_XY];
DATA_TYPE sw_a[SIZE_AB], sw_b[SIZE_AB], sw_x[SIZE_XY], sw_y[SIZE_XY], sw_temp[SIZE_XY];

void init(DATA_TYPE A[SIZE_AB], DATA_TYPE x[SIZE_XY]) {

  	int i, j;

 	for (i = 0; i < N; i++) {
    	x[i] = ((DATA_TYPE) i) / N;

		for (j = 0; j < N; j++) {
			A[i*N + j] = ((DATA_TYPE) i*j) / N;
		}
    }

}

void abstract_kernel (void * arguments) {

	init(hw_a, hw_x);

	gesummv( hw_a, hw_b, hw_x, hw_y, hw_temp);

}

void gesummv_sw(DATA_TYPE A[SIZE_AB], DATA_TYPE B[SIZE_AB], DATA_TYPE x[SIZE_XY], DATA_TYPE y[SIZE_XY], DATA_TYPE tmp[SIZE_XY]) {
	int i,j;

	for (i = 0; i < N; i++) {
		tmp[i] = 0;
		y[i] = 0;
		for (j = 0; j < N; j++) {
			tmp[i] = A[i*N + j] * x[j] + tmp[i];
			y[i] = B[i*N + j] * x[j] + y[i];
		}

		y[i] = ALPHA * tmp[i] + BETA * y[i];
	}

}

int compare_result(char * input_fileName, char * hw_result) {

	init(sw_a, sw_x);

	// Run software version
	gesummv_sw( sw_a, sw_b, sw_x, sw_y, sw_temp);

	int i;
	for (i=0; i<SIZE_XY; i++) {
		if(fabs(hw_y[i] - sw_y[i]) > 0.00001){
			return -1; // Failure
		}
	}

	return 0; // Success
}
