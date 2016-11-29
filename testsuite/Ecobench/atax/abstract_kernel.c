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
� Copyright

**********************************************************************/

#include "abstract_kernel.h"

DATA_TYPE hw_A[SIZE_A], hw_x[SIZE_X], hw_y[SIZE_Y], hw_tmp[SIZE_TMP];
DATA_TYPE sw_A[SIZE_A], sw_x[SIZE_X], sw_y[SIZE_Y], sw_tmp[SIZE_TMP];

void init(DATA_TYPE x[SIZE_X], DATA_TYPE A[SIZE_A], DATA_TYPE y[SIZE_Y], DATA_TYPE tmp[SIZE_TMP]) {

	int i, j;

	for (i = 0; i < NX; i++) {
		x[i] = i * M_PI;
		for (j = 0; j < NY; j++) {
			A[i*NY + j] = ((DATA_TYPE) i*(j)) / NX;
		}
	}
	
	for (i= 0; i < NY; i++) {
		y[i] = 0;
	}
	
	for (i = 0; i < NX; i++) {
      	tmp[i] = 0;
	}
}

void abstract_kernel (void * arguments) {

	init( hw_x, hw_A, hw_y, hw_tmp );

	atax( hw_A, hw_x, /*hw_y,*/ hw_tmp );

}

void atax_sw(DATA_TYPE A[SIZE_A], DATA_TYPE x[SIZE_X], /*DATA_TYPE y[SIZE_Y],*/ DATA_TYPE tmp[SIZE_TMP]) {
	int i,j;
/*
	for (i= 0; i < NY; i++) {
		y[i] = 0;
	}
*/
	for (i = 0; i < NX; i++) {
      	tmp[i] = 0;

      	for (j = 0; j < NY; j++) {
			tmp[i] = tmp[i] + A[i*NY + j] * x[j];
		}
/*
		for (j = 0; j < NY; j++) {
			y[j] = y[j] + A[i*NY + j] * tmp[i];
		}
		*/
    }

}

int compare_result(char * input_fileName, char * hw_result) {

	init( sw_x, sw_A, sw_y, sw_tmp );

	// Run software version
	atax_sw( sw_A, sw_x, /*sw_y,*/ sw_tmp );

	int i;
	for (i=0; i<SIZE_TMP; i++) {
		if(fabs(hw_tmp[i] - sw_tmp[i]) > 0.00001){
			return -1; // Failure
		}
	}
	/*
	for (i=0; i<SIZE_Y; i++) {
		if(fabs(hw_y[i] - sw_y[i]) > 0.00001){
			return -1; // Failure
		}
	}*/

	return 0; // Success
}
