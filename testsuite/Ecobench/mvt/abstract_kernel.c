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

DATA_TYPE hw_a[SIZE_A], hw_x1[SIZE], hw_x2[SIZE], hw_y1[SIZE], hw_y2[SIZE];
DATA_TYPE sw_a[SIZE_A], sw_x1[SIZE], sw_x2[SIZE], sw_y1[SIZE], sw_y2[SIZE];

void init(DATA_TYPE a[SIZE_A], DATA_TYPE x1[SIZE], DATA_TYPE x2[SIZE], DATA_TYPE y_1[SIZE], DATA_TYPE y_2[SIZE]) {

	int i, j;

	for (i=0; i<N; i++) {
		x1[i] = 0.0;
		x2[i] = 0.0;
		y_1[i] = 0.0;
		y_2[i] = 0.0;

		for (j=0; j<N; j++) {
			a[i*N + j] = (DATA_TYPE)(i+j+1.0)/N;
		}
	}

}

void abstract_kernel (void * arguments) {

	init(hw_a, hw_x1, hw_x2, hw_y1, hw_y2);

	mvt(hw_a, hw_x1, /*hw_x2,*/ hw_y1/*, hw_y2*/);

}

void mvt_sw(DATA_TYPE a[SIZE_A], DATA_TYPE x1[SIZE], DATA_TYPE x2[SIZE], DATA_TYPE y1[SIZE], DATA_TYPE y2[SIZE])
{

	int i, j, k, l;

	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
       			x1[i] = x1[i] + a[i*N + j] * y1[j];
        	}
    }

	for (k=0; k<N; k++) {
		for (l=0; l<N; l++) {
 		       	x2[k] = x2[k] + a[k*N + l] * y2[l];
      	}
    }

}

int compare_result(char * input_fileName, char * hw_result) {

	init(sw_a, sw_x1, sw_x2, sw_y1, sw_y2);

	// Run software version
	mvt_sw( sw_a, sw_x1, sw_x2, sw_y1, sw_y2 );

	int i;
	for (i=0; i<SIZE; i++) {
		if(fabs(hw_x1[i] - sw_x1[i]) > 0.00001){
			return -1; // Failure
		}
	}
/*
	for (i=0; i<SIZE; i++) {
		if(fabs(hw_x2[i] - sw_x2[i]) > 0.00001){
			return -1; // Failure
		}
	}
*/
	return 0; // Success
}
