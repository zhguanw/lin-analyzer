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

DATA_TYPE hw_a[SIZE_A], hw_r[SIZE_R], hw_s[SIZE_S], hw_p[SIZE_P], hw_q[SIZE_Q];
DATA_TYPE sw_a[SIZE_A], sw_r[SIZE_R], sw_s[SIZE_S], sw_p[SIZE_P], sw_q[SIZE_Q];

void init(DATA_TYPE A[SIZE_A], DATA_TYPE p[SIZE_P], DATA_TYPE r[SIZE_P]) {

	int i, j;

  	for (i = 0; i < NX; i++) {
		r[i] = i * M_PI;

		for (j = 0; j < NY; j++) {
			A[i*NY + j] = ((DATA_TYPE) i*j) / NX;
		}
 	}

	for (i = 0; i < NY; i++) {
    	p[i] = i * M_PI;
	}

}

void abstract_kernel (void * arguments) {

	init(hw_a, hw_p, hw_r);

	bicg( hw_a, hw_r, hw_s, hw_p, hw_q);

}

void bicg_sw(DATA_TYPE A[SIZE_A], DATA_TYPE r[SIZE_R], DATA_TYPE s[SIZE_S], DATA_TYPE p[SIZE_P], DATA_TYPE q[SIZE_Q]) {

	int i,j;

  	for (i = 0; i < NY; i++) {
		s[i] = 0.0;
	}

	for (i = 0; i < NX; i++) {
		q[i] = 0.0;
		for (j = 0; j < NY; j++) {
			s[j] = s[j] + r[i] * A[i*NY + j];
			q[i] = q[i] + A[i*NY + j] * p[j];
	  	}
	}

}

int compare_result(char * input_fileName, char * hw_result) {

	init(sw_a, sw_p, sw_r);

	// Run software version
	bicg_sw( sw_a, sw_r, sw_s, sw_p, sw_q);

	int i;
	for (i=0; i<SIZE_S; i++) {
		if(fabs(hw_s[i] - sw_s[i]) > 0.00001){
			return -1; // Failure
		}
	}

	for (i=0; i<SIZE_Q; i++) {
		if(fabs(hw_q[i] - sw_q[i]) > 0.00001){
			return -1; // Failure
		}
	}

	return 0; // Success
}
