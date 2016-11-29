#include "atax.h"

void atax(DATA_TYPE A[SIZE_A], DATA_TYPE x[SIZE_X], /* DATA_TYPE y[SIZE_Y],*/ DATA_TYPE tmp[SIZE_TMP]) {

	int i,j;

	for (i = 0; i < NX; i++) {
		DATA_TYPE sum_tmp = (DATA_TYPE) 0;
      	for (j = 0; j < NY; j++) {
			sum_tmp += A[i*NY + j] * x[j];
		}
		tmp[i] = sum_tmp;
    }
	/*
	for (j = 0; j < NY; j++) {
		for (i = 0; i < NX; i++) {
			y[j] = y[j] + A[i*NY + j] * tmp[i];
		}
	}*/

}
