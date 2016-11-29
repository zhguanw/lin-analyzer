#include "bicg.h"

void bicg(DATA_TYPE A[SIZE_A], DATA_TYPE r[SIZE_R], DATA_TYPE s[SIZE_S], DATA_TYPE p[SIZE_P], DATA_TYPE q[SIZE_Q]) {

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
