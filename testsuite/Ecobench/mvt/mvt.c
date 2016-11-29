#include "mvt.h"

void mvt(DATA_TYPE a[SIZE_A], DATA_TYPE x1[SIZE], /*DATA_TYPE x2[SIZE],*/ DATA_TYPE y1[SIZE]/*, DATA_TYPE y2[SIZE]*/) {

	int i, j, k, l;

	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
       			x1[i] = x1[i] + a[i*N + j] * y1[j];
        	}
    }
/*
	for (k=0; k<N; k++) {
		for (l=0; l<N; l++) {
 		       	x2[k] = x2[k] + a[k*N + l] * y2[l];
      	}
    }
*/
}
