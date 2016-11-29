#include "gesummv.h"

void gesummv(DATA_TYPE A[SIZE_AB], DATA_TYPE B[SIZE_AB], DATA_TYPE x[SIZE_XY], DATA_TYPE y[SIZE_XY], DATA_TYPE tmp[SIZE_XY])
{
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
