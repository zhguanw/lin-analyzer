#ifndef COMMON_H
#define COMMON_H

typedef float DATA_TYPE;

/* Problem size */
#define NX 256
#define NY 256

#define SIZE_A NX*NY
#define SIZE_R NX
#define SIZE_S NY
#define SIZE_P NY
#define SIZE_Q NX

#ifndef M_PI
#define M_PI 3.14159
#endif

/* Declared constant values for ALPHA and BETA (same as values in PolyBench 2.0) */
#define ALPHA 1
#define BETA 1

#endif
