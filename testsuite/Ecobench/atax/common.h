#ifndef COMMON_H
#define COMMON_H

typedef float DATA_TYPE;

/* Problem size */
#define NX 128
#define NY 128

#define SIZE_A NX*NY
#define SIZE_X NY
#define SIZE_Y NY
#define SIZE_TMP NX

#ifndef M_PI
#define M_PI 3.14159
#endif

/* Declared constant values for ALPHA and BETA (same as values in PolyBench 2.0) */
#define ALPHA 32412
#define BETA 2123

#endif
