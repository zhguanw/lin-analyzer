#ifndef COMMON_H
#define COMMON_H

typedef float DATA_TYPE;

/* Problem size */
#define NI 128
#define NJ 128
#define NK 128
#define SIZE_A NI*NK
#define SIZE_B NK*NJ
#define SIZE_C NI*NJ

/* Declared constant values for ALPHA and BETA (same as values in PolyBench 2.0) */
#define ALPHA 32412.0f
#define BETA 2123.0f

#endif
