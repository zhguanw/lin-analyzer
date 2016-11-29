/**********************************************************************
 * Author: Guanwen Zhong
 * Associated Filename: abstract_kernel.h
 * Purpose: The abstract_kernel file of 2D Convolution kernel
 * Origin: Ecobench benchmarch v1.0

This benchmark is modified from Polybench_GPU v1.0. For more information, please check:
http://xxxxx

**********************************************************************
Reserved.
© Copyright

**********************************************************************/
#ifndef ABSTRACT_KERNEL_H
#define ABSTRACT_KERNEL_H

/*********************************************************************
 * The following part is an implemetation of an abstraction function.
 *
**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "mm.h"

void init(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B]);

void abstract_kernel (void * arguments);

void mm2_sw(DATA_TYPE A[SIZE_A], DATA_TYPE B[SIZE_B], DATA_TYPE C[SIZE_C]);

int compare_result(char * input_fileName, char * hw_result);

#endif // End of ABSTRACT_KERNEL_H
