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
#include "convolution2d.h"

void init(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE]);

void abstract_kernel (void * arguments);

void conv2D_sw(DATA_TYPE A[SIZE], DATA_TYPE B[SIZE]);

int compare_result(char * input_fileName, char * hw_result);

#endif // End of ABSTRACT_KERNEL_H
