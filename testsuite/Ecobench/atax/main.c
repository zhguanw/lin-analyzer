/**********************************************************************
Author: Guanwen Zhong
Associated Filename: main.c
Purpose: Testbench file for stencil 3D
Origin: Stencil 3D benchmarch from Parboil 2.5

This benchmark is modified from cpu-based stencil application in
Parboil 2.5. For more information, please check:
http://impact.crhc.illinois.edu/Parboil/parboil_download_page.aspx

**********************************************************************
Reserved.
© Copyright

**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "abstract_kernel.h"

#define SUCCESS 0
#define FAILURE -1

//void abstract_kernel (void * arguments);

int main(int argc, char** argv) {

	printf("CPU-based 7 points stencil codes****\n");
	printf("Original version is from Parboil Benchmark version 2.5\n");
	printf("This version maintained by Guanwen Zhong <guanwen@comp.nus.edu.sg>  ****\n");
	printf("It is better to provide input path as an argument\n");
	
	/// Read input data from input.data
	char * inputData;
	char * goldenData;

	char inputData_name[256];
	char goldenData_name[256];
#ifdef BIG_DATASET
	if (argc == 2) {
		strcpy(inputData_name, argv[1]);
		//strcpy(goldenData_name, argv[1]);
		strcat(inputData_name, "biginput.data");
		//strcat(goldenData_name, "biggolden.data");
	}else {
		strcpy(inputData_name, "biginput.data");
		//strcpy(goldenData_name, "biggolden.data");
	}

#else

#ifdef MEDIUM_DATASET
	//char * inputData_name = "mediuminput.data";
	//char * goldenData_name = "mediumgolden.data";
	if (argc == 2) {
		strcpy(inputData_name, argv[1]);
		//strcpy(goldenData_name, argv[1]);
		strcat(inputData_name, "mediuminput.data");
		//strcat(goldenData_name, "mediumgolden.data");
	}else {
		strcpy(inputData_name, "mediuminput.data");
		//strcpy(goldenData_name, "mediumgolden.data");
	}
#else// Small data set
	//char * inputData_name = "smallinput.data";
	//char * goldenData_name = "smallgolden.data";
	if (argc == 2) {
		strcpy(inputData_name, argv[1]);
		//strcpy(goldenData_name, argv[1]);
		strcat(inputData_name, "smallinput.data");
		//strcat(goldenData_name, "smallgolden.data");
	}else {
		strcpy(inputData_name, "smallinput.data");
		//strcpy(goldenData_name, "smallgolden.data");
	}
#endif // End of MEDIUM_DATASET

#endif // End of BIG_DATASET

	//printf("InputData Name = %s\n", inputData_name);
	//printf("GoldenData Name = %s\n", goldenData_name);

	FILE *fp;
	int size_file;
	int rt_fread = 0;
	fp = fopen(inputData_name, "rb");
	if (fp == NULL) {
		printf("InputData File Open Error!\n");
		return FAILURE;
	}else {
		// Set the current position of fp to SEEK_END,
		// that is SEEK_CUR = SEEK_END
		fseek(fp, 0L, SEEK_END);
		size_file = ftell(fp);
		// Set SEEK_CUR back to the beginning of fp (SEEK_SET)
		fseek(fp, 0L, SEEK_SET);
		inputData = (char *) malloc(size_file);
		rt_fread = fread(inputData, size_file, 1, fp);
		if (rt_fread != 1) {
			printf("InputData fread return Error!\n");
			return FAILURE;
		}
		//printf("InputData File read Success!\n");
		fclose(fp);
	}

	/// Run kernel
	abstract_kernel(inputData);

	/// Check output data from reference data set
	/*
	FILE *fp_golden;
	size_file = 0;
	rt_fread = 0;
	fp_golden = fopen(goldenData_name, "rb");
	if (fp_golden == NULL) {
		printf("GoldenData File Open Error!\n");
		return FAILURE;
	}else {
		// Set the current position of fp to SEEK_END,
		// that is SEEK_CUR = SEEK_END
		fseek(fp_golden, 0L, SEEK_END);
		size_file = ftell(fp_golden);
		// Set SEEK_CUR back to the beginning of fp (SEEK_SET)
		fseek(fp_golden, 0L, SEEK_SET);
		goldenData = (char *) malloc(size_file);
		rt_fread = fread(goldenData, size_file, 1, fp_golden);
		if (rt_fread != 1) {
			printf("GoldenData fread return Error!\n");
			return FAILURE;
		}
		printf("GoldenData File read Success!\n");
		fclose(fp_golden);
	}
	*/

	//assert(!memcmp(inputData, goldenData, size_file) && "Output result is wrong!\n");
	/*
	struct args_list_t *inputD = (struct args_list_t*) inputData;
	struct args_list_t *golden = (struct args_list_t*) goldenData;
	int i = 0;
	for (i=0; i<nx*ny*nz; i++) {
		if (inputD->Anext[i] != golden->Anext[i]) {
			printf("inputD->Anext[%d] = %f\n", i, inputD->Anext[i]);
			printf("golden->Anext[%d] = %f\n", i, golden->Anext[i]);
			printf("Output result Error!\n");
			return FAILURE;
		}
	}
	*/
	int rt_val = compare_result(inputData_name, inputData);
	if (rt_val == 0) {
		free(inputData);
		//free(goldenData);
		printf("\n\n\tSuccess!\n\n");
		return SUCCESS;
	}else {
		printf("\n\n\tFailure!\n\n");
		return FAILURE;
	}
}

