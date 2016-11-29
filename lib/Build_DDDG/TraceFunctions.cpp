#include "profile_h/TraceFunctions.h"

//FILE *full_trace_file;
gzFile full_trace_file;

int initp = 0;

int inst_count = 0;

void trace_logger_init() {
	std::string trace_file_name = inputPath + "dynamic_trace.gz";
	//std::string trace_file_name = "dynamic_trace.gz";
	//full_trace_file = fopen("dynamic_trace", "w");
	//full_trace_file = fopen(trace_file_name.c_str(), "w");
	full_trace_file = gzopen(trace_file_name.c_str(), "w");

	//if (full_trace_file == NULL) {
	if (full_trace_file == Z_NULL) {
		//perror("Failed to open logfile \"dynamic_trace\"");
		//exit(-1);
		std::string err_str = "Error, can not open file " + trace_file_name;
		assert(false && err_str.c_str() );
	}
	//atexit(trace_logger_fin);
}

void trace_logger_fin() {
	//fclose(full_trace_file);
	std::cout << "DEBUG-INFO: [profiling_trace-generation] Dynamic trace generated successfully\n";
	gzclose(full_trace_file);
}

void trace_logger_log0(int line_number, char *name, char *bbid, char *instid, int opcode) {
	if (!initp) {
		trace_logger_init();
		initp = 1;
	}
	//fprintf(full_trace_file, "\n0,%d,%s,%s,%s,%d,%d\n", line_number, name, bbid, instid, opcode, inst_count);
	gzprintf(full_trace_file, "\n0,%d,%s,%s,%s,%d,%d\n", line_number, name, bbid, instid, opcode, inst_count);
	inst_count++;
}

void trace_logger_log_int(int line, int size, int64_t value, int is_reg, char *label) {
	assert(initp == 1 && "initp is not equal to 1");
	if (line == RESULT_LINE){
		//fprintf(full_trace_file, "r,%d,%ld,%d,%s\n", size, value, is_reg, label);
		gzprintf(full_trace_file, "r,%d,%ld,%d,%s\n", size, value, is_reg, label);
	}
	else if (line == FORWARD_LINE){
		//fprintf(full_trace_file, "f,%d,%ld,%d,%s\n", size, value, is_reg, label);
		gzprintf(full_trace_file, "f,%d,%ld,%d,%s\n", size, value, is_reg, label);
	}
	else{
		//fprintf(full_trace_file, "%d,%d,%ld,%d,%s\n", line, size, value, is_reg, label);
		gzprintf(full_trace_file, "%d,%d,%ld,%d,%s\n", line, size, value, is_reg, label);
	}
}

void trace_logger_log_double(int line, int size, double value, int is_reg, char *label) {
	assert(initp == 1 && "initp is not equal to 1");
	if (line == RESULT_LINE){
		//fprintf(full_trace_file, "r,%d,%f,%d,%s\n", size, value, is_reg, label);
		gzprintf(full_trace_file, "r,%d,%f,%d,%s\n", size, value, is_reg, label);
	}
	else if (line == FORWARD_LINE){
		//fprintf(full_trace_file, "f,%d,%f,%d,%s\n", size, value, is_reg, label);
		gzprintf(full_trace_file, "f,%d,%f,%d,%s\n", size, value, is_reg, label);
	}
	else{
		//fprintf(full_trace_file, "%d,%d,%f,%d,%s\n", line, size, value, is_reg, label);
		gzprintf(full_trace_file, "%d,%d,%f,%d,%s\n", line, size, value, is_reg, label);
	}
}

void trace_logger_log_int_noreg(int line, int size, int64_t value, int is_reg) {
	assert(initp == 1 && "initp is not equal to 1");
	if (line == RESULT_LINE){
		//fprintf(full_trace_file, "r,%d,%ld,%d\n", size, value, is_reg);
		gzprintf(full_trace_file, "r,%d,%ld,%d\n", size, value, is_reg);
	}
	else if (line == FORWARD_LINE) {
		//fprintf(full_trace_file, "f,%d,%ld,%d\n", size, value, is_reg);
		gzprintf(full_trace_file, "f,%d,%ld,%d\n", size, value, is_reg);
	}
	else{
		//fprintf(full_trace_file, "%d,%d,%ld,%d\n", line, size, value, is_reg);
		gzprintf(full_trace_file, "%d,%d,%ld,%d\n", line, size, value, is_reg);
	}
}

void trace_logger_log_double_noreg(int line, int size, double value, int is_reg) {
	assert(initp == 1 && "initp is not equal to 1");
	if (line == RESULT_LINE) {
		//fprintf(full_trace_file, "r,%d,%f,%d\n", size, value, is_reg);
		gzprintf(full_trace_file, "r,%d,%f,%d\n", size, value, is_reg);
	}
	else if (line == FORWARD_LINE) {
		//fprintf(full_trace_file, "f,%d,%f,%d\n", size, value, is_reg);
		gzprintf(full_trace_file, "f,%d,%f,%d\n", size, value, is_reg);
	}
	else {
		//fprintf(full_trace_file, "%d,%d,%f,%d\n", line, size, value, is_reg);
		gzprintf(full_trace_file, "%d,%d,%f,%d\n", line, size, value, is_reg);
	}
}