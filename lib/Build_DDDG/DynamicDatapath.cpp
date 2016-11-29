#include "profile_h/DynamicDatapath.h"

DynamicDatapath::DynamicDatapath(std::string kernel_name, std::string trace_file_name, std::string input_path, std::string lp_name, unsigned lp_level, unsigned lp_unroll_factor, bool enable_pipelining, unsigned IL_asap) : BaseDatapath(kernel_name, trace_file_name, "", input_path, lp_name, lp_level, lp_unroll_factor, IL_asap) {
	std::cout << "===========================\n";
	std::cout << "DEBUG-INFO: [trace-analysis_dynamic-trace] Analyzing DDDG for loop " << lp_name << "\n";
	std::cout << "===========================\n";

	/// Initialization
	// initialize_graph();
	
	// FPGA flow
	//removeInductionDependence();
	//removePhiNodes();
	initBaseAddress();

	/// completePartition() and scratchpadPartition() should be called when we enable optimization, otherwise it will
	/// make calculateFPGAResRequired(...) invoked by asap_scheduling(...) or alap_scheduling(...) failed because of
	/// modified array name
	//completePartition();
	//scratchpadPartition();
	//loopUnrolling();

	//loopFlatten();
	//calculateInstructionDistribution();
	//calculateTimestampDDDG();
	//parallelismProfileDDDG();

	/// Schedule and Calculate execution cycles
	//uint64_t exe_cycles = run_fpga_simulation();
	if (enable_pipelining == false) {

		if (show_dddg_bf_opt == true) {
			/// Output the Graph
			dumpGraph();
		}

		uint64_t exe_cycles = fpga_estimation();
		std::cout << "DEBUG-INFO: [trace-analysis_dynamic-trace] FPGA execution time = " << exe_cycles << "(cycles)" << std::endl;
	}
	else {
		VERBOSE_PRINT(std::cout << "DEBUG-INFO: [trace-analysis_dynamic-trace] Analyzing recurrence initiation interval\n");
		IL_asap_ii = fpga_estimation_one_more_subtrace_for_recII_calculation();
	}
	
	//uint64_t exe_cycles = run_parallel_simulation();
	//cout << "Parallel Execution cycles = " << num_cycles << endl;

	// GPU flow
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [trace-analysis_dynamic-trace] Finished\n");
	VERBOSE_PRINT(std::cout << "-------------------\n");
}

DynamicDatapath::~DynamicDatapath() {}

unsigned DynamicDatapath::getIL_asap_ii() const {
	return IL_asap_ii;
}