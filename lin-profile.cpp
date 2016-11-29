#include "llvm/IR/LLVMContext.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/CommandLine.h"
//#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include "llvm/PassManager.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/IRReader/IRReader.h"
#include "llvm/Bitcode/BitcodeWriterPass.h"
#include "llvm/Support/ToolOutputFile.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/Debug.h"
#include "llvm/ADT/Triple.h"

#include "profile_h/lin-profile.h"

#define DEBUG_TYPE "lin-analyzer"

using namespace llvm;

std::string inputFileName;
std::string inputPath;
std::string outputPath;
std::string config_filename;
std::vector<std::string> kernel_names;

bool enable_profiling_time_only;
bool memory_trace_gen;
bool enable_no_trace;
bool show_cfg_detailed;
bool show_cfg_only;
bool show_dddg_bf_opt;
bool show_dddg_af_opt;
bool verbose_print;
bool enable_store_buffer;
bool enable_shared_load_removal;
bool disable_shared_load_removal;
bool enable_repeated_store_removal;
bool enable_tree_height_reduction_float;
bool enable_tree_height_reduction_integer;
bool enable_memory_disambiguation;
bool increase_load_latency;
bool disable_fp_unit_threshold;
bool enable_extra_scalar;
bool enable_rw_rw_memory;
bool target_vc707;
std::vector<std::string> target_loops;

/// Global map
BBidFreqMap BBidFreq;
BranchidFreqMap BranchIDFreq;

const std::string help_message = "----------------------------------------------------------------------\n"
																 "    Lin-analyzer: A High Level Analysis Tool for FPGA Accelerators\n"
																 "----------------------------------------------------------------------\n"
																 "Usage:	lin-analyzer [file.bc] -Ipath=[path] -config=[filename] [kernel-name] -Opath=[path] -TargetLoops=[index] [options]\n"
																 "Options:\n"
																 "	-h, --help               Help information.\n"
																 "	-profiling-time          Profile kernels without FPGA estimation.\n"
																 "	-mem-trace               Obtain memory trace for access pattern analysis without FPGA estimation. This\n"
																 "	                         option should be used with -profiling-time together.\n"
																 "	-no-trace                Disable dynamic trace generation.\n"
																 "	-TargetLoops             Specify target loops focused. Eg., -TargetLoops=2,3: only analyse loop 2 and 3.\n"
																 "	-cfg-detailed            Show CFG with detailed instructions.\n"
																 "	-cfg-only                Show CFG only with basic blocks.\n"
																 "	-dddg-bf-opt             Show DDDG before optimization. May slow down program if input size is large\n"
																 "	-dddg-af-opt             Show DDDG after optimization. May slow down program if input size is large\n"
																 "	-verbose                 Verbose mode, print more information.\n"
																 "	-dis-store-buffer        Disable store-buffer optimization.\n"
																 "	-shared-load-removal     Enable shared-load-removal optimization.\n"
																 "	-dis-shared-load-removal Disable shared-load-removal opt., even for completely unrolling config.\n"
																 "	-dis-rp-store-removal    Disable repeated-store-removal optimization.\n"
																 "	-THR-float               Enable tree height reduction optimization for floating point operations.\n"
																 "	-THR-integer             Enable tree height reduction optimization for integer operations.\n"
																 "	-memory-disambig         Enable memory disambiguation optimization.\n"
																 "	-dis-fp-threshold        Disable floating point unit threshold and area budget is unlimited.\n"
																 "	-en-extra-scalar         Sometimes, result might be shifted, this option is used to improve prediction.\n"
																 "	-en-rw-rw-memory         Use memory with two ports and each supports read and write operation. Default is\n"
																 "	                         read-only and read-write.\n"
																 "	-vc707                   Target for Xilinx Virtex7 VC707. Default is Xilinx Zedboard or ZC702\n"
																 "\n"
																 "Report bugs to guanwen<guanwen@comp.nus.edu.sg>\n"
																 "----------------------------------------------------------------------\n";

int main(int argc, char **argv) {

	parse_input_arguments(argc, argv);

	errs() << "********************************************************\n";
	errs() << "********************************************************\n\n";
	errs() << "     Lin-analyzer: An FPGA High-Level Analysis Tool\n\n";
	errs() << "********************************************************\n";
	errs() << "********************************************************\n";

	LLVMContext &Context = getGlobalContext();
	SMDiagnostic Err;

	std::unique_ptr<Module> M;

	// Note: The path should include \5C in Unix or / in Windows in the end!

	// Store the location of InputFilename
	inputFileName = InputFilename.c_str();

	M.reset(ParseIRFile(InputFilename, Err, Context));
	//M.reset(ParseIRFile("E:/llvm/home/henry/llvm_3.5/input/test1.bc", Err, Context));

	if (!M.get()) {
		Err.print(argv[0], errs());
		return 1;
	}

	Module &mod = *M.get();
	Triple TheTriple(mod.getTargetTriple());

	if (OutputFilename.empty()) {
		if (InputFilename == "-") {
			OutputFilename = "-";
		}
		else {
			const std::string &IFN = InputFilename;
			int Len = IFN.length();
			// If the source ends in .bc, strip it off
			if (IFN[Len - 3] == '.' && IFN[Len - 2] == 'b' && IFN[Len - 1] == 'c') {
				OutputFilename = std::string(IFN.begin(), IFN.end() - 3) + "_trace.bc";
			}
		}
	}

	// Figure out what stream we are supposed to write to...
	std::unique_ptr<tool_output_file> Out;
	std::string ErrorInfo;
	Out.reset(new tool_output_file(OutputFilename.c_str(), ErrorInfo, sys::fs::F_None));
	if (!ErrorInfo.empty()) {
		errs() << ErrorInfo << "\n";
		return 1;
	}

	PassManager Passes;
	Passes.add(createLoopNumberPass());
	Passes.add(createAssignBasicBlockIDPass());
	Passes.add(createAssignLoadStoreIDPass());
	Passes.add(createExtractLoopInfoPass());
#ifndef BUILD_DDDG_H
	Passes.add(createQueryBasicBlockIDPass());
	// The query functions of Load/Store ID are implemented in AssignLoadStoreIDPass.
	Passes.add(createQueryLoadStoreIDPass()); 
#if defined(TWO_TIME_PROFILING) || defined(ENABLE_INSTDISTRIBUTION)
	Passes.add(createGetLoopBoundPass());
#endif // End of TWO_TIME_PROFILING or ENABLE_INSTDISTRIBUTION

#ifdef ENABLE_INSTDISTRIBUTION
	Passes.add(createGetInstDistributionPass());
#else // define ENABLE_INSTDISTRIBUTION
	Passes.add(createCodeInstrumentPass());
	Passes.add(createAnalysisProfilingPass());
#endif // End of ENABLE_INSTDISTRIBUTION

#else // define BUILD_DDDG_H
	Passes.add(createInstrumentForDDDGPass());
#endif // End of ifndef BUILD_DDDG_H
	Passes.add(createBitcodeWriterPass(Out->os()));

	// Before executing passes, print the final values of the LLVM options.
	//cl::PrintOptionValues();

	Passes.run(mod);

	// Declare success.
	Out->keep();

	return 0;
}

void parse_input_arguments(int argc, char **argv) {
	//assert(argc > 3 && "Note: Please provide path of input data, bitcode and names of kernel functions\n");
	if (argc < 2) {
		// Print message
		errs() << "Missing input arguments (\"lin-analyzer --help\" or \"lin-analyzer -h\" for help)\n\n";
		exit(-1);
	} else if (argc < 3) {
		std::string first_argv = argv[1];
		if (!first_argv.compare("-h") || !first_argv.compare("--help")) {
			//Print Help message
			errs() << help_message;
			exit(-1);
		}
		else {
			errs() << "Missing input arguments (\"lin-analyzer --help\" or \"lin-analyzer -h\" for help)\n\n";
			exit(-1);
		}
	}
	InputFilename = argv[1];
	size_t pos = InputFilename.find(".bc");
	if (pos == std::string::npos) {
		errs() << "Wrong input argument, input bitcode file should be *.bc (\"lin-analyzer --help\" or \"lin-analyzer -h\" for help)\n\n";
		exit(-1);
	}
	errs() << "Input bitcode file: " << InputFilename << "\n";
	std::string input_path = argv[2];
	size_t pos_path = input_path.find("-Ipath=");
	if (pos_path == std::string::npos) {
		errs() << "Wrong input argument, path should be valid (\"lin-analyzer --help\" or \"lin-analyzer -h\" for help)\n\n";
		exit(-1);
	}
	pos_path = input_path.find("=");
	inputPath = input_path.substr(pos_path + 1);
#ifdef _MSC_VER
	inputPath += "\\";
#else
	inputPath += "/";
#endif // End of _MSC_VER
	errs() << "Input path: " << inputPath << "\n";

	outputPath = inputPath;

	std::string config_name = argv[3];
	size_t pos_conf = config_name.find("-config=");
	if (pos_conf == std::string::npos) {
		errs() << "Wrong input argument, missing configuration file name (\"lin-analyzer --help\" or \"lin-analyzer -h\" for help)\n\n";
		exit(-1);
	}
	pos_conf = config_name.find("=");
	config_filename = config_name.substr(pos_conf + 1);

	kernel_names.clear();
	kernel_names.push_back(argv[4]);

	enable_profiling_time_only = false;
	memory_trace_gen = false;
	enable_no_trace = false;
	show_cfg_detailed = false;
	show_cfg_only = false;
	show_dddg_bf_opt = false;
	show_dddg_af_opt = false;
	verbose_print = false;
	enable_store_buffer = true;
	enable_shared_load_removal = false;
	disable_shared_load_removal = false;
	enable_repeated_store_removal = true;
	enable_tree_height_reduction_float = false;
	enable_tree_height_reduction_integer = false;
	enable_memory_disambiguation = false;
	disable_fp_unit_threshold = false;
	enable_extra_scalar = false;
	enable_rw_rw_memory = false;
	target_vc707 = false;
	target_loops.clear();

	for (int i = 5; i < argc; i++) {
		std::string tmp_str = argv[i];
		if (!tmp_str.compare("-profiling-time")) {
			enable_profiling_time_only = true;
		}
		else if (!tmp_str.compare("-mem-trace")) {
			memory_trace_gen = true;
		}
		else if (!tmp_str.compare("-no-trace")) {
			enable_no_trace = true;
		}
		else if (tmp_str.find("-TargetLoops")!=std::string::npos) {
			// Split strings
			size_t pos_name = tmp_str.find("=");
			std::string loop_indexes = tmp_str.substr(pos_name+1);
			if (loop_indexes == "") {
				errs() << "Please specify indexes of target loops, starting from 0. Initialization loops will be counted.\n\n";
				exit(-1);
			}
			else {
				size_t has_comma = loop_indexes.find(",");
				if (has_comma == std::string::npos) {
					target_loops.push_back(loop_indexes);
				}
				else {
					while (loop_indexes.find(",") != std::string::npos) {
						size_t pos_comma = loop_indexes.find(",");
						std::string lp_index = loop_indexes.substr(0, pos_comma);
						loop_indexes.erase(0, pos_comma+1);
						if (lp_index != "") {
							target_loops.push_back(lp_index);
						}
					}
					// Get the last element
					if (loop_indexes != "") {
						target_loops.push_back(loop_indexes);
					}
				}
			}
		}
		else if (!tmp_str.compare("-cfg-detailed")) {
			show_cfg_detailed = true;
		}
		else if (!tmp_str.compare("-cfg-only")) {
			show_cfg_only = true;
		}
		else if (!tmp_str.compare("-dddg-bf-opt")) {
			show_dddg_bf_opt = true;
		}
		else if (!tmp_str.compare("-dddg-af-opt")) {
			show_dddg_af_opt = true;
		}
		else if (!tmp_str.compare("-verbose")) {
			verbose_print = true;
		}
		else if (!tmp_str.compare("-dis-store-buffer")) {
			enable_store_buffer = false;
		}
		else if (!tmp_str.compare("-shared-load-removal")) {
			enable_shared_load_removal = true;
		}
		else if (!tmp_str.compare("-dis-shared-load-removal")) {
			disable_shared_load_removal = true;
		}
		else if (!tmp_str.compare("-dis-rp-store-removal")) {
			enable_repeated_store_removal = false;
		}
		else if (!tmp_str.compare("-THR-float")) {
			enable_tree_height_reduction_float = true;
		}
		else if (!tmp_str.compare("-THR-integer")) {
			enable_tree_height_reduction_integer = true;
		}
		else if (!tmp_str.compare("-memory-disambiguation")) {
			enable_memory_disambiguation = true;
		}
		else if (!tmp_str.compare("-dis-fp-threshold")) {
			disable_fp_unit_threshold = true;
		}
		else if (!tmp_str.compare("-en-extra-scalar")) {
			enable_extra_scalar = true;
		} 
		else if (!tmp_str.compare("-en-rw-rw-memory")) {
			enable_rw_rw_memory = true;
		}
		else if (!tmp_str.compare("-vc707")) {
			target_vc707 = true;
		}
		else if (!tmp_str.compare("-h") || !tmp_str.compare("--help")) {
			//Print Help message
			errs() << help_message;
			exit(-1);
		}
		else if (tmp_str.find("-Opath=") != std::string::npos) {
			size_t pos_oPath = tmp_str.find("=");
			outputPath = tmp_str.substr(pos_path + 1);
#ifdef _MSC_VER
			outputPath += "\\";
#else
			outputPath += "/";
#endif // End of _MSC_VER
		}
		else {
			errs() << "Wrong input argument (\"lin-analyzer --help\" or \"lin-analyzer -h\" for help)\n\n";
			exit(-1);
		}

	}

	// Loading kernel names into kernel_names vector
	DEBUG(dbgs() << "We only focus on the kernel for this application: \n");
	DEBUG(dbgs() << "\tKernel name: " << kernel_names[0] << "\n");

	DEBUG(dbgs() << "Please make sure all functions within a kernel function are included.");
	DEBUG(dbgs() << "We also need to consider these functions. Otherwise, the tool will ");
	DEBUG(dbgs() << "ignore all these functions and cause problems!\n");

}