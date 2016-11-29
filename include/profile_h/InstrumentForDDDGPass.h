#ifndef INSTRUMENT_FOR_DDDG_PASS_H
#define INSTRUMENT_FOR_DDDG_PASS_H

#include <vector>
#include <cmath>
#include "llvm/Pass.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Metadata.h"
#include "llvm/IR/DebugInfo.h"
#include "llvm/IR/Verifier.h"
#include "llvm/ExecutionEngine/GenericValue.h"
#include "llvm/ExecutionEngine/JIT.h"
#include "llvm/Support/raw_os_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/ManagedStatic.h"
#include "llvm/Support/TargetSelect.h"

#include "llvm/Transforms/Utils/Cloning.h"
#include "llvm/Transforms/Utils/ValueMapper.h"

#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"

#include "profile_h/SlotTracker.h"
#include "profile_h/TraceFunctions.h"
#include "profile_h/generic_func.h"
#include "profile_h/BaseDatapath.h"
#include "profile_h/DDDG.h"
#include "profile_h/DynamicDatapath.h"

#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>

#define NUM_OF_INTRINSICS 35
#define NUM_OF_LLVM_INTRINSICS 33
#define RESULT_LINE 19134
#define FORWARD_LINE 24601
#define DMA_STORE 98
#define DMA_LOAD 99

using namespace std;

namespace llvm {

	class InstrumentForDDDG : public ModulePass {

	public:
		static char ID;
		InstrumentForDDDG();
		void getAnalysisUsage(AnalysisUsage &AU) const;
		bool doInitialization(Module &M);
		bool runOnModule(Module &M);

		char ** str_split(char *a_str, const char a_delim, unsigned *size);
		int trace_or_not(char* func);
		bool is_tracking_function(string func);
		int getMemSize(Type *T);

		/// Function used to instrument LLVM-IR
		void print_line(BasicBlock::iterator itr, int line, int line_number, char *func_or_reg_id,
										char *bbID, char *instID, int opty,	int datasize = 0, Value *value = NULL, 
										bool is_reg = 0);

		void insertInstid(std::string inst_id, unsigned op_code);

		void insertInstid2bbName(std::string inst_id, std::string bbName);

		bool getInstId(Instruction *itr, char* bbid, char* instid, int &instc);

		void getBBId(Value *BB, char *bbid);

		bool performOnBasicBlock(BasicBlock &BB);

		void remove_config(std::string kernel_name, std::string input_path);
		void parse_config(std::string kernel_name, std::string input_path);
		void getUnrollingConfiguration(lpNameLevelPair2headBBnameMapTy& lpNameLvPair2headerBBMap);
		bool readUnrollingConfig(loopName2levelUnrollVecMapTy& lpName2levelUrPairVecMap, std::unordered_map<int, int > &unrolling_config);

		void loopBasedTraceAnalysis();

		void open_summary_file(ofstream& summary_file, std::string kernel_name);
		void close_summary_file(ofstream& summary_file);

	private:
		Function *TL_log0, *TL_log_int, *TL_log_double, *TL_log_int_noreg, *TL_log_double_noreg;
		Module *curr_module;

		SlotTracker *st;
		Function *curr_function;

		void extract_memory_trace_for_access_pattern();

		char **functions;
		unsigned num_of_functions;

		//loopName2levelUnrollVecMapTy loopName2levelUnrollVecMap;
		std::unordered_map<int, int> unrollingConfig;
		std::vector<std::string> pipeline_loop_levelVec;
	};

	//Embedded Profiler Engine
	class ProfilingEngine {

	public:

		ProfilingEngine(Module &M, Function* log0Fn, Function* log_intFn, Function* log_doubleFn,
										Function* log_int_noregFn, Function* log_double_noregFn);
		/*
		explicit EmbeddedProfilerEngine(Module &M, Function* loadFn, Function* storeFn) :
		Mod(M), RecordLoadFn(loadFn), RecordStoreFn(storeFn) {}*/
		void runOnProfiler();

	private:

		Module& Mod;

		Function* log0_Fn;
		Function* log_int_Fn;
		Function* log_double_Fn;
		Function* log_int_noreg_Fn;
		Function* log_double_noreg_Fn;

	};

	struct ProfilingJITContext {
		ProfilingEngine *P;
		ProfilingJITContext();

	private:
		//uint64_t LastBBid;
	};

	struct ProfilingJITSingletonContext{
		ProfilingJITSingletonContext(ProfilingEngine *P);
		~ProfilingJITSingletonContext();
	};

	static ManagedStatic<ProfilingJITContext> GlobalContextDDDG;

} // End of llvm namespace

#endif // End of INSTRUMENT_FOR_DDDG_PASS_H
