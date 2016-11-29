add_llvm_library(LLVMLinProfiler
	InstrumentForDDDGPass.cpp
	#SlotTracker.cpp
	#TraceFunctions.cpp
	#BaseDatapath.cpp
	#DDDG.cpp
	#DynamicDatapath.cpp
	#file_func.cpp
	#opcode_func.cpp
	#generic_func.cpp
	GetInstDistributionPass.cpp
	LoopNumberPass.cpp
	AssignBasicBlockID.cpp
	AssignLoadStoreID.cpp
	QueryBasicBlockID.cpp
	QueryLoadStoreID.cpp
	ExtractLoopInfoPass.cpp
	GetLoopBoundPass.cpp
	CodeInstrumentPass.cpp
	TraceDataStructure.cpp
	AnalysisProfilingPass.cpp
	)

llvm_config(LLVMLinProfiler jit)

#MESSAGE(STATUS "add llvm library LLVMLinProfiler")
