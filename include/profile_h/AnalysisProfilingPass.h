#ifndef ANALYSISPROFILINGPASS_H
#define ANALYSISPROFILINGPASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/LoopInfo.h"

#include "profile_h/lin-profile.h"
#include "profile_h/QueryBasicBlockIDPass.h"
#include "profile_h/QueryLoadStoreIDPass.h"
#include "profile_h/ExtractLoopInfoPass.h"
#include "profile_h/TraceDataStructure.h"
#include "profile_h/CodeInstrumentPass.h"
#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"

namespace llvm{

	class AnalysisProfiling : public ModulePass {
	
	public:
		static char ID;
		AnalysisProfiling();
		void getAnalysisUsage(AnalysisUsage &AU) const;
		bool runOnModule(Module &M);
		uint64_t getBBFreq(BasicBlock* BB) const;
		uint64_t getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const;
		
		bool testbbIDpass(Module &M);

	private:
		QueryBasicBlockID* bbIDpass;
		QueryLoadStoreID* lsIDpass;

		/// Used in VerifyProfile()
		BBFreqMap BBFreq;
		BranchFreqMap BranchFreq;

		//typedef std::map<std::string, std::vector<std::pair<unsigned, unsigned> > > loop2subloopMapTy;
		/// loopName --> subloop <depth, bound> map
		/// Move to auxiliary.h
		//typedef std::map<std::string, std::map<unsigned, unsigned> > loop2subloopMapTy;
		loop2subloopMapTy loop2subloopMap;
		std::map<unsigned, unsigned> subloopDepthandLoopBound;
		
		///Move to auxiliary.h
		//typedef std::map<std::string, std::vector<std::string> > funcName2loopsMapTy;
		//funcName2loopsMapTy funcName2loopsMap;

		/// loopName2BoundList: the whole loop name -> loop bounds list
		/// Loop bounds sequence: 
		/// The top-level loop -> lower-level loop -> ... -> the innermost-level loop
		/// Move to auxiliary.h
		//typedef std::map<std::string, std::list<unsigned> > loopName2BoundListTy;
		//typedef std::map<std::string, loopName2BoundListTy> funcName2loopBoundsTy;
		//funcName2loopBoundsTy funcName2loopBounds;

		/// Move to auxiliary.h
		//typedef std::map<std::string, loop2subloopMapTy> func2loopInfoMapTy;
		//func2loopInfoMapTy func2loopInfoMap;

		/// Used to recognize the relationship between basic block ids and loop name
		/// Move to auxiliary.h
		//typedef std::map<uint64_t, std::string> BBids2loopNameMapTy;
		//BBids2loopNameMapTy BBids2loopNameMap;

		/// Move to auxiliary.h
		//typedef std::map<std::string, BBids2loopNameMapTy> funcName2loopInfoMapTy;
		/// Move to GetLoopBound.h
		//funcName2loopInfoMapTy funcName2loopInfoMap;

		/// vector used to store all loops' basic block IDs
		//std::vector<uint64_t> bbID_for_loops_vec;

		std::map<unsigned, unsigned> loopDepth2BoundMap;

		void VerifyProfile();
		void extract_loopInfo_inFunc(Module &M);

		void calculateLoopBound(Function* func);
		void modifyTrace_EvalIterationIndex();
		void trace_eval_iteration_index(Function* func);

		void analyze_Dynamic_DataDep();
		void analyze_Arithmetic_Intensity();

		// Write information of dynamic data dependence and arithmetic intensity of loops into a file
		void write_dynDataDep_ArithIntensity_file();

		// vector of (load/store address, trace_entry ID) pair
		//typedef std::vector<std::pair<uintptr_t, uint64_t> > lsAddrTraceIDpairVecTy;
		//typedef std::list<std::pair<uintptr_t, RecordTrace> > lsAddrTraceIDpairListTy;
		typedef std::map<uintptr_t, RecordTrace> lsAddrTraceIDpairMapTy;
		typedef lsAddrTraceIDpairMapTy::iterator iterator_t;
		iterator_t findEntryInlsAddrTraceIDpairVec(lsAddrTraceIDpairMapTy &lsAddrTraceIDpairVec, uintptr_t ls_addr);

		std::map<uint64_t, unsigned> loopHeaderBBid2loopDepthMap;

		//typedef std::multimap<std::string, RecordTrace> loopName2Trace_fMMapTy;

		//typedef std::multimap<std::string, loopName2Trace_fMMapTy > funcN2loopNToTrace_fMMapTy;
		//funcN2loopNToTrace_fMMapTy funcN2loopNToTrace_fMMap;
		typedef std::multimap<std::string, RecordTrace> funcName2Trace_entryMMapTy;
		funcName2Trace_entryMMapTy funcName2Trace_entryMMap;

		// loopName --> subTrace_f_entry Map
		//typedef std::map<std::string, TraceType> loopName2subTraceMapTy;
		typedef std::pair<std::string, std::string> funcLoopNamePairTy;
		// (funcionName, loopName) -->  loopName2subTraceMapTy Map
		typedef std::multimap<funcLoopNamePairTy, TraceType> funcName2loopTraceMMapTy;
		funcName2loopTraceMMapTy funcName2loopTraceMMap;

		// (functionName, loopName) -> number of binary operators executed
		typedef std::map<std::pair<std::string, std::string>, uint64_t> fnlpName2BinaryOpNumMapTy;
		fnlpName2BinaryOpNumMapTy fnlpName2BinaryOpNumMap;

		// (functionName, loopName) -> number of memory bytes transferred
		typedef std::map<std::pair<std::string, std::string>, uint64_t> fnlpName2MemByteMapTy;
		fnlpName2MemByteMapTy fnlpName2MemByteMap;

		// (functionName, loopName) -> arithmetic intensity
		typedef std::map<std::pair<std::string, std::string>, float> fnlpName2ArithIntensityMapTy;
		fnlpName2ArithIntensityMapTy fnlpName2ArithIntensityMap;

		// Loop Basic Block ID -> the corresponding loop depth
		typedef std::map<uint64_t, unsigned> loopBBid2loopDepthMapTy;
		loopBBid2loopDepthMapTy loopBBid2loopDepthMap;
		void getloopBBid2loopDepthMap(Function* func);
	};

	/// Write the trace into a file
	void write_trace_file(TraceType& trace, std::string fname);

} // End namespace llvm

#endif // End of ANALYSISPROFILINGPASS_H