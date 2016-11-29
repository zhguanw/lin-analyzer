#ifndef GETLOOPBOUNDPASS_H
#define GETLOOPBOUNDPASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/InstVisitor.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/IR/Verifier.h"
#include "llvm/ADT/SmallString.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/ExecutionEngine/GenericValue.h"
#include "llvm/Transforms/Utils/Cloning.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/ExecutionEngine/JIT.h"
#include "llvm/Support/ManagedStatic.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Transforms/Utils/ValueMapper.h"
#include "llvm/Bitcode/ReaderWriter.h"
#include "llvm/Support/ToolOutputFile.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/CommandLine.h"
//#include "llvm/IR/BasicBlock.h"

#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"
#include "profile_h/ExtractLoopInfoPass.h"
#include "profile_h/QueryBasicBlockIDPass.h"
#include "profile_h/QueryLoadStoreIDPass.h"
#include "profile_h/lin-profile.h"
//#include "profile_h/CodeInstrumentPass.h"

#include "llvm/Transforms/Scalar.h"
#include "llvm/Analysis/LoopInfo.h"

#include <list>

static llvm::cl::opt<std::string>
OutputModifiedLBBCFilename("o", llvm::cl::desc("<output modified bitcode file for getting loop bounds>"), llvm::cl::value_desc("filename"));

namespace llvm {

	class GetLoopBound : public ModulePass {
		
	public:
		static char ID;
		ValueToValueMapTy VMap;

		GetLoopBound();
		void getAnalysisUsage(AnalysisUsage &AU) const;
		bool doInitialization(Module &M);
		bool runOnModule(Module &M);

		void getLoopBBID2HeaderAndloopDepthPairMap(Module &M);
		bool isLoopHeader(uint64_t bbid);
		void getHeaderID_list_vec(Module &M);
		void extract_loopInfo_inFunc(Module &M);
		void calculateLoopBound(Function* func);

	private:
		const DataLayout* DL;
		QueryBasicBlockID* bbIDpass;
		QueryLoadStoreID* lsIDpass;

		Type* Int64Type;
		Type* VoidType;
		Type* VoidPtrType;

		Type* VoidPtrBBType;

		Function* BBFreqIncFn;
		Function* RecordLoadFn;
		Function* RecordStoreFn;
		BasicBlock* LastBB;

		Module* NewM;

		//typedef std::map<std::string, std::map<unsigned, unsigned> > loop2subloopMapTy;
		loop2subloopMapTy loop2subloopMap;
		std::map<unsigned, unsigned> subloopDepthandLoopBound;

		/// Loop Basic Block ID --> (the corresponding Loop Header BB id, loop depth)
		//typedef std::map<uint64_t, std::pair<uint64_t, unsigned> > loopBBid2HeaderAndloopDepthPairMapTy;
		//loopBBid2HeaderAndloopDepthPairMapTy loopBBid2HeaderAndloopDepthPairMap;

		/// headerID_list_vec contains a vector of headerID_list
		/// headerID_list:  The header BB id of the innermost loop --> The header BB id of the outermost loop
		//typedef std::vector<std::list<uint64_t> > headerID_list_vecTy;
		//headerID_list_vecTy headerID_list_vec;

		//typedef std::map<std::string, headerID_list_vecTy> funcName2headIDlistVecMapTy;
		//funcName2headIDlistVecMapTy funcName2headIDlistVecMap;

		std::vector<std::string> functionNameVec;

		void instrumentFunction(Function* F);

		void instrumentBasicBlock(BasicBlock* BB);

		void insertCounterBefore(Instruction* Inst, bool Increase, bool edgeIncrease);

		void write_Modified_bitcode_to_File();
	};

	//Embedded Profiler Engine
	class EmbeddedProfilerEngine_lb {

	public:

		//BBidFreqMap BBFreq;
		//BranchidFreqMap BranchFreq;

		EmbeddedProfilerEngine_lb(Module &M, Function* bbFreqfn);
		/*
		explicit EmbeddedProfilerEngine(Module &M, Function* loadFn, Function* storeFn) :
		Mod(M), RecordLoadFn(loadFn), RecordStoreFn(storeFn) {}*/

		void runOnProfiler();

		//void increaseEdgeCounter(BasicBlock* src, BasicBlock* dest);
		void increaseEdgeCounter(uint64_t src_id, uint64_t dest_id, uint64_t increase_edge);

		uint64_t getBBidFreq(uint64_t BBid) const;

		uint64_t getBranchFreq(uint64_t Src_id, uint64_t Dest_id) const;

	private:

		Module& Mod;
		Function* BBFreqFn;

		//ValueToValueMapTy V2VMap;
		//typedef std::map<uint64_t, uint64_t> BBFreqMap;
		//typedef std::map<llvm::BasicBlock*, BBFreqMap> BranchFreqMap;

		/// Used in VerifyProfile()
		//BBFreqMap BBFreq;
		//BranchFreqMap BranchFreq;

		void VerifyProfile();
	};

	struct ProfilerJITContext_lb {
		EmbeddedProfilerEngine_lb *P;
		ProfilerJITContext_lb();
		uint64_t getLastBB() const {
			return LastBBid;
		}

		void setLastBB(uint64_t lastbbid) {
			LastBBid = lastbbid;
		}
	private:
		uint64_t LastBBid;
	};

	struct ProfilerJITSingletonContext_lb{
		ProfilerJITSingletonContext_lb(EmbeddedProfilerEngine_lb *P);
		~ProfilerJITSingletonContext_lb();
	};

	static ManagedStatic<ProfilerJITContext_lb> GlobalContext_lb;

} // End namespace LLVM

namespace {
	std::vector<uint64_t> EntryBBid_vec;
}

#endif // End of GETLOOPBOUNDPASS_H