#ifndef CODEINSTRUMENTPASS_H
#define CODEINSTRUMENTPASS_H

#include "llvm/Pass.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/IR/Verifier.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/CallSite.h"
#include "llvm/ExecutionEngine/GenericValue.h"
#include "llvm/Transforms/Utils/Cloning.h"
#include "llvm/IR/DataLayout.h"

// Query id of BB and load/store instructions
#include "profile_h/QueryBasicBlockIDPass.h"
#include "profile_h/QueryLoadStoreIDPass.h"
#include "profile_h/AssignLoadStoreIDPass.h"
//#include "profile_h/LoopNumberPass.h"
#include "profile_h/TraceDataStructure.h"

#include "profile_h/Passes.h"

#include "llvm/ExecutionEngine/JIT.h"
#include "llvm/Support/ManagedStatic.h"
#include "llvm/Support/TargetSelect.h"

#include "llvm/IR/InstVisitor.h"

#include "llvm/Transforms/Utils/ValueMapper.h"
//#include "llvm/Transforms/Utils/Cloning.h"

#include "llvm/Bitcode/ReaderWriter.h"
#include "llvm/Support/ToolOutputFile.h"
#include "llvm/Support/FileSystem.h"

#include "profile_h/lin-profile.h"
#include "profile_h/auxiliary.h"
#include <fstream>
#include <list>

#define VERIFY_PROFILE

static llvm::cl::opt<std::string>
OutputModifiedBCFilename("o", llvm::cl::desc("<output modified bitcode file>"), llvm::cl::value_desc("filename"));


namespace llvm {

	class CodeInstrument : public ModulePass,
						   public InstVisitor<CodeInstrument>{

	public:
		static char ID;
		ValueToValueMapTy VMap;

		CodeInstrument();

		void getAnalysisUsage(AnalysisUsage &AU) const;

		bool doInitialization(Module &M);

		bool runOnModule(Module& M);

		/// runOnFunction is just used to show the CFG graph
		bool runOnFunction(Function& Fn);

		void visitLoadInst(LoadInst& inst);

		void visitStoreInst(StoreInst& inst);

		bool isLoopHeader(uint64_t bbid);

		void getLoopBBID2HeaderAndloopDepthPairMap(Module &M);

		void getHeaderID_list_vec(Module &M);

		//bool runOnBasicBlock(BasicBlock &BB);

		//uint64_t getBBFreq(BasicBlock* BB) const;

		//uint64_t getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const;

		bool ignoreFunc(std::string func_name);
#ifdef VERIFY_PROFILE
		void VerifyProfile();

		/// Functions for accessing Basic Block frequency
		uint64_t getBBFreq(BasicBlock* BB) const;

		uint64_t getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const;
#endif // End of VERIFY_PROFILE	

	private:
		const DataLayout* DL;
		QueryBasicBlockID* bbIDpass;
		QueryLoadStoreID* lsIDpass;

		Type* Int8Type;
		Type* Int32Type;
		Type* Int64Type;
		Type* VoidType;
		Type* VoidPtrType;

		Type* VoidPtrBBType;

		Function* BBFreqIncFn;
		Function* RecordLoadFn;
		Function* RecordStoreFn;
		BasicBlock* LastBB;

		Module* NewM;

		std::vector<std::string> functionNameVec;
		//btree::btree_map<uint64_t, std::pair<BasicBlock*, uint64_t> > lsID2BBandBBidPairMap;
		std::map<uint64_t, std::pair<BasicBlock*, uint64_t> > lsID2BBandBBidPairMap;
		void getlsID2BBandBBidPairMap();

#ifdef VERIFY_PROFILE
		/// Used in VerifyProfile()
		BBFreqMap BBFreq;
		BranchFreqMap BranchFreq;
		/// End of Used in VerifyProfile()
#endif // End of VERIFY_PROFILE	
		//TraceType Trace;

		std::vector<std::string> FuncName_Array;

		void instrumentFunction(Function* F);

		void instrumentBasicBlock(BasicBlock* BB);

		void insertCounterBefore(Instruction* Inst, bool Increase, bool edgeIncrease);

		void instrumentLoadInst(Instruction* inst);

		void instrumentStoreInst(Instruction* inst);

		//void recordBB(uint64_t BBPtr, uint64_t inc);

		//void increaseEdgeCounter(BasicBlock* src, BasicBlock* dest);

		void write_Modified_bitcode_to_File();

		void calculate_exitBB(Module &M);
	};

	//Embedded Profiler Engine
	class EmbeddedProfilerEngine {
	
	public:

		//BBidFreqMap BBFreq;
		//BranchidFreqMap BranchFreq;

		EmbeddedProfilerEngine(Module &M, Function* bbFreqfn, Function* loadFn, Function* storeFn);
		/*
		explicit EmbeddedProfilerEngine(Module &M, Function* loadFn, Function* storeFn) :
			Mod(M), RecordLoadFn(loadFn), RecordStoreFn(storeFn) {}*/

		void runOnProfiler();

		//void increaseEdgeCounter(BasicBlock* src, BasicBlock* dest);
		void increaseEdgeCounter(BasicBlock* src, uint64_t src_id, BasicBlock* dest, uint64_t dest_id, uint64_t increase_edge);

		//uint64_t getBBFreq(BasicBlock* BB) const;

		//uint64_t getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const;

	private:

		Module& Mod;
		Function* BBFreqFn;
		Function* RecordLoadFn;
		Function* RecordStoreFn;
		//ValueToValueMapTy V2VMap;
		//typedef std::map<llvm::BasicBlock*, uint64_t> BBFreqMap;
		//typedef std::map<llvm::BasicBlock*, BBFreqMap> BranchFreqMap;

		//void VerifyProfile() const;
	};

	struct ProfilerJITContext {
		EmbeddedProfilerEngine *P;
		ProfilerJITContext();
		uint64_t getLastBB() const {
			return LastBBid;
		}

		void setLastBB(uint64_t lastbbid) {
			LastBBid = lastbbid;
		}
	private:
		uint64_t LastBBid;
	};

	struct ProfilerJITSingletonContext{
		ProfilerJITSingletonContext(EmbeddedProfilerEngine *P);
		~ProfilerJITSingletonContext();
	};

	/// Utility function: Use to record address of a load/store instruction
	/// Given an LLVM value, insert a cast expressions or cast instructions as
	/// necessary to make the value the specified type.
	///
	/// \param V - The value which needs to be of the specified type.
	/// \param Ty - The type to which V should be casted (if necessary).
	/// \param Name - The name to assign the new casted value (if one is created).
	/// \param InsertPt - The point where a cast instruction should be inserted
	/// \return An LLVM value of the desired type, which may be the original value
	///         passed into the function, a constant cast expression of the passed
	///         in value, or an LLVM cast instruction.
	static inline Value *castTo(Value *V,
		Type *Ty,
		Twine Name,
		Instruction *InsertPt) {
		// Assert that we're not trying to cast a NULL value.
		assert(V && "castTo: trying to cast a NULL Value!\n");

		// Don't bother creating a cast if it's already the correct type.
		if (V->getType() == Ty)
			return V;

		// If it's a constant, just create a constant expression.
		if (Constant *C = dyn_cast<Constant>(V)) {
			Constant *CE = ConstantExpr::getZExtOrBitCast(C, Ty);
			return CE;
		}

		// Otherwise, insert a cast instruction.
		return CastInst::CreateZExtOrBitCast(V, Ty, Name, InsertPt);
	} // End of Utility function

	/// Helper function for printing trace information
	void print_trace_info(TraceType& trace);
	//void write_trace_file(TraceType& trace, std::string fname);
	std::list<uint64_t> getHeaderID_list(uint64_t Aheader_id);

	static ManagedStatic<ProfilerJITContext> GlobalContext;

	// Vector used to record load/store information at runtime
	//std::vector<RecordLoadStore> recordls;

} // End namespace LLVM

namespace {
	/// Load/Store ID to Load/Store Instruction Map, used for testing purpose
	std::map<unsigned, llvm::Instruction*> testing_lsID2lsInstMap;
	std::map<unsigned, llvm::BasicBlock*> testing_bbID2bbMap;

	/// Loop Basic Block ID --> (the corresponding Loop Header BB id, loop depth)
	//typedef std::map<uint64_t, std::pair<uint64_t, unsigned> > loopBBid2HeaderAndloopDepthPairMapTy;
	//loopBBid2HeaderAndloopDepthPairMapTy loopBBid2HeaderAndloopDepthPairMap;

	/*
	/// Loop Basic Block ID --> Vector<Evaluation Iteration Index list>
	//typedef std::map<uint64_t, std::vector<std::list<unsigned> > > lpBBid2EvalIterIndxMapTy;
	/// Evaluation Iteration Index sequence:  The top-level -> lower-level --> the innermost-level loop
	typedef std::map<uint64_t, std::list<unsigned> > lpBBid2EvalIterIndxMapTy;
	lpBBid2EvalIterIndxMapTy lpBBid2EvalIterIndxMap;
	*/

	/// headerID_list_vec contains a vector of headerID_list
	/// headerID_list:  The header BB id of the innermost loop --> The header BB id of the outermost loop
	//typedef std::vector<std::list<uint64_t> > headerID_list_vecTy;
	//headerID_list_vecTy headerID_list_vec;
	//llvm::BBidFreqMap BBidFreq;
	//llvm::BranchidFreqMap BranchIDFreq;

	//typedef std::map<std::string, headerID_list_vecTy> funcName2headIDlistVecMapTy;
	//funcName2headIDlistVecMapTy funcName2headIDlistVecMap;

	typedef std::map<std::string, std::vector<uint64_t> > funcName2exitBBidMapTy;
	funcName2exitBBidMapTy funcName2exitBBidMap;
	std::vector<uint64_t> whole_exitBBid_vec;
	
	typedef std::map<uint64_t, uint64_t> exitBBid2headerBBidMapTy;
	exitBBid2headerBBidMapTy exitBBid2headerBBidMap;

}

#endif // End CODEINSTRUMENTPASS_H