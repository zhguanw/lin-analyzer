#ifndef EXTRACTLOOPINFOPASS_H
#define EXTRACTLOOPINFOPASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Support/raw_os_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/IR/Metadata.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/IR/Dominators.h"
#include "llvm/Analysis/CFGPrinter.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/Analysis/ScalarEvolution.h"

#include "profile_h/QueryBasicBlockIDPass.h"
#include "profile_h/QueryLoadStoreIDPass.h"
#include "profile_h/LoopNumberPass.h"
#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"
#include "profile_h/DDDG.h"

namespace llvm {

	/// Basic Block --> Loop Level Map
	//typedef std::map<BasicBlock*, unsigned> BB2LoopLevelMapTy;
	/// Load/Store ID --> (Basic Block, Loop Level Map) Map
	typedef std::map<unsigned, std::pair<BasicBlock*, unsigned> > LSID2BBLoopLevelPairMapTy;
	/// Loop Name (NOT loop level) --> LSID2BBLoopLevelPairMapTy Map
	typedef std::map<std::string, LSID2BBLoopLevelPairMapTy > LoopID2LSInfoMapTy;
	/// Function Name --> vector<LoopID2LSInfoMapTy> Map
	typedef std::map<std::string, LoopID2LSInfoMapTy > Func2LoopInfoMapTy;

	class ExtractLoopInfo : public LoopPass {
	public:
		static char ID;

		ExtractLoopInfo();
		bool doInitialization(Loop *lp, LPPassManager &LPM);
		void getAnalysisUsage(AnalysisUsage &AU) const;
		bool runOnLoop(Loop* lp, LPPassManager &LPM);
		bool isPerfectNest(Loop *L, LoopInfo* li);
		bool hasNoMemoryOps(BasicBlock *b);
		LoopID2LSInfoMapTy getFunc2LoopInfomap(std::string func_name) const;

	private:
		//LoopNumber* loopNumberpass;
		//QueryLoadStoreID* lsIDpass;

		//BB2LoopLevelMapTy BB2LoopLevelmap;
		LSID2BBLoopLevelPairMapTy LSID2_BB2LoopLevel_map;
		LoopID2LSInfoMapTy LoopID2LSInfomap;
		Func2LoopInfoMapTy Func2LoopInfomap;

		//std::vector<LSID2BBLoopLevelPairMapTy> LSID2BBLoopLevelMap_vec;
		//std::vector<LoopID2LSInfoMapTy> LoopID2LSInfoMap_vec;
		std::vector<Function* > exploredFunc;
		unsigned loopID;
		unsigned depth;
		NamedMDNode* NMD;
		NamedMDNode* loopNumNMD;
		unsigned countNumLoopInaFunc;
		std::vector<Value*> LoopsMetadataNode;
		//std::map<std::string, unsigned> funcName2loopNumMap;

		/// Load/Store Instruction -> Load/StoreID
		std::map<Instruction*, std::pair<unsigned, unsigned> > LoadStoreIDMap;
		std::map<unsigned, Instruction*> LineID2LSMap;
		NamedMDNode* recordLSMDNode;
		bool alreadyCheck;

		typedef std::map<BasicBlock*, uint64_t> BB2BBidMapTy;
		BB2BBidMapTy BB2BBidMap;

		unsigned depth_record;

		std::vector<Loop *> exploredLoop;
	};

} // End namespace LLVM

#endif // End of definition of EXTRACTLOOPINFOPASS_H