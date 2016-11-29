/// This pass is used to calculate the number of 
/// loops inside functions.

#ifndef LOOPNUMBERPASS_H
#define LOOPNUMBERPASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Support/raw_os_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/Transforms/Scalar.h"

#include "profile_h/Passes.h"
#include "profile_h/lin-profile.h"

namespace llvm {

	class LoopNumber : public LoopPass {
		
	public:
		static char ID;
		LoopNumber();
		bool doInitialization(Loop *lp, LPPassManager &LPM);
		void getAnalysisUsage(AnalysisUsage &AU) const;
		bool runOnLoop(Loop* lp, LPPassManager &LPM);
		/*
		unsigned getNumberLoopOfFunc(std::string func_name) const;
		unsigned getNumberFuncHasLoop() const;
		*/

	private:
		//unsigned func_hasLoop_number;
		unsigned loop_counter;
		//std::map<std::string, unsigned> loopNumInaFunc;
		NamedMDNode* NMD;
		std::vector<Value*> LoopsMetadataNode;
		std::vector<Function* > exploredFunc;
	};

} // End of llvm namespace

#endif // End of definition of LOOPNUMBERPASS_H