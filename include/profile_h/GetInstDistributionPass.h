#ifndef GETINSTDISTRIBUTION_H
#define GETINSTDISTRIBUTION_H

#include <fstream>

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/InstVisitor.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/Verifier.h"
//#include "llvm/IR/BasicBlock.h"

#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/Casting.h"

#include "llvm/Transforms/Scalar.h"
#include "llvm/ADT/Statistic.h"

#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"
#include "profile_h/QueryBasicBlockIDPass.h"

namespace llvm{
	
	struct instDis {
		uint64_t meminst;
		uint64_t ldinst;
		uint64_t stinst;
		uint64_t brinst;
		uint64_t compinst;
		uint64_t others;
	};

	typedef std::map<uint64_t, instDis> bbID2instDisMapTy;

	class GetInstDistribution : public ModulePass{

	public:
		static char ID;
		GetInstDistribution();

		void getAnalysisUsage(AnalysisUsage &AU) const;

		bool runOnModule(Module &M);

	private:
		bbID2instDisMapTy bbID2instDisMap;
		QueryBasicBlockID* bbIDpass;

	};

} // End namespace LLVM

#endif // End GETINSTDISTRIBUTION_H