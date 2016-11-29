#ifndef QueryBasicBlockID_PASS_H
#define QueryBasicBlockID_PASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
//#include "llvm/IR/InstVisitor.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/IR/Verifier.h"
#include "llvm/IR/Constants.h"
//#include "llvm/IR/BasicBlock.h"

#include "profile_h/AssignBasicBlockIDPass.h"
#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"

#include <map>

namespace llvm{

	class QueryBasicBlockID : public ModulePass{

	public:
		static char ID;
		QueryBasicBlockID();

		void getAnalysisUsage(AnalysisUsage &AU) const;

		unsigned getBBid(BasicBlock* BB) const;

		BasicBlock* getBBbyID(unsigned id);

		uint64_t getNumOfBB() const;

		bool runOnModule(Module &M);

	private:
		uint64_t numOfBB;
		unsigned int counter;
		std::map<std::string, unsigned> BBstr2IDMap;
		std::map<BasicBlock*, unsigned> BB2IDMap;
		std::map<unsigned, BasicBlock*> ID2BBMap;

	};
} // End namespace LLVM

#endif // End QueryBasicBlockID_PASS_H