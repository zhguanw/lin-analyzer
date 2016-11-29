#ifndef ASSIGNLOADSTOREID_PASS_H
#define ASSIGNLOADSTOREID_PASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/InstVisitor.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/IR/Verifier.h"
#include "llvm/ADT/SmallString.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/ADT/Statistic.h"

#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"
//#include "llvm/IR/BasicBlock.h"

#include "llvm/Transforms/Scalar.h"
#include "llvm/Analysis/LoopInfo.h"

namespace llvm{

	class AssignLoadStoreID : public ModulePass,
		public InstVisitor<AssignLoadStoreID>{

	public:
		static char ID;
		AssignLoadStoreID();

		void getAnalysisUsage(AnalysisUsage &AU) const;

		bool runOnModule(Module &M);

		void visitLoad(LoadInst& I);

		void visitStore(StoreInst& I);

		/*
		unsigned getlsID(Instruction* I) const;

		Instruction* getlsbylsID(unsigned lsid) const;

		unsigned getlsLineID(Instruction* I) const;

		unsigned getInstLineID(Instruction* I) const;

		Instruction* getInstbyInstLineID(unsigned instLineid) const;
		*/

	private:
		unsigned counter;
		unsigned instID;
		
		
		// Query ID of Load/Store instructions
		typedef std::map<Instruction*, unsigned> InstToIDMapTy;
		InstToIDMapTy Inst2IDmap;
		typedef std::map<unsigned, Instruction*> IDToInstMapTy;
		IDToInstMapTy ID2Instmap;
		/*
		typedef std::map<unsigned, unsigned> InstLineToLoadStoreIDMapTy;
		InstLineToLoadStoreIDMapTy LineID2lsIDmap;
		std::map<Instruction*, std::pair<unsigned, unsigned>> lsInst2lsIDMap;
		typedef std::map<unsigned, Instruction*> LSIDToLSinstMapTy;
		LSIDToLSinstMapTy lsID2lsInstMap;
		*/

		NamedMDNode *NMD;

		MDNode* assignID(Instruction* I, unsigned id);
	};

} // End namespace LLVM

#endif // End ASSIGNLOADSTOREID_PASS_H