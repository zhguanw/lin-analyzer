#ifndef QUERYLOADSTOREID_PASS_H
#define QUERYLOADSTOREID_PASS_H

#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
//#include "llvm/IR/InstVisitor.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/IR/Verifier.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/ADT/Statistic.h"
//#include "llvm/IR/BasicBlock.h"

#include "profile_h/AssignLoadStoreIDPass.h"
#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"

#include <vector>

#include <map>

namespace llvm{
	class QueryLoadStoreID : public ModulePass{

	public:
		static char ID;
		QueryLoadStoreID();

		void getAnalysisUsage(AnalysisUsage &AU) const;

		//unsigned getLoadStoreid(Instruction* I) const;

		//Instruction* getLoadStorebyID(unsigned id) const;

		bool runOnModule(Module &M);

		// Query function members
		unsigned getlsID(Instruction* I) const;

		Instruction* getlsbylsID(unsigned lsid) const;

		unsigned getlsLineID(Instruction* I) const;

		unsigned getInstLineID(Instruction* I) const;

		Instruction* getInstbyInstLineID(unsigned instLineid) const;

		BasicBlock* getBBbylsID(unsigned lsid) const;

	private:

		/*
		std::map<std::string, unsigned> INSTstr2IDMap;
		std::map<Instruction*, unsigned> INST2IDMap;
		std::map<unsigned, Instruction*> ID2INSTMap;
		*/

		// Query ID of Load/Store instructions
		typedef std::map<Instruction*, unsigned> InstToIDMapTy;
		//InstToIDMapTy Inst2IDmap;
		typedef std::map<unsigned, Instruction*> IDToInstMapTy;
		//IDToInstMapTy ID2Instmap;

		std::vector<Instruction*> worklist;

		typedef std::map<unsigned, unsigned> InstLineToLoadStoreIDMapTy;
		InstLineToLoadStoreIDMapTy LineID2lsIDmap;
		std::map<Instruction*, std::pair<unsigned, unsigned>> lsInst2lsIDMap;
		typedef std::map<unsigned, Instruction*> LSIDToLSinstMapTy;
		LSIDToLSinstMapTy lsID2lsInstMap;

	};

} // End namespace LLVM

#endif // End QUERYLOADSTOREID_PASS_H