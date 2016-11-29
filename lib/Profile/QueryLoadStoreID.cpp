#include "profile_h/QueryLoadStoreIDPass.h"

#define DEBUG_TYPE "query-laodstore-id"
using namespace llvm;

static const char *mdKindName = "recordLoadStoreID"; // Record Load/Store ID

STATISTIC(NumInst, "The number of total instructions: #");
STATISTIC(NumLoadInst, "The number of load instructions: #");
STATISTIC(NumStoreInst, "The number of store instructions: #");

QueryLoadStoreID::QueryLoadStoreID() : ModulePass(ID) {
	DEBUG(dbgs() << "Initialize QueryLoadStoreID pass\n");
	initializeQueryLoadStoreIDPass(*PassRegistry::getPassRegistry());
}

void QueryLoadStoreID::getAnalysisUsage(AnalysisUsage &AU) const{
	//AU.addRequired<AssignLoadStoreID>();
	//AU.setPreservesCFG();
	AU.setPreservesAll();
}

/*
unsigned QueryLoadStoreID::getLoadStoreid(Instruction* I) const{
	std::map<Instruction*, unsigned>::const_iterator it;
	it = INST2IDMap.find(I);
	if (it != INST2IDMap.end()) {
		return INST2IDMap.at(I);
	}
	else {
		errs() << "Basic Block is not in the BB2IDMap, ERROR\n";
		return false;
	}
}

Instruction* QueryLoadStoreID::getLoadStorebyID(unsigned id) const{
	std::map<unsigned, Instruction*>::const_iterator it;
	it = ID2INSTMap.find(id);
	if (it != ID2INSTMap.end()) {
		return ID2INSTMap.at(id);
	}
	else {
		errs() << "Can not find ID of the Basic Block, error!\n";
		return NULL;
	}
}*/

bool QueryLoadStoreID::runOnModule(Module &M) {
	DEBUG(dbgs() << "\n\nBegin QueryLoadStoreID Pass :\n");
	// Store all instructions in a worklist;
	Module::iterator FI, FE;
	std::vector<std::string>::iterator fn_it;

	unsigned instid = 0;
	unsigned counter = 0;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		if (FI->isDeclaration()) {
			continue;
		}
		fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (fn_it == kernel_names.end()) {
			// Ignore uninteresting functions
			continue;
		}

		for (inst_iterator I = inst_begin(FI), E = inst_end(FI); I != E; ++I) {
			worklist.push_back(&*I);
			if (isa<LoadInst>(&*I) || isa<StoreInst>(&*I)) {
				counter++;
			}
			instid++;
		}
	}

	// Print number of instructions in the module:
	NumInst = worklist.size();
	DEBUG(dbgs() << "\tNumber of the total instructions in the kernal functions: " << NumInst << "\n\n");

	// Create a named Metadata node:
	NamedMDNode *NMD = M.getNamedMetadata(mdKindName);

	// Scan the NamedMDNode we get and extract its BB-ID pairs
	for (unsigned i = 0; i < NMD->getNumOperands(); i++) {
		MDNode* node = NMD->getOperand(i);
		MDString* Instmdstr = dyn_cast<MDString>(node->getOperand(0));
		std::string Inststr = Instmdstr->getString();

		ConstantInt* lineID = dyn_cast<ConstantInt>(node->getOperand(1));
		unsigned lineid = static_cast<unsigned> (lineID->getZExtValue());

		ConstantInt* lsID = dyn_cast<ConstantInt>(node->getOperand(2));
		unsigned lsid = static_cast<unsigned> (lsID->getZExtValue());

		Instruction* inst = worklist[lineid];
		//errs() << "Instruction = " << Instmdstr->getString() << "  -> ID = " \
		//					<< *ID << "\n";
		if (Inststr == "LoadInst") {
			NumLoadInst++;
		}
		else if (Inststr == "StoreInst") {
			NumStoreInst++;
		}
		else {
			// Do nothing
		}

		lsID2lsInstMap.insert(std::make_pair(lsid, inst));
		lsInst2lsIDMap.insert(std::make_pair(inst, std::make_pair(lineid, lsid)));
	}

	// Check how many load/store instructions
	DEBUG(dbgs() << "\tLoad/Store instruction number: " << counter << "\n");
	DEBUG(dbgs() << "\tLoad instruction number: " << NumLoadInst << "\n");
	DEBUG(dbgs() << "\tStore instruction number: " << NumStoreInst << "\n");

	unsigned num_LoadStore = lsID2lsInstMap.size();
	if (num_LoadStore != counter) {
		assert(false && "Error: Number of Load/Store instructions is wrong\n");
	}

	//Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		return false;
	}

	DEBUG(dbgs() << "\nEnd QueryLoadStoreID Pass\n\n");
	return false;
}

unsigned QueryLoadStoreID::getlsID(Instruction* I) const {
	if (isa<LoadInst>(I) || isa<StoreInst>(I)) {
		std::map<Instruction*, std::pair<unsigned, unsigned>>::const_iterator it;
		it = lsInst2lsIDMap.find(I);
		if (it != lsInst2lsIDMap.end()) {
			return lsInst2lsIDMap.at(I).second;
		}
		else {
			DEBUG(dbgs() << "The instruction caused the problem: " << *I << "\n");
			assert(false && "Error: The instruction is not a Load or store instruction\n");
			return 1;
		}
	}
	else {
		DEBUG(dbgs() << "The instruction caused the problem: " << *I << "\n");
		assert(false && "Error: The instruction is not a Load or store instruction\n");
		return 1;
	}
}

Instruction* QueryLoadStoreID::getlsbylsID(unsigned lsid) const {
	LSIDToLSinstMapTy::const_iterator it;
	it = lsID2lsInstMap.find(lsid);
	if (it != lsID2lsInstMap.end()) {
		return lsID2lsInstMap.at(lsid);
	}
	else {
		assert(false && "Error: Could not find the lsID\n");
		return NULL;
	}
}

unsigned QueryLoadStoreID::getlsLineID(Instruction* I) const {
	if (isa<LoadInst>(I) || isa<StoreInst>(I)) {
		std::map<Instruction*, std::pair<unsigned, unsigned>>::const_iterator it;
		it = lsInst2lsIDMap.find(I);
		if (it != lsInst2lsIDMap.end()) {
			return lsInst2lsIDMap.at(I).first;
		}
		else {
			DEBUG(dbgs() << "The instruction caused the problem: " << *I << "\n");
			assert(false && "Error: The instruction is not a Load or store instruction\n");
			return 1;
		}
	}
	else {
		assert(false && "Error: The instruction is not a Load or store instruction\n");
		return 1;
	}
}

unsigned QueryLoadStoreID::getInstLineID(Instruction* I) const {
	if (isa<Instruction>(I)) {
		std::map<Instruction*, std::pair<unsigned, unsigned>>::const_iterator it;
		it = lsInst2lsIDMap.find(I);
		if (it != lsInst2lsIDMap.end()) {
			return lsInst2lsIDMap.at(I).first;
		}
		else {
			DEBUG(dbgs() << "The instruction caused the problem: " << *I << "\n");
			assert(false && "Error: The instruction is not a Load or store instruction\n");
			return 1;
		}
	}
	else {
		assert(false && "Error: The instruction is not a Load or store instruction\n");
		return 1;
	}
}

Instruction* QueryLoadStoreID::getInstbyInstLineID(unsigned instLineid) const {
	
	if (instLineid < worklist.size()) {
		return worklist.at(instLineid);
	}
	else {
		assert(false && "Error: Could not find the line ID\n");
		return NULL;
	}
}

BasicBlock* QueryLoadStoreID::getBBbylsID(unsigned lsid) const {
	Instruction* inst = getlsbylsID(lsid);
	return inst->getParent();
}

char QueryLoadStoreID::ID = 0;
INITIALIZE_PASS(QueryLoadStoreID, "queryLoadStoreid",
				"This pass is used to get the unique id of specific Load/Store",
				false,
				false
				)

	ModulePass *llvm::createQueryLoadStoreIDPass() {
	return new 	QueryLoadStoreID();
}