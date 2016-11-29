#include "profile_h/QueryBasicBlockIDPass.h"

#define DEBUG_TYPE "query-bb-id"
using namespace llvm;

static const char *mdKindName = "recordBBid"; // Record Basic Block ID

QueryBasicBlockID::QueryBasicBlockID() : ModulePass(ID) {
	numOfBB = 0;
	counter = 0;
	numOfBB = 0;
	DEBUG(dbgs() << "Initialize QueryBasicBlockID pass\n");
	initializeQueryBasicBlockIDPass(*PassRegistry::getPassRegistry());
}

void QueryBasicBlockID::getAnalysisUsage(AnalysisUsage &AU) const {
	//AU.addRequired<AssignBasicBlockID>();
	//AU.setPreservesCFG();
	AU.setPreservesAll();
}

unsigned QueryBasicBlockID::getBBid(BasicBlock* BB) const {
	std::map<BasicBlock*, unsigned>::const_iterator it;
	it = BB2IDMap.find(BB);
	if (it != BB2IDMap.end()) {
		return BB2IDMap.at(BB);
	}
	else {
		assert(false && "Error: Basic Block is not in the BB2IDMap\n");
		return 1;
	}
}

BasicBlock* QueryBasicBlockID::getBBbyID(unsigned id) {
	std::map<unsigned, BasicBlock*>::const_iterator it;
	it = ID2BBMap.find(id);
	if (it != ID2BBMap.end()) {
		return ID2BBMap.at(id);
	}
	else {
		assert(false && "Error: Cannot find ID of the Basic Block!\n");
		return NULL;
	}
}

bool QueryBasicBlockID::runOnModule(Module &M) {
	DEBUG(dbgs() << "\n\nBegin QueryBasicBlockID Pass :\n");
	// Create a named Metadata node:
	
	NamedMDNode *NMD = M.getNamedMetadata(mdKindName);

	/// Scan the NamedMDNode we get and extract its BB-ID pairs
	/// Actually, we can remove the metadata for BB and function name
	for (unsigned i = 0; i < NMD->getNumOperands(); i++) {
		MDNode* node = NMD->getOperand(i);
		MDString* Funcmdstr = dyn_cast<MDString>(node->getOperand(0));
		std::string Func = Funcmdstr->getString();
		MDString* BBmdstr = dyn_cast<MDString>(node->getOperand(1));
		std::string BB = BBmdstr->getString();
		ConstantInt* ID = dyn_cast<ConstantInt>(node->getOperand(2));
		unsigned id = static_cast<unsigned> (ID->getZExtValue());
		//errs() << "BB = " << BBmdstr->getString() << "  -> ID = " \
		//					<< *ID << "\n";
		BBstr2IDMap.insert(std::make_pair(Func+BB, id));
	}

	// Convert BBstr to Basic Block
	std::vector<std::string>::iterator fn_it;
	unsigned counter = 0;
	for (Module::iterator FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		if (FI->isDeclaration()) {
			continue;
		}
		fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (fn_it == kernel_names.end()) {
			// Ignore uninteresting functions
			continue;
		}
		for (Function::iterator BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
			// Check if we have BB in BBstr2IDMap, then we store the BB with id into 
			// BB2IDMap
			/*
			std::map<std::string, unsigned>::iterator it = BBstr2IDMap.find(BI->getName());
			if (BBstr2IDMap.end() != it){
				BB2IDMap.insert(std::make_pair(BI, it->second));
				ID2BBMap.insert(std::make_pair(it->second, BI));
			}*/
			BB2IDMap.insert(std::make_pair(BI, counter));
			ID2BBMap.insert(std::make_pair(counter, BI));
			counter++;
		}
	}
	numOfBB = BB2IDMap.size();
	if (numOfBB != counter) {
		assert(false && "Error: Number of Basic Block is wrong\n");
	}

	//Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		return false;
	}

	DEBUG(dbgs() << "\nEnd QueryBasicBlockID Pass\n\n");
	return false;
}

uint64_t QueryBasicBlockID::getNumOfBB() const {
	return numOfBB;
}

char QueryBasicBlockID::ID = 0;
INITIALIZE_PASS(QueryBasicBlockID, "queryBBid",
	"This pass is used to get the unique id of specific basic block",
	false,
	false
	)

ModulePass *llvm::createQueryBasicBlockIDPass() {
	return new 	QueryBasicBlockID();
}