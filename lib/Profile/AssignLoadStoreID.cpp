#include "profile_h/AssignLoadStoreIDPass.h"

#define DEBUG_TYPE "assign-loadstore-id"
using namespace llvm;

static const char *mdKindName = "recordLoadStoreID"; 
// Record Load/Store ID and its memory address

//STATISTIC(NumInst, "The number of total instructions: #");
//STATISTIC(NumLoadInst, "The number of load instructions: #");
//STATISTIC(NumStoreInst, "The number of store instructions: #");

/*
static inline void intToStr(intptr_t V, SmallString<36> &S) {
	static char encoding_table[] = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
		'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
		'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
		'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
		'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
		'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
		'w', 'x', 'y', 'z', '0', '1', '2', '3',
		'4', '5', '6', '7', '8', '9' }; //, '+', '/'};
	const unsigned table_size = array_lengthof(encoding_table);

	assert(V && "Cannot convert 0 yet!");
	while (V) {
		unsigned char Digit = V % table_size;
		S += encoding_table[Digit];
		V /= table_size;
	}
}

static StringRef translatePtr2Str(void *V, SmallString<36> &S) {
	S.push_back('_');
	intToStr(intptr_t(V), S);
	S.push_back('_');
	return S.str();
}*/

AssignLoadStoreID::AssignLoadStoreID() : ModulePass(ID) {
	DEBUG(dbgs() << "Initialize AssignLoadStoreID pass\n");
	counter = 0;
	initializeAssignLoadStoreIDPass(*PassRegistry::getPassRegistry());
}

void AssignLoadStoreID::getAnalysisUsage(AnalysisUsage &AU) const {
	//AU.addRequiredID(LoopSimplifyID);
	//AU.addPreservedID(LoopSimplifyID);

	//AU.addRequiredID(LCSSAID);
	//AU.addPreservedID(LCSSAID);

	AU.setPreservesCFG();
}

bool AssignLoadStoreID::runOnModule(Module &M) {
	DEBUG(dbgs() << "\n\nBegin AssignLoadStoreID Pass :\n");
	// Create a named Metadata node:
	NMD = M.getOrInsertNamedMetadata(mdKindName);
			
	// Initialize Load/Store ID
	counter = 0;
	instID = 0;

	//errs() << "\nCheck Inst2ID map \n";
	Module::iterator FI, FE;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {

		if (FI->isDeclaration()) {
			continue;
		}

		std::vector<std::string>::iterator it;
		it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (it == kernel_names.end()) {
			// We are not interested in this function
			continue;
		}

		for (inst_iterator I = inst_begin(FI), E = inst_end(FI); I != E; ++I) {
			Instruction* inst = &*I;
			Inst2IDmap.insert(std::make_pair(inst, instID));
			ID2Instmap.insert(std::make_pair(instID, inst));
			// Check Inst2IDmap:
			//errs() << "Instruction: " << *inst << "  ID: " << instID << "\n";
			instID++;
		}
	}

	// Print number of instructions in the module:
	//errs() << "\nNumber of instructions in the module: " << instID << "\n\n";

	/*
	Value *Elts[] = {
		MDString::get(M.getContext(), "test1")
	};
	MDNode *Node = MDNode::get(M.getContext(), Elts);
	NMD->addOperand(Node);
	NMD->dump(); */

	// Scan through all instructions in the module and assign unique ID
	visit(M);

	// Check how many basic blocks
	DEBUG(dbgs() << "\tLoad/Store instruction number: " << counter << "\n");
	//errs() << "Load instruction number: " << NumLoadInst << "\n";
	//errs() << "Store instruction number: " << NumStoreInst << "\n";

	//Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		return false;
	}

	DEBUG(dbgs() << "End AssignLoadStoreID Pass\n\n");
	return true;
}

void AssignLoadStoreID::visitLoad(LoadInst& I) {
	Function* fn = I.getParent()->getParent();
	std::vector<std::string>::iterator it;
	it = std::find(kernel_names.begin(), kernel_names.end(), fn->getName());
	if (it != kernel_names.end()) {
		// We only care interesting kernel functions
		NMD->addOperand(assignID(&I, counter));
		counter++;
		//NumLoadInst++;
	}
}

void AssignLoadStoreID::visitStore(StoreInst& I) {
	Function* fn = I.getParent()->getParent();
	std::vector<std::string>::iterator it;
	it = std::find(kernel_names.begin(), kernel_names.end(), fn->getName());
	if (it != kernel_names.end()) {
		// We only care interesting kernel functions
		NMD->addOperand(assignID(&I, counter));
		counter++;
		//NumStoreInst++;
	}
}

/*
unsigned AssignLoadStoreID::getlsID(Instruction* I) const {
	if (isa<LoadInst>(I) && isa<StoreInst>(I)) {
		std::map<Instruction*, std::pair<unsigned, unsigned>>::const_iterator it;
		it = lsInst2lsIDMap.find(I);
		if (it != lsInst2lsIDMap.end()) {
			return lsInst2lsIDMap.at(I).second;
		}
		else {
			errs() << "Error! The instruction is not a Load or store instruction\n";
			return -1;
		}
	}
	else {
		errs() << "Error! The instruction is not a Load or store instruction\n";
		return -1;
	}
}

Instruction* AssignLoadStoreID::getlsbylsID(unsigned lsid) const {
	LSIDToLSinstMapTy::const_iterator it;
	it = lsID2lsInstMap.find(lsid);
	if (it != lsID2lsInstMap.end()) {
		return lsID2lsInstMap.at(lsid);
	}
	else {
		errs() << "Error! Could not find the lsID\n";
		return NULL;
	}
}

unsigned AssignLoadStoreID::getlsLineID(Instruction* I) const {
	if (isa<LoadInst>(I) && isa<StoreInst>(I)) {
		std::map<Instruction*, std::pair<unsigned, unsigned>>::const_iterator it;
		it = lsInst2lsIDMap.find(I);
		if (it != lsInst2lsIDMap.end()) {
			return lsInst2lsIDMap.at(I).first;
		}
		else {
			errs() << "Error! The instruction is not a Load or store instruction\n";
			return -1;
		}
	}
	else {
		errs() << "Error! The instruction is not a Load or store instruction\n";
		return -1;
	}
}

unsigned AssignLoadStoreID::getInstLineID(Instruction* I) const {
	if (isa<Instruction>(I)) {
		InstToIDMapTy::const_iterator it;
		it = Inst2IDmap.find(I);
		if (it != Inst2IDmap.end()) {
			return Inst2IDmap.at(I);
		}
		else {
			errs() << "Error! The instruction is not a Load or store instruction\n";
			return -1;
		}
	}
	else {
		errs() << "Error! The instruction is not a Load or store instruction\n";
		return -1;
	}
}


Instruction* AssignLoadStoreID::getInstbyInstLineID(unsigned instLineid) const {
	IDToInstMapTy::const_iterator it;
	it = ID2Instmap.find(instLineid);
	if (it != ID2Instmap.end()) {
		return ID2Instmap.at(instLineid);
	}
	else {
		errs() << "Error! Could not find the line ID\n";
		return NULL;
	}
}
*/

MDNode* AssignLoadStoreID::assignID(Instruction* I, unsigned id) {
	// Fetch the context in which the enclosing module was defined
	LLVMContext &Context = I->getContext();

	// Get instruction's ID
	unsigned instid = Inst2IDmap.at(I);
	//errs() << "\n Load/Store Line ID: " << instid << "\n";

	if (isa<LoadInst>(I)) {
		// Create a metadata node that contains ID as a constant:
		Value* lsID[3] = {
			MDString::get(Context, "LoadInst"),
			ConstantInt::get(Type::getInt32Ty(Context), instid),
			ConstantInt::get(Type::getInt32Ty(Context), id)
		};

		//lsID2lsInstMap.insert(std::make_pair(id, I));
		//lsInst2lsIDMap.insert(std::make_pair(I, std::make_pair(instid, id)));
		return MDNode::getWhenValsUnresolved(Context, lsID, false);
	}
	else if (isa<StoreInst>(I)) {
		// Create a metadata node that contains ID as a constant:
		Value* lsID[3] = {
			MDString::get(Context, "StoreInst"),
			ConstantInt::get(Type::getInt32Ty(Context), instid),
			ConstantInt::get(Type::getInt32Ty(Context), id)
		};
		//lsID2lsInstMap.insert(std::make_pair(id, I));
		//lsInst2lsIDMap.insert(std::make_pair(I, std::make_pair(instid, id)));
		return MDNode::getWhenValsUnresolved(Context, lsID, false);
	}
	else {
		assert(false && "Error, we invoke assignID for non-load/store instructions");
		return MDNode::get(Context, MDString::get(Context, "Error, we invoke assignID for non-load/store instructions"));
	}
}

char AssignLoadStoreID::ID = 0;
INITIALIZE_PASS(AssignLoadStoreID, "assignLoadStoreID",
				"This pass is used to assign unique id for each load/store \
				instruction",
				false,
				true	/*false  We need to modify the code later. In initial step, we
						 just check the unique ID first*/
				)

ModulePass *llvm::createAssignLoadStoreIDPass() {
	return new 	AssignLoadStoreID();
}