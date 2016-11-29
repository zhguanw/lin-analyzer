#include "profile_h/AssignBasicBlockIDPass.h"

#define DEBUG_TYPE "assign-bb-id"
using namespace llvm;

static const char *mdKindName = "recordBBid"; // Record Basic Block ID

funcBBNmPair2numInstInBBMapTy funcBBNmPair2numInstInBBMap;
getElementPtrName2arrayNameMapTy getElementPtrName2arrayNameMap;

AssignBasicBlockID::AssignBasicBlockID() : ModulePass(ID) {
	errs() << "===========================\n";
	errs() << "DEBUG-INFO: [static-analysis_bb-info-extraction] Extracting basic block information\n";
	errs() << "===========================\n";
	counter = 0;
	DEBUG(dbgs() << "Initialize AssignBasicBlockID pass\n");
	initializeAssignBasicBlockIDPass(*PassRegistry::getPassRegistry());
}

void AssignBasicBlockID::getAnalysisUsage(AnalysisUsage &AU) const {
	//AU.addRequiredID(LoopSimplifyID);
	//AU.addPreservedID(LoopSimplifyID);

	//AU.addRequiredID(LCSSAID);
	//AU.addPreservedID(LCSSAID);

	AU.setPreservesCFG();
}

bool AssignBasicBlockID::runOnModule(Module &M) {
	DEBUG(dbgs() << "\n\nBegin AssignBasicBlockID Pass :\n");
	// Create a named Metadata node:
	NamedMDNode *NMD = M.getOrInsertNamedMetadata(mdKindName);
			
	// Initialize BB ID
	counter = 0;

	// Scan through all basic block in the module and assign unique ID
	Module::iterator FI, FE;
	Function::iterator BI, BE;
	BasicBlock::iterator bbit, bbie;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		//runOnFunction(*FI);
		if (FI->isDeclaration()) {
			continue;
		}

		std::vector<std::string>::iterator it;
		it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (it == kernel_names.end()) {
			// We are not interested in this function
			continue;
		}

		std::string func_name = FI->getName();
		DEBUG(dbgs() << "Function name = " << func_name << "\n");
		for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
			NMD->addOperand(assignID(BI, counter));
			counter++;
			
			// Count instructions inside a basic block
			unsigned num_inst = 0;
			std::string bb_name = BI->getName();
			DEBUG(dbgs() << "BB name = " << bb_name << "\n");
			for (bbit = BI->begin(), bbie = BI->end(); bbit != bbie; ++bbit) {
				if (CallInst* CI = dyn_cast<CallInst>(bbit)) {
					if (CI->isTailCall() == false) {
						num_inst++;
					}
					else {
						// Do not count for tail call instructions.
						//errs() << "This call instruction is a tail call\n";
					}
				}
				else {
					//errs() << *bbit << "\n";
					num_inst++;
				}

				if (GetElementPtrInst* getElePtrInst = dyn_cast<GetElementPtrInst>(bbit)) {
					std::string getElementPtrName = getElePtrInst->getName().str();
					std::string arrayName = getElePtrInst->getOperand(0)->getName().str();
					getElementPtrName2arrayNameMapTy::iterator it_arrayName = getElementPtrName2arrayNameMap.find(getElementPtrName);
					if (it_arrayName == getElementPtrName2arrayNameMap.end()) {
						getElementPtrName2arrayNameMap.insert(std::make_pair(getElementPtrName, arrayName));
					}
				}
			}
			//errs() << "num inst = " << num_inst << "\n";
			funcBBNmPair2numInstInBBMap.insert(std::make_pair(std::make_pair(func_name, bb_name), num_inst));
		}
	}

	// Check how many basic blocks
	DEBUG(dbgs() << "\tBasic Block number: " << counter << "\n");

	//Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		return false;
	}

	DEBUG(dbgs() << "End AssignBasicBlockID Pass\n\n");
	return true;
}

MDNode* AssignBasicBlockID::assignID(BasicBlock *BB, unsigned id) {
	// Fetch the context in which the enclosing module was defined
	LLVMContext &Context = BB->getContext();
	std::string BBname = BB->getName();
	std::string func_name = BB->getParent()->getName();

	// Create a metadata node that contains ID as a constant:
	Value* ID[]{
		MDString::get(Context, func_name),
		MDString::get(Context, BBname),
		ConstantInt::get(Type::getInt32Ty(Context), id)
	};

	return MDNode::getWhenValsUnresolved(Context, ArrayRef<Value*>(ID, 3), false);
}
/*
MDNode* assignID (Instruction *I, unsigned id) {
	// Fetch the context in which the enclosing module was defined
	LLVMContext &Context = I->getContext();

	// Create a metadata node that contains ID as a constant:
	Value* ID[]{
		MDString::get(Context, I->getParent()->getName()),
		ConstantInt::get(Type::getInt32Ty(Context), id)
	};
			
	return MDNode::getWhenValsUnresolved(Context, ArrayRef<Value*>(ID, 2), true);
}*/

char AssignBasicBlockID::ID = 0;

INITIALIZE_PASS(AssignBasicBlockID, "assignBBid",
				"This pass is used to assign unique id for each basic block",
				false,
				true
				)

/*
INITIALIZE_PASS_BEGIN(AssignBasicBlockID, "assignBBid", 
				"This pass is used to assign unique id for each basic block",
				false,
				true	
				)
INITIALIZE_PASS_DEPENDENCY(LoopInfo)
INITIALIZE_PASS_DEPENDENCY(LoopSimplify)
INITIALIZE_PASS_END(AssignBasicBlockID, "assignBBid",
				"This pass is used to assign unique id for each basic block",
				false,
				true	
				)
*/

ModulePass *llvm::createAssignBasicBlockIDPass() {
	return new 	AssignBasicBlockID();
}