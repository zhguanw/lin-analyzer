/// Format of the NamedMDNode:
/// !Loop_Number_of_kernel = !{!1, !2, !3, ...}
/// !1 = metadata !{metadata !"matrixmul", metadata !"loop1"}
/// !2 = metadata !{metadata !"matrixmul", metadata !"loop2"}
/// !3 = metadata !{metadata !"main", metadata !"loop1"}
/// ...

/// Therefore, to get loop number in a function, we need to check the number of 
/// metadata node that contains the function name.
///

#include "profile_h/LoopNumberPass.h"

#define DEBUG_TYPE "loop-number"

static const char *mdKindName = "Loop_Number_of_kernel";

using namespace llvm;

LoopNumber::LoopNumber() : LoopPass(ID) {
	//func_hasLoop_number = 0;
	//loopNumInaFunc.clear();
	loop_counter = 0;
	LoopsMetadataNode.clear();
	DEBUG(dbgs() << "Initialize LoopNumber pass\n");
	initializeLoopNumberPass(*PassRegistry::getPassRegistry());
}

bool LoopNumber::doInitialization(Loop* lp, LPPassManager &LPM) {
	Module* M = lp->getHeader()->getParent()->getParent();
	if (!M->getNamedMetadata(mdKindName)) {
		NMD = M->getOrInsertNamedMetadata(mdKindName);
	}

	return true;
}

void LoopNumber::getAnalysisUsage(AnalysisUsage &AU) const {
	AU.addRequired<LoopInfo>();
	AU.addPreserved<LoopInfo>();

	AU.addRequiredID(LoopSimplifyID);
	AU.addPreservedID(LoopSimplifyID);

	//AU.addRequiredID(LCSSAID);
	//AU.addPreservedID(LCSSAID);

	AU.setPreservesCFG();
	AU.setPreservesAll();
	// FIXME: Loop unroll requires LCSSA. And LCSSA requires dom info.
	// If loop unroll does not preserve dom info then LCSSA pass on next
	// loop will receive invalid dom info.
	// For now, recreate dom info, if loop is unrolled.
	//AU.addPreserved<DominatorTreeWrapperPass>();
}

bool LoopNumber::runOnLoop(Loop* lp, LPPassManager &LPM) {

	/// Check whether loop is in a simplify form
	if (!lp->isLoopSimplifyForm()) {
		assert(false && "Loop is not in a Simplify Form!\n");
	}

	Function* func = lp->getHeader()->getParent();
	/// Ignore uninterested functions
	std::vector<std::string>::iterator fnName_it;
	fnName_it = std::find(kernel_names.begin(), kernel_names.end(), func->getName());
	if (fnName_it == kernel_names.end()) {
		// This function is uninteresting, just ignore it.
		return false;
	}

	LLVMContext &Context = func->getContext();

	std::vector<Function*>::iterator it;
	it = std::find(exploredFunc.begin(), exploredFunc.end(), func);
	if (it == exploredFunc.end()) {
		// Reset loop counter for a new function.
		loop_counter = 0;
		exploredFunc.push_back(func);
	}

	/*
	std::string func_name = func->getName();
	std::map<std::string, unsigned>::iterator it;
	it = loopNumInaFunc.find(func_name);
	if (it == loopNumInaFunc.end()) {
		// Create an entry of the map for this function and initialize the number
		// of loops in the function as 0
		loopNumInaFunc.insert(std::make_pair(func_name, 0));
		MDString* func_MDname = MDString::get(Context, func->getName());
		LoopsMetadataNode.push_back(func_MDname);
		// Initialize loop_counter
		loop_counter = 0;
	}*/

	unsigned depth = lp->getLoopDepth();
	/// Trace the top-level loop of a loop
	if (depth == 1) {
		loop_counter++;

		MDString* func_name = MDString::get(Context, func->getName());
		LoopsMetadataNode.push_back(func_name);
		std::string loop_name = "loop" + std::to_string(loop_counter);
		LoopsMetadataNode.push_back(MDString::get(Context, loop_name));
		MDNode* MD = MDNode::getWhenValsUnresolved(Context, ArrayRef<Value*>(LoopsMetadataNode), false);
		NMD->addOperand(MD);
		LoopsMetadataNode.clear();
	}

	/// Update loop number in a function
	//loopNumInaFunc[func_name] = loop_counter;

	return true;
}

/*
unsigned LoopNumber::getNumberLoopOfFunc(std::string func_name) const{
	return loopNumInaFunc.at(func_name);
}

unsigned LoopNumber::getNumberFuncHasLoop() const {
	return loopNumInaFunc.size();
}*/

char LoopNumber::ID = 0;
INITIALIZE_PASS_BEGIN(LoopNumber, "LoopNumber",
					  "This pass is used to calculate the number of loops inside functions.",
					  false,
					  true)
INITIALIZE_PASS_DEPENDENCY(LoopInfo)
INITIALIZE_PASS_DEPENDENCY(LoopSimplify)
INITIALIZE_PASS_END(LoopNumber, "LoopNumber",
					"This pass is used to calculate the number of loops inside functions.",
					false,
					true)

Pass* llvm::createLoopNumberPass() {
	return new LoopNumber();
}