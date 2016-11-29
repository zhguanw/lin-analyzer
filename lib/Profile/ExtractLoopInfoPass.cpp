#include "profile_h/ExtractLoopInfoPass.h"
#include "profile_h/auxiliary.h"

#define DEBUG_TYPE "extract-loopinfo"

namespace llvm {
	fnlpNamePair2BinaryOpBBidMapTy fnlpNamePair2BinaryOpBBidMap;
	funcName2loopNumMapTy funcName2loopNumMap;
	BB2loopNameMapTy BB2loopNameMap;
#ifdef BUILD_DDDG_H
	bbFuncNamePair2lpNameLevelPairMapTy bbFuncNamePair2lpNameLevelPairMap;
	bbFuncNamePair2lpNameLevelPairMapTy headerBBFuncnamePair2lpNameLevelPairMap;
	bbFuncNamePair2lpNameLevelPairMapTy exitBBFuncnamePair2lpNameLevelPairMap;
#endif // End of BUILD_DDDG_H
	LpName2numLevelMapTy LpName2numLevelMap;
}

static const char *mdKindName = "Loop_Info_of_this_kernel";
static const char *loopNumberMDNodeName = "Loop_Number_of_kernel";
static const char *recordLoadStoreIDMDNodeName = "recordLoadStoreID";
lpNameLevelPair2headBBnameMapTy lpNameLevelPair2headBBnameMap;
lpNameLevelPair2headBBnameMapTy lpNameLevelPair2exitingBBnameMap;
wholeloopName2loopBoundMapTy wholeloopName2loopBoundMap;
wholeloopName2perfectOrNotMapTy wholeloopName2perfectOrNotMap;

using namespace llvm;

ExtractLoopInfo::ExtractLoopInfo() : LoopPass(ID) {
	errs() << "DEBUG-INFO: [static-analysis_loop-info-extraction] Extracting loop information\n";
	DEBUG(dbgs() << "Initialize ExtractLoopInfo pass\n");
	loopID = 0;
	depth = 0;
	depth_record = 0;
	countNumLoopInaFunc = 0;
	alreadyCheck = 0;
	LoopsMetadataNode.clear();
	DEBUG(dbgs() << "Initialize ExtractLoopInfo pass\n");
	initializeExtractLoopInfoPass(*PassRegistry::getPassRegistry());
}

bool ExtractLoopInfo::doInitialization(Loop *lp, LPPassManager &LPM) {
	Function* fn = lp->getHeader()->getParent();
	std::vector<std::string>::iterator fn_it;
	fn_it = std::find(kernel_names.begin(), kernel_names.end(), fn->getName());
	if (fn_it == kernel_names.end()) {
		// We do not care uninteresting functions
		return false;
	}

	/*
	/// Initialize FuncName2NumLoopsMap with (functionName, 0)
	for (unsigned i = 0; i < kernel_names.size(); i++) {
		FuncName2NumLoopsMap.insert(std::make_pair(kernel_names.at(i), 0));
	}
	*/

	Module* M = lp->getHeader()->getParent()->getParent();
	if (!M->getNamedMetadata(mdKindName)) {
		NMD = M->getOrInsertNamedMetadata(mdKindName);
	}

	/// Check function that contains this loop
	//errs() << "Function in ExtractLoopInfo::doInitialization: ";
#ifndef ENABLE_INSTDISTRIBUTION
	//errs() << lp->getHeader()->getParent()->getName() << "\n";
#endif // End of ENABLE_INSTDISTRIBUTION

	/// Clear all contents in maps
	//BB2LoopLevelmap.clear();
	//LSID2_BB2LoopLevel_map.clear();
	//LoopID2LSInfomap.clear();
	//Func2LoopInfomap.clear();

	/// Check Load/Store LineID and LSID
	if (alreadyCheck == 0) {
		unsigned LineID = 0;
		unsigned LoadStoreID = 0;
		Module::iterator FI, FE;

		for (FI = M->begin(), FE = M->end(); FI != FE; ++FI) {
			if (FI->isDeclaration()) {
				continue;
			}
			fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
			if (fn_it == kernel_names.end()) {
				// Ignore uninteresting functions
				continue;
			}
			funcName2loopNumMap.insert(std::make_pair(FI->getName(), 0));
			// Check Load/Store Instruction's ID
			for (inst_iterator I = inst_begin(FI), E = inst_end(FI); I != E; ++I){
				Instruction* inst = &*I;
				if (isa<LoadInst>(inst) || isa<StoreInst>(inst)) {
					LoadStoreIDMap.insert(std::make_pair(inst, std::make_pair(LineID, LoadStoreID)));
					LineID2LSMap.insert(std::make_pair(LineID, inst));
					LoadStoreID++;
				}

				LineID++;
			}
		}
		/// Check whether the load/store ID is correct:
		if (M->getNamedMetadata(recordLoadStoreIDMDNodeName)) {
			recordLSMDNode = M->getNamedMetadata(recordLoadStoreIDMDNodeName);
			unsigned size_LSMDNode = recordLSMDNode->getNumOperands();
			assert("IR has been changed, Load/Store LineID Mismatch!\n" && size_LSMDNode == LoadStoreIDMap.size());

			for (unsigned i = 0; i < size_LSMDNode; i++) {
				MDNode* tempMDNode = recordLSMDNode->getOperand(i);
				ConstantInt* line_id_ptr = dyn_cast<ConstantInt>(tempMDNode->getOperand(1));
				ConstantInt* ls_id_ptr = dyn_cast<ConstantInt>(tempMDNode->getOperand(2));
				uint64_t line_id = line_id_ptr->getZExtValue();
				uint64_t ls_id = ls_id_ptr->getZExtValue();
				Instruction* Inst = LineID2LSMap.at(line_id);
				if (ls_id != LoadStoreIDMap.at(Inst).second) {
					assert("IR has been changed, Load/Store ID Mismatch!\n" && size_LSMDNode == LoadStoreIDMap.size());
				}
			}
		}
		// clear LineID2LSMap
		LineID2LSMap.clear();
		/// End checking

		/// Get the basic block ids of the whole program. This part is used to trace the header ID of all loops
		uint64_t bb_counter = 0;
		Function::iterator BI, BE;
		for (FI = M->begin(), FE = M->end(); FI != FE; ++FI) {
			if (FI->isDeclaration()) {
				continue;
			}
			fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
			if (fn_it == kernel_names.end()) {
				// Ignore uninteresting functions
				continue;
			}
			for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
				BB2BBidMap.insert(std::make_pair(BI, bb_counter));
				bb_counter++;
			}
		}

		/// Get the number of functions have loops in the module and 
		/// the number of loops in a function here

		if (M->getNamedMetadata(loopNumberMDNodeName)) {
			loopNumNMD = M->getNamedMetadata(loopNumberMDNodeName);
			unsigned size_loopNumNMD = loopNumNMD->getNumOperands();
			for (unsigned i = 0; i < size_loopNumNMD; i++) {
				MDNode* tempMDNode = loopNumNMD->getOperand(i);
				MDString* func_name_mdstr = dyn_cast<MDString>(tempMDNode->getOperand(0));
				std::string func_name_str = func_name_mdstr->getString();
				funcName2loopNumMapTy::iterator it;
				it = funcName2loopNumMap.find(func_name_str);
				if (it != funcName2loopNumMap.end()) {
					it->second++;
				}
			}
		}

		// Do not need to do the contents in this if statement again
		alreadyCheck = 1; 
	}

	/// We have changed IR
	return true;
}

void ExtractLoopInfo::getAnalysisUsage(AnalysisUsage &AU) const{

	AU.addRequired<LoopInfo>();
	AU.addPreserved<LoopInfo>();
	AU.addRequiredID(LoopSimplifyID);
	AU.addPreservedID(LoopSimplifyID);
	AU.addRequiredID(LCSSAID);
	AU.addPreservedID(LCSSAID);
	AU.addRequired<ScalarEvolution>();
	AU.addPreserved<ScalarEvolution>();
	// FIXME: Loop unroll requires LCSSA. And LCSSA requires dom info.
	// If loop unroll does not preserve dom info then LCSSA pass on next
	// loop will receive invalid dom info.
	// For now, recreate dom info, if loop is unrolled.
	AU.addPreserved<DominatorTreeWrapperPass>();

	//AU.setPreservesCFG();
}

bool ExtractLoopInfo::runOnLoop(Loop* lp, LPPassManager &LPM) {
	DEBUG(dbgs() << "\n\nBegin ExtractLoopInfo Pass :\n");
	LoopInfo* LI = &getAnalysis<LoopInfo>();
	ScalarEvolution *SE = &getAnalysis<ScalarEvolution>();
	//lsIDpass = &getAnalysis<QueryLoadStoreID>();
	//loopNumberpass = &getAnalysis<LoopNumber>();

	//DominatorTree& DT = getAnalysis<DominatorTreeWrapperPass>().getDomTree();

	/// Check whether loop is in a simplify form
	if (!lp->isLoopSimplifyForm()) {
		assert(false && "Error: Loop is not in a Simplify Form!\n");
	}

	Function* func = lp->getHeader()->getParent();
	std::string funcName = func->getName();
	std::vector<std::string>::iterator fn_it;
	fn_it = std::find(kernel_names.begin(), kernel_names.end(), func->getName());
	if (fn_it == kernel_names.end()) {
		// We do not consider uninteresting functions
		return false;
	}

	LLVMContext &Context = func->getContext();
	
	std::vector<Function*>::iterator it;
	it = std::find(exploredFunc.begin(), exploredFunc.end(), func);
	if (it == exploredFunc.end()) {
		// We only print the CFG of a function unexplored
		if (show_cfg_detailed == true) {
			func->viewCFG();
		}
		
		if (show_cfg_only == true) {
			func->viewCFGOnly();
		}
		// New function, we need to add function name to LoopsMetadataNode as 
		// the first operand and calculate the number of loops inside it again
		MDString* func_name = MDString::get(Context, funcName);
		LoopsMetadataNode.push_back(func_name);
		countNumLoopInaFunc = 0;

		exploredFunc.push_back(func);
	}

	unsigned NumLoopInaFunc = funcName2loopNumMap.find(funcName)->second;

	/// We only record the top-level loop.
	/// Approach: Insert a global named metadata node and attach metadata node to it
	/// Format:
	/// "Loop Info of this kernel", MetadataNode1, MetadataNode2, ...
	/// MetadataNode1: Function Name1, Loop1MetadataNode, Loop1MetadataNode, ...
	/// Loop1MetadataNode: "loop1"+ "_" + (std::string) depth, depth, loop bound (if applicable) and loop header BB id
	/// Loop2MetadataNode: "loop2"+ "_" + (std::string) depth, depth, loop bound (if applicable) and loop header BB id
	///
	/// MetadataNode2: Function Name2, Loop1MetadataNode, Loop1MetadataNode, ...
	/// ...
	depth = lp->getLoopDepth();
	depth_record++;

	// Add LoopiMetadataNode (i = 1,2,...,n) to LoopsMetadataNode. The first entry is the 
	// innermost loop
	//std::string loop_name = "loop" + std::to_string(countNumLoopInaFunc+1) + "_" + std::to_string(depth);
	//std::string whole_lp_name = "loop" + std::to_string(countNumLoopInaFunc + 1);
	std::string loop_name = funcName + "_loop" + "-" + std::to_string(NumLoopInaFunc - countNumLoopInaFunc - 1) + "_" + std::to_string(depth);
	std::string whole_lp_name = funcName + "_loop" + "-" + std::to_string(NumLoopInaFunc - countNumLoopInaFunc - 1);
	LoopsMetadataNode.push_back(MDString::get(Context, loop_name));
	LoopsMetadataNode.push_back(ConstantInt::get(Type::getInt32Ty(Context), depth));
	// Detect loop bound of this loop level if available. Otherwise, we set 0 at his field
	// FIXME: Current implementation, we just assign 0 to this field
	LoopsMetadataNode.push_back(ConstantInt::get(Type::getInt32Ty(Context), 0));
	// Add loop header id in this named_metadatanode
	BasicBlock* headerBB = lp->getHeader();
	//errs() << "Header BB name: " << headerBB->getName() << "\n";
	//errs() << "Exit BB name: " << lp->getExitBlock()->getName() << " of loop: " << loop_name << "\n";
	LoopsMetadataNode.push_back(ConstantInt::get(Type::getInt32Ty(Context), BB2BBidMap.at(headerBB)));

#ifndef ENABLE_INSTDISTRIBUTION
	DEBUG(dbgs() << "Anylyzing Loop " << whole_lp_name << " and its loop depth is: " << depth << "\n");
#endif // End of ENABLE_INSTDISTRIBUTION

	/// Try to get the loop bound of this loop
	BasicBlock* latchBB = lp->getLoopLatch();
	unsigned tripCount = 0;
	if (latchBB) {
		tripCount = SE->getSmallConstantTripCount(lp, latchBB);
		wholeloopName2loopBoundMap.insert(std::make_pair(loop_name, tripCount));
		//errs() << "Latch BB = " << latchBB->getName() << "\n";
		//errs() << "loop bound = " << tripCount << "\n";
	}

	if (isPerfectNest(lp, LI) == false) {
		wholeloopName2perfectOrNotMap.insert(std::make_pair(loop_name, false));
		//errs() << "Non perfect loop name (whole): " << loop_name << "\n";
	}
	else {
		wholeloopName2perfectOrNotMap.insert(std::make_pair(loop_name, true));
		//errs() << "Perfect loop name (whole): " << loop_name << "\n";
	}

	//BasicBlock* bb_exit = lp->getExitBlock();
	//BasicBlock* bb_exiting = lp->getExitingBlock();
	//BasicBlock* bb_latch = lp->getLoopLatch();
	//errs() << "bb_exit = " << bb_exit->getName() << "\n";
	//errs() << "bb_exiting = " << bb_exiting->getName() << "\n";
	//errs() << "bb_latch = " << bb_latch->getName() << "\n";
	
	//exitBBFuncnamePair2lpNameLevelPairMap
	assert(lp->getExitingBlock() != NULL && "Can not find exiting basic block\n");
	bbFuncNamePairTy bbFuncPair = std::make_pair(lp->getExitingBlock()->getName(), funcName);
	lpNameLevelPairTy lpNameLevelPair = std::make_pair(whole_lp_name, depth);
	std::pair<std::string, std::string> lpNameLevelStrPair = std::make_pair(whole_lp_name, std::to_string(depth));
	exitBBFuncnamePair2lpNameLevelPairMap.insert(std::make_pair(bbFuncPair, lpNameLevelPair));
	lpNameLevelPair2exitingBBnameMap.insert(std::make_pair(lpNameLevelStrPair, bbFuncPair.first));

	if (depth == 1) { // The top-level loop
		// Count loops inside the function
		countNumLoopInaFunc++;

		/// This loop level is the top level, we need to explore all its load/store instructions
		/// and construct a LSID2_BB2LoopLevel_map.
		std::vector<BasicBlock*> BBvec_temp;
		BBvec_temp = lp->getBlocks();
		for (unsigned i = 0; i < BBvec_temp.size(); i++) {
			BasicBlock* BB = BBvec_temp.at(i);
			unsigned loop_level = LI->getLoopFor(BB)->getLoopDepth();
			//BB2LoopLevelmap.insert(std::make_pair(BB, loop_level));
			BB2loopNameMap.insert(std::make_pair(BB, whole_lp_name));

			// Explore all load/store instructions inside Basic Blocks of a loop
#ifndef ENABLE_INSTDISTRIBUTION
			//errs() << "Explore all load/store instructions inside Basic Blocks of a loop \n";
#endif // End of ENABLE_INSTDISTRIBUTION
			BasicBlock::iterator it, ie;
			for (it = BB->begin(), ie = BB->end(); it != ie; ++it) {
				if (isa<LoadInst>(it) || isa<StoreInst>(it)) {
					Instruction* inst_temp = &*it;
					unsigned lsid = LoadStoreIDMap.at(inst_temp).second;
					LSID2_BB2LoopLevel_map.insert(std::make_pair(lsid, std::make_pair(BB, loop_level)));
#ifndef ENABLE_INSTDISTRIBUTION
					//errs() << "Instruction: " << *it << "\t Belonged BB = " << BB->getName() << "\t Loop_level = " << loop_level << "\n";
#endif // End of ENABLE_INSTDISTRIBUTION
				}
			}
		}  // Now we have explored all load/store instructions inside a loop (whole loop) of a function

		/// Insert the Load/Store information (LSID2_BB2LoopLevel_map) into LoopID2LSInfomap
		//std::string lp_name = "loop" + std::to_string(countNumLoopInaFunc);
		std::string function_name = func->getName();
		std::string lp_name = function_name + "_loop" + "-" + std::to_string(NumLoopInaFunc - countNumLoopInaFunc);
		LoopID2LSInfomap.insert(std::make_pair(lp_name, LSID2_BB2LoopLevel_map));
		//LoopID2LSInfomap.insert(std::make_pair(whole_lp_name, LSID2_BB2LoopLevel_map));
		/// Clear LSID2_BB2LoopLevel_map in order to store new load/store info of other loops in the same function
		LSID2_BB2LoopLevel_map.clear();
		
		/// Store LoopID2LSInfomap into a Func2LoopInfomap
		/// FIXME: Need a pass to calculate number of loops inside a function and number of functions inside a module,
		/// then we can use
		///		if (countNumLoopInaFunc == NumLoopInaFunc) {
		///			Func2LoopInfomap.insert(...);
		///		}

		
		if (countNumLoopInaFunc == NumLoopInaFunc) {
			Func2LoopInfomap.insert(std::make_pair(function_name, LoopID2LSInfomap));
			LoopID2LSInfomap.clear();

			// Add Metadata into the global NamedMetadata Node
			MDNode* MD = MDNode::getWhenValsUnresolved(Context, ArrayRef<Value*>(LoopsMetadataNode), false);
			NMD->addOperand(MD);

			// We finish record one (or one nested) loop, need to clear contents in LoopsMetadataNode
			LoopsMetadataNode.clear();
		}

		// Add Metadata into the global NamedMetadata Node
		//MDNode* MD = MDNode::getWhenValsUnresolved(Context, ArrayRef<Value*>(LoopsMetadataNode), false);
		//NMD->addOperand(MD);

		// We finish record one (or one nested) loop, need to clear contents in LoopsMetadataNode
		//LoopsMetadataNode.clear();

		// Calculate number of arithmatic operations for loops inside functions except for the "main" function
		// The map used to trace binary operators: (functionName, loopName) -> std::vector<BB id> > fnlpNamePair2BinaryOpBBidMap.
		// After we get the basic block frequency, we can get the total number of binary operations executed by
		// multiplication.
		std::vector<uint64_t> BBidVec;
		if (func->getName() != "main") {
			std::vector<BasicBlock*> BBvec = lp->getBlocks();
			BasicBlock::iterator IT, IE;
			unsigned size_BBvec = BBvec.size();
			unsigned arithmetic_counter = 0;
			for (unsigned i = 0; i < size_BBvec; i++) {
				BasicBlock* BBtemp = BBvec.at(i);
				for (IT = BBtemp->begin(), IE = BBtemp->end(); IT != IE; ++IT) {
					// Find arithmetic instructions except for induction variables
					// Arithmetic instructions are called "Binary Operations" in LLVM

					// Should we include "Bitwise Binary Operations" as parts of 
					// arithmetic instructions?
					if (isa<BinaryOperator>(IT)) {
						// Ignore binary operations for induction variables
						//errs() << "BB name = " << BBtemp->getName() << ";  IT->getName() = " << IT->getName() << "\n";
						if (IT->getName().find("indvars") == std::string::npos) {
							//errs() << "BB name = " << BBtemp->getName() << ";  IT = " << *IT << "\n";
							BBidVec.push_back(BB2BBidMap.at(BBtemp));
						}
					}
				}
			}

			// So far, we record all binary operators of a loop. We need to store it into fnlpNamePair2BinaryOpBBidMap
			std::pair<std::string, std::string> fnlpName = std::make_pair(function_name, lp_name);
			//std::pair<std::string, std::string> fnlpName = std::make_pair(func->getName(), whole_lp_name);
			fnlpNamePair2BinaryOpBBidMap.insert(std::make_pair(fnlpName, BBidVec));
			BBidVec.clear();

			#ifdef BUILD_DDDG_H
			for (unsigned i = 0; i < size_BBvec; i++) {
				BasicBlock* BBtemp = BBvec.at(i);
				unsigned lplevel = LI->getLoopFor(BBtemp)->getLoopDepth();
				std::string lplevelStr = std::to_string(lplevel);
				std::string bbName = BBtemp->getName();
				bbFuncNamePairTy bbfnName = std::make_pair(bbName, function_name);
				lpNameLevelPairTy lpNameLevel = std::make_pair(lp_name, lplevel);
				
				if (LI->isLoopHeader(BBtemp)) {
					lpNameLevelPair2headBBnameMap.insert(std::make_pair(std::make_pair(lp_name, lplevelStr), bbName));
					headerBBFuncnamePair2lpNameLevelPairMap.insert(std::make_pair(bbfnName, lpNameLevel));
				}

				bbFuncNamePair2lpNameLevelPairMap.insert(std::make_pair(bbfnName, lpNameLevel));
			}

			LpName2numLevelMap.insert(std::make_pair(lp_name, depth_record));
			depth_record = 0;
			#endif // End of BUILD_DDDG_H
		}
	}

	std::vector<BasicBlock*> BBvec;
	BBvec = lp->getBlocks();
#ifndef ENABLE_INSTDISTRIBUTION
	/*
	errs() << "Check basic blocks inside the " << depth << "-level loop of Loop " << countNumLoopInaFunc << "\n";
	for (unsigned i = 0; i < BBvec.size(); i++) {
		errs() << "BasicBlock Name: " << BBvec[i]->getName() << "\n";
		errs() << "\t Loop level: " << LI->getLoopFor(BBvec[i])->getLoopDepth() << "\n";
	}
	errs() << "\n";
	*/
#endif // End of ENABLE_INSTDISTRIBUTION
	
	//headerBB = lp->getHeader();
	//BasicBlock* latchBB = lp->getLoopLatch();
	//BasicBlock* exitingBB = lp->getExitingBlock();
	//BasicBlock* exitBB = lp->getExitBlock();
	//errs() << "headerBB = " << headerBB->getName() << "\n";
	//errs() << "latchBB = " << latchBB->getName() << "\n";
	//errs() << "exitingBB = " << exitingBB->getName() << "\n";
	//errs() << "exitBB = " << exitBB->getName() << "\n\n";
	
	//}

	DEBUG(dbgs() << "End ExtractLoopInfo Pass\n");
	return true;
}

bool ExtractLoopInfo::isPerfectNest(Loop *L, LoopInfo* li) {

	//check to see if we're the innermost nest
	if (L->getBlocks().size() == 1) {
		//body = *L->block_begin();
		std::vector<Loop *>::iterator found = std::find(exploredLoop.begin(), exploredLoop.end(), L);
		if (found == exploredLoop.end()) {
			exploredLoop.push_back(L);
		}
		else {
			exploredLoop.push_back(L->getParentLoop());
		}
		return true;
	}
	else {
		//do we have a single subloop?
		if (L->getSubLoops().size() != 1){
			exploredLoop.push_back(L);
			return false;
		}
		//make sure all our non-nested loop blocks are innocuous

		std::vector<Loop *>::iterator found = std::find(exploredLoop.begin(), exploredLoop.end(), L);
		if (found != exploredLoop.end()) {
			// We have checked this recursive loop level already, no need to check again. The upper loop
			// level is a perfect loop.
			exploredLoop.push_back(L->getParentLoop());
			return true;
		}

		for (Loop::block_iterator b = L->block_begin(), e = L->block_end(); b
			!= e; b++)
		{
			BasicBlock *block = *b;
			if (li->getLoopFor(block) == L)
			{
				if (!hasNoMemoryOps(block)) {
					exploredLoop.push_back(L);
					return false;
				}
			}
		}

		//recursively check subloops
		return isPerfectNest(*L->begin(), li);
	}
}

bool ExtractLoopInfo::hasNoMemoryOps(BasicBlock *b) {
	//errs() << "Basic Block name: " << b->getName() << "\n";
	for (BasicBlock::iterator I = b->begin(), E = b->end(); I != E; I++)
	{
		switch (I->getOpcode())
		{
		//case Instruction::Free:
		case Instruction::FAdd:
		case Instruction::FSub:
		case Instruction::FMul:
		case Instruction::FDiv:
		case Instruction::FCmp:
			//std::cerr << *I << "\n";
			return false;
		}
	}
	return true;
}

LoopID2LSInfoMapTy ExtractLoopInfo::getFunc2LoopInfomap(std::string func_name) const {
	return Func2LoopInfomap.at(func_name);
}

char ExtractLoopInfo::ID = 0;
/*
INITIALIZE_PASS(ExtractLoopInfo, "extractLoopInfo",
	"This pass is used to extract basic loop info eg loop bound & loop levels and insert \
					 a function to trace loop iteration at runtime",
					 false,
					 true
					 ) */

INITIALIZE_PASS_BEGIN(ExtractLoopInfo, "extractLoopInfo",
					  "This pass is used to extract basic loop info eg loop bound & loop levels and insert \
					  a function to trace loop iteration at runtime",
					  false,
					  true)
INITIALIZE_PASS_DEPENDENCY(LoopNumber)
INITIALIZE_PASS_DEPENDENCY(LoopInfo)
INITIALIZE_PASS_DEPENDENCY(LoopSimplify)
INITIALIZE_PASS_DEPENDENCY(LCSSA)
INITIALIZE_PASS_DEPENDENCY(ScalarEvolution)
INITIALIZE_PASS_DEPENDENCY(DominatorTreeWrapperPass)
INITIALIZE_PASS_END(ExtractLoopInfo, "extractLoopInfo",
					"This pass is used to extract basic loop info eg loop bound & loop levels and insert \
					a function to trace loop iteration at runtime",
					false,
					true)


Pass* llvm::createExtractLoopInfoPass() {
	return new ExtractLoopInfo();
}