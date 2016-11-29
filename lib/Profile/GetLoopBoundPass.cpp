#include "profile_h/GetLoopBoundPass.h"

#define DEBUG_TYPE "get-loop-bound"

namespace llvm {
	funcName2loopsMapTy funcName2loopsMap;
	func2loopInfoMapTy func2loopInfoMap;
	funcName2loopBoundsTy funcName2loopBounds;
	BBids2loopNameMapTy BBids2loopNameMap;
	funcName2loopInfoMapTy funcName2loopInfoMap;

	//typedef std::map<uint64_t, uint64_t> BBidFreqMap;
	BBidFreqMap BBidFreq;
	BranchidFreqMap BranchIDFreq;

	BBheaderID2loopBoundMapTy BBheaderID2loopBoundMap;
}

using namespace llvm;

static const char *LoopInfoMDNodeName = "Loop_Info_of_this_kernel";
static const char *loopNumberMDNodeName = "Loop_Number_of_kernel";

GetLoopBound::GetLoopBound() : ModulePass(ID) {
	errs() << "Initialize GetLoopBound pass\n";
}

void GetLoopBound::getAnalysisUsage(AnalysisUsage &AU) const {
	AU.addRequired<QueryBasicBlockID>();
	AU.addRequired<QueryLoadStoreID>();
	AU.addPreserved<QueryBasicBlockID>();
	AU.addPreserved<QueryLoadStoreID>();

	AU.addRequired<LoopInfo>();
	AU.addPreserved<LoopInfo>();
	AU.setPreservesCFG();
}

bool GetLoopBound::doInitialization(Module &M) {

	// Obtain different Types
	Int64Type = IntegerType::getInt64Ty(M.getContext());
	VoidType = Type::getVoidTy(M.getContext());

	LastBB = nullptr;

	/// Clone Module
	NewM = CloneModule(&M, VMap);

	/// recordBB function arguments' type: 
	/// Return type:	VoidType
	/// First argument:  Int64Type    -- Used to record the ID of a basic block
	/// Second argument: Int64Type  -- Used to trace whether there are call instructions from a basic block. If yes,
	///								   we do not increase the frequency of the basic block; otherwise, we increase it.
	BBFreqIncFn = cast<Function>(NewM->getOrInsertFunction("recordBB_loopBound", VoidType, Int64Type, Int64Type, Int64Type, nullptr));

	return false;
}

bool GetLoopBound::runOnModule(Module &M){
	errs() << "\n\nBegin GetLoopBound Pass :\n";

	bbIDpass = &getAnalysis<QueryBasicBlockID>();
	lsIDpass = &getAnalysis<QueryLoadStoreID>();

#ifdef ENABLE_INSTDISTRIBUTION

#else
	/// Record the IDs of basic blocks that are loop headers
	getLoopBBID2HeaderAndloopDepthPairMap(M);

	/// Record the IDs of header Basic Block of loops
	getHeaderID_list_vec(M);
#endif // End of ENABLE_INSTDISTRIBUTION

	Module::iterator FI, FE;
	Function::iterator BI, BE;
	BasicBlock::iterator I, E;

	DL = NewM->getDataLayout();

	/// Instrument Function for recording basic block frequency:
	std::vector<std::string>::iterator fn_it;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		if (FI->isDeclaration()) {
			continue;
		}
		fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (fn_it == kernel_names.end()) {
			// Ignore uninteresting functions
			continue;
		}

		instrumentFunction(FI);
	}

	/// Write the new module bitcode into a file
	write_Modified_bitcode_to_File();

	/// Before run on the modified bitcode, initialize BBidFreq first, otherwise
	/// it will just ignore the basic block that will not be executed and cause 
	/// program crashed.
	uint64_t sizeOfBB = bbIDpass->getNumOfBB();
	for (unsigned i = 0; i < sizeOfBB; i++) {
		BBidFreq.insert(std::make_pair(i, 0));
	}

	/// Integrate JIT profiling engine and run the embedded profiler
	EmbeddedProfilerEngine_lb P(*NewM, BBFreqIncFn);
	//EmbeddedProfilerEngine P(M, BBFreqIncFn, RecordLoadFn, RecordStoreFn);
	//EmbeddedProfilerEngine P(M, RecordLoadFn, RecordStoreFn);
	P.runOnProfiler();

#ifdef ENABLE_INSTDISTRIBUTION

#else
	/// Calculate loop bounds for all loops in the kernel
	/// 1. Get the names of functions that have loops
	extract_loopInfo_inFunc(M);

	/// 2. Calculate the loop bounds of loops
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {

		if (FI->isDeclaration()) {
			errs() << "Function name: " << FI->getName() << "\n";
			continue;
		}
		fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (fn_it == kernel_names.end()) {
			// We do not consider uninteresting functions
			continue;
		}

		funcName2loopsMapTy::iterator it;
		it = funcName2loopsMap.find(FI->getName());
		if (it != funcName2loopsMap.end()) {
			calculateLoopBound(FI);
			/*
			/// Record the IDs of basic blocks that are loop headers
			LoopInfo* LI = &getAnalysis<LoopInfo>(*FI);
			for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
			if (LI->isLoopHeader(BI)) {
			// The basic block is a loop header
			uint64_t bbid = bbIDpass->getBBid(BI);
			unsigned lp_depth = LI->getLoopDepth(BI);
			loopHeaderBBid2loopDepthMap.insert(std::make_pair(bbid, lp_depth));
			}
			} /// Now we have already store loopHeaderBBid2loopDepthMap
			*/
		} /// We have finished calculating loop bounds of this function
	}
#endif // End of ENABLE_INSTDISTRIBUTION

	// Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		assert(false && "Module is broken!\n");
	}

	// Clear BBidFreq, BranchIDFreq
#ifndef ENABLE_INSTDISTRIBUTION
	BBidFreq.clear();
	BranchIDFreq.clear();
#endif // End of ENABLE_INSTDISTRIBUTION

	delete NewM;

	errs() << "\nEnd GetLoopBound pass\n\n";
	return false;
}

/// FIXME: This code will fail if more than one loop inside a kernel. 
/// Bug fixed
void GetLoopBound::getLoopBBID2HeaderAndloopDepthPairMap(Module &M) {
	/// Get function name we focus and store into functionNameVec
	NamedMDNode* NMD_funcName;
	assert(M.getNamedMetadata(loopNumberMDNodeName) && "No metadata node named Loop_Number_of_kernel, error!\n");
	NMD_funcName = M.getNamedMetadata(loopNumberMDNodeName);
	unsigned size_NMDfuncName = NMD_funcName->getNumOperands();
	std::vector<std::string>::iterator it_fnName;
	for (unsigned i = 0; i < size_NMDfuncName; i++) {
		MDNode* node = NMD_funcName->getOperand(i);
		unsigned size_node = node->getNumOperands();
		MDString* func_name = dyn_cast<MDString>(node->getOperand(0));
		it_fnName = std::find(functionNameVec.begin(), functionNameVec.end(), func_name->getString());
		if (it_fnName == functionNameVec.end()) {
			functionNameVec.push_back(func_name->getString());
		}
	} // We get the function names we care about.

	/// Record the IDs of basic blocks that are loop headers
	Module::iterator FI, FE;
	Function::iterator BI, BE;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		std::vector<std::string>::iterator it;
		it = std::find(functionNameVec.begin(), functionNameVec.end(), FI->getName());
		// We only care about the specific functions
		if (it != functionNameVec.end()) {
			LoopInfo* LI = &getAnalysis<LoopInfo>(*FI);
			for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
				errs() << "BB name = " << BI->getName() << "\n";
				// We only care about the basic blocks in loops
				if (LI->getLoopFor(BI) != NULL) {
					// The basic block is a loop header
					uint64_t bbid = bbIDpass->getBBid(BI);
					unsigned lp_depth = LI->getLoopDepth(BI);
					if (LI->isLoopHeader(BI)) {
						loopBBid2HeaderAndloopDepthPairMap.insert(std::make_pair(bbid, std::make_pair(bbid, lp_depth)));
					}
					else {
						BasicBlock* headerBB = LI->getLoopFor(BI)->getHeader();
						uint64_t headerBBid = bbIDpass->getBBid(headerBB);
						loopBBid2HeaderAndloopDepthPairMap.insert(std::make_pair(bbid, std::make_pair(headerBBid, lp_depth)));
					}
				}
			}
		}
	} // We have got loopHeaderBBid2loopDepthMap.
}

bool GetLoopBound::isLoopHeader(uint64_t bbid) {
	BasicBlock* BB = bbIDpass->getBBbyID(bbid);
	Function* func = BB->getParent();
	LoopInfo* LI = &getAnalysis<LoopInfo>(*func);
	return LI->isLoopHeader(BB);
}

void GetLoopBound::getHeaderID_list_vec(Module &M){
	/// headerID_list_vec contains a vector of headerID_list
	/// headerID_list:  The header BB id of the innermost loop --> The header BB id of the outermost loop

	/// Get the loop header list for all loops in the module
	NamedMDNode* NMD_loopInfo;
	assert(M.getNamedMetadata(LoopInfoMDNodeName) && "CodeInstrumentPass No metadata node named Loop_Info_of_this_kernel, error!\n");
	NMD_loopInfo = M.getNamedMetadata(LoopInfoMDNodeName);
	unsigned size_loopInfo = NMD_loopInfo->getNumOperands();
	for (unsigned i = 0; i < size_loopInfo; i++) {
		MDNode* node_loopInfo = NMD_loopInfo->getOperand(i);
		unsigned size_NumOp = node_loopInfo->getNumOperands();
		MDString* func_name;
		//MDString* subloop_name;
		ConstantInt* subloop_depth;
		//ConstantInt* subloop_bound;
		// Add loop header bb id previous
		ConstantInt* subloop_headerID;

		std::list<uint64_t> headerID_list;

		func_name = dyn_cast<MDString>(node_loopInfo->getOperand(0));
		unsigned loop_counter = 0;
		for (unsigned j = 1; j < size_NumOp; j += 4) {
			subloop_depth = dyn_cast<ConstantInt>(node_loopInfo->getOperand(j + 1));

			// Add loop header bb id previous
			subloop_headerID = dyn_cast<ConstantInt>(node_loopInfo->getOperand(j + 3));

			unsigned depth = static_cast<unsigned> (subloop_depth->getZExtValue());

			// Add loop header bb id previous
			unsigned headerID = static_cast<uint64_t> (subloop_headerID->getZExtValue());

			// The list is starting from the header ID of the innermost loop
			headerID_list.push_back(headerID);

			if (depth == 1) { // The top level of a loop
				headerID_list_vec.push_back(headerID_list);
				headerID_list.clear();
			}
		}
		funcName2headIDlistVecMap.insert(std::make_pair(func_name->getString(), headerID_list_vec));
	}
}

/// FIXME: This code has a bug
/// Bugs fixed, but we need to be aware that we count loop N is the first loop in the kernel and
/// loop 1 is the last loop in the kernel, maybe we can change this in 'Loop_Info_of_this_kernel' metadata
void GetLoopBound::extract_loopInfo_inFunc(Module &M) {
	NamedMDNode* NMD_loopNumber;
	assert(M.getNamedMetadata(loopNumberMDNodeName) && "No metadata node named Loop_Number_of_kernel, error!\n");
	NMD_loopNumber = M.getNamedMetadata(loopNumberMDNodeName);
	unsigned size_loopNumber = NMD_loopNumber->getNumOperands();
	//funcName2loopsMapTy::iterator it_fnName2lp;
	std::string funcName;
	std::vector<std::string > loopsName;
	uint64_t lp_Num_counter = 0;
	uint64_t lp_Num_of_a_func;
	for (unsigned i = 0; i < size_loopNumber; i++) {
		MDNode* node = NMD_loopNumber->getOperand(i);
		unsigned size_node = node->getNumOperands();
		MDString* func_name;
		for (unsigned j = 0; j < size_node; j++) {
			if (j == 0) {
				func_name = dyn_cast<MDString>(node->getOperand(j));
				funcName = func_name->getString();
			}
			else {
				MDString* mdstr_loop = dyn_cast<MDString>(node->getOperand(j));
				loopsName.push_back(mdstr_loop->getString());
				lp_Num_counter++;
			}
		}
		lp_Num_of_a_func = funcName2loopNumMap.at(funcName);

		if (lp_Num_counter == lp_Num_of_a_func) {
			funcName2loopsMap.insert(std::make_pair(funcName, loopsName));
			loopsName.clear();
		}
	}

	NamedMDNode* NMD_loopInfo;
	assert(M.getNamedMetadata(LoopInfoMDNodeName) && "AnalysisProfilingPass No metadata node named Loop_Info_of_this_kernel, error!\n");
	NMD_loopInfo = M.getNamedMetadata(LoopInfoMDNodeName);
	unsigned size_loopInfo = NMD_loopInfo->getNumOperands();
	for (unsigned i = 0; i < size_loopInfo; i++) {
		MDNode* node_loopInfo = NMD_loopInfo->getOperand(i);
		unsigned size_NumOp = node_loopInfo->getNumOperands();
		MDString* func_name;
		MDString* subloop_name;
		ConstantInt* subloop_depth;
		ConstantInt* subloop_bound;
		// Add loop header bb id previous
		ConstantInt* subloop_headerID;

		func_name = dyn_cast<MDString>(node_loopInfo->getOperand(0));

		//BBids2loopNameMapTy BBids2loopNameMap;

		unsigned loop_counter = 0;
		for (unsigned j = 1; j < size_NumOp; j += 4) {
			subloop_name = dyn_cast<MDString>(node_loopInfo->getOperand(j));
			subloop_depth = dyn_cast<ConstantInt>(node_loopInfo->getOperand(j + 1));
			subloop_bound = dyn_cast<ConstantInt>(node_loopInfo->getOperand(j + 2));
			// Add loop header bb id previous
			subloop_headerID = dyn_cast<ConstantInt>(node_loopInfo->getOperand(j + 3));

			unsigned depth = static_cast<unsigned> (subloop_depth->getZExtValue());
			unsigned bound = static_cast<unsigned> (subloop_bound->getZExtValue());
			// Add loop header bb id previous
			unsigned headerID = static_cast<uint64_t> (subloop_headerID->getZExtValue());

			subloopDepthandLoopBound.insert(std::make_pair(depth, bound));

			if (depth == 1) { // The top level of a loop
				std::string loop_name = funcName2loopsMap.at(func_name->getString()).at(loop_counter);
				loop2subloopMap.insert(std::make_pair(loop_name, subloopDepthandLoopBound));
				subloopDepthandLoopBound.clear();

				/// Get all basic blocks inside a loop and store them into loopName2BBidsMap
				BasicBlock* HeaderTopLevel = bbIDpass->getBBbyID(headerID);
				Function* Func = HeaderTopLevel->getParent();
				LoopInfo* LI = &getAnalysis<LoopInfo>(*Func);
				std::vector<BasicBlock*> BBvec_temp = LI->getLoopFor(HeaderTopLevel)->getBlocks();
				//std::vector<uint64_t> BBidvec_temp;
				unsigned size_vec = BBvec_temp.size();
				for (unsigned k = 0; k < size_vec; k++) {
					BasicBlock* BB_temp = BBvec_temp.at(k);
					uint64_t BBid_temp = bbIDpass->getBBid(BB_temp);
					//BBidvec_temp.push_back(BBid_temp);
					BBids2loopNameMap.insert(std::make_pair(BBid_temp, loop_name));
					//bbID_for_loops_vec.push_back(BBid_temp);
				}
				//loopName2BBidsMap.insert(std::make_pair(loop_name, BBidvec_temp));
				BBvec_temp.clear();
				//BBidvec_temp.clear();

				loop_counter++;
			}
		}

		/// Used to recognize the relationship between basic block ids and loop name
		funcName2loopInfoMap.insert(std::make_pair(func_name->getString(), BBids2loopNameMap));

		func2loopInfoMap.insert(std::make_pair(func_name->getString(), loop2subloopMap));
		loop2subloopMap.clear();
	}

} // End of GetLoopBound::extract_loopInfo_inFunc(Module &M)

void GetLoopBound::calculateLoopBound(Function* func) {
	/// We need to record (HeaderID -> Loop Bound) Map
	errs() << "function: " << func->getName() << "\n";
	LoopInfo* LI = &getAnalysis<LoopInfo>(*func);

	std::map<unsigned, uint64_t> depth2BBfreqMap;
	std::map<unsigned, uint64_t> depth2BBHeaderIDMap;
	std::map<unsigned, unsigned> depth2LoopBoundMap;
	std::map<std::string, std::map<unsigned, uint64_t> > loopName2depthBBfreqMap;
	std::map<std::string, std::map<unsigned, uint64_t> > loopName2depthHeaderIDMap;

	unsigned depth = 0;
	uint64_t bbfreq = 0;
	uint64_t headerID = 0;
	unsigned loop_counter = 1;
	unsigned num_nestedLevel = 0;
	Function::iterator BI, BE;
	for (BI = func->begin(), BE = func->end(); BI != BE; ++BI) {
		errs() << "Basic Block = " << BI->getName() << "\n";
		// Do not consider basic blocks outside loops
		if (!LI->getLoopFor(BI)) {
			continue;
		}

		// Do not consider basic blocks that are not a header of loops
		if (!LI->isLoopHeader(BI)) {
			continue;
		}

		depth = LI->getLoopDepth(BI);
		headerID = bbIDpass->getBBid(BI);
		bbfreq = BBidFreq.at(headerID);
		depth2BBfreqMap.insert(std::make_pair(depth, bbfreq));
		depth2BBHeaderIDMap.insert(std::make_pair(depth, headerID));

		std::string loop_name = "loop" + std::to_string(loop_counter);
		num_nestedLevel = func2loopInfoMap.at(func->getName()).at(loop_name).size();
		if (depth == num_nestedLevel) {
			// We have explore all header basic blocks of the loop, "loop_name", store 
			// its depth and bbfreq information into loopName2depthBBfreqMap
			loopName2depthBBfreqMap.insert(std::make_pair(loop_name, depth2BBfreqMap));
			loopName2depthHeaderIDMap.insert(std::make_pair(loop_name, depth2BBHeaderIDMap));
			depth2BBfreqMap.clear();
			depth2BBHeaderIDMap.clear();
			loop_counter++;
		}

		errs() << "Loop Name =" << loop_name << "  loop header = " << BI->getName() \
			<< "  Depth = " << depth << "  BBfreq = " << bbfreq << "\n";
	}

	/// Calculate loop bounds for all loops in a function
	/// For the top level,    bound1 = header1BB_freq - 1
	/// For the second level, bound2 = header2BB_freq / header1BB_freq
	/// For the third level,  bound3 = header3BB_freq / header2BB_freq
	/// For the i-th level,   bound#i# = header#i#BB_freq / header#i-1#_freq

	// Loop bounds sequence: 
	// The top-level loop -> lower-level loop -> ... -> the innermost-level loop
	unsigned num_loops = func2loopInfoMap.at(func->getName()).size();
	loopName2BoundListTy loopName2BoundList;
	for (unsigned i = 0; i < num_loops; i++) {
		std::string lp_name = "loop" + std::to_string(i + 1);
		num_nestedLevel = func2loopInfoMap.at(func->getName()).at(lp_name).size();
		unsigned bound = 0;
		unsigned higher_bound = 0;
		std::list<unsigned> loopBoundList;
		for (unsigned j = 1; j < num_nestedLevel + 1; j++) {
			if (j == 1) {
				bound = loopName2depthBBfreqMap.at(lp_name).at(j) - 1;
			}
			else {
				bound = loopName2depthBBfreqMap.at(lp_name).at(j) / higher_bound;
			}

			higher_bound = loopName2depthBBfreqMap.at(lp_name).at(j);
			headerID = loopName2depthHeaderIDMap.at(lp_name).at(j);
			// Store the loop bound into func2loopInfoMap
			func2loopInfoMapTy::iterator it = func2loopInfoMap.find(func->getName());
			if (it != func2loopInfoMap.end()) {
				it->second.at(lp_name).find(j)->second = bound;
			}
			loopBoundList.push_back(bound);
			BBheaderID2loopBoundMap.insert(std::make_pair(headerID, bound));
		}
		assert(loopBoundList.size() == num_nestedLevel && "Size of loopBoundList is NOT equal to num_nestedLevel!\n");
		loopName2BoundList.insert(std::make_pair("loop" + std::to_string(i + 1), loopBoundList));
		loopBoundList.clear();
	} /// End of calculating loop bounds for all loops in a function
	funcName2loopBounds.insert(std::make_pair(func->getName(), loopName2BoundList));
	loopName2BoundList.clear();

	/*
	LoopInfo::iterator li, le;
	for (li = LI->begin(), le = LI->end(); li != le; ++li) {
	std::vector<Loop*> Subloops = (*li)->getSubLoops();
	unsigned size_subloops = Subloops.size();
	for (unsigned i = 0; i < size_subloops; i++) {
	Loop* lp = Subloops.at(i);
	std::vector<BasicBlock*> BBvec_temp;
	BBvec_temp = lp->getBlocks();
	for (unsigned j = 0; j < BBvec_temp.size(); j++) {
	errs() << "Basic Block = " << BBvec_temp.at(j)->getName() << "\n";
	}
	}
	}
	*/
} // End of GetLoopBound::calculateLoopBound(Function* func)

void GetLoopBound::instrumentFunction(Function* F) {
	for (Function::iterator BI = F->begin(), BE = F->end(); BI != BE; ++BI) {
		instrumentBasicBlock(BI);
	}
}

void GetLoopBound::instrumentBasicBlock(BasicBlock* BB) {
	bool edge_increase;
	if (BB == &BB->getParent()->getEntryBlock()) {
		edge_increase = false;
		EntryBBid_vec.push_back(bbIDpass->getBBid(BB));
	}
	else {
		edge_increase = true;
	}
	for (BasicBlock::iterator I = BB->begin(), E = BB->end(); I != E; ++I) {
		Instruction* Inst = I;

		//errs() << "Instruction: " << *Inst << "\n";
		if (isa<CallInst>(Inst)) {
			//errs() << "Instruction in if(... CallInst ...): " << *Inst << "\n";
			CallSite CS(Inst);
			Function* callFunc = CS.getCalledFunction();

			BasicBlock::iterator InsertPos = I;
#ifndef ENABLE_INSTDISTRIBUTION
			errs() << "Before insertCounterBefore Instr: " << *I << "\n";
#endif // End of ENABLE_INSTDISTRIBUTION
			insertCounterBefore(++InsertPos, false, edge_increase);

#ifndef ENABLE_INSTDISTRIBUTION
			errs() << "Instruction in BB " << *I << "\n";
#endif // End of ENABLE_INSTDISTRIBUTION
			continue;
		}
	}

	insertCounterBefore(BB->getTerminator(), true, edge_increase);
}

void GetLoopBound::insertCounterBefore(Instruction* Inst, bool increase, bool edgeIncrease) {
	/*
	errs() << "insertCounterBefore Instr: " << *Inst << "\n";

	BasicBlock *OldBB = Inst->getParent();
	errs() << "insertCounterBefore OldBB: " << OldBB->getName() << "\n";
	Value *Args[] = {
	ConstantInt::get(Int64Type, uint64_t(OldBB)),
	ConstantInt::get(Int64Type, uint64_t(increase ? 1 : 0))
	};

	CallInst::Create(BBFreqIncFn, Args, "", Inst);
	*/
#ifndef ENABLE_INSTDISTRIBUTION
	errs() << "insertCounterBefore Instr: " << *Inst << "\n";
#endif // End of ENABLE_INSTDISTRIBUTION
	BasicBlock *OldBB = Inst->getParent();
	unsigned bb_id = bbIDpass->getBBid(OldBB);
	Instruction *MappedInst = cast<Instruction>((Value*)VMap.lookup(Inst));

	// Get the entry instruction of a Basic Block
	Instruction* entryNonPhiInst = OldBB->getFirstNonPHI();
	Instruction* MappedFirstInst = cast<Instruction>((Value*)VMap.lookup(entryNonPhiInst));
#ifndef ENABLE_INSTDISTRIBUTION
	errs() << "insertCounterBefore OldBB: " << OldBB->getName() << "\n";
#endif // End of ENABLE_INSTDISTRIBUTION
	Value *Args[] = {
		ConstantInt::get(Int64Type, uint64_t(bb_id)),
		ConstantInt::get(Int64Type, uint64_t(increase ? 1 : 0)),
		ConstantInt::get(Int64Type, uint64_t(edgeIncrease ? 1 : 0))
	};

	//CallInst::Create(BBFreqIncFn, Args, "", Inst);
	if (increase == false) {
		CallInst::Create(BBFreqIncFn, Args, "", MappedInst);
	}
	else {
		// We move the recordBB CallInst into the beginning of a basic block, because we need to trace the 
		// evaluation iteration index for all load/store instructions within loops.
		CallInst::Create(BBFreqIncFn, Args, "", MappedFirstInst);
	}

}

/// Write the new module bitcode into a file
void GetLoopBound::write_Modified_bitcode_to_File() {

	if (inputFileName == "-") {
		OutputModifiedLBBCFilename = "-";
	}
	else {
		const std::string &IFN = inputFileName;
		int Len = IFN.length();
		// If the source ends in .bc, strip it off
		if (IFN[Len - 3] == '.' && IFN[Len - 2] == 'b' && IFN[Len - 1] == 'c') {
			OutputModifiedLBBCFilename = std::string(IFN.begin(), IFN.end() - 3) + "_Modified_LoopBound.bc";
		}
	}

	std::unique_ptr<tool_output_file> Out;
	std::string ErrorInfo;
	Out.reset(new tool_output_file(OutputModifiedLBBCFilename.c_str(), ErrorInfo, sys::fs::F_None));
	if (!ErrorInfo.empty()) {
		errs() << ErrorInfo << "\n";
		assert(false && "ErrorInfo is not empty in OutputModifiedBCFilename\n");
	}

	WriteBitcodeToFile(NewM, Out->os());
	Out->os().close();
	if (Out->os().has_error()) {
		errs() << "could not write bitcode file: " << OutputModifiedLBBCFilename.c_str() << "\n";
		Out->os().clear_error();
		assert(false && "could not write bitcode file OutputModifiedBCFilename\n");
	}

	Out->keep();
} /// We have successfully written the modified bitcode into a new file

void recordBB_loopBound(uint64_t BBid, uint64_t inc, uint64_t edge_increase) {

	//errs() << "BBid = " << BBid << "\n";
	//BasicBlock* BB = reinterpret_cast<BasicBlock*>(BBPtr);
	uint64_t lastBBid = GlobalContext_lb->getLastBB();

	if (inc) {
		GlobalContext_lb->P->increaseEdgeCounter(lastBBid, BBid, edge_increase);
	}
	GlobalContext_lb->setLastBB(BBid);

}

EmbeddedProfilerEngine_lb::EmbeddedProfilerEngine_lb(Module &M, Function* bbFreqfn) :
	Mod(M), BBFreqFn(bbFreqfn) {
}

void EmbeddedProfilerEngine_lb::runOnProfiler(){
	ProfilerJITSingletonContext_lb JTSC(this);
	//std::unique_ptr<Module> M(generateInstrumentedModule());

	// Before running profiling, we need to clear the vector recordls
	//recordls.clear();

	// TODO: Why I always fail when using CloneModule?
	//Module* M = CloneModule(&Mod, V2VMap);
	//std::unique_ptr<Module> M(newM);
	Module *M = (&Mod);

	// TODO: Insert instrumental code to extract the trace.
	EngineBuilder builder(M);
	ExecutionEngine *EE = builder.setEngineKind(EngineKind::JIT).create();
	// Where is getBBFreqInc()?
	//EE->addGlobalMapping(getBBFreqInc(), reinterpret_cast<void*>(IncreaseBBCounter));
	EE->addGlobalMapping(BBFreqFn, reinterpret_cast<void*>(recordBB_loopBound));

	if (!EE) {
		errs() << ": Failed to construct ExecutionEngine: " << "\n";
	}

	Function *EntryFn = M->getFunction("main");

	// Nothing we can do if we cannot find the entry function of the module.
	if (EntryFn == nullptr)
		return;

	// Here we go. Run the module.
	EE->runStaticConstructorsDestructors(false);

	// Trigger compilation separately so code regions that need to be
	// invalidated will be known.
	(void)EE->getPointerToFunction(EntryFn);

	// TODO: Accept argv from user script?
	// FIXME: Support passing arguments
	std::vector<std::string> Argv;
	Argv.push_back(M->getModuleIdentifier());
	Argv.push_back(inputPath);

	// Run main.
	int Result = EE->runFunctionAsMain(EntryFn, Argv, nullptr);

	if (Result != 0) {
		errs() << "Module return nonzero result during trace generation: "
			<< Result << '\n';
	}

	// Run static destructors.
	EE->runStaticConstructorsDestructors(true);

	VerifyProfile();

	//delete newM;
}

void EmbeddedProfilerEngine_lb::increaseEdgeCounter(uint64_t src_id, uint64_t dest_id, uint64_t increase_edge){
	//++BBFreq[dest];
	++BBidFreq[dest_id];

	// We introduce increase_edge to ignore the condition that src_id = INITIAL_BIGVALUE
	if (increase_edge) {
		//++BranchFreq[dest][src];
		++BranchIDFreq[dest_id][src_id];
	}

} // End of EmbeddedProfilerEngine_lb::increaseEdgeCounter()

/// Functions for accessing Basic Block frequency
uint64_t EmbeddedProfilerEngine_lb::getBBidFreq(uint64_t BBid) const {
	BBidFreqMap::const_iterator I = BBidFreq.find(BBid);
	return (I == BBidFreq.end() ? 0 : I->second);
} // End of EmbeddedProfilerEngine_lb::getBBidFreq


uint64_t EmbeddedProfilerEngine_lb::getBranchFreq(uint64_t Src_id, uint64_t Dest_id) const {
	BranchidFreqMap::const_iterator I = BranchIDFreq.find(Dest_id);
	if (I == BranchIDFreq.end()) {
		return 0;
	}

	const BBidFreqMap& CurFreqs = I->second;
	BBidFreqMap::const_iterator J = CurFreqs.find(Src_id);

	return (J == CurFreqs.end() ? 0 : J->second);
} // End of EmbeddedProfilerEngine_lb::getBranchFreq()

void EmbeddedProfilerEngine_lb::VerifyProfile() {
	typedef BBidFreqMap::const_iterator bbid_iterator;
	typedef BranchidFreqMap::const_iterator brid_iterator;
	typedef BBidFreqMap::const_iterator bb_iterator;

	errs() << "\nVerify the profile data: \n";
	/// Checking the Basic Block Frequency and Edge Frequency
	for (bbid_iterator I = BBidFreq.begin(), E = BBidFreq.end(); I != E; ++I) {
		uint64_t BBid = I->first;
		// We do not count the call edge frequency
		std::vector<uint64_t>::iterator EntryBBid_it;
		EntryBBid_it = std::find(EntryBBid_vec.begin(), EntryBBid_vec.end(), BBid);
		if (EntryBBid_it != EntryBBid_vec.end())
			continue;

		uint64_t Freq = I->second;

		// If a basic block is not executed at all, its frequency will be equal to 0,
		// we just ignore it.
		if (Freq == 0) {
			continue;
		}

		BranchidFreqMap::const_iterator BrFreqsAt = BranchIDFreq.find(BBid);
		assert(BrFreqsAt != BranchIDFreq.end() && "Freq info not found!");
		const BBidFreqMap &BrFreqs = BrFreqsAt->second;
		for (bb_iterator I = BrFreqs.begin(), E = BrFreqs.end(); I != E; ++I)
			Freq -= I->second;

		if (Freq != 0) {
			llvm_unreachable("Bad profile infomration!");
		}
	}
} // End of EmbeddedProfilerEngine_lb::VerifyProfile()

ProfilerJITContext_lb::ProfilerJITContext_lb() : P(nullptr), LastBBid(INITIAL_BIGVALUE) {
	// If we have a native target, initialize it to ensure it is linked in and
	// usable by the JIT.
	InitializeNativeTarget();
	InitializeNativeTargetAsmPrinter();
	InitializeNativeTargetAsmParser();
}

ProfilerJITSingletonContext_lb::ProfilerJITSingletonContext_lb(EmbeddedProfilerEngine_lb *P) {
	GlobalContext_lb->P = P;
}

ProfilerJITSingletonContext_lb::~ProfilerJITSingletonContext_lb() {
	GlobalContext_lb->P = nullptr;
	//TODO: Other clearup?
}

char GetLoopBound::ID = 0;
INITIALIZE_PASS(GetLoopBound, "Analyze loop bounds of kernel at runtime",
								"This pass is used to analyze loop bounds at runtime",
								false,
								false
								)

ModulePass *llvm::createGetLoopBoundPass() {
	return new 	GetLoopBound();
}