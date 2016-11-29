#include "profile_h/CodeInstrumentPass.h"

#define DEBUG_TYPE "code-Instrument"

//#define IGNORE_RECORDBB

namespace llvm {
	//typedef std::map<uint64_t, uint64_t> BBidFreqMap;
	//BBidFreqMap BBidFreq;
	//BranchidFreqMap BranchIDFreq;
#ifndef USE_LIST_FOR_RECORDLS
	std::vector<RecordLoadStore> recordls;
#else
	std::list<RecordLoadStore> recordls;
#endif // End of USE_LIST_FOR_RECORDLS

	TraceType Trace_f;

	/// Loop Basic Block ID --> Vector<Evaluation Iteration Index list>
	//typedef std::map<uint64_t, std::vector<std::list<unsigned> > > lpBBid2EvalIterIndxMapTy;
	/// Evaluation Iteration Index sequence:  The top-level -> lower-level --> the innermost-level loop
	lpBBid2EvalIterIndxMapTy lpBBid2EvalIterIndxMap;

	//typedef std::vector<uintptr_t> lsAddressVecTy;
	//lsAddressVecTy lsAddressVec;

	headerID_list_vecTy headerID_list_vec;
	loopBBid2HeaderAndloopDepthPairMapTy loopBBid2HeaderAndloopDepthPairMap;
	funcName2headIDlistVecMapTy funcName2headIDlistVecMap;

	//lpBBid2EvalIterIndxMapTy lpBBid2EvalIterIndxMap;
}

using namespace llvm;

static const char *LoopInfoMDNodeName = "Loop_Info_of_this_kernel";
static const char *loopNumberMDNodeName = "Loop_Number_of_kernel";

CodeInstrument::CodeInstrument() : ModulePass(ID) {
	errs() << "Initialize CodeInstrument pass\n";
	FuncName_Array = { "recordBB", "recordLoad", "recordStore", "printf" };
	initializeCodeInstrumentPass(*PassRegistry::getPassRegistry());
}

void CodeInstrument::getAnalysisUsage(AnalysisUsage &AU) const{

	AU.addRequired<QueryBasicBlockID>();
	AU.addRequired<QueryLoadStoreID>();
	AU.addPreserved<QueryBasicBlockID>();
	AU.addPreserved<QueryLoadStoreID>();


	AU.addRequired<LoopInfo>();
	AU.addPreserved<LoopInfo>();
	//AU.addRequired<LoopNumber>();
	//AU.addPreserved<LoopNumber>();
	//AU.addRequired<DataLayout>();

	AU.setPreservesCFG();
}

bool CodeInstrument::doInitialization(Module &M) {

	// Obtain different Types
	Int8Type = IntegerType::getInt8Ty(M.getContext());
	Int32Type = IntegerType::getInt32Ty(M.getContext());
	Int64Type = IntegerType::getInt64Ty(M.getContext());
	VoidType = Type::getVoidTy(M.getContext());
	//VoidPtrType = PointerType::getUnqual(Int32Type);
	VoidPtrType = PointerType::getUnqual(Int8Type);

	LastBB = nullptr;

	/// Clone Module
	NewM = CloneModule(&M, VMap);
	
	/// recordBB function arguments' type: 
	/// Return type:	VoidType
	/// First argument:  Int64Type    -- Used to record the ID of a basic block
	/// Second argument: Int64Type  -- Used to trace whether there are call instructions from a basic block. If yes,
	///								   we do not increase the frequency of the basic block; otherwise, we increase it.
	//BBFreqIncFn = cast<Function>(M.getOrInsertFunction("recordBB", VoidType, Int64Type, Int64Type, nullptr));
	//BBFreqIncFn = cast<Function>(M.getOrInsertFunction("recordBB", VoidType, Int64Type, Int64Type, Int64Type, nullptr));
	BBFreqIncFn = cast<Function>(NewM->getOrInsertFunction("recordBB", VoidType, Int64Type, Int64Type, Int64Type, nullptr));

	/// recordLoad function arguments' type: 
	/// Return type:	 VoidType
	/// First argument:  Int64Type	 -- Used to record the ID of a basic block that contains this load/store instruction
	/// Second argument: Int64Type   -- Used to record the ID of a load instruction
	/// Third argument:  VoidPtrType -- Used to record memory address of a load instruction
	/// Fourth argument: Int64Type	 -- Used to record memory size of a load instuction
	/// Fifth argument:	 nullptr	 -- Not sure its functionnality
	//RecordLoadFn = cast<Function>(M.getOrInsertFunction("recordLoad", VoidType, Int64Type, VoidPtrType, Int64Type, nullptr));
	RecordLoadFn = cast<Function>(NewM->getOrInsertFunction("recordLoad", VoidType, Int64Type, Int64Type, VoidPtrType, Int64Type, nullptr));

	/// recordStore function arguments' type:
	/// Return type:     VoidType
	/// First argument:  Int64Type	 -- Used to record the ID of a basic block that contains this load/store instruction
	/// Second argument: Int64Type	 -- Used to record the ID of a store instruction
	/// Third argument:  VoidPtrType -- Used to record memory address of a store instruction
	/// Fourth argument: Int64Type	 -- Used to record memory size of a store instuction
	/// Fifth argument:  nullptr	 -- Not sure its functionnality
	//RecordStoreFn = cast<Function>(M.getOrInsertFunction("recordStore", VoidType, Int64Type, VoidPtrType, Int64Type, nullptr));
	RecordStoreFn = cast<Function>(NewM->getOrInsertFunction("recordStore", VoidType, Int64Type, Int64Type, VoidPtrType, Int64Type, nullptr));
	//visit(worklist.begin(), worklist.end());

	return false;
}

bool CodeInstrument::runOnModule(Module& M) {
	errs() << "\n\nBegin CodeInstrument Pass :\n";
	bbIDpass = &getAnalysis<QueryBasicBlockID>();
	lsIDpass = &getAnalysis<QueryLoadStoreID>();

	Module::iterator FI, FE;
	Function::iterator BI, BE;
	BasicBlock::iterator I, E;

#ifndef GETLOOPBOUNDPASS_H
	/// Record the IDs of basic blocks that are loop headers
	/// The invocation is moved to GetLoopBound Pass
	getLoopBBID2HeaderAndloopDepthPairMap(M);

	/// Record the IDs of header Basic Block of loops
	/// The invocation is moved to GetLoopBound Pass
	getHeaderID_list_vec(M);
#endif // End of GETLOOPBOUNDPASS_H

	// Obtain Datalayout
	//DL = M.getDataLayout();
	DL = NewM->getDataLayout();
	//DataLayoutPass = &getAnalysis<DataLayout>();
	
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

	/// Instrument Load/Store instructions
	std::vector<Instruction*> worklist;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		if (FI->isDeclaration()) {
			continue;
		}
		fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (fn_it == kernel_names.end()) {
			// Ignore uninteresting functions
			continue;
		}
		for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
			for (I = BI->begin(), E = BI->end(); I != E; ++I) {
				worklist.push_back(I);
			}
		}
	}

	visit(worklist.begin(), worklist.end());
	/// End of instrumentation of Load/Store instructions

	// Check bbIDpass & lsIDpass
	errs() << "Test bbIDpass & lsIDpass: \n";
	unsigned check_BBcounter = 0;
	unsigned check_INSTcounter = 0;
	unsigned check_LScounter = 0;

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

			errs() << "Name of BB: " << BI->getName() << "  ID: " << bbIDpass->getBBid(BI) << "\n";
			//errs() << "Name of BB: " << BI->getName() << "  ID: " << lsIDpass->getBBid(BI) << "\n";
			assert(check_BBcounter == bbIDpass->getBBid(BI) && "\nError: BB ID mismatches!\n");
			/*
			if (check_BBcounter != bbIDpass->getBBid(BI)) {
				errs() << "\nError: BB ID mismatches!\n";
				return true;
			}*/
			check_BBcounter++;
			
			/// Insert Basic Block's ID to Basic Block Map, used for recording basic block frequency
			testing_bbID2bbMap.insert(std::make_pair(bbIDpass->getBBid(BI), BI));

			// Check load/store ID
			BasicBlock::iterator I, E;
			for (I = BI->begin(), E = BI->end(); I != E; ++I) {
				if (isa<StoreInst>(I) || isa<LoadInst>(I)) {
					/*
					/// No need to check line ID, because the module has been changed and insert function entry, 
					/// thus the line IDs are already changed! However, to check load/store IDs are still meaningful.
					if (check_INSTcounter != lsIDpass->getInstLineID(I)) {
						errs() << "\nLoad/Store Instruction Line ID mismatches!\n";
					}*/
					assert(check_LScounter == lsIDpass->getlsID(I) && "\nLoad/Store Instruction ID mismatches!\n");
					/*
					if (check_LScounter != lsIDpass->getlsID(I)) {
						errs() << "\nLoad/Store Instruction ID mismatches!\n";
						return true;
					}*/
					check_LScounter++;

					/// Insert Load/Store ID to Load/Store Instruction Map, used for testing purpose
					testing_lsID2lsInstMap.insert(std::make_pair(lsIDpass->getlsID(I), I));
				}
				check_INSTcounter++;
			}
		}
	} 
	errs() << "Pass!\n"; // End testing bbIDpass & lsIDpass
	
	/// Pre-calculate lsID2BBandBBidPairMap for copying trace later
	getlsID2BBandBBidPairMap();

	/// Write the new module bitcode into a file
	write_Modified_bitcode_to_File();
	
	/// Calculate the exit basic blocks for all loops in kernel functions ** except 
	/// for the "main" function **. This vector is used to disable increasing "lastElement"
	/// when a loop exits from a header basic block to an exit basic block.
	calculate_exitBB(M);

	/// Before run on the modified bitcode, initialize BBidFreq first, otherwise
	/// it will just ignore the basic block that will not be executed and cause 
	/// program crashed.
	uint64_t sizeOfBB = bbIDpass->getNumOfBB();
	for (unsigned i = 0; i < sizeOfBB; i++) {
		BBidFreq.insert(std::make_pair(i, 0));
	}

	/// Integrate JIT profiling engine and run the embedded profiler
	EmbeddedProfilerEngine P(*NewM, BBFreqIncFn, RecordLoadFn, RecordStoreFn);
	//EmbeddedProfilerEngine P(M, BBFreqIncFn, RecordLoadFn, RecordStoreFn);
	//EmbeddedProfilerEngine P(M, RecordLoadFn, RecordStoreFn);

#ifdef ENABLE_TIMER
	Timer BaseLine("runProfiling: ");
	BaseLine.startTimer();
#endif // End of ENABLE_TIMER
	P.runOnProfiler();
#ifdef ENABLE_TIMER
	BaseLine.stopTimer();
#endif // End of ENABLE_TIMER

#ifdef VERIFY_PROFILE
	/// Verify the basic block frequency profiling result
	VerifyProfile();
#endif // End of VERIFY_PROFILE
	/// Finish profiling
	
	/// Store all load/store records into our Trace and print the trace
	Trace_f.clear(); // Clear vector Trace first
	unsigned sizeOfTrace = recordls.size();
	RecordTrace trace_entry;
#ifndef USE_LIST_FOR_RECORDLS
	for (unsigned i = 0; i < sizeOfTrace; i++) {
		RecordLoadStore LoadStoreTrace = recordls.at(i);
		//BasicBlock* BBrecord = lsIDpass->getBBbylsID(LoadStoreTrace.getlsID());
		uint64_t lsid_tr = LoadStoreTrace.getlsID();
		BasicBlock* BBrecord = lsID2BBandBBidPairMap.at(lsid_tr).first;
		uint64_t BBrecordID = lsID2BBandBBidPairMap.at(lsid_tr).second;
		//uint64_t BBrecordID = bbIDpass->getBBid(BBrecord);
		//Instruction* inst = lsIDpass->getlsbylsID(lsid_tr);
		Instruction* inst = testing_lsID2lsInstMap.at(lsid_tr);
		unsigned lineid = lsIDpass->getInstLineID(inst);
		/// FIXME: ignore lineid first, since we do not need to use it so far
		//unsigned lineid = 0;
		uint64_t traceid = i;
		trace_entry.setRecordTrace(LoadStoreTrace, lineid, BBrecord, BBrecordID, traceid);
		Trace_f.push_back(trace_entry);
	}
#else

#ifdef ENABLE_TIMER
	Timer TraceTime("getTraceTime: ");
	TraceTime.startTimer();
#endif // End of ENABLE_TIMER
	unsigned i = 0;
	std::list<RecordLoadStore>::iterator it, ie;
	for (it = recordls.begin(), ie = recordls.end(); it != ie; ++it) {
		RecordLoadStore LoadStoreTrace = *it;
		//BasicBlock* BBrecord = lsIDpass->getBBbylsID(LoadStoreTrace.getlsID());
		uint64_t lsid_tr = LoadStoreTrace.getlsID();
		BasicBlock* BBrecord = lsID2BBandBBidPairMap.at(lsid_tr).first;
		uint64_t BBrecordID = lsID2BBandBBidPairMap.at(lsid_tr).second;
		//uint64_t BBrecordID = bbIDpass->getBBid(BBrecord);
		//Instruction* inst = lsIDpass->getlsbylsID(lsid_tr);
		Instruction* inst = testing_lsID2lsInstMap.at(lsid_tr);
		unsigned lineid = lsIDpass->getInstLineID(inst);
		/// FIXME: ignore lineid first, since we do not need to use it so far
		//unsigned lineid = 0;
		uint64_t traceid = i;
		trace_entry.setRecordTrace(LoadStoreTrace, lineid, BBrecord, BBrecordID, traceid);
		Trace_f.push_back(trace_entry);
		i++;
	}
#ifdef ENABLE_TIMER
	TraceTime.stopTimer();
#endif // End of ENABLE_TIMER

#endif // End of USE_LIST_FOR_RECORDLS

	//print_trace_info(Trace_f);

	/*
	/// Create output trace file name and write to the file
	std::string TraceFilename;
	if (TraceFilename.empty()) {
		if (inputFileName == "-") {
			TraceFilename = "-";
		}
		else {
			const std::string &IFN = inputFileName;
			int Len = IFN.length();
			// If the source ends in .bc, strip it off
			if (IFN[Len - 3] == '.' && IFN[Len - 2] == 'b' && IFN[Len - 1] == 'c') {
				TraceFilename = std::string(IFN.begin(), IFN.end() - 3) + ".trace";
			}
		}
	}
	write_trace_file(Trace_f, TraceFilename);
	*/

	// Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		assert(false && "Module is broken!\n");
	}
	

	/*
	/// Check CFG graph
	Module::iterator MI, ME;
	for (MI = M.begin(), ME = M.end(); MI != ME; ++MI) {
		runOnFunction(*MI);
	}*/
	delete NewM;
	errs() << "\nEnd CodeInstrument Pass\n\n";
	return false;
	//return true;
}

/// Calculate the exit basic blocks for all loops in kernel functions ** except 
/// for the "main" function **. This vector is used to disable increasing "lastElement"
/// when a loop exits from a header basic block to an exit basic block.
void CodeInstrument::calculate_exitBB(Module &M) {
	Module::iterator FI, FE;
	std::vector<std::string>::iterator fn_it;

	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		// We only focus on analyze kernel functions instead of recordBB, recordLoad/Store, 
		// printf, main, etc., functions
		if (FI->isDeclaration()) {
			errs() << "Function name: " << FI->getName() << "\n";
			continue;
		}
		fn_it = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (fn_it == kernel_names.end()) {
			// We do not consider uninteresting functions
			continue;
		}

		if (ignoreFunc(FI->getName()) == false) {
			if (FI->getName() != "main") {
				// FIXME: If there is a loop that is not canonical loop, then the function
				//		  "getExitBlock()" may be failed. Please be aware of this problem
				std::vector<uint64_t> exitBBidVec;
				LoopInfo* LI = &getAnalysis<LoopInfo>(*FI);
				errs() << "Function name = " << FI->getName() << "\n";
				assert(funcName2headIDlistVecMap.size() != 0 && "Error, funcName2headIDlistVecMap has size 0!\n");
				headerID_list_vecTy headerID_list_vec_temp = funcName2headIDlistVecMap.at(FI->getName());
				unsigned size_vec = headerID_list_vec_temp.size();
				for (unsigned i = 0; i < size_vec; i++) {
					std::list<uint64_t> headerID_list_temp = headerID_list_vec_temp.at(i);
					std::list<uint64_t>::iterator it_hid, ie_hid;
					for (it_hid = headerID_list_temp.begin(), ie_hid = headerID_list_temp.end(); it_hid != ie_hid; ++it_hid) {
						BasicBlock* headerBB_temp = bbIDpass->getBBbyID(*it_hid);
						Loop* lp_temp = LI->getLoopFor(headerBB_temp);
						BasicBlock* exitBB_temp = lp_temp->getExitBlock();
						if (exitBB_temp != NULL) {
							uint64_t bbid = bbIDpass->getBBid(exitBB_temp);
							exitBBidVec.push_back(bbid);
							whole_exitBBid_vec.push_back(bbid);
							/// FIXME: If there are multiple exit basic blocks for a header, then it might be a problem here.
							///		   Currently, we only consider one header is associated with one exit basic block
							exitBBid2headerBBidMap.insert(std::make_pair(bbid, *it_hid));
						}
					}
				}
				funcName2exitBBidMap.insert(std::make_pair(FI->getName(), exitBBidVec));
				exitBBidVec.clear();
			}
		}
	}
} /// End calculate_exitBB

/// Write the new module bitcode into a file
void CodeInstrument::write_Modified_bitcode_to_File() {

	if (inputFileName == "-") {
		OutputModifiedBCFilename = "-";
	}
	else {
		const std::string &IFN = inputFileName;
		int Len = IFN.length();
		// If the source ends in .bc, strip it off
		if (IFN[Len - 3] == '.' && IFN[Len - 2] == 'b' && IFN[Len - 1] == 'c') {
			OutputModifiedBCFilename = std::string(IFN.begin(), IFN.end() - 3) + "_Modified.bc";
		}
	}

	std::unique_ptr<tool_output_file> Out;
	std::string ErrorInfo;
	Out.reset(new tool_output_file(OutputModifiedBCFilename.c_str(), ErrorInfo, sys::fs::F_None));
	if (!ErrorInfo.empty()) {
		errs() << ErrorInfo << "\n";
		assert(false && "ErrorInfo is not empty in OutputModifiedBCFilename\n");
	}

	WriteBitcodeToFile(NewM, Out->os());
	Out->os().close();
	if (Out->os().has_error()) {
		errs() << "could not write bitcode file: " << OutputModifiedBCFilename.c_str() << "\n";
		Out->os().clear_error();
		assert(false && "could not write bitcode file OutputModifiedBCFilename\n");
	}

	Out->keep();
} /// We have successfully written the modified bitcode into a new file

/// FIXME: This code will fail if more than one loop inside a kernel. 
/// Bug fixed
void CodeInstrument::getLoopBBID2HeaderAndloopDepthPairMap(Module &M){
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

void CodeInstrument::getHeaderID_list_vec(Module &M){
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

/// Pre-calculate lsID2BBandBBidPairMap for copying trace later
void CodeInstrument::getlsID2BBandBBidPairMap() {
	std::map<unsigned, llvm::Instruction*>::iterator it, ie;
	for (it = testing_lsID2lsInstMap.begin(), ie = testing_lsID2lsInstMap.end(); it != ie; ++it) {
		uint64_t lsid = it->first;
		Instruction* inst = it->second;
		BasicBlock* bb = inst->getParent();
		uint64_t bbid = bbIDpass->getBBid(bb);
		lsID2BBandBBidPairMap.insert(std::make_pair(lsid, std::make_pair(bb, bbid)));
	}
}

bool CodeInstrument::runOnFunction(Function& Fn){
	Fn.viewCFGOnly();
	return false;
}

bool CodeInstrument::isLoopHeader (uint64_t bbid) {
	BasicBlock* BB = bbIDpass->getBBbyID(bbid);
	Function* func = BB->getParent();
	LoopInfo* LI = &getAnalysis<LoopInfo>(*func);
	return LI->isLoopHeader(BB);
}

/*
bool CodeInstrument::runOnBasicBlock(BasicBlock& BB) {	
	return true;
}*/

void CodeInstrument::instrumentFunction(Function* F) {
	for (Function::iterator BI = F->begin(), BE = F->end(); BI != BE; ++BI) {
		instrumentBasicBlock(BI);
	}
}

void CodeInstrument::instrumentBasicBlock(BasicBlock* BB) {
	bool edge_increase;
	if (BB == &BB->getParent()->getEntryBlock()) {
		edge_increase = false;
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
			std::string Funname = callFunc->stripPointerCasts()->getName().str();
			
			if (ignoreFunc(Funname)) {
				errs() << "Instruction in Funname: " << Funname << "\n";
			}
			else {
				BasicBlock::iterator InsertPos = I;
				errs() << "Before insertCounterBefore Instr: " << *I << "\n";
				insertCounterBefore(++InsertPos, false, edge_increase);
			}
			errs() << "Instruction in BB " << *I << "\n";
			continue;
		}
	}

	insertCounterBefore(BB->getTerminator(), true, edge_increase);
}

bool CodeInstrument::ignoreFunc(std::string func_name) {
	std::vector<std::string>::iterator it;
	it = std::find(FuncName_Array.begin(), FuncName_Array.end(), func_name);
	if (it != FuncName_Array.end()) {
		return true;
	}
	if (func_name.find("llvm.memcpy") != std::string::npos) {
		return true;
	}
	
	return false;
}

void CodeInstrument::insertCounterBefore(Instruction* Inst, bool increase, bool edgeIncrease) {
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

	errs() << "insertCounterBefore Instr: " << *Inst << "\n";
	
	BasicBlock *OldBB = Inst->getParent();
	unsigned bb_id = bbIDpass->getBBid(OldBB);
	Instruction *MappedInst = cast<Instruction>((Value*)VMap.lookup(Inst));

	// Get the entry instruction of a Basic Block
	Instruction* entryNonPhiInst = OldBB->getFirstNonPHI();
	Instruction* MappedFirstInst = cast<Instruction>((Value*)VMap.lookup(entryNonPhiInst));

	errs() << "insertCounterBefore OldBB: " << OldBB->getName() << "\n";
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

void CodeInstrument::visitLoadInst (LoadInst& inst) {
	/*
	Function* func = inst.getParent()->getParent();
	/// We do not consider loops inside the 'main' function
	if (func->getName() != "main") {
		instrumentLoadInst(&inst);
	}*/

	Function* func = inst.getParent()->getParent();
	/// We do not consider loops inside the 'main' function
	if (func->getName() != "main") {
		/// we only consider load/store instructions inside loops
		BasicBlock* bb = inst.getParent();
		uint64_t bbid = bbIDpass->getBBid(bb);

		loopBBid2HeaderAndloopDepthPairMapTy::iterator it;
		it = loopBBid2HeaderAndloopDepthPairMap.find(bbid);
		if (it != loopBBid2HeaderAndloopDepthPairMap.end()) {
			// This basic block is inside a loop
			instrumentLoadInst(&inst);
		}
		// we ignore load/store instructions outside a loop
	}
}

void CodeInstrument::visitStoreInst (StoreInst& inst) {
	/*
	Function* func = inst.getParent()->getParent();
	/// We do not consider loops inside the 'main' function
	if (func->getName() != "main") {
		instrumentStoreInst(&inst);
	}*/

	Function* func = inst.getParent()->getParent();
	/// We do not consider loops inside the 'main' function
	if (func->getName() != "main") {
		/// we only consider load/store instructions inside loops
		BasicBlock* bb = inst.getParent();
		uint64_t bbid = bbIDpass->getBBid(bb);

		loopBBid2HeaderAndloopDepthPairMapTy::iterator it;
		it = loopBBid2HeaderAndloopDepthPairMap.find(bbid);
		if (it != loopBBid2HeaderAndloopDepthPairMap.end()) {
			// This basic block is inside a loop
			instrumentStoreInst(&inst);
		}
		// we ignore load/store instructions outside a loop
	}
}

void CodeInstrument::instrumentLoadInst (Instruction* inst) {
	// For Debug purpose:
	errs() << "Instrumented Load instruction:  " << *inst << "\n";
	/*
	LoadInst* ldInst = dyn_cast<LoadInst>(inst);
	Value* LDPointer = ldInst->getPointerOperand();
	Type* VTy = LDPointer->getType()->getPointerElementType();

	// Get ID of a store instruction
	Value* loadid = ConstantInt::get(Int64Type, lsIDpass->getlsID(inst));

	// Cast the pointer to a void pointer type
	Value* loadPtr = ldInst->getPointerOperand();
	loadPtr = castTo(loadPtr, VoidPtrType, loadPtr->getName(), ldInst);

	// Get the size of stored data
	uint64_t data_size = DL->getTypeStoreSize(VTy);
	Value* load_size = ConstantInt::get(Int64Type, data_size);

	std::vector< Value* > Args;
	Args.push_back(loadid);
	Args.push_back(loadPtr);
	Args.push_back(load_size);
	*/
	// For Debug purpose:
	//errs() << (int64_t) DL->getTypeStoreSize(VTy) << "\n";

	Instruction *MappedInst = cast<Instruction>((Value*)VMap.lookup(inst));
	LoadInst* ldInst = dyn_cast<LoadInst>(MappedInst);
	Value* LDPointer = ldInst->getPointerOperand();
	Type* VTy = LDPointer->getType()->getPointerElementType();

	// Get ID of a basic block that contains this instruction
	BasicBlock* BB = inst->getParent();
	Value* bbid = ConstantInt::get(Int64Type, bbIDpass->getBBid(BB));

	// Get ID of a store instruction
	Value* loadid = ConstantInt::get(Int64Type, lsIDpass->getlsID(inst));

	// Cast the pointer to a void pointer type
	Value* loadPtr = ldInst->getPointerOperand();
	loadPtr = castTo(loadPtr, VoidPtrType, loadPtr->getName(), ldInst);

	// Get the size of stored data
	// FIXME: Do I need to use the DL of cloned Module?
	uint64_t data_size = DL->getTypeStoreSize(VTy);
	Value* load_size = ConstantInt::get(Int64Type, data_size);

	std::vector< Value* > Args;
	Args.push_back(bbid);
	Args.push_back(loadid);
	Args.push_back(loadPtr);
	Args.push_back(load_size);
	//CallInst::Create(RecordLoadFn, Args, "", inst);
	CallInst::Create(RecordLoadFn, Args, "", MappedInst);
}

void CodeInstrument::instrumentStoreInst (Instruction* inst) {
	// For Debug purpose:
	errs() << "Instrumented Store instruction:  " << *inst << "\n";
	/*
	StoreInst* stInst = dyn_cast<StoreInst>(inst);
	Value* STPointer = stInst->getPointerOperand();
	Type* VTy = STPointer->getType()->getPointerElementType();

	// Get ID of a store instruction
	Value* storeid = ConstantInt::get(Int64Type, lsIDpass->getlsID(inst));
	
	// Cast the pointer to a void pointer type
	Value* storePtr = stInst->getPointerOperand();
	storePtr = castTo(storePtr, VoidPtrType, storePtr->getName(), stInst);

	// Get the size of stored data
	uint64_t data_size = DL->getTypeStoreSize(VTy);
	Value* store_size = ConstantInt::get(Int64Type, data_size);

	std::vector< Value* > Args;
	Args.push_back(storeid);
	Args.push_back(storePtr);
	Args.push_back(store_size);
	*/
	Instruction *MappedInst = cast<Instruction>((Value*)VMap.lookup(inst));

	StoreInst* stInst = dyn_cast<StoreInst>(MappedInst);
	Value* STPointer = stInst->getPointerOperand();
	Type* VTy = STPointer->getType()->getPointerElementType();

	// Get ID of a basic block that contains this instruction
	BasicBlock* BB = inst->getParent();
	Value* bbid = ConstantInt::get(Int64Type, bbIDpass->getBBid(BB));

	// Get ID of a store instruction
	Value* storeid = ConstantInt::get(Int64Type, lsIDpass->getlsID(inst));

	// Cast the pointer to a void pointer type
	Value* storePtr = stInst->getPointerOperand();
	storePtr = castTo(storePtr, VoidPtrType, storePtr->getName(), stInst);

	// Get the size of stored data
	// FIXME: Do I need to use the DL of cloned Module?
	uint64_t data_size = DL->getTypeStoreSize(VTy);
	Value* store_size = ConstantInt::get(Int64Type, data_size);

	std::vector< Value* > Args;
	Args.push_back(bbid);
	Args.push_back(storeid);
	Args.push_back(storePtr);
	Args.push_back(store_size);
	//CallInst::Create(RecordStoreFn, Args, "", inst);
	CallInst::Create(RecordStoreFn, Args, "", MappedInst);
}

#ifdef VERIFY_PROFILE
void CodeInstrument::VerifyProfile() {
	typedef BBidFreqMap::const_iterator bbid_iterator;
	typedef BranchidFreqMap::const_iterator brid_iterator;
	typedef BBFreqMap::const_iterator bb_iterator;

	errs() << "\nVerify the profile data: \n";
	errs() << "Number of Load/Store instructions: " << recordls.size() << "\n";
	unsigned load_counter = 0;
	unsigned store_counter = 0;

#ifndef USE_LIST_FOR_RECORDLS
	for (unsigned i = 0; i < recordls.size(); i++) {
		switch (recordls.at(i).getlsRecordType()) {
		case RecordType::LOADType:
			load_counter++;
			break;
		case RecordType::STOREType:
			store_counter++;
			break;
		default:
			errs() << "Default\n";// Do nothing
		}
	}
#else
	std::list<RecordLoadStore>::iterator it, ie;
	for (it = recordls.begin(), ie = recordls.end(); it != ie; ++it) {
		switch ( it->getlsRecordType()) {
		case RecordType::LOADType:
			load_counter++;
			break;
		case RecordType::STOREType:
			store_counter++;
			break;
		default:
			errs() << "Default\n";// Do nothing
		}
	}
#endif // End of USE_LIST_FOR_RECORDLS

	errs() << "Number of Load instructions: " << load_counter << "\n";
	errs() << "Number of Store instructions: " << store_counter << "\n";
	
	/// Convert Basic Block IDs in EmbeddedProfilerEngine::BBidFreqMap & BranchidFreqMap into Basic 
	/// Block * in CodeInstrument::BBFreqMap & BranchFreqMap
	/// BBFreq
	errs() << "BBFreq: \n";
	//for (bbid_iterator I = GlobalContext->P->BBFreq.begin(), E = GlobalContext->P->BBFreq.end(); I != E; ++I) {
	for (bbid_iterator I = BBidFreq.begin(), E = BBidFreq.end(); I != E; ++I) {
		unsigned bbid = I->first;
		unsigned freq = I->second;
		errs() << "(" << bbid << ", " << freq << ")\n";
		BasicBlock *BB = bbIDpass->getBBbyID(bbid);
		BBFreq.insert(std::make_pair(BB, freq));
	}
	/// BranchFreq
	errs() << "BranchFreq: \n";
	BBFreqMap bbFreq_temp;
	//for (brid_iterator I = GlobalContext->P->BranchFreq.begin(), E = GlobalContext->P->BranchFreq.end(); I != E; ++I) {
	for (brid_iterator I = BranchIDFreq.begin(), E = BranchIDFreq.end(); I != E; ++I) {
		unsigned bbid_dest = I->first;
		BasicBlock* BB_dest = bbIDpass->getBBbyID(bbid_dest);
		for (bbid_iterator it = I->second.begin(), ie = I->second.end(); it != ie; ++it) {
			unsigned bbid_src = it->first;
			unsigned br_freq = it->second;
			BasicBlock* BB_src = bbIDpass->getBBbyID(bbid_src);
			bbFreq_temp.insert(std::make_pair(BB_src, br_freq));
			errs() << bbid_dest << " -> (" << bbid_src << ", " << br_freq << ")\n";
		}
		
		BranchFreq.insert(std::make_pair(BB_dest, bbFreq_temp));
		bbFreq_temp.clear();
	}
	
	/// Checking the Basic Block Frequency and Edge Frequency
	for (bb_iterator I = BBFreq.begin(), E = BBFreq.end(); I != E; ++I) {
		BasicBlock *BB = I->first;
		// We do not count the call edge frequency
		if (BB == &BB->getParent()->getEntryBlock())
			continue;

		uint64_t Freq = I->second;

		// If a basic block is not executed at all, its frequency will be equal to 0,
		// we just ignore it.
		if (Freq == 0) {
			continue;
		}

		errs() << "BasicBlock:  " << BB->getName() \
			   << "  Frequency: " << Freq << "\n";

		BranchFreqMap::const_iterator BrFreqsAt = BranchFreq.find(BB);
		assert(BrFreqsAt != BranchFreq.end() && "Freq info not found!");
		const BBFreqMap &BrFreqs = BrFreqsAt->second;
		for (bb_iterator I = BrFreqs.begin(), E = BrFreqs.end(); I != E; ++I)
			Freq -= I->second;

		if (Freq != 0) {
			llvm_unreachable("Bad profile infomration!");
		}
	}
	
} // End of CodeInstrument::VerifyProfile()


/// Functions for accessing Basic Block frequency
uint64_t CodeInstrument::getBBFreq(BasicBlock* BB) const {
	BBFreqMap::const_iterator I = BBFreq.find(BB);
	return (I == BBFreq.end() ? 0 : I->second);
} // End of CodeInstrument::getBBFreq


uint64_t CodeInstrument::getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const {
	BranchFreqMap::const_iterator I = BranchFreq.find(Dest);
	if (I == BranchFreq.end()) {
		return 0;
	}

	const BBFreqMap& CurFreqs = I->second;
	BBFreqMap::const_iterator J = CurFreqs.find(Src);

	return (J == CurFreqs.end() ? 0 : J->second);
} // End of CodeInstrument::getBranchFreq()
#endif // End of VERIFY_PROFILE

/// Insert record functions below
void recordBB(uint64_t BBid, uint64_t inc, uint64_t edge_increase) {
#ifndef IGNORE_RECORDBB
	//errs() << "BBid = " << BBid << "\n";
	//BasicBlock* BB = reinterpret_cast<BasicBlock*>(BBPtr);
	std::map<unsigned, BasicBlock*>::iterator it;
	it = testing_bbID2bbMap.find(BBid);
	assert(it != testing_bbID2bbMap.end() && "Can not find BBid in testing_bbID2bbMap! Error\n");
	BasicBlock* BB = testing_bbID2bbMap.at(BBid);
	uint64_t lastBBid = GlobalContext->getLastBB();
	BasicBlock* lastBB;
	if (lastBBid == INITIAL_BIGVALUE) {
		/// Initial lastBBid in GlobalContext, just ignore it
		lastBB = nullptr;
	}
	else {
		lastBB = testing_bbID2bbMap.at(lastBBid);
	}
	
	if (inc) {
		GlobalContext->P->increaseEdgeCounter(lastBB, lastBBid, BB, BBid, edge_increase);
		
		/// Record evaluation iteration index here, we only need to consider basic blocks inside loops
		loopBBid2HeaderAndloopDepthPairMapTy::iterator it;
		it = loopBBid2HeaderAndloopDepthPairMap.find(BBid);
		if (it != loopBBid2HeaderAndloopDepthPairMap.end()) {
			// This basic block is inside a loop
			// Check whether this basic block is a loop header
			uint64_t headerid = it->second.first;
			//std::vector<std::list<unsigned> > vecEvalIterIndex;
			std::list<unsigned> eval_iteration_index;

			/// The Basic Block is a loop header
			if (BBid == headerid) {
				// This basic block is a loop header
				unsigned loop_depth = it->second.second;
				lpBBid2EvalIterIndxMapTy::iterator it_eval;

				it_eval = lpBBid2EvalIterIndxMap.find(BBid);
				if (it_eval == lpBBid2EvalIterIndxMap.end()) {
					// We are going to record the evaluation iteration index for this basic block ID here
					// For the top level loop header, we only need one iteration variable.
					if (loop_depth == 1) {
						eval_iteration_index.push_back(0);
					}
					else {
						/// Find the target headerID_list from headerID_list_vec
						std::list<uint64_t> headerIDList = getHeaderID_list(BBid);

						for (unsigned i = 0; i < loop_depth; i++) {
							// Initialize the eval_iteration_index
							if (i == loop_depth-1) {
								eval_iteration_index.push_back(0);
							}
							else {
								// Trace the basic block id of different level loop, starting from the 
								// top-level loop. Based on the basic block id, we search the 
								// lpBBid2EvalIterIndxMap to get the value of the corresponding list and
								// assign it to current basic block's eval_iteration_index

								// Solution: trace the header basic block using Metadata Node, then we 
								//           pre-get a list that contains header IDs for each level loop

								// Get header id of the upper-level loop, otherwise delete it from the tail of headerIDList
								// The reason is that we only need to update the evaluation iteration index based on the
								// header basic block of the parent-level loop. For example, the evaluation iteration index
								// for header bb id of the top-level loop is 9, then the evaluation interation index of header
								// basic block has also "9" in the field for the top-level loop
								if (i == loop_depth-2) {
									uint64_t upper_HeaderID = headerIDList.back();
									lpBBid2EvalIterIndxMapTy::iterator it_header_eval;
									it_header_eval = lpBBid2EvalIterIndxMap.find(upper_HeaderID);
									assert(it_header_eval != lpBBid2EvalIterIndxMap.end() && "can not find header in lpBBid2EvalIterIndxMap. Error!\n");
									eval_iteration_index = lpBBid2EvalIterIndxMap.at(upper_HeaderID);
								}
								else {
									headerIDList.pop_back();
								}

							}
						}
					}

					//vecEvalIterIndex.push_back(eval_iteration_index);
					//lpBBid2EvalIterIndxMap.insert(std::make_pair(BBid, vecEvalIterIndex));
					lpBBid2EvalIterIndxMap.insert(std::make_pair(BBid, eval_iteration_index));
					//vecEvalIterIndex.clear();
					/// RecordBB function should be inserted into the beginning of the basic block??
					/// Because we can use the evaluation iteration index when invoking recordLoad 
					/// or recordStore functions.
				}
				else {
					// We have recorded previous, just need to update the index. 
					// The last element in eval_iteration_index is the innermost loop we explore so far
					eval_iteration_index.assign(it_eval->second.begin(), it_eval->second.end());
					unsigned lastElement = eval_iteration_index.back();

					/// Find the target headerID_list from headerID_list_vec
					/// headerIDList : the innermost-level --> the top-level loop
					std::list<uint64_t> headerIDList_temp = getHeaderID_list(BBid);
					if (BBid == headerIDList_temp.back()) {
						// If the BBid is the header of the top-level loop, 
						// then we just update the eval-iteration.
						eval_iteration_index.pop_back();
						//eval_iteration_index.push_back(lastElement + 1);
						//it_eval->second.assign(eval_iteration_index.begin(), eval_iteration_index.end());
					}
					else {
						std::list<uint64_t>::iterator it_headerID;
						it_headerID = std::find(headerIDList_temp.begin(), headerIDList_temp.end(), BBid);
						assert(it_headerID != headerIDList_temp.end() && "Error, we can not find BBid in headerIDList_temp\n");
						it_headerID++;
						uint64_t upper_headerid = *it_headerID;
						std::list<unsigned> upperHeader_EvalList = lpBBid2EvalIterIndxMap.at(upper_headerid);
						eval_iteration_index.assign(upperHeader_EvalList.begin(), upperHeader_EvalList.end());
						//eval_iteration_index.push_back(lastElement + 1);
						//it_eval->second.assign(eval_iteration_index.begin(), eval_iteration_index.end());
					}
					/// FIXME: 
					/// If we perform two times profiling, first we get loop bounds of each loop levels in loops,
					/// then we need to make sure the counter will not exceed the corresponding loop bound. To do that,
					/// We just need to insert if condition below to wrap 'eval_iteration_index.push_back(lastElement + 1);',
					/// and remove 'modifyTrace_EvalIterationIndex()' in AnalysisProfiling Pass
					///
					/// If we perform one time profiling, just leave it and reserve 'modifyTrace_EvalIterationIndex()' in 
					/// AnalysisProfiling Pass
#ifndef GETLOOPBOUNDPASS_H // One-time profiling
					eval_iteration_index.push_back(lastElement + 1);
#else // Two-time profiling
					uint64_t lp_bound = BBheaderID2loopBoundMap.at(BBid);
					eval_iteration_index.push_back( (lastElement + 1) % lp_bound );
#endif // End of GETLOOPBOUNDPASS_H
					it_eval->second.assign(eval_iteration_index.begin(), eval_iteration_index.end());
					/*
					std::vector<uint64_t>::iterator it_exitBBidvec;
					it_exitBBidvec = std::find(whole_exitBBid_vec.begin(), whole_exitBBid_vec.end(), BBid);
					if (it_exitBBidvec != whole_exitBBid_vec.end()) {
						// Do not increase lastElement if we exit from a header basic block to a exit basic block
						eval_iteration_index.push_back(lastElement);
					}
					else {
						eval_iteration_index.push_back(lastElement + 1);
					}
					it_eval->second.assign(eval_iteration_index.begin(), eval_iteration_index.end());
					*/
				}
				eval_iteration_index.clear();
			} /// The Basic Block is NOT a loop header
			else {
				/// We need to find out this basic block is belonged to which loop header, and then we 
				/// copy the loop header's eval_iteration_index into this basic block's eval_iteration_index
				loopBBid2HeaderAndloopDepthPairMapTy::iterator it_loopBBid2Pair;
				it_loopBBid2Pair = loopBBid2HeaderAndloopDepthPairMap.find(BBid);
				assert(it_loopBBid2Pair != loopBBid2HeaderAndloopDepthPairMap.end() && "Can not find BB id in loopBBid2HeaderAndloopDepthPairMap, ERROR!\n");
				uint64_t header_of_BBid = it_loopBBid2Pair->second.first;
				eval_iteration_index.assign(lpBBid2EvalIterIndxMap.at(header_of_BBid).begin(), lpBBid2EvalIterIndxMap.at(header_of_BBid).end());

				std::vector<uint64_t>::iterator it_exitBBidvec;
				it_exitBBidvec = std::find(whole_exitBBid_vec.begin(), whole_exitBBid_vec.end(), BBid);
				if (it_exitBBidvec != whole_exitBBid_vec.end()) {
					// Need to decrease lastElement if we exit from a header basic block to a exit basic block. Because we have increase it when we 
					// visit the header before this exit basic block.
					// Before doing that, we need to find the corresponding header basic block id first, then we need to substract the lastElement from
					// lpBBid2EvalIterIndxMap
					uint64_t relatedHeaderid = exitBBid2headerBBidMap.at(BBid);
					assert(lpBBid2EvalIterIndxMap.find(relatedHeaderid)!=lpBBid2EvalIterIndxMap.end() && "Error, header ID should be in lpBBid2EvalIterIndxMap\n");
					std::list<unsigned> eval_list_temp = lpBBid2EvalIterIndxMap.at(relatedHeaderid);
					unsigned lastIter = eval_list_temp.back() - 1;
					eval_list_temp.pop_back();
					eval_list_temp.push_back(lastIter);
					lpBBid2EvalIterIndxMap.at(relatedHeaderid).assign(eval_list_temp.begin(), eval_list_temp.end());
				}
				
				lpBBid2EvalIterIndxMapTy::iterator it_lpBBid2Eval;
				it_lpBBid2Eval = lpBBid2EvalIterIndxMap.find(BBid);
				if (it_lpBBid2Eval != lpBBid2EvalIterIndxMap.end()) {
					// If we already record this basic block into this lpBBid2EvalIterIndxMap map, then we just need to update it
					lpBBid2EvalIterIndxMap.at(BBid).assign(eval_iteration_index.begin(), eval_iteration_index.end());
				}
				else {
					// If we do not record this basic block into this lpBBid2EvalIterIndxMap map, then we just need to record it
					lpBBid2EvalIterIndxMap.insert(std::make_pair(BBid, eval_iteration_index));
				}
				eval_iteration_index.clear();
			}
		}
		else {
			/// For the basic blcoks that are not belonged to any loops, we just do nothing here.
		}
		
	}
	GlobalContext->setLastBB(BBid);
#endif // End of IGNORE_RECORDBB
}

void recordLoad(uint64_t BBid, uint64_t id, unsigned char* ptr, uint64_t load_size) {
	/// Check the BB and Function
	Instruction* inst = testing_lsID2lsInstMap.at(id);
	BasicBlock* BB = inst->getParent();
	Function* func = BB->getParent();
	//errs() << "Instrumented Instruction = " << *inst << "  BB = " << BB->getName() \
	//	   << "BBid = " << BBid << " Function = " << func->getName() << "\n";

	RecordLoadStore load_record(RecordType::LOADType, id, ptr, load_size);
	lpBBid2EvalIterIndxMapTy::iterator it;
	it = lpBBid2EvalIterIndxMap.find(BBid);
	if (it != lpBBid2EvalIterIndxMap.end()) {
		std::list<unsigned> eval_iter_index_temp = lpBBid2EvalIterIndxMap.at(BBid);
		load_record.setEval_Iteration_Index(eval_iter_index_temp);
	}

	recordls.push_back(load_record);

	/// lsAddressVec is used to analyze the dynamic data dependence
	//lsAddressVec.push_back(load_record.getls_address());
}

void recordStore(uint64_t BBid, uint64_t id, unsigned char* ptr, uint64_t store_size) {
	/// Check the BB and Function
	Instruction* inst = testing_lsID2lsInstMap.at(id);
	BasicBlock* BB = inst->getParent();
	Function* func = BB->getParent();
	//errs() << "Instrumented Instruction = " << *inst << "  BB = " << BB->getName() \
	//	   << "BBid = " << BBid << " Function = " << func->getName() << "\n";

	RecordLoadStore store_record(RecordType::STOREType, id, ptr, store_size);
	lpBBid2EvalIterIndxMapTy::iterator it;
	it = lpBBid2EvalIterIndxMap.find(BBid);
	if (it != lpBBid2EvalIterIndxMap.end()) {
		std::list<unsigned> eval_iter_index_temp = lpBBid2EvalIterIndxMap.at(BBid);
		store_record.setEval_Iteration_Index(eval_iter_index_temp);
	}

	recordls.push_back(store_record);

	/// lsAddressVec is used to analyze the dynamic data dependence
	//lsAddressVec.push_back(store_record.getls_address());
}

/// The above functions are used for recording BB, load/store instructions

std::list<uint64_t> llvm::getHeaderID_list(uint64_t Aheader_id) {
	unsigned size_vec = headerID_list_vec.size();
	for (unsigned i = 0; i < size_vec; i++) {
		std::list<uint64_t> headerList = headerID_list_vec.at(i);
		std::list<uint64_t>::iterator it;
		it = std::find(headerList.begin(), headerList.end(), Aheader_id);
		if (it != headerList.end()) {
			// We find the headerList, just need to return it
			return headerList;
		}
	}

	assert(false && "Can not find any headerList in headerID_list_vec, ERROR!\n");
	// It is impossible to execute this return, just return random stuff
	return headerID_list_vec.at(0);
}

EmbeddedProfilerEngine::EmbeddedProfilerEngine(Module &M, Function* bbFreqfn, Function* loadFn, Function* storeFn) :
Mod(M), BBFreqFn(bbFreqfn), RecordLoadFn(loadFn), RecordStoreFn(storeFn){
}

void EmbeddedProfilerEngine::increaseEdgeCounter(BasicBlock* src, uint64_t src_id, BasicBlock* dest, uint64_t dest_id, uint64_t increase_edge){
	//++BBFreq[dest];
	++BBidFreq[dest_id];
	if (src == nullptr) {
		//errs() << "Initial value, src is a nullptr\n";
	}
	else {
		//errs() << "uint64_t(src):  " << uint64_t(src);
		//errs() << "    BB src = " << src->getName() << "\n";
		//errs() << "    BB src id = " << src_id << "\n";
		//errs() << "uint64_t(dest):  " << uint64_t(dest);
		//errs() << "    BB dest id = " << dest_id << "\n";
		//std::string dest_name = dest->getName();
		//errs() << "    BB dest = " << dest_name << "\n";
		//errs() << "        Belonged to Function = " << dest->getParent()->getName() << "\n";
	}
	
	/*
	// Ignore the call edges.  Why?
	if (dest == &dest->getParent()->getEntryBlock())
		return;
	*/
	if (increase_edge) {
		//++BranchFreq[dest][src];
		++BranchIDFreq[dest_id][src_id];
	}
}

/*
void EmbeddedProfilerEngine::VerifyProfile() const {
	typedef BBFreqMap::const_iterator bb_iterator;
	errs() << "\nVerify the profile data: \n";
	errs() << "Number of Load/Store instructions: " << recordls.size() << "\n";
	unsigned load_counter = 0;
	unsigned store_counter = 0;
	for (unsigned i = 0; i < recordls.size(); i++) {
		switch (recordls.at(i).getlsRecordType()) {
			case RecordType::LOADType:
				load_counter++;
				break;
			case RecordType::STOREType:
				store_counter++;
				break;
			default:
				errs() << "Default\n";// Do nothing
		}
	}
	errs() << "Number of Load instructions: " << load_counter << "\n";
	errs() << "Number of Store instructions: " << store_counter << "\n";

	
	for (bb_iterator I = BBFreq.begin(), E = BBFreq.end(); I != E; ++I) {
		BasicBlock *BB = I->first;
		// We do not count the call edge frequency
		if (BB == &BB->getParent()->getEntryBlock())
			continue;

		uint64_t Freq = I->second;

		errs() << "BasicBlock:  " << BB->getName() \
			   << "  Frequency: " << Freq << "\n";

		BranchFreqMap::const_iterator BrFreqsAt = BranchFreq.find(BB);
		assert(BrFreqsAt != BranchFreq.end() && "Freq info not found!");
		const BBFreqMap &BrFreqs = BrFreqsAt->second;
		for (bb_iterator I = BrFreqs.begin(), E = BrFreqs.end(); I != E; ++I)
			Freq -= I->second;

		if (Freq != 0) {
			llvm_unreachable("Bad profile infomration!");
		}
	}
	
} // End of EmbeddedProfilerEngine::VerifyProfile()
*/

/*
uint64_t EmbeddedProfilerEngine::getBBFreq(BasicBlock* BB) const {
	BBFreqMap::const_iterator I = BBFreq.find(BB);
	return (I == BBFreq.end() ? 0 : I->second);
} // End of EmbeddedProfilerEngine::getBBFreq
*/

/*
uint64_t EmbeddedProfilerEngine::getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const {
	BranchFreqMap::const_iterator I = BranchFreq.find(Dest);
	if (I == BranchFreq.end()) {
		return 0;
	}

	const BBFreqMap& CurFreqs = I->second;
	BBFreqMap::const_iterator J = CurFreqs.find(Src);

	return (J == CurFreqs.end() ? 0 : J->second);
} // End of EmbeddedProfilerEngine::getBranchFreq
*/


void EmbeddedProfilerEngine::runOnProfiler(){
	ProfilerJITSingletonContext JTSC(this);
	//std::unique_ptr<Module> M(generateInstrumentedModule());

	// Before running profiling, we need to clear the vector recordls
	recordls.clear();

	// TODO: Why I always fail when using CloneModule?
	//Module* M = CloneModule(&Mod, V2VMap);
	//std::unique_ptr<Module> M(newM);
	Module *M = (&Mod);

	// TODO: Insert instrumental code to extract the trace.
	EngineBuilder builder(M);
	ExecutionEngine *EE = builder.setEngineKind(EngineKind::JIT).create();
	// Where is getBBFreqInc()?
	//EE->addGlobalMapping(getBBFreqInc(), reinterpret_cast<void*>(IncreaseBBCounter));
	EE->addGlobalMapping(BBFreqFn, reinterpret_cast<void*>(recordBB));
	EE->addGlobalMapping(RecordLoadFn, reinterpret_cast<void*>(recordLoad));
	EE->addGlobalMapping(RecordStoreFn, reinterpret_cast<void*>(recordStore));

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

	//VerifyProfile();

	//delete newM;
}

ProfilerJITContext::ProfilerJITContext() : P(nullptr), LastBBid(INITIAL_BIGVALUE) {
	// If we have a native target, initialize it to ensure it is linked in and
	// usable by the JIT.
	InitializeNativeTarget();
	InitializeNativeTargetAsmPrinter();
	InitializeNativeTargetAsmParser();
}


ProfilerJITSingletonContext::ProfilerJITSingletonContext(EmbeddedProfilerEngine *P) {
	GlobalContext->P = P;
}

ProfilerJITSingletonContext::~ProfilerJITSingletonContext() {
	GlobalContext->P = nullptr;
	//TODO: Other clearup?
}

void llvm::print_trace_info(TraceType& trace) {
	errs() << "\n\nPrinting Trace Information: \n";
	errs() << "RecordType\t\tlsID\t\tlineID\t\tBB\t\tFunc\t\tAddress\t\tls_Size\n";
	unsigned size = trace.size();
#ifndef USE_LIST_FOR_TRACE_F
	for (unsigned i = 0; i < size; i++) {
		RecordTrace trace_entry = trace.at(i);
		errs() << trace_entry.getRecordType() << "\t";
		errs() << trace_entry.getlsID() << "\t";
		errs() << trace_entry.getlineID() << "\t";
		errs() << trace_entry.getBelongedBB()->getName() << "\t";
		errs() << trace_entry.getBelongedFunc()->getName() << "\t";
		errs() << trace_entry.getls_address() << "\t";
		errs() << trace_entry.getls_size() << "\n";
	}
#else
	TraceType::iterator it, ie;
	for (it = trace.begin(), ie = trace.end(); it != ie; ++it) {
		RecordTrace trace_entry = *it;
		errs() << trace_entry.getRecordType() << "\t";
		errs() << trace_entry.getlsID() << "\t";
		errs() << trace_entry.getlineID() << "\t";
		errs() << trace_entry.getBelongedBB()->getName() << "\t";
		errs() << trace_entry.getBelongedFunc()->getName() << "\t";
		errs() << trace_entry.getls_address() << "\t";
		errs() << trace_entry.getls_size() << "\n";
	}
#endif // End of USE_LIST_FOR_TRACE_F
	errs() << "Finish Printing Trace Information!\n\n";
}

/*
void llvm::write_trace_file(TraceType& trace, std::string fname) {

	std::ofstream trace_file;
	trace_file.open(fname);
	trace_file << "Eval_Iteration_Index contains iteration number starting from the top-level loop to the innermost-level loop\n";
	trace_file << "\n\nPrinting Trace Information: \n";
	trace_file << "Eval_Iteration_Index\tInstruction\tRecordType\tlsID\tlineID\tBB\tBBid\tFunc\tAddress\tls_Size\n";
	unsigned size = trace.size();
	for (unsigned i = 0; i < size; i++) {
		RecordTrace trace_entry = trace.at(i);

		uint64_t BBid = trace_entry.getBelongedBBid();
		lpBBid2EvalIterIndxMapTy::iterator it;
		it = lpBBid2EvalIterIndxMap.find(BBid);
		if (it != lpBBid2EvalIterIndxMap.end()) {
			/// Print the Evaluation Iteration Index
			std::list<unsigned> eval_iter_list = trace_entry.getls_record().getEval_Iteration_Index();
			std::list<unsigned>::iterator itr_i, itr_e;
			for (itr_i = eval_iter_list.begin(), itr_e = eval_iter_list.end(); itr_i != itr_e; itr_i++) {
				trace_file << *itr_i << "\t";
			}
		}
		else {
			trace_file << "\t\t\t";
		}
		trace_file << "\t";

		trace_file << trace_entry.getRecordType() << "\t\t";
		trace_file << trace_entry.getlsID() << "\t\t";
		trace_file << trace_entry.getlineID() << "\t\t";
		trace_file << trace_entry.getBelongedBB()->getName().str() << "\t\t";
		trace_file << BBid << "\t\t";
		trace_file << trace_entry.getBelongedFunc()->getName().str() << "\t\t";
		trace_file << trace_entry.getls_address() << "\t\t";
		trace_file << trace_entry.getls_size() << "\t\t";
		
		trace_file << "\n";
		
	}
	trace_file << "Finish Printing Trace Information!\n\n";
	trace_file.close();
}
*/

char CodeInstrument::ID = 0;
INITIALIZE_PASS(CodeInstrument, "codeInstrument",
				"This pass is used to instrument the source code",
				false,
				false
				//true // We need to modify the source code
				)

ModulePass* llvm::createCodeInstrumentPass() {
	return new CodeInstrument();
}