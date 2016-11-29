#include "profile_h/AnalysisProfilingPass.h"

#define DEBUG_TYPE "analysis-profiling"

namespace llvm {
	DataDepType DynDataDepVec;
	DataDepType RefinedDynDataDepVec;
	DataDepType RefinedInterDynDataDepVec;
}

using namespace llvm;

static const char *LoopInfoMDNodeName = "Loop_Info_of_this_kernel";
static const char *loopNumberMDNodeName = "Loop_Number_of_kernel";


AnalysisProfiling::AnalysisProfiling() : ModulePass(ID) {
	errs() << "Initialize AnalysisProfiling pass\n";
}

void AnalysisProfiling::getAnalysisUsage(AnalysisUsage &AU) const{
	AU.addRequired<QueryBasicBlockID>();
	AU.addPreserved<QueryBasicBlockID>();

	AU.addRequired<QueryLoadStoreID>();
	AU.addPreserved<QueryLoadStoreID>();

	AU.addRequired<LoopInfo>();
	AU.addPreserved<LoopInfo>();

	AU.setPreservesCFG();
}

bool AnalysisProfiling::runOnModule(Module &M){
	/// Get the QueryBasicBlockID & QueryLoadStoreID passes
	bbIDpass = &getAnalysis<QueryBasicBlockID>();
	lsIDpass = &getAnalysis<QueryLoadStoreID>();

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

	/// Get the names of functions that have loops
	/// This following function is moved to GetLoopBound Pass
#ifndef GETLOOPBOUNDPASS_H
	extract_loopInfo_inFunc(M);
#endif // End of GETLOOPBOUNDPASS_H

	Module::iterator FI, FE;
	std::vector<std::string>::iterator fn_it;
	
	/// Calculate the loop bounds of loops and get loopBBid2loopDepthMap
	Function::iterator BI, BE;
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
#ifndef GETLOOPBOUNDPASS_H
			calculateLoopBound(FI);
#endif // End of GETLOOPBOUNDPASS_H
			getloopBBid2loopDepthMap(FI);
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

#ifndef GETLOOPBOUNDPASS_H
#ifdef ENABLE_TIMER
	Timer modifyTrace("modifyTrace_EvalIterationIndex()");
	modifyTrace.startTimer();
#endif // End of ENABLE_TIMER
	/// Modify the evaluation iteration index in the trace, so that it can be aware of 
	/// loop bounds &&
	/// Refine the Trace_f, so that we can only need to focus on trace_entries of loops
	modifyTrace_EvalIterationIndex();
#ifdef ENABLE_TIMER
	modifyTrace.stopTimer();
#endif // End of ENABLE_TIMER
#endif // End of GETLOOPBOUNDPASS_H

	/*
	/// Trace the evaluation iteration index of loops in functions except 
	/// for the "main" function
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		// We only focus on analyze kernel functions instead of main function
		funcName2loopsMapTy::iterator it;
		it = funcName2loopsMap.find(FI->getName());
		if (it != funcName2loopsMap.end()) {
			if (FI->getName() != "main") {
				trace_eval_iteration_index(FI);
			}
		}
	}*/

#ifdef ENABLE_TIMER
	Timer analyze_Dynamic("analyze_Dynamic_DataDep()");
	analyze_Dynamic.startTimer();
#endif // End of ENABLE_TIMER
	/// Analyze the dynamic data dependence here
	analyze_Dynamic_DataDep();
#ifdef ENABLE_TIMER
	analyze_Dynamic.stopTimer();
#endif // End of ENABLE_TIMER

	/// Analyze the memory intensity here
	analyze_Arithmetic_Intensity();

	// Write information of dynamic data dependence and arithmetic intensity of loops 
	// into a file
	write_dynDataDep_ArithIntensity_file();

	//bool IDpass = testbbIDpass(M);

	//VerifyProfile();

	// Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		return false;
	}

	return false;
}

/// FIXME: This code has a bug
/// Bugs fixed, but we need to be aware that we count loop N is the first loop in the kernel and
/// loop 1 is the last loop in the kernel, maybe we can change this in 'Loop_Info_of_this_kernel' metadata
void AnalysisProfiling::extract_loopInfo_inFunc(Module &M) {
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
		unsigned NumLoopInaFunc = funcName2loopNumMap.find(func_name->getString())->second;

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
				//std::string loop_name = funcName2loopsMap.at(func_name->getString()).at(loop_counter);
				std::string loop_name = funcName2loopsMap.at(func_name->getString()).at(NumLoopInaFunc - loop_counter - 1);
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

}

void AnalysisProfiling::calculateLoopBound(Function* func) {
	errs() << "function: " << func->getName() << "\n";
	LoopInfo* LI = &getAnalysis<LoopInfo>(*func);
	
	std::map<unsigned, uint64_t> depth2BBfreqMap;
	std::map<unsigned, unsigned> depth2LoopBoundMap;
	std::map<std::string, std::map<unsigned, uint64_t> > loopName2depthBBfreqMap;

	unsigned depth = 0;
	uint64_t loopbbid = 0;
	uint64_t bbfreq = 0;
	unsigned loop_counter = 1;
	unsigned num_nestedLevel = 0;
	Function::iterator BI, BE;
	for (BI = func->begin(), BE = func->end(); BI != BE; ++BI) {
		errs() << "Basic Block = " << BI->getName() << "\n";
		// Do not consider basic blocks outside loops
		if (!LI->getLoopFor(BI)) {
			continue;
		}

		// Get loopBBid2loopDepthMap
		depth = LI->getLoopDepth(BI);
		loopbbid = bbIDpass->getBBid(BI);
		loopBBid2loopDepthMap.insert(std::make_pair(loopbbid, depth));

		// Do not consider basic blocks that are not a header of loops
		if (!LI->isLoopHeader(BI)) {
			continue;
		}

		bbfreq = BBidFreq.at(loopbbid);
		depth2BBfreqMap.insert(std::make_pair(depth, bbfreq));

		std::string loop_name = "loop" + std::to_string(loop_counter);
		num_nestedLevel = func2loopInfoMap.at(func->getName()).at(loop_name).size();
		if (depth == num_nestedLevel) {
			// We have explore all header basic blocks of the loop, "loop_name", store 
			// its depth and bbfreq information into loopName2depthBBfreqMap
			loopName2depthBBfreqMap.insert(std::make_pair(loop_name, depth2BBfreqMap));
			depth2BBfreqMap.clear();
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
		for (unsigned j = 1; j < num_nestedLevel+1; j++) {
			if (j==1) {
				bound = loopName2depthBBfreqMap.at(lp_name).at(j) - 1;
			}
			else {
				bound = loopName2depthBBfreqMap.at(lp_name).at(j) / higher_bound;
			}

			higher_bound = loopName2depthBBfreqMap.at(lp_name).at(j);
			// Store the loop bound into func2loopInfoMap
			func2loopInfoMapTy::iterator it = func2loopInfoMap.find(func->getName());
			if (it != func2loopInfoMap.end()) {
				it->second.at(lp_name).find(j)->second = bound;
			}
			loopBoundList.push_back(bound);
		}
		assert(loopBoundList.size() == num_nestedLevel && "Size of loopBoundList is NOT equal to num_nestedLevel!\n");
		loopName2BoundList.insert(std::make_pair("loop"+std::to_string(i+1),loopBoundList));
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
	
}

/// Get loopBBid2loopDepthMap
void AnalysisProfiling::getloopBBid2loopDepthMap(Function* func){
	errs() << "function: " << func->getName() << "\n";
	LoopInfo* LI = &getAnalysis<LoopInfo>(*func);

	unsigned depth = 0;
	uint64_t loopbbid = 0;

	Function::iterator BI, BE;
	for (BI = func->begin(), BE = func->end(); BI != BE; ++BI) {
		errs() << "Basic Block = " << BI->getName() << "\n";
		// Do not consider basic blocks outside loops
		if (!LI->getLoopFor(BI)) {
			continue;
		}

		// Get loopBBid2loopDepthMap
		depth = LI->getLoopDepth(BI);
		loopbbid = bbIDpass->getBBid(BI);
		loopBBid2loopDepthMap.insert(std::make_pair(loopbbid, depth));
	}
}

/// Modify the evaluation iteration index in the trace, so that it can be aware of 
/// loop bounds &&
/// Refine the Trace_f, so that we can only need to focus on trace_entries of loops
void AnalysisProfiling::modifyTrace_EvalIterationIndex() {
	/// Evaluation Iteration Index sequence:  
	/// The top-level -> lower-level -> ... -> the innermost-level loop
	std::string funcName_old = "";
	//loopName2Trace_fMMapTy loopName2Trace_fMMap;

	std::string loopName_old = "";
	TraceType subTrace_f;

	funcLoopNamePairTy funcLoopNamePair;

	TraceType::iterator it_tr, ie_tr;
	for (it_tr = Trace_f.begin(), ie_tr = Trace_f.end(); it_tr != ie_tr; ++it_tr) {
		//RecordTrace trace_entry = *it_tr;
		std::list<unsigned> eval_iter_list = it_tr->getls_record().getEval_Iteration_Index();

		// Figure out this trace_entry belongs to which loop
		Function* func = it_tr->getBelongedFunc();
//#ifdef IGNORE_FIRST
		//LoopInfo* LI = &getAnalysis<LoopInfo>(*func);
//#endif // End of IGNORE_FIRST
		uint64_t BBid = it_tr->getBelongedBBid();

		std::string funcName_new = func->getName();

		if (funcName_old == "") {
			funcName_old = funcName_new;
		}

		if (funcName_old != funcName_new) {
			funcLoopNamePair = std::make_pair(funcName_old, loopName_old);
			funcName2loopTraceMMap.insert(std::make_pair(funcLoopNamePair, subTrace_f));
			subTrace_f.clear();
			loopName_old = "";

			funcName_old = funcName_new;
		}

		// We need to bypass the basic block ids that are not belonged to any loops
		//std::vector<uint64_t>::iterator it_bbid2lpMap;
		//it_bbid2lpMap = std::find(bbID_for_loops_vec.begin(), bbID_for_loops_vec.end(), BBid);
		BBids2loopNameMapTy::iterator it_bbid2lpMap;
		it_bbid2lpMap = BBids2loopNameMap.find(BBid);
		//if (it_bbid2lpMap != bbID_for_loops_vec.end()) {
		
		/// FIXME: This if used to check whether the trace_entry is belonged to a loop can be removed, because
		/// we have already modified visitLoad/StoreInst functions in CodeInstrumentPass to only instrument
		/// load/store instructions inside loops. But we leave it and when we have time, we need to remove
		/// this part
		if (it_bbid2lpMap != BBids2loopNameMap.end()) {
//#ifdef IGNORE_FIRST
			// We don't care about the trace_entries of main function
			//if (funcName_new != "main") {

			std::string loop_name = funcName2loopInfoMap.at(funcName_new).at(BBid);

			if (loopName_old == "") {
				loopName_old = loop_name;
			}

			// Figure out this trace_entry is in which loop level
			//unsigned depth = LI->getLoopFor(it_tr->getBelongedBB())->getLoopDepth();
			//unsigned depth = loopBBid2loopDepthMap.find(BBid)->second;

			// Get the loop bound for this loop level and 
			// Modify the Evaluation Iteration Index in the trace_entry and store it back
			std::list<unsigned> boundList = funcName2loopBounds.at(funcName_new).at(loop_name);
			std::list<unsigned> new_eval_iter_list;
			std::list<unsigned>::iterator it, ie;

			for (it = eval_iter_list.begin(), ie = eval_iter_list.end(); it != ie; ++it) {
				unsigned bound_val = boundList.front();
				unsigned iter_val = *it;
				new_eval_iter_list.push_back(iter_val % bound_val);

				boundList.pop_front();
			}

			// Store the new eval_iter_list into the trace file
			it_tr->setEval_Iteration_Index(new_eval_iter_list);
			new_eval_iter_list.clear();

			/// Refine the Trace_f, so that we can only need to focus on trace_entries of loops
			//funcName2Trace_entryMMap.insert(std::pair<std::string, RecordTrace>(funcName_new, *it_tr));
			//}
			if (loopName_old != loop_name && funcName_old == funcName_new) {
				funcLoopNamePairTy funcLoopNamePair = std::make_pair(funcName_old, loopName_old);
				funcName2loopTraceMMap.insert(std::make_pair(funcLoopNamePair, subTrace_f));
				subTrace_f.clear();
				loopName_old = loop_name;
			}
			subTrace_f.push_back(*it_tr);

//#endif // End of IGNORE_FIRST
		}
		eval_iter_list.clear();
	}
	/// Store the last subTrace_f
	funcLoopNamePair = std::make_pair(funcName_old, loopName_old);
	funcName2loopTraceMMap.insert(std::make_pair(funcLoopNamePair, subTrace_f));
	subTrace_f.clear();
} // End of modifyTrace_EvalIterationIndex()

void AnalysisProfiling::trace_eval_iteration_index(Function* func) {
	/*
	LoopInfo* LI = &getAnalysis<LoopInfo>(*func);

	unsigned size_trace = Trace_f.size();
	for (unsigned i = 0; i < size_trace; i++) {
		RecordTrace trace_entry = Trace_f.at(i);
		BasicBlock* BB = trace_entry.getBelongedBB();
		unsigned lsid = trace_entry.getlsID();
		Instruction* inst = lsIDpass->getlsbylsID(lsid);


	}
	*/
}

void AnalysisProfiling::analyze_Dynamic_DataDep() {
	// recordls    : contains all information of interested load/store instrutions, while
	// Trace_f     : not only contains all information inside recordls, but also includes 
	//               corresponding basic block information, such as pointers and IDs of 
	//               basic blocks.
	// lsAddressVec: only contains load/store addresses, which are information of memory accesses

	// The approach we use to extract dynamic data dependence is called 'Pairwise' Method.
	// This method is not so efficient, but currently we do not care about efficiency.
	// We can improve this part if we want in the future.

	/// FIXME: Assume we have two load/store instructions with the same address in the same iteration,
	///        if there is a load/store instruction in other iteration, when we calculate the dynamic
	///        data dependence, we will have 2 entries. We need to combine these two entries.

	/// We need to check whether the location is intra-loop or inter-loops. If it is inter-loops, then
	/// we just ignore it at this moment, because we only focus on dynamic data dependence within loops
	/// so far. But it will be interesting if program has features like dataflow.

	/// 1. Seperate the trace into a funcName2loopTraceMap
	// This part is done in AnalysisProfiling::modifyTrace_EvalIterationIndex function

	/// 2. Analyze data dependence 'DynDataDepVec' of subTrace of each loop inside interested function
	/*
	unsigned size_func = funcName2loopInfoMap.size();
	std::vector<std::string> funcNameVec;
	for (unsigned i = 0; i < size_func; i++) {
	
	}*/
#ifndef GETLOOPBOUNDPASS_H
	funcName2loopTraceMMapTy::iterator it, ie;
	lsAddrTraceIDpairMapTy lsAddrTraceIDpairVec;
	for (it = funcName2loopTraceMMap.begin(), ie = funcName2loopTraceMMap.end(); it != ie; ++it) {
		std::pair<std::string, std::string> funcloopNamePair = it->first;
		TraceType::iterator it_tr, ie_tr;
		//lsAddrTraceIDpairMapTy lsAddrTraceIDpairVec;
		uint64_t i = 0;
		for (it_tr = it->second.begin(), ie_tr = it->second.end(); it_tr != ie_tr; ++it_tr) {
			i++;
			//std::pair<uintptr_t, uint64_t> temp_pair = std::make_pair(it_tr->getls_address(), it_tr->getTraceID());
			//std::pair<uintptr_t, RecordTrace> temp_pair = std::make_pair(it_tr->getls_address(), *it_tr);
			uintptr_t ls_addr = it_tr->getls_address();
			RecordTrace rt = *it_tr;
			lsAddrTraceIDpairMapTy::iterator temp_it;
			lsAddrTraceIDpairMapTy::iterator temp_ie = lsAddrTraceIDpairVec.end();
			temp_it = findEntryInlsAddrTraceIDpairVec(lsAddrTraceIDpairVec, it_tr->getls_address());

			if (temp_it == lsAddrTraceIDpairVec.end()) {
				//lsAddrTraceIDpairVec.push_back(temp_pair); // Move this statement to the place outside if statement
				//lsAddrTraceIDpairVec.push_back(std::make_pair(ls_addr, rt));
				lsAddrTraceIDpairVec.insert(std::make_pair(ls_addr, rt)); // Move this statement to the place outside if statement
			}
			else
			{
				// We have found an entry with the same address in the history vector, lsAddrTraceIDpairVec.
				DepType depty = DepType::INITType;
				//RecordTrace srcTrace = Trace_f.at(temp_it->second);
				RecordTrace srcTrace = temp_it->second;
				uint64_t srclsid = srcTrace.getlsID();
				uint64_t dstlsid = it_tr->getlsID();
				bool ignore = false;

				if (srcTrace.getRecordType() == "STOREType" && it_tr->getRecordType() == "LOADType") {
					depty = DepType::RAWType; // Flow dependence
				}
				else if (srcTrace.getRecordType() == "LOADType" && it_tr->getRecordType() == "STOREType") {
					depty = DepType::WARType; // Anti dependence
				}
				else if (srcTrace.getRecordType() == "STOREType" && it_tr->getRecordType() == "STOREType") {
					depty = DepType::WAWType; // Output dependence
				}
				else {
					// For Read after Read dependence, we don't care about it.
					ignore = true;
				}
				
				if (ignore == false) {
					fn_lpNamePairTy srcfn_lpNamePair = srcTrace.getfn_lpNamePair();
					bool inter_loop, inter_function;
					if (srcfn_lpNamePair.first == funcloopNamePair.first) {
						inter_function = false;
						if (srcfn_lpNamePair.second == funcloopNamePair.second) {
							inter_loop = false;
						}
						else {
							inter_loop = true;
						}
					}
					else {
						inter_function = true;
						inter_loop = true;
					}

					std::list<unsigned> src_eval_iter_list = srcTrace.getls_record().getEval_Iteration_Index();
					std::list<unsigned> dst_eval_iter_list = it_tr->getls_record().getEval_Iteration_Index();
					unsigned size_src = src_eval_iter_list.size();
					unsigned size_dst = dst_eval_iter_list.size();
					std::vector<int64_t> depDistance;

					if (inter_loop || inter_function) {
						// Do not update depDistance
					}
					else {
						// Update depDistance
						if (size_src > size_dst) {
							for (unsigned i = 0; i < (size_src - size_dst); ++i) {
								dst_eval_iter_list.push_back(INITIAL_BIGVALUE);
							}
						}
						else if (size_src < size_dst) {
							for (unsigned i = 0; i < (size_dst - size_src); ++i) {
								src_eval_iter_list.push_back(INITIAL_BIGVALUE);
							}
						}

						std::list<unsigned>::iterator it_src, ie_src;
						std::list<unsigned>::iterator it_dst = dst_eval_iter_list.begin();
						for (it_src = src_eval_iter_list.begin(), ie_src = src_eval_iter_list.end(); it_src != ie_src; ++it_src) {
							if ((*it_src) == INITIAL_BIGVALUE || (*it_dst == INITIAL_BIGVALUE)) {
								depDistance.push_back(INITIAL_BIGVALUE);
							}
							else {
								depDistance.push_back((*it_dst) - (*it_src));
							}

							it_dst++;
						}

					} // Already updated depDistance
					
					DynDataDep temp_dataDep(depty, srclsid, dstlsid, depDistance, srcfn_lpNamePair, funcloopNamePair, inter_loop, inter_function);
					DynDataDepVec.push_back(temp_dataDep);

					// Need to erase the searched content inside lsAddrTraceIDpairVec
					lsAddrTraceIDpairVec.erase(temp_it);
					//lsAddrTraceIDpairVec.push_back(temp_pair);
					//lsAddrTraceIDpairVec.push_back(std::make_pair(ls_addr, rt));
					lsAddrTraceIDpairVec.insert(std::make_pair(ls_addr, rt));
				}
				else { // ignore = true, RAR dependence
					// Do not erase the searched content inside lsAddrTraceIDpairVec
					// Do not add temp_pair into lsAddrTraceIDpairVec;
				}
			}
			//lsAddrTraceIDpairVec.push_back(temp_pair);
		}
		/// If we want to analyze dependency across loops or functions, we need to comment the following 
		/// statement
		//lsAddrTraceIDpairVec.clear();
	}
#else // Enable TWO_TIME_PROFILING
	lsAddrTraceIDpairMapTy lsAddrTraceIDpairVec;
	TraceType::iterator it_tr, ie_tr;
	uint64_t i = 0;
	for (it_tr = Trace_f.begin(), ie_tr = Trace_f.end(); it_tr != ie_tr; ++it_tr) {
		i++;
		//std::pair<uintptr_t, uint64_t> temp_pair = std::make_pair(it_tr->getls_address(), it_tr->getTraceID());
		//std::pair<uintptr_t, RecordTrace> temp_pair = std::make_pair(it_tr->getls_address(), *it_tr);
		fn_lpNamePairTy funcloopNamePair = it_tr->getfn_lpNamePair();
		uintptr_t ls_addr = it_tr->getls_address();
		RecordTrace rt = *it_tr;
		lsAddrTraceIDpairMapTy::iterator temp_it;
		lsAddrTraceIDpairMapTy::iterator temp_ie = lsAddrTraceIDpairVec.end();
		temp_it = findEntryInlsAddrTraceIDpairVec(lsAddrTraceIDpairVec, it_tr->getls_address());

		if (temp_it == lsAddrTraceIDpairVec.end()) {
			//lsAddrTraceIDpairVec.push_back(temp_pair); // Move this statement to the place outside if statement
			//lsAddrTraceIDpairVec.push_back(std::make_pair(ls_addr, rt));
			lsAddrTraceIDpairVec.insert(std::make_pair(ls_addr, rt)); // Move this statement to the place outside if statement
		}
		else
		{
			// We have found an entry with the same address in the history vector, lsAddrTraceIDpairVec.
			DepType depty = DepType::INITType;
			//RecordTrace srcTrace = Trace_f.at(temp_it->second);
			RecordTrace srcTrace = temp_it->second;
			uint64_t srclsid = srcTrace.getlsID();
			uint64_t dstlsid = it_tr->getlsID();
			bool ignore = false;

			if (srcTrace.getRecordType() == "STOREType" && it_tr->getRecordType() == "LOADType") {
				depty = DepType::RAWType; // Flow dependence
			}
			else if (srcTrace.getRecordType() == "LOADType" && it_tr->getRecordType() == "STOREType") {
				depty = DepType::WARType; // Anti dependence
			}
			else if (srcTrace.getRecordType() == "STOREType" && it_tr->getRecordType() == "STOREType") {
				depty = DepType::WAWType; // Output dependence
			}
			else {
				// For Read after Read dependence, we don't care about it.
				ignore = true;
			}

			if (ignore == false) {
				fn_lpNamePairTy srcfn_lpNamePair = srcTrace.getfn_lpNamePair();
				bool inter_loop, inter_function;
				if (srcfn_lpNamePair.first == funcloopNamePair.first) {
					inter_function = false;
					if (srcfn_lpNamePair.second == funcloopNamePair.second) {
						inter_loop = false;
					}
					else {
						inter_loop = true;
					}
				}
				else {
					inter_function = true;
					inter_loop = true;
				}

				std::list<unsigned> src_eval_iter_list = srcTrace.getls_record().getEval_Iteration_Index();
				std::list<unsigned> dst_eval_iter_list = it_tr->getls_record().getEval_Iteration_Index();
				unsigned size_src = src_eval_iter_list.size();
				unsigned size_dst = dst_eval_iter_list.size();
				std::vector<int64_t> depDistance;

				if (inter_loop || inter_function) {
					// Do not update depDistance
				}
				else {
					// Update depDistance
					if (size_src > size_dst) {
						for (unsigned i = 0; i < (size_src - size_dst); ++i) {
							dst_eval_iter_list.push_back(INITIAL_BIGVALUE);
						}
					}
					else if (size_src < size_dst) {
						for (unsigned i = 0; i < (size_dst - size_src); ++i) {
							src_eval_iter_list.push_back(INITIAL_BIGVALUE);
						}
					}

					std::list<unsigned>::iterator it_src, ie_src;
					std::list<unsigned>::iterator it_dst = dst_eval_iter_list.begin();
					for (it_src = src_eval_iter_list.begin(), ie_src = src_eval_iter_list.end(); it_src != ie_src; ++it_src) {
						if ((*it_src) == INITIAL_BIGVALUE || (*it_dst == INITIAL_BIGVALUE)) {
							depDistance.push_back(INITIAL_BIGVALUE);
						}
						else {
							depDistance.push_back((*it_dst) - (*it_src));
						}

						it_dst++;
					}

				} // Already updated depDistance

				DynDataDep temp_dataDep(depty, srclsid, dstlsid, depDistance, srcfn_lpNamePair, funcloopNamePair, inter_loop, inter_function);
				DynDataDepVec.push_back(temp_dataDep);

				// Need to erase the searched content inside lsAddrTraceIDpairVec
				lsAddrTraceIDpairVec.erase(temp_it);
				//lsAddrTraceIDpairVec.push_back(temp_pair);
				//lsAddrTraceIDpairVec.push_back(std::make_pair(ls_addr, rt));
				lsAddrTraceIDpairVec.insert(std::make_pair(ls_addr, rt));
			}
			else { // ignore = true, RAR dependence
				// Do not erase the searched content inside lsAddrTraceIDpairVec
				// Do not add temp_pair into lsAddrTraceIDpairVec;
			}
		}
		//lsAddrTraceIDpairVec.push_back(temp_pair);
	}
	/// If we want to analyze dependency across loops or functions, we need to comment the following 
	/// statement
	//lsAddrTraceIDpairVec.clear();
#endif // End of #ifndef GETLOOPBOUNDPASS_H

	/*
	unsigned counter = 0;
	RecordTrace trace_entry;
	lsAddressVecTy::iterator it, ie, temp_it;
	for (it = lsAddressVec.begin(), ie = lsAddressVec.end(); it != ie; ++it) {
		uintptr_t temp_addr = *it;
		temp_it = std::find(it+1, lsAddressVec.end(), temp_addr);
		if (temp_it != lsAddressVec.end()) {
			
		}

		counter++;
	}*/

	/// 3. RefinedDynDataDepVec:
	///				Refine DynDataDepVec so that we can only focus on loop-carried and loop-independent
	///				data dependence. 
	///
	///		 RefinedInterDynDataDepVec:
	///				Refine DynDataDepVec so that we can only focus on inter-loop and inter-function data
	///				dependence.

#ifndef USE_LIST_FOR_DATADEPTYPE
	unsigned size_DepVec = DynDataDepVec.size();
	/*
	for (unsigned i = 0; i < size_DepVec; i++) {
		bool isloopIndependent = DynDataDepVec.at(i).isLoopIndependent();
		if (isloopIndependent == false) {
			RefinedDynDataDepVec.push_back(DynDataDepVec.at(i));
		}
	}
	*/
	for (unsigned i = 0; i < size_DepVec; i++) {
		bool isinterloopdependent = DynDataDepVec.at(i).isinterloopdependent();
		bool isinterfuncdependent = DynDataDepVec.at(i).isinterfuncdependent();
		if ( (isinterloopdependent == false) && (isinterfuncdependent == false) ) {
			RefinedDynDataDepVec.push_back(DynDataDepVec.at(i));
		}
		if ( (isinterloopdependent == true) || (isinterfuncdependent == true) ) {
			RefinedInterDynDataDepVec.push_back(DynDataDepVec.at(i));
		}
	}

#else
	DataDepType::iterator it_DDD, ie_DDD;
	/*
	for (it_DDD = DynDataDepVec.begin(), ie_DDD = DynDataDepVec.end(); it_DDD != ie_DDD; it_DDD++) {
		bool isloopIndependent = it_DDD->isLoopIndependent();
		if (isloopIndependent == false) {
			RefinedDynDataDepVec.push_back(*it_DDD);
		}
	}
	*/
	for (it_DDD = DynDataDepVec.begin(), ie_DDD = DynDataDepVec.end(); it_DDD != ie_DDD; it_DDD++) {
		bool isinterloopdependent = it_DDD->isinterloopdependent();
		bool isinterfuncdependent = it_DDD->isinterfuncdependent();
		if ((isinterloopdependent == false) && (isinterfuncdependent == false)) {
			RefinedDynDataDepVec.push_back(*it_DDD);
		}
		if ((isinterloopdependent == true) || (isinterfuncdependent == true)) {
			RefinedInterDynDataDepVec.push_back(*it_DDD);
		}
	}
#endif // End of USE_LIST_FOR_DATADEPTYPE


}

void AnalysisProfiling::analyze_Arithmetic_Intensity() {
	// 1. Calculate number of arithmatic operations for loops inside functions except for the "main" function
	//    The map used to trace binary operators: (functionName, loopName) -> std::vector<BB id> > fnlpNamePair2BinaryOpBBidMap.
	//    After we get the basic block frequency, we can get the total number of binary operations executed by
	//    multiplication.
	//    The map, "fnlpNamePair2BinaryOpBBidMap", is calculated in ExtractLoopInfoPass
	//    BBid frequency information can be obtained from BBidFreqMap, which has the type 'typedef std::map<uint64_t, uint64_t> BBidFreqMap'
	uint64_t binaryOpNum = 0;
	fnlpNamePair2BinaryOpBBidMapTy::iterator it, ie;
	for (it = fnlpNamePair2BinaryOpBBidMap.begin(), ie = fnlpNamePair2BinaryOpBBidMap.end(); it != ie; ++it) {
		std::vector<uint64_t> BinaryOpBBid = it->second;
		unsigned size_BinaryOpBBid = BinaryOpBBid.size();
		for (unsigned i = 0; i < size_BinaryOpBBid; i++) {
			binaryOpNum += BBidFreq.at( BinaryOpBBid.at(i) );
		}

		fnlpName2BinaryOpNumMap.insert(std::make_pair(it->first, binaryOpNum));
		// reset binaryOpNum
		binaryOpNum = 0;
	}

	// 2. Calculate number of transferring bytes for loops inside functions && arithmetic intensity
	//    In this part, we need to use the multimap, 'funcName2loopTraceMMap'.
	//    The data type of 'funcName2loopTraceMMap' is shown below:
	//           (funcionName, loopName) -->  loopName2subTraceMapTy Map
	//           typedef std::multimap<funcLoopNamePairTy, TraceType> funcName2loopTraceMMapTy;
	//           funcName2loopTraceMMapTy funcName2loopTraceMMap;
	//    Since 'funcName2loopTraceMMap' is a multiple map, it might record more than once execution
	//    of a loop in a function. But we only need to calculate trace for only one execution and ignore 
	//    the rest.

	//    Arithmetic intensity: is defined as the number of operations performed per word of memory transferred
	//    Here arithmeticIntensity = number of binary operators / number of memory bytes transferred
	uint64_t memory_size_byte = 0;
	fnlpName2BinaryOpNumMapTy::iterator it_op, ie_op;
	for (it_op = fnlpName2BinaryOpNumMap.begin(), ie_op = fnlpName2BinaryOpNumMap.end(); it_op != ie_op; ++it_op) {
		fn_lpNamePairTy fnlpName = it_op->first;

#ifndef GETLOOPBOUNDPASS_H
		std::pair<funcName2loopTraceMMapTy::iterator, funcName2loopTraceMMapTy::iterator> begin_end;
		begin_end = funcName2loopTraceMMap.equal_range(fnlpName);

		funcName2loopTraceMMapTy::iterator itr = begin_end.first;
		unsigned size_subTrace = itr->second.size();
#ifndef USE_LIST_FOR_TRACE_F
		for (unsigned i = 0; i < size_subTrace; i++) {
			RecordTrace trace_entry = itr->second.at(i);
			memory_size_byte += trace_entry.getls_size();
		}
#else
		TraceType::iterator it_st, ie_st;
		for (it_st = itr->second.begin(), ie_st = itr->second.end(); it_st != ie_st; ++it_st) {
			RecordTrace trace_entry = *it_st;
			memory_size_byte += trace_entry.getls_size();
		}
#endif // End of USE_LIST_FOR_TRACE_F

#else // Enable TWO_TIME_PROFILING

		TraceType::iterator it_tr, ie_tr;
		for (it_tr = Trace_f.begin(), ie_tr = Trace_f.end(); it_tr != ie_tr; ++it_tr) {
			fn_lpNamePairTy fn_lpNamePair = it_tr->getfn_lpNamePair();
			if (fn_lpNamePair == fnlpName) {
				memory_size_byte += it_tr->getls_size();
			}
		}

#endif // End of #ifndef GETLOOPBOUNDPASS_H

		float arithmeticIntensity = (float) it_op->second / (float) memory_size_byte;

		fnlpName2MemByteMap.insert(std::make_pair(fnlpName, memory_size_byte));
		fnlpName2ArithIntensityMap.insert(std::make_pair(fnlpName, arithmeticIntensity));
		// reset memory_size_byte
		memory_size_byte = 0;
	}

}

// Write information of dynamic data dependence and arithmetic intensity of loops into a file
void AnalysisProfiling::write_dynDataDep_ArithIntensity_file() {
	/// We need two map here
	/// 1. For arithmetic intensity:	fnlpName2ArithIntensityMap
	/// 2. For dynamic data dependence: RefinedDynDataDepVec

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
				TraceFilename = std::string(IFN.begin(), IFN.end() - 3) + ".DDDAI";
			}
		}
	}

	std::ofstream trace_file;
	trace_file.open(TraceFilename);
	trace_file << "Dynamic Data Dependence and Arithmetic Intensity report\n\n";
	
	unsigned sizeAI = fnlpName2ArithIntensityMap.size();
	unsigned sizeDataDep = RefinedDynDataDepVec.size();
	unsigned sizeInterDataDep = RefinedInterDynDataDepVec.size();
	fnlpName2ArithIntensityMapTy::iterator it, ie;
	for (it = fnlpName2ArithIntensityMap.begin(), ie = fnlpName2ArithIntensityMap.end(); it != ie; ++it) {
		std::string funcName = it->first.first;
		std::string loopName = it->first.second;
		trace_file << "Function Name:\t" << funcName << "\t Loop Name:\t" << it->first.second << "\n";
		trace_file << "------------------------\n";
		trace_file << "Arithmetic Intensity:\t" << it->second << "\n\n";
		trace_file << "-----------\n";
		trace_file << "Dynamic Data Dependence Information: \n";
		trace_file << "(srcLSid -> dstLSid)\tDepType\t\tDepDistance\tCrossLevel\tLoopCarried\tInterLoop\tInterFunction\n";
#ifndef USE_LIST_FOR_DATADEPTYPE
		/// FIXME: We need to check both srcFunctionName and dstFunctionName (the same with loop name),
		/// But currently we just ignore it, since it will not affect the correctness of our program
		for (unsigned i = 0; i < sizeDataDep; i++) {
			if (funcName != RefinedDynDataDepVec.at(i).getsrcFuncName()) {
				continue;
			}

			if (loopName != RefinedDynDataDepVec.at(i).getsrcLoopName()) {
				continue;
			}
			std::pair<uint64_t, uint64_t> depPair = RefinedDynDataDepVec.at(i).getDepPair();
			trace_file << "(" << depPair.first << " -> " << depPair.second << ")\t\t\t\t";
			trace_file << RefinedDynDataDepVec.at(i).getDepType() << "\t\t";
			trace_file << "( ";
			for (unsigned j = 0; j < RefinedDynDataDepVec.at(i).getDepDistance().size(); j++) {
				int64_t temp = RefinedDynDataDepVec.at(i).getDepDistance().at(j);
				if (temp != INITIAL_BIGVALUE) {
					trace_file << temp << " ";
				}
				else {
					trace_file << " ";
				}
			}
			trace_file << ")\t";

			if (RefinedDynDataDepVec.at(i).isCrossLevel() == true) {
				trace_file << "true\t\t";
			}
			else {
				trace_file << "false\t\t";
			}

			if (RefinedDynDataDepVec.at(i).isloopcarrieddependent() == true) {
				trace_file << "true\t\t";
			}
			else {
				trace_file << "false\t\t";
			}
			
			if (RefinedDynDataDepVec.at(i).isinterloopdependent() == true) {
				trace_file << "true\t\t";
			}
			else {
				trace_file << "false\t\t";
			}

			if (RefinedDynDataDepVec.at(i).isinterfuncdependent() == true) {
				trace_file << "true\n";
			}
			else {
				trace_file << "false\n";
			}

		}
#else // Define USE_LIST_FOR_DATADEPTYPE
		DataDepType::iterator it_d, ie_d;
		for (it_d = RefinedDynDataDepVec.begin(), ie_d = RefinedDynDataDepVec.end(); it_d != ie_d; it_d++) {
			if (funcName != it_d->getsrcFuncName()) {
				continue;
			}

			if (loopName != it_d->getsrcLoopName()) {
				continue;
			}
			std::pair<uint64_t, uint64_t> depPair = it_d->getDepPair();
			trace_file << "(" << depPair.first << " -> " << depPair.second << ")\t\t\t\t";
			trace_file << it_d->getDepType() << "\t\t";
			trace_file << "( ";
			for (unsigned j = 0; j < it_d->getDepDistance().size(); j++) {
				int64_t temp = it_d->getDepDistance().at(j);
				if (temp != INITIAL_BIGVALUE) {
					trace_file << temp << " ";
				}
				else {
					trace_file << " ";
				}
			}
			trace_file << ")\t";

			if (it_d->isCrossLevel() == true) {
				trace_file << "true\t\t";
			}
			else {
				trace_file << "false\t\t";
			}

			if (it_d->isloopcarrieddependent() == true) {
				trace_file << "true\t\t";
			}
			else {
				trace_file << "false\t\t";
			}

			if (it_d->isinterloopdependent() == true) {
				trace_file << "true\t\t";
			}
			else {
				trace_file << "false\t\t";
			}

			if (it_d->isinterfuncdependent() == true) {
				trace_file << "true\n";
			}
			else {
				trace_file << "false\n";
			}
		}
#endif // End of USE_LIST_FOR_DATADEPTYPE
		trace_file << "-----------\n";
		trace_file << "------------------------\n\n";
	}

	trace_file << "Dynamic Data Dependence Information across loops or functions: \n";
	trace_file << "(srcLSid -> dstLSid)\tDepType\t\tSourcePair To DestinationPair\tInterLoop\tInterFunction\n";
#ifndef USE_LIST_FOR_DATADEPTYPE
	for (unsigned i = 0; i < sizeInterDataDep; i++) {
		std::pair<uint64_t, uint64_t> depPair = RefinedInterDynDataDepVec.at(i).getDepPair();
		trace_file << "(" << depPair.first << " -> " << depPair.second << ")\t\t\t\t";
		trace_file << RefinedInterDynDataDepVec.at(i).getDepType() << "\t\t";
		std::string srcfuncName = RefinedInterDynDataDepVec.at(i).getsrcFuncName();
		std::string srcloopName = RefinedInterDynDataDepVec.at(i).getsrcLoopName();
		std::string dstfuncName = RefinedInterDynDataDepVec.at(i).getdstFuncName();
		std::string dstloopName = RefinedInterDynDataDepVec.at(i).getdstLoopName();
		trace_file << "(" << srcfuncName << ", " << srcloopName << ") -> ";
		trace_file << "(" << dstfuncName << ", " << dstloopName << ")\t\t";

		if (RefinedInterDynDataDepVec.at(i).isinterloopdependent() == true) {
			trace_file << "true\t\t";
		}
		else {
			trace_file << "false\t\t";
		}

		if (RefinedInterDynDataDepVec.at(i).isinterfuncdependent() == true) {
			trace_file << "true\n";
		}
		else {
			trace_file << "false\n";
		}
	}
#else // define USE_LIST_FOR_DATADEPTYPE
	DataDepType::iterator it_d, ie_d;
	for (it_d = RefinedInterDynDataDepVec.begin(), ie_d = RefinedInterDynDataDepVec.end(); it_d != ie_d; it_d++) {
		std::pair<uint64_t, uint64_t> depPair = it_d->getDepPair();
		trace_file << "(" << depPair.first << " -> " << depPair.second << ")\t\t\t\t";
		trace_file << it_d->getDepType() << "\t\t";
		std::string srcfuncName = it_d->getsrcFuncName();
		std::string srcloopName = it_d->getsrcLoopName();
		std::string dstfuncName = it_d->getdstFuncName();
		std::string dstloopName = it_d->getdstLoopName();
		trace_file << "(" << srcfuncName << ", " << srcloopName << ") -> ";
		trace_file << "(" << dstfuncName << ", " << dstloopName << ")\t\t";

		if (it_d->isinterloopdependent() == true) {
			trace_file << "true\t\t";
		}
		else {
			trace_file << "false\t\t";
		}

		if (it_d->isinterfuncdependent() == true) {
			trace_file << "true\n";
		}
		else {
			trace_file << "false\n";
		}
	}
#endif // End of USE_LIST_FOR_DATADEPTYPE
	trace_file << "-----------\n";
	trace_file << "------------------------\n\n";

	trace_file << "Finish printing information of dynamic data dependence and arithmetic intensity!\n\n";
	trace_file.close();
}

AnalysisProfiling::iterator_t AnalysisProfiling::findEntryInlsAddrTraceIDpairVec(lsAddrTraceIDpairMapTy &lsAddrTraceIDpairVec, uintptr_t ls_addr) {
	return lsAddrTraceIDpairVec.find(ls_addr);
	/*
	lsAddrTraceIDpairListTy::iterator it, ie;
	for (it = lsAddrTraceIDpairVec.begin(), ie = lsAddrTraceIDpairVec.end(); it != ie; ++it) {
		if (it->first == ls_addr) {
			return it;
		}
	}
	return lsAddrTraceIDpairVec.end();
	*/
}

bool AnalysisProfiling::testbbIDpass(Module &M) {
	// Check bbIDpass & lsIDpass
	errs() << "Test bbIDpass & lsIDpass: \n";
	unsigned check_BBcounter = 0;
	unsigned check_INSTcounter = 0;
	unsigned check_LScounter = 0;

	for (Module::iterator FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		for (Function::iterator BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {

			errs() << "Name of BB: " << BI->getName() << "  ID: " << bbIDpass->getBBid(BI) << "\n";
			//errs() << "Name of BB: " << BI->getName() << "  ID: " << lsIDpass->getBBid(BI) << "\n";
			/*
			if (check_BBcounter != bbIDpass->getBBid(BI)) {
				errs() << "\nError: BB ID mismatches!\n";
				return false;
			}*/
			check_BBcounter++;
		}
	}
	errs() << "Pass!\n"; // End testing bbIDpass & lsIDpass
	return true;
}

void AnalysisProfiling::VerifyProfile() {
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
		switch (it->getlsRecordType()) {
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
		errs() << "BasicBlock:  " << BB->getName();

		if (BB == &BB->getParent()->getEntryBlock())
			continue;

		uint64_t Freq = I->second;

		errs() << "  Frequency: " << Freq << "\n";

		// If a basic block is not executed at all, its frequency will be equal to 0,
		// we just ignore it.
		if (Freq == 0) {
			continue;
		}

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
uint64_t AnalysisProfiling::getBBFreq(BasicBlock* BB) const {
	BBFreqMap::const_iterator I = BBFreq.find(BB);
	return (I == BBFreq.end() ? 0 : I->second);
} // End of CodeInstrument::getBBFreq


uint64_t AnalysisProfiling::getBranchFreq(BasicBlock* Src, BasicBlock* Dest) const {
	BranchFreqMap::const_iterator I = BranchFreq.find(Dest);
	if (I == BranchFreq.end()) {
		return 0;
	}

	const BBFreqMap& CurFreqs = I->second;
	BBFreqMap::const_iterator J = CurFreqs.find(Src);

	return (J == CurFreqs.end() ? 0 : J->second);
} // End of CodeInstrument::getBranchFreq()

/// Write the trace into a file
void llvm::write_trace_file(TraceType& trace, std::string fname) {

	std::ofstream trace_file;
	trace_file.open(fname);
	trace_file << "Eval_Iteration_Index contains iteration number starting from the top-level loop to the innermost-level loop\n";
	trace_file << "\n\nPrinting Trace Information: \n";
	trace_file << "Eval_Iteration_Index\tRecordType\tlsID\tlineID\tBB\tBBid\tFunc\tAddress\tls_Size\n";
	unsigned size = trace.size();
#ifndef USE_LIST_FOR_TRACE_F
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
#else
	TraceType::iterator it, ie;
	for (it = trace.begin(), ie = trace.end(); it != ie; ++it) {
		RecordTrace trace_entry = *it;

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
#endif // End of USE_LIST_FOR_TRACE_F
	trace_file << "Finish Printing Trace Information!\n\n";
	trace_file.close();
}

char AnalysisProfiling::ID = 0;
INITIALIZE_PASS(AnalysisProfiling, "Analyze the profiling trace",
				"This pass is used to analyze the profiling trace",
				false,
				false
				)

ModulePass *llvm::createAnalysisProfilingPass() {
	return new 	AnalysisProfiling();
}