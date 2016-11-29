#include "profile_h/InstrumentForDDDGPass.h"
#include "profile_h/auxiliary.h"

#define DEBUG_TYPE "instrument-code-for-building-dddg"

staticInstID2OpcodeMapTy staticInstID2OpcodeMap;
instName2bbNameMapTy instName2bbNameMap;
headerBBFuncNamePair2lastInstMapTy headerBBFuncNamePair2lastInstMap;
headerBBFuncNamePair2lastInstMapTy exitingBBFuncNamePair2lastInstMap;
loopName2levelUnrollVecMapTy loopName2levelUnrollVecMap;
std::ofstream summary;

using namespace llvm;

char list_of_intrinsics[NUM_OF_INTRINSICS][25] = {
	"llvm.memcpy",	// standard C lib
	"llvm.memmove",
	"llvm.memset",
	"llvm.sqrt",
	"llvm.powi",
	"llvm.sin",
	"llvm.cos",
	"llvm.pow",
	"llvm.exp",
	"llvm.exp2",
	"llvm.log",
	"llvm.log10",
	"llvm.log2",
	"llvm.fma",
	"llvm.fabs",
	"llvm.copysign",
	"llvm.floor",
	"llvm.ceil",
	"llvm.trunc",
	"llvm.rint",
	"llvm.nearbyint",
	"llvm.round",
	"llvm.bswap",	//bit manipulation
	"llvm.ctpop",
	"llvm.ctlz",
	"llvm.cttz",
	"llvm.sadd.with.overflow",	//arithmetic with overflow
	"llvm.uadd.with.overflow",
	"llvm.ssub.with.overflow",
	"llvm.usub.with.overflow",
	"llvm.smul.with.overflow",
	"llvm.umul.with.overflow",
	"llvm.fmuladd",		//specialised arithmetic
	"dmaLoad",
	"dmaStore",
};

InstrumentForDDDG::InstrumentForDDDG() : ModulePass(ID) {
	DEBUG(dbgs() << "-------------------\n");
	DEBUG(dbgs() << "\n\tInitialize InstrumentForDDDG pass\n");
	staticInstID2OpcodeMap.clear();
	instName2bbNameMap.clear();
	initializeInstrumentForDDDGPass(*PassRegistry::getPassRegistry());
}

void InstrumentForDDDG::getAnalysisUsage(AnalysisUsage &AU) const {
	AU.setPreservesCFG();
}

bool InstrumentForDDDG::doInitialization(Module &M) {

	// Add external trace_logger function declaratio
	TL_log0 = cast<Function>(M.getOrInsertFunction("trace_logger_log0", Type::getVoidTy(M.getContext()),
																									Type::getInt64Ty(M.getContext()),
																									Type::getInt8PtrTy((M.getContext())),
																									Type::getInt8PtrTy((M.getContext())),
																									Type::getInt8PtrTy((M.getContext())),
																									Type::getInt64Ty(M.getContext()),
																									NULL));

	TL_log_int = cast<Function>(M.getOrInsertFunction("trace_logger_log_int", Type::getVoidTy(M.getContext()),
																										Type::getInt64Ty(M.getContext()),
																										Type::getInt64Ty(M.getContext()),
																										Type::getInt64Ty(M.getContext()),
																										Type::getInt64Ty(M.getContext()),
																										Type::getInt8PtrTy((M.getContext())),
																										NULL));

	TL_log_double = cast<Function>(M.getOrInsertFunction("trace_logger_log_double", Type::getVoidTy(M.getContext()),
																												Type::getInt64Ty(M.getContext()),
																												Type::getInt64Ty(M.getContext()),
																												Type::getDoubleTy(M.getContext()),
																												Type::getInt64Ty(M.getContext()),
																												Type::getInt8PtrTy((M.getContext())),
																												NULL));

	TL_log_int_noreg = cast<Function>(M.getOrInsertFunction("trace_logger_log_int_noreg", Type::getVoidTy(M.getContext()),
																													Type::getInt64Ty(M.getContext()),
																													Type::getInt64Ty(M.getContext()),
																													Type::getInt64Ty(M.getContext()),
																													Type::getInt64Ty(M.getContext()),
																													NULL));

	TL_log_double_noreg = cast<Function>(M.getOrInsertFunction("trace_logger_log_double_noreg", Type::getVoidTy(M.getContext()),
																															Type::getInt64Ty(M.getContext()),
																															Type::getInt64Ty(M.getContext()),
																															Type::getDoubleTy(M.getContext()),
																															Type::getInt64Ty(M.getContext()),
																															NULL));

	char func_string[256] = "";
	//func_string = getenv("WORKLOAD");
	if (kernel_names.empty()) {
		assert(false && "ERROR: Please set kernel_names!\n");
	}

	char temp_str[128] = "";
	unsigned size_ker_name = kernel_names.size();
	for (unsigned i = 0; i < size_ker_name; i++) {
		strcpy(temp_str, kernel_names.at(i).c_str());
		if (i == size_ker_name-1) {
			strcat( func_string, temp_str );
		}
		else {
			strcat( func_string, strcat(temp_str, ",") );
		}
	}

	if (func_string == NULL)
	{
		//errs() << "Please set WORKLOAD as an environment variable!\n";
		return false;
	}
	functions = str_split(func_string, ',', &num_of_functions);

	for (unsigned i = 0; i < num_of_functions; i++) {
		DEBUG(dbgs() << "Function " << i << " = " << *(functions + i) << "\n");
	}

	st = createSlotTracker(&M);
	st->initialize();
	curr_module = &M;
	curr_function = NULL;
	return false;
}

char ** InstrumentForDDDG::str_split(char *a_str, const char a_delim, unsigned *size) {
	int count = 0;
	char *tmp = a_str;
	char *last_comma = 0;
	char delim[2];
	delim[0] = a_delim;
	delim[1] = 0;

	while (*tmp) {
		if (a_delim == *tmp) {
			count++;
			last_comma = tmp;
		}
		tmp++;
	}
	count++;

	char **result;
	result = (char **)malloc(sizeof(char *)* count);
	if (result) {
		int idx = 0;
		char * token = strtok(a_str, delim);
		while (token) {
			assert(idx < count);
			*(result + idx) = strdup(token);
			idx++;
			token = strtok(0, delim);
		}
	}
	*size = count;
	return result;
}

int InstrumentForDDDG::trace_or_not(char* func) {
	if (is_tracking_function(func)) {
		return 1;
	}
		
	for (int i = 0; i < NUM_OF_INTRINSICS; i++) {
		if (strstr(func, list_of_intrinsics[i]) == func) {
			// TODO: Super hacky way of ensuring that dmaLoad and dmaStore always
			// get tracked (by adding them as llvm intrinsics). We should come up
			// with a better way of doing this...
			if (i < NUM_OF_LLVM_INTRINSICS) {
				return i + 2;
			}
			else {
				return 1;
			}
		}
	}

	return -1;
}

bool InstrumentForDDDG::is_tracking_function(string func) {
	for (unsigned i = 0; i < num_of_functions; i++) {
		if (strcmp(*(functions + i), func.c_str()) == 0) {
			return true;
		}
	}
	return false;
}

int InstrumentForDDDG::getMemSize(Type *T) {
	int size = 0;
	if (T->isPointerTy()) {
		return 8 * 8;
	}
	//return getMemSize(T->getPointerElementType());
	else if (T->isFunctionTy()) {
		size = 0;
	}
	else if (T->isLabelTy()) {
		size = 0;
	}
	else if (T->isStructTy())	{
		StructType *S = dyn_cast<StructType>(T);
		for (unsigned i = 0; i != S->getNumElements(); i++)	{
			Type *t = S->getElementType(i);
			size += getMemSize(t);
		}
	}
	else if (T->isFloatingPointTy()) {
		switch (T->getTypeID())	{
		case llvm::Type::HalfTyID:        ///<  1: 16-bit floating point typ
			size = 16; break;
		case llvm::Type::FloatTyID:       ///<  2: 32-bit floating point type
			size = 4 * 8; break;
		case llvm::Type::DoubleTyID:      ///<  3: 64-bit floating point type
			size = 8 * 8; break;
		case llvm::Type::X86_FP80TyID:    ///<  4: 80-bit floating point type (X87)
			size = 10 * 8; break;
		case llvm::Type::FP128TyID:       ///<  5: 128-bit floating point type (112-bit mantissa)
			size = 16 * 8; break;
		case llvm::Type::PPC_FP128TyID:   ///<  6: 128-bit floating point type (two 64-bits, PowerPC)
			size = 16 * 8; break;
		default:
			fprintf(stderr, "!!Unknown floating point type size\n");
			assert(false && "Unknown floating point type size");
		}
	}
	else if (T->isIntegerTy()) {
		size = cast<IntegerType>(T)->getBitWidth();
	}
	else if (T->isVectorTy()) {
		size = cast<VectorType>(T)->getBitWidth();
	}
	else if (T->isArrayTy()) {
		ArrayType *A = dyn_cast<ArrayType>(T);
		size = (int)A->getNumElements()* A->getElementType()->getPrimitiveSizeInBits();
	}
	else {
		fprintf(stderr, "!!Unknown data type: %d\n", T->getTypeID());
		assert(false && "Unknown data type");
	}

	return size;
}

/// Function used to instrument LLVM-IR
void InstrumentForDDDG::print_line(BasicBlock::iterator itr, int line, int line_number, char *func_or_reg_id,
																	 char *bbID, char *instID, int opty, int datasize, Value *value,
																	 bool is_reg) {
	CallInst *tl_call;
	IRBuilder<> IRB(itr);

	Value *v_line, *v_opty, *v_value, *v_linenumber;
	v_line = ConstantInt::get(IRB.getInt64Ty(), line);
	v_opty = ConstantInt::get(IRB.getInt64Ty(), opty);
	v_linenumber = ConstantInt::get(IRB.getInt64Ty(), line_number);

	//print line 0
	if (line == 0) {
		Constant *v_func_id = ConstantDataArray::getString(curr_module->getContext(), func_or_reg_id, true);
		ArrayType* ArrayTy_0 = ArrayType::get(IntegerType::get(curr_module->getContext(), 8), (strlen(func_or_reg_id) + 1));
		GlobalVariable *gvar_array = new GlobalVariable(*curr_module, ArrayTy_0,
			true, GlobalValue::PrivateLinkage, 0, ".str");
		gvar_array->setInitializer(v_func_id);
		std::vector<Constant*> indices;
		ConstantInt *zero = ConstantInt::get(curr_module->getContext(), APInt(32, StringRef("0"), 10));
		indices.push_back(zero);
		indices.push_back(zero);
		Constant * vv_func_id = ConstantExpr::getGetElementPtr(gvar_array, indices);

		Constant *v_bb = ConstantDataArray::getString(curr_module->getContext(), bbID, true);
		ArrayType * ArrayTy_bb = ArrayType::get(IntegerType::get(curr_module->getContext(), 8), (strlen(bbID) + 1));
		GlobalVariable *gvar_array_bb = new GlobalVariable(*curr_module, ArrayTy_bb,
																											  true, GlobalValue::PrivateLinkage, 0, ".str");
		gvar_array_bb->setInitializer(v_bb);
		ConstantInt *zero_bb = ConstantInt::get(curr_module->getContext(), APInt(32, StringRef("0"), 10));
		std::vector<Constant*> indices_bb;
		indices_bb.push_back(zero_bb);
		indices_bb.push_back(zero_bb);
		Constant * vv_bb = ConstantExpr::getGetElementPtr(gvar_array_bb, indices_bb);

		Constant *v_inst = ConstantDataArray::getString(curr_module->getContext(), instID, true);
		ArrayType *ArrayTy_instid = ArrayType::get(IntegerType::get(curr_module->getContext(), 8), (strlen(instID) + 1));
		GlobalVariable *gvar_array_instid = new GlobalVariable(*curr_module, ArrayTy_instid,
																														true, GlobalValue::PrivateLinkage, 0, ".str");
		gvar_array_instid->setInitializer(v_inst);
		std::vector<Constant*> indices_instid;
		ConstantInt *zero_instid = ConstantInt::get(curr_module->getContext(), APInt(32, StringRef("0"), 10));
		indices_instid.push_back(zero_instid);
		indices_instid.push_back(zero_instid);
		Constant * vv_inst = ConstantExpr::getGetElementPtr(gvar_array_instid, indices_instid);
		tl_call = IRB.CreateCall5(TL_log0, v_linenumber, vv_func_id, vv_bb, vv_inst, v_opty);

		// record instruction id into staticInstID2OpcodeMap
		insertInstid(instID, (unsigned) opty);
		insertInstid2bbName(instID, bbID);
		if ((unsigned) opty == LLVM_IR_Br) {
			bbFuncNamePairTy bbfnName = std::make_pair(bbID, func_or_reg_id);
			bbFuncNamePair2lpNameLevelPairMapTy::iterator it_header = headerBBFuncnamePair2lpNameLevelPairMap.find(bbfnName);
			headerBBFuncNamePair2lastInstMapTy::iterator it_lastInst = headerBBFuncNamePair2lastInstMap.find(bbfnName);
			bool foundInheaderBBFuncnamePair2lpNameLevelPairMap = it_header != headerBBFuncnamePair2lpNameLevelPairMap.end();
			bool foundInheaderBBFuncNamePair2firstInstMap = it_lastInst != headerBBFuncNamePair2lastInstMap.end();
			if (foundInheaderBBFuncnamePair2lpNameLevelPairMap && !foundInheaderBBFuncNamePair2firstInstMap) {
				headerBBFuncNamePair2lastInstMap.insert(std::make_pair(bbfnName, instID));
			}

			bbFuncNamePair2lpNameLevelPairMapTy::iterator it_exiting = exitBBFuncnamePair2lpNameLevelPairMap.find(bbfnName);
			headerBBFuncNamePair2lastInstMapTy::iterator it_ltIt_exiting = exitingBBFuncNamePair2lastInstMap.find(bbfnName);
			bool foundInexitBBFuncnamePair2lpNameLevelPairMap = it_exiting != exitBBFuncnamePair2lpNameLevelPairMap.end();
			bool foundInexitingBBFuncNamePair2lastInstMap = it_ltIt_exiting != exitingBBFuncNamePair2lastInstMap.end();
			if (foundInexitBBFuncnamePair2lpNameLevelPairMap && !foundInexitingBBFuncNamePair2lastInstMap) {
				exitingBBFuncNamePair2lastInstMap.insert(std::make_pair(bbfnName, instID));
			}
		}

	}
	//print line with reg
	else {
		Value *v_size;
		v_size = ConstantInt::get(IRB.getInt64Ty(), datasize);
		Value *v_is_reg;
		v_is_reg = ConstantInt::get(IRB.getInt64Ty(), is_reg);

		//if (func_or_reg_id != NULL)
		if (is_reg) {
			assert(func_or_reg_id != NULL);
			Constant *v_reg_id = ConstantDataArray::getString(curr_module->getContext(), func_or_reg_id, true);
			ArrayType* ArrayTy_0 = ArrayType::get(IntegerType::get(curr_module->getContext(), 8),
																						(strlen(func_or_reg_id) + 1));
			GlobalVariable *gvar_array = new GlobalVariable(*curr_module, ArrayTy_0,
																											true, GlobalValue::PrivateLinkage, 0, ".str");
			gvar_array->setInitializer(v_reg_id);
			std::vector<Constant*> indices;
			ConstantInt *zero = ConstantInt::get(curr_module->getContext(), APInt(32, StringRef("0"), 10));
			indices.push_back(zero);
			indices.push_back(zero);
			Constant * vv_reg_id = ConstantExpr::getGetElementPtr(gvar_array, indices);

			if (value != NULL) {
				if (opty == llvm::Type::IntegerTyID) {
					v_value = IRB.CreateZExt(value, IRB.getInt64Ty());
					tl_call = IRB.CreateCall5(TL_log_int, v_line, v_size, v_value, v_is_reg, vv_reg_id);
				}
				else if (opty >= llvm::Type::HalfTyID &&opty <= llvm::Type::PPC_FP128TyID) {
					v_value = IRB.CreateFPExt(value, IRB.getDoubleTy());
					tl_call = IRB.CreateCall5(TL_log_double, v_line, v_size, v_value, v_is_reg, vv_reg_id);
				}
				// deal with functions individually
				else if (opty == llvm::Type::PointerTyID) {
					v_value = IRB.CreatePtrToInt(value, IRB.getInt64Ty());
					tl_call = IRB.CreateCall5(TL_log_int, v_line, v_size, v_value, v_is_reg, vv_reg_id);
				}
				else {
					fprintf(stderr, "normal data else: %d, %s\n", opty, func_or_reg_id);
				}
			}
			//else if (value == NULL &&  bbID != NULL && strcmp(bbID, "phi") == 0 )
			//{
			//v_value = ConstantInt::get(IRB.getInt64Ty(), 999);
			//tl_call = IRB.CreateCall5(TL_log_int, v_line, v_size, v_value, v_is_reg, vv_reg_id);
			//}
			else {
				v_value = ConstantInt::get(IRB.getInt64Ty(), 0);
				tl_call = IRB.CreateCall5(TL_log_int, v_line, v_size, v_value, v_is_reg, vv_reg_id);
			}
			//else
			//fprintf(stderr, "normal data else: %d, %s\n",opty, func_or_reg_id);
		} 
		// is_reg = 0
		else {
			if (value != NULL) {
				if (opty == llvm::Type::IntegerTyID) {
					v_value = IRB.CreateZExt(value, IRB.getInt64Ty());
					tl_call = IRB.CreateCall4(TL_log_int_noreg, v_line, v_size, v_value, v_is_reg);
				}
				else if (opty >= llvm::Type::HalfTyID &&opty <= llvm::Type::PPC_FP128TyID) {
					v_value = IRB.CreateFPExt(value, IRB.getDoubleTy());
					tl_call = IRB.CreateCall4(TL_log_double_noreg, v_line, v_size, v_value, v_is_reg);
				}
				// deal with functions individually
				else if (opty == llvm::Type::PointerTyID) {
					v_value = IRB.CreatePtrToInt(value, IRB.getInt64Ty());
					tl_call = IRB.CreateCall4(TL_log_int_noreg, v_line, v_size, v_value, v_is_reg);
				}
				else {
					fprintf(stderr, "value not empty, normal data else: %d\n", opty);
				}
			}
			//else if (value == NULL &&  bbID != NULL && strcmp(bbID, "phi") == 0 )
			//{
			//v_value = ConstantInt::get(IRB.getInt64Ty(), 999);
			//tl_call = IRB.CreateCall4(TL_log_int_noreg, v_line, v_size, v_value, v_is_reg);
			//}
			else {
				v_value = ConstantInt::get(IRB.getInt64Ty(), 0);
				tl_call = IRB.CreateCall4(TL_log_int_noreg, v_line, v_size, v_value, v_is_reg);
			}
			//fprintf(stderr, "normal data else: %d\n",opty);
		}
	}
} // End of print_line(...)

void InstrumentForDDDG::insertInstid(std::string inst_id, unsigned op_code) {
	staticInstID2OpcodeMap.insert(std::make_pair(inst_id, op_code));
}

void InstrumentForDDDG::insertInstid2bbName(std::string inst_id, std::string bbName){
	instName2bbNameMap.insert(std::make_pair(inst_id, bbName));
}

bool InstrumentForDDDG::getInstId(Instruction *itr, char* bbid, char* instid, int &instc) {
	int id = st->getLocalSlot(itr);
	bool f = itr->hasName();
	if (f) {
		strcpy(instid, (char*)itr->getName().str().c_str());
		return true;
	}

	if (!f && id >= 0) {
		sprintf(instid, "%d", id);
		return true;
	}
	else if (!f && id == -1) {
		char tmp[10];
		char dash[5] = "-";
		sprintf(tmp, "%d", instc);
		if (bbid != NULL) {
			strcpy(instid, bbid);
		}
		strcat(instid, dash);
		strcat(instid, tmp);
		instc++;
		return true;
	}

	return false;
}

void InstrumentForDDDG::getBBId(Value *BB, char *bbid) {
	int id;
	id = st->getLocalSlot(BB);
	bool hasName = BB->hasName();
	if (hasName) {
		strcpy(bbid, (char*)BB->getName().str().c_str());
	}

	if (!hasName && id >= 0) {
		sprintf(bbid, "%d", id);
	}
	else if (!hasName && id == -1) {
		fprintf(stderr, "!!This block does not have a name or a ID!\n");
	}
}

void InstrumentForDDDG::extract_memory_trace_for_access_pattern() {

	std::string file_name = inputPath + "mem_trace.txt";
	std::ofstream mem_trace;
	mem_trace.open(file_name);
	if (mem_trace.is_open()) {
		std::cout << "DEBUG-INFO: [Mem-trace] Start" << std::endl;
	}
	else {
		assert(false && "Error: Cannot open mem_trace file!\n");
	}

	//FILE *tracefile;
	gzFile tracefile;
	//tracefile = fopen(trace_name.c_str(), "r");
	std::string tracefile_name = inputPath + "dynamic_trace.gz";
	tracefile = gzopen(tracefile_name.c_str(), "r");
	if (tracefile == Z_NULL) {
		std::string err_str = "Error! gzfile " + tracefile_name + " can not open!";
		assert(false && err_str.c_str());
	}

	char buffer[256];
	char curr_static_function[256];
	char instid[256], bblockid[256];
	int microop;
	int line_num;
	int count;
	bool trace_entry = false;
	while (!gzeof(tracefile)) {
		if (gzgets(tracefile, buffer, sizeof(buffer)) == Z_NULL)
			continue;
		std::string wholeline(buffer);
		size_t pos_end_tag = wholeline.find(",");

		if (pos_end_tag == std::string::npos) {
			continue;
		}

		std::string tag = wholeline.substr(0, pos_end_tag);
		std::string line_left = wholeline.substr(pos_end_tag + 1);
		if ((trace_entry == false) && (tag.compare("0") == 0)){
			//parse_instruction_line(line_left);
			sscanf(line_left.c_str(), "%d,%[^,],%[^,],%[^,],%d,%d\n", &line_num, curr_static_function, bblockid, instid, &microop, &count);
			std::string func_name(curr_static_function);
			std::string bb_name(bblockid);
			std::string inst_name(instid);
			// Get loop name
			if (microop == 27 || microop == 28) {
				// Only consider load and store
				trace_entry = true;
				bbFuncNamePair2lpNameLevelPairMapTy::iterator it_found = bbFuncNamePair2lpNameLevelPairMap.find(std::make_pair(bb_name, func_name));
				if (it_found != bbFuncNamePair2lpNameLevelPairMap.end()) {
					llvm::lpNameLevelPairTy lpNamelpLevelPair = it_found->second;
					std::string loop_name = lpNamelpLevelPair.first;
					std::string whole_loop_name = loop_name + "-" + to_string(lpNamelpLevelPair.second);
					unsigned int num_levels = LpName2numLevelMap.at(loop_name);
					mem_trace << whole_loop_name << "," << num_levels << "," << instid << ",";
					if (microop == 27) {
						// Load instruction
						mem_trace << "load,";
					}
					else if (microop == 28) {
						mem_trace << "store,";
					}
					else {
						// Do nothing
					}
					mem_trace << count << ",";
				}
			}

		}
		else if ((trace_entry == true) && (((tag.compare("1") == 0) && microop == 27) || ((tag.compare("2") == 0) && microop == 28))) {
			// Load instruction for obtaining load/store addresses
			unsigned int tmp1_v = 0;
			unsigned int tmp2_v = 0;
			unsigned long int addr_v = 0;
			char predecessor_name[256];
			sscanf(line_left.c_str(), "%d,%ld,%d,%s\n", &tmp1_v, &addr_v, &tmp2_v, predecessor_name);
			mem_trace << addr_v << ",";
		}
		else if ((trace_entry == true) && (((tag.compare("r") == 0) && microop == 27) || ((tag.compare("1") == 0) && microop == 28))) {
			// Load instruction for obtaining load/store values
			unsigned int tmp1_v = 0;
			unsigned int tmp2_v = 0;
			float its_value = 0.0f;
			char its_name[256];
			sscanf(line_left.c_str(), "%d,%f,%d,%s\n", &tmp1_v, &its_value, &tmp2_v, its_name);
			mem_trace << its_value << std::endl;
			trace_entry = false;
		}
		else {
			// Do nothing here
		}

	}

	gzclose(tracefile);
	mem_trace.close();
	std::cout << "DEBUG-INFO: [Mem-trace] Finished" << std::endl;
}

bool InstrumentForDDDG::runOnModule(Module &M) {
	errs() << "DEBUG-INFO: [instrumentation_instrument-pass] Start instrumentation\n";
	bool result;
	Module::iterator FI, FE;
	Function::iterator BI, BE;
	std::vector<std::string>::iterator it_name;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		it_name = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		if (it_name != kernel_names.end()) {
			for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
				result = performOnBasicBlock(*BI);
			}
		}
	}
	if (enable_no_trace == false) {
		VERBOSE_PRINT(errs() << "DEBUG-INFO: [profiling_profiling-engine] Start profiling engine\n");
		/// Integrate JIT profiling engine and run the embedded profiler
		ProfilingEngine P(M, TL_log0, TL_log_int, TL_log_double, TL_log_int_noreg, TL_log_double_noreg);
		P.runOnProfiler();
		/// Finished Profiling
		VERBOSE_PRINT(errs() << "DEBUG-INFO: [profiling_profiling-engine] Finished\n");
		if (memory_trace_gen == true) {
			extract_memory_trace_for_access_pattern();
		}
	}
	else {
		errs() << "DEBUG-INFO: [profiling_profiling-engine] We already have generated dynamic trace for this application, no need to generate again.\n";
	}

	if (enable_profiling_time_only == true) {
		errs() << "DEBUG-INFO: [profiling_profiling-engine] Only get profiling time without estimating FPGA execution time\n";
		return result;
	}

	// Verify the module
	std::string ErrorStr;
	raw_string_ostream OS(ErrorStr);
	if (verifyModule(M, &OS)){
		errs() << OS.str() << "\n";
		assert(false && "DEBUG-INFO: [instrumentation_instrument-pass] Module is broken!\n");
	}

	loopBasedTraceAnalysis();

	VERBOSE_PRINT(errs() << "DEBUG-INFO: [instrumentation_instrument-pass] Finished\n");
	return result;
}

bool InstrumentForDDDG::performOnBasicBlock(BasicBlock &BB) {
	Function *F = BB.getParent();
	int instc = 0;
	char funcName[256];

	if (curr_function != F) {
		st->purgeFunction();
		st->incorporateFunction(F);
		curr_function = F;
	}
	strcpy(funcName, curr_function->getName().str().c_str());
	if (!is_tracking_function(funcName)) {
		return false;
	}

	//cout << "Tracking function: " << funcName << endl;
	
	//deal with phi nodes
	BasicBlock::iterator insertp = BB.getFirstInsertionPt();
	BasicBlock::iterator itr = BB.begin();
	if (dyn_cast<PHINode>(itr)) {
		for (; PHINode *N = dyn_cast<PHINode>(itr); itr++) {
			//errs() << "\tTracking instruction: " << *itr << "\n";
			//errs() << "\t\tin function: " << funcName << "\n";
			Value *curr_operand = NULL;
			bool is_reg = 0;
			int size = 0, opcode;
			char bbid[256], instid[256];
			char operR[256];
			//int DataSize, value;

			int line_number = -1;

			getBBId(&BB, bbid);
			getInstId(itr, bbid, instid, instc);
			opcode = itr->getOpcode();

			if (MDNode *N = itr->getMetadata("dbg")) {
				DILocation Loc(N);                      // DILocation is in DebugInfo.h
				line_number = Loc.getLineNumber();
			}
			print_line(insertp, 0, line_number, funcName, bbid, instid, opcode);

			//for instructions using registers
			int i, num_of_operands = itr->getNumOperands();

			if (num_of_operands > 0) {
				char phi[5] = "phi";
				for (i = num_of_operands - 1; i >= 0; i--) {
					curr_operand = itr->getOperand(i);
					//is_reg = 0;
					is_reg = curr_operand->hasName();
					if (Instruction *I = dyn_cast<Instruction>(curr_operand)) {
						int flag = 0;
						is_reg = getInstId(I, NULL, operR, flag);
						assert(flag == 0);
						if (curr_operand->getType()->isVectorTy()) {
							print_line(insertp, i + 1, -1, operR, phi, NULL,
												 curr_operand->getType()->getTypeID(),
												 getMemSize(curr_operand->getType()),
												 NULL,
												 is_reg);
						}
						else {
							print_line(insertp, i + 1, -1, operR, phi, NULL,
												 I->getType()->getTypeID(),
												 getMemSize(I->getType()),
												 NULL,
												 is_reg);
						}
					}
					else if (curr_operand->getType()->isVectorTy()) {
						char operand_id[256];
						strcpy(operand_id, curr_operand->getName().str().c_str());
						print_line(insertp, i + 1, -1, operand_id, phi, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 is_reg);
					}
					else{
						char operand_id[256];
						strcpy(operand_id, curr_operand->getName().str().c_str());
						print_line(insertp, i + 1, -1, operand_id, phi, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 curr_operand,
											 is_reg);
					}
				}
			}

			if (!itr->getType()->isVoidTy()) {
				is_reg = 1;
				if (itr->getType()->isVectorTy()) {
					print_line(insertp, RESULT_LINE, -1, instid, NULL, NULL,
										 itr->getType()->getTypeID(),
										 getMemSize(itr->getType()),
											NULL, is_reg);
				}
				else if (itr->isTerminator()) {
					fprintf(stderr, "It is terminator...\n");
				}
				else {
					print_line(insertp, RESULT_LINE, -1, instid, NULL, NULL,
										 itr->getType()->getTypeID(),
										 getMemSize(itr->getType()),
										 itr, is_reg);
				}
			}
		}
	}

	//for ALL instructions
	BasicBlock::iterator nextitr;
	for (BasicBlock::iterator itr = insertp; itr != BB.end(); itr = nextitr) {
		//errs() << "\tTracking instruction: " << *itr << "\n";
		//errs() << "\t\tin function: " << funcName << "\n";
		Value *curr_operand = NULL;
		bool is_reg = 0;
		int size = 0, opcode;
		char bbid[256], instid[256];
		char operR[256];
		//int DataSize, value;
		int line_number = -1;

		nextitr = itr;
		nextitr++;

		//Get static BasicBlock ID: produce bbid
		getBBId(&BB, bbid);
		//Get static instruction ID: produce instid
		getInstId(itr, bbid, instid, instc);

		//Get opcode: produce opcode
		opcode = itr->getOpcode();

		if (MDNode *N = itr->getMetadata("dbg")) {
			DILocation Loc(N);                      // DILocation is in DebugInfo.h
			line_number = Loc.getLineNumber();
		}
		int callType = -1;
		if (CallInst *I = dyn_cast<CallInst>(itr)) {
			char callfunc[256];
			Function *fun = I->getCalledFunction();
			if (fun) {
				strcpy(callfunc, fun->getName().str().c_str());
			}
			callType = trace_or_not(callfunc);
			if (callType == -1) {
				continue;
			}
		}

		int i, num_of_operands = itr->getNumOperands();
		if (itr->getOpcode() == Instruction::Call && callType == 1) {

			CallInst *CI = dyn_cast<CallInst>(itr);
			Function *fun = CI->getCalledFunction();
			strcpy(operR, (char*)fun->getName().str().c_str());
			if (fun->getName().str().find("dmaLoad") != std::string::npos) {
				print_line(itr, 0, line_number, funcName, bbid, instid, DMA_LOAD);
			}
			else if (fun->getName().str().find("dmaStore") != std::string::npos) {
				print_line(itr, 0, line_number, funcName, bbid, instid, DMA_STORE);
			}
			else {
				print_line(itr, 0, line_number, funcName, bbid, instid, opcode);
			}

			curr_operand = itr->getOperand(num_of_operands - 1);
			is_reg = curr_operand->hasName();
			assert(is_reg);
			print_line(itr, num_of_operands, -1, operR, NULL, NULL,
								 curr_operand->getType()->getTypeID(),
								 getMemSize(curr_operand->getType()),
								 curr_operand,
								 is_reg);

			const Function::ArgumentListType &Args(fun->getArgumentList());
			int num_of_call_operands = CI->getNumArgOperands();
			int call_id = 0;
			for (Function::ArgumentListType::const_iterator arg_it = Args.begin(), arg_end = Args.end(); arg_it != arg_end; ++arg_it) {
				char curr_arg_name[256];
				strcpy(curr_arg_name, (char *)arg_it->getName().str().c_str());

				curr_operand = itr->getOperand(call_id);
				is_reg = curr_operand->hasName();
				if (Instruction *I = dyn_cast<Instruction>(curr_operand)) {
					int flag = 0;
					is_reg = getInstId(I, NULL, operR, flag);
					assert(flag == 0);
					if (curr_operand->getType()->isVectorTy()) {
						print_line(itr, call_id + 1, -1, operR, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 is_reg);
						print_line(itr, FORWARD_LINE, -1, curr_arg_name, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 true);
					}
					else {
						print_line(itr, call_id + 1, -1, operR, NULL, NULL,
											 I->getType()->getTypeID(),
											 getMemSize(I->getType()),
											 curr_operand,
											 is_reg);
						print_line(itr, FORWARD_LINE, -1, curr_arg_name, NULL, NULL,
											 I->getType()->getTypeID(),
										 	 getMemSize(I->getType()),
											 curr_operand,
											 true);
					}
				}
				else {
					if (curr_operand->getType()->isVectorTy()) {
						char operand_id[256];
						strcpy(operand_id, curr_operand->getName().str().c_str());
						print_line(itr, call_id + 1, -1, operand_id, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 is_reg);
						print_line(itr, FORWARD_LINE, -1, curr_arg_name, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 true);
					}

					else if (curr_operand->getType()->isLabelTy()) {
						char label_id[256];
						getBBId(curr_operand, label_id);
						print_line(itr, call_id + 1, -1, label_id, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 is_reg);
						print_line(itr, FORWARD_LINE, -1, curr_arg_name, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 true);
					}
					// is function
					else if (curr_operand->getValueID() == 2) {
						char func_id[256];
						strcpy(func_id, curr_operand->getName().str().c_str());
						print_line(itr, call_id + 1, -1, func_id, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 is_reg);
						print_line(itr, FORWARD_LINE, -1, curr_arg_name, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 NULL,
											 true);
					}
					else {
						char operand_id[256];
						strcpy(operand_id, curr_operand->getName().str().c_str());
						print_line(itr, call_id + 1, -1, operand_id, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 curr_operand,
											 is_reg);
						print_line(itr, FORWARD_LINE, -1, curr_arg_name, NULL, NULL,
											 curr_operand->getType()->getTypeID(),
											 getMemSize(curr_operand->getType()),
											 curr_operand,
											 true);
					}
				}
				call_id++;
			}
		}
		/// If this instruction is not a Call instruction, then execute the following code.
		else
		{
			print_line(itr, 0, line_number, funcName, bbid, instid, opcode);
			if (num_of_operands > 0) {
				for (i = num_of_operands - 1; i >= 0; i--) {
					curr_operand = itr->getOperand(i);
					is_reg = curr_operand->hasName();
					//char arg_label_in_callee[256];

					//for instructions using registers
					if (Instruction *I = dyn_cast<Instruction>(curr_operand)) {
						int flag = 0;
						is_reg = getInstId(I, NULL, operR, flag);
						assert(flag == 0);
						if (curr_operand->getType()->isVectorTy()) {
							print_line(itr, i + 1, -1, operR, NULL, NULL,
												 curr_operand->getType()->getTypeID(),
												 getMemSize(curr_operand->getType()),
												 NULL,
												 is_reg);
						}
						else {
							print_line(itr, i + 1, -1, operR, NULL, NULL,
												 I->getType()->getTypeID(),
												 getMemSize(I->getType()),
												 curr_operand,
												 is_reg);
						}
					} /// If this operand is not an instruction (maybe a label), then execute the following code
					else {
						if (curr_operand->getType()->isVectorTy()) {
							char operand_id[256];
							strcpy(operand_id, curr_operand->getName().str().c_str());
							print_line(itr, i + 1, -1, operand_id, NULL, NULL,
												 curr_operand->getType()->getTypeID(),
												 getMemSize(curr_operand->getType()),
												 NULL,
												 is_reg);
						}

						else if (curr_operand->getType()->isLabelTy()) {
							char label_id[256];
							getBBId(curr_operand, label_id);
							print_line(itr, i + 1, -1, label_id, NULL, NULL,
												 curr_operand->getType()->getTypeID(),
												 getMemSize(curr_operand->getType()),
												 NULL,
												 is_reg);
						}
						// is function
						else if (curr_operand->getValueID() == 2) {
							char func_id[256];
							strcpy(func_id, curr_operand->getName().str().c_str());
							print_line(itr, i + 1, -1, func_id, NULL, NULL,
												 curr_operand->getType()->getTypeID(),
												 getMemSize(curr_operand->getType()),
												 NULL,
												 is_reg);
						}
						else {
							char operand_id[256];
							strcpy(operand_id, curr_operand->getName().str().c_str());
							print_line(itr, i + 1, -1, operand_id, NULL, NULL,
												 curr_operand->getType()->getTypeID(),
												 getMemSize(curr_operand->getType()),
												 curr_operand,
												 is_reg);
						}
					}
				}
			}
		}

		//for call instruction

		//handle function result
		if (!itr->getType()->isVoidTy()) {
			is_reg = 1;
			if (itr->getType()->isVectorTy()) {
				print_line(nextitr, RESULT_LINE, -1, instid, NULL, NULL,
									 itr->getType()->getTypeID(),
									 getMemSize(itr->getType()),
									 NULL, is_reg);
			}
			else if (itr->isTerminator()) {
				DEBUG(dbgs() << "It is terminator...\n");
			}
			else {
				print_line(nextitr, RESULT_LINE, -1, instid, NULL, NULL,
									 itr->getType()->getTypeID(),
									 getMemSize(itr->getType()),
									 itr, is_reg);
			}
		}
	}

	return true;
}

void InstrumentForDDDG::remove_config(std::string kernel_name, std::string input_path) {
	//std::string input_kernel = input_path + kernel_name;
	std::string input_kernel = outputPath + kernel_name;
	std::string pipelining(input_kernel + "_pipelining_config");
	std::string unrolling(input_kernel + "_unrolling_config");
	std::string array_info(input_kernel + "_array_info");
	std::string partition(input_kernel + "_partition_config");
	std::string comp_partition(input_kernel + "_complete_partition_config");

	ifstream pipe(pipelining);
	if (pipe.is_open()) {
		pipe.close();
		if ( remove(pipelining.c_str()) != 0) {
			assert(false && "Error: deleting loop pipelining configuration file!\n");
		}
	}

	ifstream unroll(unrolling);
	if (unroll.is_open()) {
		unroll.close();
		if (remove(unrolling.c_str()) != 0) {
			assert(false && "Error: deleting loop unrolling configuration file!\n");
		}
	}

	ifstream arrayInfo(array_info);
	if (arrayInfo.is_open()) {
		arrayInfo.close();
		if (remove(array_info.c_str()) != 0) {
			assert(false && "Error: deleting array information file!\n");
		}
	}

	ifstream part(partition);
	if (part.is_open()) {
		part.close();
		if (remove(partition.c_str()) != 0) {
			assert(false && "Error: deleting array partitioning configuration file!\n");
		}
	}

	ifstream comp_part(comp_partition);
	if (comp_part.is_open()) {
		comp_part.close();
		if (remove(comp_partition.c_str()) != 0) {
			assert(false && "Error: deleting completely array partitioning configuration file!\n");
		}
	}

}

void InstrumentForDDDG::parse_config(std::string kernel_name, std::string input_path)
{
	ifstream config_file;
	//std::string input_kernel = input_path + kernel_name;
	std::string input_kernel = outputPath + kernel_name;
	//std::string config_file_name = input_kernel + "_configuration";
	//std::string config_file_name = input_path + config_filename;
	std::string config_file_name = config_filename;
	config_file.open(config_file_name);
	if (!config_file.is_open()) {
		assert(false && "Error: missing configuration file!\n");
	}
	std::string wholeline;

	std::vector<std::string> pipelining_config;
	std::vector<std::string> unrolling_config;
	std::vector<std::string> partition_config;
	std::vector<std::string> comp_partition_config;
	std::vector<std::string> array_info;

	pipeline_loop_levelVec.clear();

	while (!config_file.eof())
	{
		wholeline.clear();
		getline(config_file, wholeline);
		if (wholeline.size() == 0)
			break;
		string type, rest_line;
		int pos_end_tag = wholeline.find(",");
		if (pos_end_tag == -1)
			break;
		type = wholeline.substr(0, pos_end_tag);
		//rest_line = wholeline.substr(pos_end_tag + 1);

		if (!type.compare("pipeline")) {
			pipelining_config.push_back(wholeline);
		}
		else if (!type.compare("unrolling")) {
			unrolling_config.push_back(wholeline);
		}
		else if (!type.compare("array")) {
			array_info.push_back(wholeline);
		}
		else if (!type.compare("partition")) {
			if (wholeline.find("complete") == std::string::npos)
				partition_config.push_back(wholeline);
			else
				comp_partition_config.push_back(wholeline);
		}
		else
		{
			// Can not detect type of this line, ignore it and do
			// nothing here.
		}
	}
	config_file.close();
	//std::string test = "_test.txt";
	std::map<std::string, unsigned> whole_lpName2comp_unroll_factorMap;
	if (pipelining_config.size() != 0)
	{
		string pipelining(input_kernel);
		pipelining += "_pipelining_config";
		//pipelining += test;
		ofstream pipe_config;
		pipe_config.open(pipelining);
		for (unsigned i = 0; i < pipelining_config.size(); ++i) {
			std::string pipelining_conf = pipelining_config.at(i);
			pipe_config << pipelining_conf << endl;

			char type[256];
			char funcName_char[256];
			unsigned loop_num, loop_level;
			sscanf(pipelining_conf.c_str(), "%[^,],%[^,],%d,%d\n", type, funcName_char, &loop_num, &loop_level);
			std::string func_name(funcName_char);
			std::string lp_name = func_name + "_loop-" + std::to_string(loop_num);
			LpName2numLevelMapTy::iterator found_lpLev = LpName2numLevelMap.find(lp_name);
			assert(found_lpLev != LpName2numLevelMap.end() && "Error: Cannot find loop name in LpName2numLevelMap!\n");
			unsigned num_level = LpName2numLevelMap[lp_name];
			std::string whole_loop_name = lp_name + "_" + std::to_string(loop_level);
			pipeline_loop_levelVec.push_back(whole_loop_name);

			assert( (loop_level<=num_level) && "Error: loop_level is larger than num_level!\n" );
			if (loop_level == num_level) {
				// Apply loop pipelining to the innermost level loop
				continue;
			}
			for (unsigned i = loop_level+1; i < num_level+1; i++) {
				std::string whole_lp_name = lp_name + "_" + std::to_string(i);
				wholeloopName2loopBoundMapTy::iterator it_whole = wholeloopName2loopBoundMap.find(whole_lp_name);
				assert((it_whole!=wholeloopName2loopBoundMap.end()) && "Error: Can not find loop name in wholeloopName2loopBoundMap!\n");
				unsigned lp_bound = wholeloopName2loopBoundMap[whole_lp_name];
				if (lp_bound == 0) {
					VERBOSE_PRINT(std::cout << "DEBUG-INFO: [parsing_configuration-extraction] loop " << lp_name << " level " << i);
					VERBOSE_PRINT(std::cout << " has a variable loop bound, can not support in current version!" << std::endl);
				}
				whole_lpName2comp_unroll_factorMap.insert(std::make_pair(whole_lp_name, lp_bound));
			}
		}
		pipe_config.close();
	}

	if (unrolling_config.size() != 0)
	{
		string file_name(input_kernel);
		file_name += "_unrolling_config";
		//file_name += test;
		ofstream output;
		output.open(file_name);
		std::vector<std::string> unroll_wholelpName_str;
		for (unsigned i = 0; i < unrolling_config.size(); ++i) {
			std::string unrolling_conf = unrolling_config.at(i);
			char type[256];
			char funcName_char[256];
			unsigned loop_num, loop_level;
			unsigned line_num, unroll_factor;
			sscanf(unrolling_conf.c_str(), "%[^,],%[^,],%d,%d,%d,%d\n", type, funcName_char, &loop_num, &loop_level, &line_num, &unroll_factor);
			std::string func_name(funcName_char);
			std::string lp_name = func_name + "_loop-" + std::to_string(loop_num);
			std::string whole_lp_name = lp_name + "_" + std::to_string(loop_level);
			unroll_wholelpName_str.push_back(whole_lp_name);
			std::map<std::string, unsigned>::iterator it_whole = whole_lpName2comp_unroll_factorMap.find(whole_lp_name);
			if (it_whole == whole_lpName2comp_unroll_factorMap.end()) {
				output << unrolling_config.at(i) << endl;
			}
			else{
				unsigned static_bound = it_whole->second;
				if (static_bound == 0) {
					output << unrolling_config.at(i) << endl;
				}
				else {
					std::string new_unroll_conf = "unrolling," + func_name + "," + std::to_string(loop_num) + ",";
					new_unroll_conf += std::to_string(loop_level) + "," + std::to_string(line_num) + "," + std::to_string(static_bound);
					output << new_unroll_conf << endl;
				}
			}
		}

		// If we apply loop pipelining at the upper loop level, but we donot specify completely unrolling pragma at the inner loop
		// levels, we need to add it to unrolling configuration file.
		std::map<std::string, unsigned>::iterator it_pipe = whole_lpName2comp_unroll_factorMap.begin();
		std::map<std::string, unsigned>::iterator ie_pipe = whole_lpName2comp_unroll_factorMap.end();
		for (; it_pipe != ie_pipe; ++it_pipe) {
			std::string whole_lp_name = it_pipe->first;
			unsigned lp_bound = it_pipe->second;
			std::vector<std::string>::iterator it_unr = unroll_wholelpName_str.begin();
			std::vector<std::string>::iterator ie_unr = unroll_wholelpName_str.end();
			std::vector<std::string>::iterator not_found = std::find(it_unr, ie_unr, whole_lp_name);
			if (not_found == unroll_wholelpName_str.end()) {
				std::size_t pos = whole_lp_name.find("_loop-");
				std::string func_name = whole_lp_name.substr(0, pos);
				std::string rest_str = whole_lp_name.substr(pos+6);
				pos = rest_str.find("_");
				unsigned loop_num = std::stoi(rest_str.substr(0, pos));
				unsigned loop_level = std::stoi(rest_str.substr(pos+1));
				int line_num = -1;
				
				std::string new_unroll_conf = "unrolling," + func_name + "," + std::to_string(loop_num) + ",";
				new_unroll_conf += std::to_string(loop_level) + "," + std::to_string(line_num) + "," + std::to_string(lp_bound);
				output << new_unroll_conf << endl;
			}
		}

		output.close();
	}

	if (array_info.size() != 0) {
		string file_name(input_kernel);
		file_name += "_array_info";
		//file_name += test;
		ofstream output;
		output.open(file_name);
		for (unsigned i = 0; i < array_info.size(); ++i)
			output << array_info.at(i) << endl;
		output.close();
	}
	else {
		assert(false && "Error: please provide array information for this kernel!\n");
	}

	if (partition_config.size() != 0)
	{
		string partition(input_kernel);
		partition += "_partition_config";
		//partition += test;
		ofstream part_config;
		part_config.open(partition);
		for (unsigned i = 0; i < partition_config.size(); ++i)
			part_config << partition_config.at(i) << endl;
		part_config.close();
	}

	if (comp_partition_config.size() != 0)
	{
		string complete_partition(input_kernel);
		complete_partition += "_complete_partition_config";
		//complete_partition += test;
		ofstream comp_config;
		comp_config.open(complete_partition);
		for (unsigned i = 0; i < comp_partition_config.size(); ++i)
			comp_config << comp_partition_config.at(i) << endl;
		comp_config.close();
	}

	increase_load_latency = false;
	if (pipelining_config.size() != 0) {
		increase_load_latency = false;
	}
	else {
		if (partition_config.size() != 0) {
			for (int i = 0; i < partition_config.size(); i++) {
				std::string  part_str = partition_config.at(i);
				unsigned size, p_factor, wordsize;
				char config_type[256];
				char type[256];
				char base_addr[256];
				sscanf(part_str.c_str(), "%[^,],%[^,],%[^,],%d,%d,%d\n", config_type, type, base_addr, &size, &wordsize, &p_factor);
				if (p_factor > 1) {
					increase_load_latency = true;
					break;
				}
			}
		}
		else {
			increase_load_latency = false;
		}
	}
}

void InstrumentForDDDG::getUnrollingConfiguration(lpNameLevelPair2headBBnameMapTy& lpNameLvPair2headerBBMap) {
	// Initialize loopName2levelUnrollPairListMap, all unrolling factors are set to 1 as default.
	loopName2levelUnrollVecMap.clear();
	lpNameLevelPair2headBBnameMapTy::iterator it = lpNameLvPair2headerBBMap.begin();
	lpNameLevelPair2headBBnameMapTy::iterator ie = lpNameLvPair2headerBBMap.end();
	for (; it != ie; ++it) {
		std::string loopName = it->first.first;
		unsigned loop_level = std::stoul(it->first.second);
		loopName2levelUnrollVecMapTy::iterator itVec = loopName2levelUnrollVecMap.find(loopName);
		if (itVec != loopName2levelUnrollVecMap.end()) {
			itVec->second.push_back(1);
		}
		else {
			std::vector<unsigned> levelUnrollPairVec;
			levelUnrollPairVec.push_back(1);
			loopName2levelUnrollVecMap.insert(std::make_pair(loopName, levelUnrollPairVec));
		}
	}

	// Read unrolling configuration and update loopName2levelUnrollPairListMap and 
	unrollingConfig.clear();
	bool succeed_or_not = readUnrollingConfig(loopName2levelUnrollVecMap, unrollingConfig);
}

bool InstrumentForDDDG::readUnrollingConfig(loopName2levelUnrollVecMapTy& lpName2levelUrPairVecMap, std::unordered_map<int, int > &unrolling_config) {
	//loopName2levelUnrollPairListMapTy vector is better
	std::string kernel_name = kernel_names.at(0);
	ifstream config_file;
	//std::string file_name(inputPath + kernel_name);
	std::string file_name(outputPath + kernel_name);
	file_name += "_unrolling_config";
	config_file.open(file_name.c_str());
	if (!config_file.is_open())
		return 0;
	while (!config_file.eof())
	{
		std::string wholeline;
		getline(config_file, wholeline);
		if (wholeline.size() == 0)
			break;
		char func[256];
		char config_type[256];
		unsigned ith_loop;
		unsigned loop_level;
		int line_num, factor;
		sscanf(wholeline.c_str(), "%[^,],%[^,],%d, %d, %d,%d\n", config_type, func, &ith_loop, &loop_level, &line_num, &factor);
		unrolling_config[line_num] = factor;
		std::string loop_name = std::string(func) + "_loop-" + std::to_string(ith_loop);
		loopName2levelUnrollVecMapTy::iterator it = lpName2levelUrPairVecMap.find(loop_name);
		if (it != lpName2levelUrPairVecMap.end()) {
			unsigned innermost_level = it->second.size();
			if ((loop_level != innermost_level) && factor > 1) {
				//std::cout << "DEBUG-INFO: [parsing_read-unrolling-configuration] For FPGA implementation, the unrolling for the innermost loop level is more interesting, we only consider this level loop." << std::endl;
				assert(loop_level<=innermost_level && "Error: Loop level information inside the unrolling_config file exceeds number of levels in this loop!\n");
				it->second[loop_level-1] = factor;
			}
			else {
				std::string whole_loop_name = loop_name + "_" + std::to_string(loop_level);
				unsigned lp_bound = wholeloopName2loopBoundMap.at(whole_loop_name);
				if (lp_bound == 0) {
					// Need to get loop bound at runtime
					it->second.at(loop_level - 1) = factor;
				}
				else {
					if ((unsigned)factor > lp_bound) {
						it->second.at(loop_level - 1) = lp_bound;
						//assert((factor<lp_bound) && "Error: Loop unrolling factor is larger than loop bound! Please use smaller loop bound\n");
					}
					else {
						it->second.at(loop_level - 1) = factor;
					}
				}

			}
		}
		else {
			std::cout << "DEBUG-INFO: [parsing_read-unrolling-configuration] Warning: Can not find the loop name inside configuration files" << std::endl;
			assert(false && "Please check whether configuration file is written in a correct way!\n");
		}
	}
	config_file.close();
	return 1;
}

void InstrumentForDDDG::loopBasedTraceAnalysis() {
	errs() << "DEBUG-INFO: [trace-analysis_loop-based-trace-analysis] Analysis loops' IL and II inside the kernel\n";
	/// Create Dynamic Data Dependence Graph
	std::string trace_file_name = inputPath + "dynamic_trace.gz";
	//std::string trace_file_name = "dynamic_trace";
	std::string kernel_name = kernel_names.at(0);

	/// Open summary file
	open_summary_file(summary, kernel_name);

	/// Remove previous configuration files
	remove_config(kernel_name, inputPath);

	/// Generate configuration files
	parse_config(kernel_name, inputPath);

	// Get unrolling configuration
	getUnrollingConfiguration(lpNameLevelPair2headBBnameMap);

	loopName2levelUnrollVecMapTy::iterator it = loopName2levelUnrollVecMap.begin();
	loopName2levelUnrollVecMapTy::iterator ie = loopName2levelUnrollVecMap.end();
	for (; it != ie; ++it) {
		bool enable_pipelining = false;
		bool skip_loop_analysis = false;
		std::string loop_name = it->first;

		if (target_loops.size() != 0) {
			skip_loop_analysis = true;
			size_t pos_index = loop_name.find("-");
			if (pos_index != std::string::npos) {
				std::string loop_index = loop_name.substr(pos_index + 1);
				std::vector<std::string>::iterator found_lp_index = std::find(target_loops.begin(), target_loops.end(), loop_index);
				if (found_lp_index != target_loops.end()) {
					skip_loop_analysis = false;
				}
			}

			if (skip_loop_analysis == true) {
				// Do not analyze this loop
				continue;
			}
		}

		std::vector<unsigned> levelUnrollVec = it->second;
		unsigned level_size = levelUnrollVec.size();

		unsigned target_loop_level = 1;
		unsigned target_unroll_factor = 1;

		for (int i = (int) level_size-1; i >= 0; i--) {
			unsigned unroll_factor = levelUnrollVec.at(i);
			std::string whole_lp_name = loop_name + "_" + std::to_string(i+1);
			wholeloopName2loopBoundMapTy::iterator it_found = wholeloopName2loopBoundMap.find(whole_lp_name);
			if (it_found != wholeloopName2loopBoundMap.end()) {
				unsigned level_bound = it_found->second;
				bool exit_flag = false;
				if (unroll_factor != level_bound) {
					target_loop_level = i + 1;
					exit_flag = true;
					//break;
				}
				target_unroll_factor = unroll_factor;

				std::vector<std::string>::iterator found_pipe = std::find(pipeline_loop_levelVec.begin(), pipeline_loop_levelVec.end(), whole_lp_name);
				if (found_pipe != pipeline_loop_levelVec.end()) {
					enable_pipelining = true;
				}
				else {
					enable_pipelining = false;
				}

				if (exit_flag) {
					break;
				}
			}
			else {
				assert(false && "Error: Can not find loop name in wholeloopName2loopBoundMap!\n");
			}
		}

		unsigned target_loop_bound;
		std::string whole_target_name = loop_name + "_" + std::to_string(target_loop_level);
		wholeloopName2loopBoundMapTy::iterator found_it = wholeloopName2loopBoundMap.find(whole_target_name);
		if (found_it != wholeloopName2loopBoundMap.end()) {
			target_loop_bound = found_it->second;
		}
		else {
			assert(false && "Error: Cannot find loop name in wholeloopName2loopBoundMap!\n");
		}

		unsigned unroll_factor = 1;
		/// Used to get recII
		unsigned IL_asap_rec = 0;
		if (enable_pipelining == true) {
			DynamicDatapath* datapath_tmp;
			unsigned target_factor = target_unroll_factor * 2;
			unroll_factor = ((target_loop_bound < target_factor) && (target_loop_bound != 0)) ? target_loop_bound : target_factor;
			datapath_tmp = new DynamicDatapath(kernel_name, trace_file_name, inputPath, loop_name, target_loop_level, unroll_factor , enable_pipelining, 0);
			IL_asap_rec = datapath_tmp->getIL_asap_ii();
			delete datapath_tmp;
		}

		errs() << "DEBUG-INFO: [trace-analysis_loop-based-trace-analysis] Building Dynamic Datapath for loop " << loop_name <<"\n";
		unroll_factor = ((target_loop_bound < target_unroll_factor) && (target_loop_bound != 0)) ? target_loop_bound : target_unroll_factor;
		DynamicDatapath DynDatapath(kernel_name, trace_file_name, inputPath, loop_name, target_loop_level, unroll_factor, false, IL_asap_rec);
		VERBOSE_PRINT(errs() << "DEBUG-INFO: [trace-analysis_loop-based-trace-analysis] Finished dynamic trace analysis for loop " << loop_name << "\n");
		VERBOSE_PRINT(errs() << "-------------------\n");
	}
	VERBOSE_PRINT(errs() << "DEBUG-INFO: [trace-analysis_loop-based-trace-analysis] Finished\n");
	close_summary_file(summary);
	errs() << "-------------------\n";
}

void InstrumentForDDDG::open_summary_file(ofstream& summary_file, std::string kernel_name) {
	errs() << "DEBUG-INFO: [Lin-Analyzer summary] Writing summary into log file\n";
	//std::string file_name(inputPath+kernel_name+"_summary.log");
	std::string file_name(outputPath + kernel_name + "_summary.log");
	summary_file.open(file_name);
	if (summary_file.is_open()) {
		summary_file << "==========================" << std::endl;
		summary_file << "   Lin-analyzer summary" << std::endl;
		summary_file << "==========================" << std::endl;
		summary_file << "function name: " << kernel_name << std::endl;
	}
	else {
		assert(false && "Error: Cannot open summary file!\n");
	}
}

void InstrumentForDDDG::close_summary_file(ofstream& summary_file) {
	summary_file.close();
	VERBOSE_PRINT(errs() << "DEBUG-INFO: [Lin-Analyzer summary] Summary file generated\n");
}

ProfilingEngine::ProfilingEngine(Module &M, Function* log0Fn, Function* log_intFn, Function* log_doubleFn,
																 Function* log_int_noregFn, Function* log_double_noregFn) : Mod(M), log0_Fn(log0Fn), 
																 log_int_Fn(log_intFn), log_double_Fn(log_doubleFn), 
																 log_int_noreg_Fn(log_int_noregFn), log_double_noreg_Fn(log_double_noregFn) {

}

void ProfilingEngine::runOnProfiler(){
	errs() << "DEBUG-INFO: [profiling_profiling-engine] Running on Profiler\n";
	ProfilingJITSingletonContext JTSC(this);
	//std::unique_ptr<Module> M(generateInstrumentedModule());
	//ValueToValueMapTy V2VMap;
	//Module* M = CloneModule(&Mod, V2VMap);
	//std::unique_ptr<Module> M(newM);
	Module *M = (&Mod);

	// TODO: Insert instrumental code to extract the trace.
	EngineBuilder builder(M);
	ExecutionEngine *EE = builder.setEngineKind(EngineKind::JIT).create();
	// Where is getBBFreqInc()?
	//EE->addGlobalMapping(getBBFreqInc(), reinterpret_cast<void*>(IncreaseBBCounter));
	EE->addGlobalMapping(log0_Fn, reinterpret_cast<void*>(trace_logger_log0));
	EE->addGlobalMapping(log_int_Fn, reinterpret_cast<void*>(trace_logger_log_int));
	EE->addGlobalMapping(log_double_Fn, reinterpret_cast<void*>(trace_logger_log_double));
	EE->addGlobalMapping(log_int_noreg_Fn, reinterpret_cast<void*>(trace_logger_log_int_noreg));
	EE->addGlobalMapping(log_double_noreg_Fn, reinterpret_cast<void*>(trace_logger_log_double_noreg));

	if (!EE) {
		assert(false && "Error: Failed to construct ExecutionEngine\n");
	}

	Function *EntryFn = M->getFunction("main");

	// Nothing we can do if we cannot find the entry function of the module.
	if (EntryFn == nullptr){
		assert(false && "Error: EntryFn is equal to nullptr!\n");
	}


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
	VERBOSE_PRINT(errs() << "DEBUG-INFO: [profiling_profiling-engine] Running main()\n");
	int Result = EE->runFunctionAsMain(EntryFn, Argv, nullptr);

	if (Result != 0) {
		errs() << "DEBUG-INFO: [profiling_profiling-engine] Module return nonzero result during trace generation: "<< Result << '\n';
		assert(false && "Error: Trace generation failed!\n");
	}
	
	trace_logger_fin();

	// Run static destructors.
	EE->runStaticConstructorsDestructors(true);

	//VerifyProfile();

	//delete M
	errs() << "DEBUG-INFO: [profiling_profiling-engine] Finished profiling: Status = " << Result << "\n";
}

ProfilingJITContext::ProfilingJITContext() : P(nullptr) {
	// If we have a native target, initialize it to ensure it is linked in and
	// usable by the JIT.
	InitializeNativeTarget();
	InitializeNativeTargetAsmPrinter();
	InitializeNativeTargetAsmParser();
}


ProfilingJITSingletonContext::ProfilingJITSingletonContext(ProfilingEngine *P) {
	GlobalContextDDDG->P = P;
}

ProfilingJITSingletonContext::~ProfilingJITSingletonContext() {
	GlobalContextDDDG->P = nullptr;
	//TODO: Other clearup?
}


char InstrumentForDDDG::ID = 0;
INITIALIZE_PASS(InstrumentForDDDG, "instrumentCodeforDDDG",
								"This pass is used to instrument code for building DDDG",
								false,
								true	/*false  We need to modify the code later. In initial step, we
											just check the unique ID first*/
								)

ModulePass *llvm::createInstrumentForDDDGPass() {
	return new 	InstrumentForDDDG();
}