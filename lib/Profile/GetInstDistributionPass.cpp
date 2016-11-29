#include "profile_h/GetInstDistributionPass.h"

#define DEBUG_TYPE "get-inst-distribution"
using namespace llvm;

STATISTIC(NumTotalInst, "The number of total instructions: #");
/// StoreInst and LoadInst
STATISTIC(NumMemInst, "The number of load/store instructions: #");
/// TerminatorInst
STATISTIC(NumBranch, "The number of branch instructions: #");
/// CmpInst and BinaryOperator
STATISTIC(NumComputation, "The number of computing operations: #");
/// Others, eg. PHINode, CallInst, CastInst, etc.
STATISTIC(NumOthers, "The number of other instructions: #");
/// LoadInst
STATISTIC(NumLoadInst, "The number of load instructions: #");
/// StoreInst
STATISTIC(NumStoreInst, "The number of store instructions: #");

GetInstDistribution::GetInstDistribution() : ModulePass(ID) {
	errs() << "\n---------------------------------------\n";
	errs() << "\n\tInitialize GetInstDistribution pass\n";
	initializeGetInstDistributionPass(*PassRegistry::getPassRegistry());
}

void GetInstDistribution::getAnalysisUsage(AnalysisUsage &AU) const {
	//AU.addRequiredID(LoopSimplifyID);
	//AU.addPreservedID(LoopSimplifyID);

	//AU.addRequiredID(LCSSAID);
	//AU.addPreservedID(LCSSAID);

	AU.addRequired<QueryBasicBlockID>();

	AU.setPreservesCFG();
}

bool GetInstDistribution::runOnModule(Module &M) {

	/// To get the instruction distribution of kernels, there are three steps 
	/// before we get the final result.
	/// 1. Get the basic block frequency by the Embedded Execution Engine
	/// 2. Get the instruction distribution of each basic block of the kernel
	/// 3. Multiply instruction distribution per basic block by basic block frequency
	/// Done

	/// Make sure basic block frequency is available
	if (BBidFreq.size() == 0 || BranchIDFreq.size() == 0) {
		errs() << "BBidFreq and BranchIDFreq are empty! Error!\n";
		return false;
	}

	/// Collect instruction distribution information per basic block
	bbIDpass = &getAnalysis<QueryBasicBlockID>();
	Module::iterator FI, FE;
	Function::iterator BI, BE;
	BasicBlock::iterator I, E;

	std::vector<std::string>::iterator it_name;
	std::string Fun_Name;
	for (FI = M.begin(), FE = M.end(); FI != FE; ++FI) {
		it_name = std::find(kernel_names.begin(), kernel_names.end(), FI->getName());
		Fun_Name = FI->getName().str();
		if (it_name != kernel_names.end()) {
			errs() << "\tFunction Name = " << FI->getName() << "\n";
			// Explore instructions within a kernel
			for (BI = FI->begin(), BE = FI->end(); BI != BE; ++BI) {
				instDis instDisBB;
				instDisBB.meminst = 0;
				instDisBB.ldinst = 0;
				instDisBB.stinst = 0;
				instDisBB.brinst = 0;
				instDisBB.compinst = 0;
				instDisBB.others = 0;

				for (I = BI->begin(), E = BI->end(); I != E; ++I) {
					Instruction* inst = &*I;
					//NumTotalInst++;
					if (isa<LoadInst>(inst)) {
						//NumMemInst++;
						//NumLoadInst++;
						instDisBB.meminst++;
						instDisBB.ldinst++;
					}
					else if (isa<StoreInst>(inst)) {
						//NumMemInst++;
						//NumStoreInst;
						instDisBB.meminst++;
						instDisBB.stinst++;
					}
					else if (isa<TerminatorInst>(inst)) {
						//NumBranch++;
						instDisBB.brinst++;
					}
					else if (isa<CmpInst>(inst) || isa<BinaryOperator>(inst)) {
						//NumComputation++;
						instDisBB.compinst++;
					}
					else {
						//NumOthers++;
						instDisBB.others++;
					}

				} // Get the instruction distribution of a basic block now

				// Store instruction distribution per basic block into bbID2instDisMap
				uint64_t bbid = bbIDpass->getBBid(BI);
				bbID2instDisMap.insert(std::make_pair(bbid, instDisBB));
			}
		}
	}
	
	/// Calculate instruction distribution of the kernel
	NumTotalInst = 0;
	NumMemInst = 0;
	NumLoadInst = 0;
	NumStoreInst = 0;
	NumBranch = 0;
	NumComputation = 0;
	NumOthers = 0;
	bbID2instDisMapTy::iterator it, ie;
	for (it = bbID2instDisMap.begin(), ie = bbID2instDisMap.end(); it != ie; it++) {
		uint64_t bb_id = it->first;
		errs() << "bb_id = " << bb_id << "\n";
		assert( (BBidFreq.find(bb_id) != BBidFreq.end()) && "Can not find bb_id in BBidFreq! Error\n");
		uint64_t freqbbid = BBidFreq.at(bb_id);
		NumLoadInst += it->second.ldinst * freqbbid;
		NumStoreInst += it->second.stinst * freqbbid; 
		NumBranch += it->second.brinst * freqbbid;
		NumComputation += it->second.compinst * freqbbid;
		NumOthers += it->second.others * freqbbid;
	}
	NumMemInst += NumLoadInst + NumStoreInst;
	NumTotalInst += NumMemInst + NumBranch + NumComputation + NumOthers;

	/// Write the instruction distribution information into a file
	std::string instDis_fileName;
	if (instDis_fileName.empty()) {
		if (inputFileName == "-") {
			instDis_fileName = "-";
		}
		else {
			const std::string &IFN = inputFileName;
			int Len = IFN.length();
			// If the source ends in .bc, strip it off
			if (IFN[Len - 3] == '.' && IFN[Len - 2] == 'b' && IFN[Len - 1] == 'c') {
				instDis_fileName = std::string(IFN.begin(), IFN.end() - 3) + ".instdis";
			}
		}
	}

	std::ofstream instdisFile;
	instdisFile.open(instDis_fileName);
	instdisFile << "Instruction Distribution of kernel " << Fun_Name << "\n";
	instdisFile << "NumMemInst;NumComputation;NumBranch;NumOthers;";
	instdisFile << "NumTotalInst;NumLoadInst;NumStoreInst\n";
	instdisFile << NumMemInst << ";" << NumComputation << ";" << NumBranch << ";";
	instdisFile << NumOthers << ";" << NumTotalInst << ";" << NumLoadInst << ";" << NumStoreInst << "\n";
	instdisFile << "Finished!\n";
	instdisFile.close();

	/// Print the information below
	errs() << "\n\tPrinting instruction distribution below:\n";
	errs() << "\t Memory instructions: \t# " << NumMemInst << "\n";
	errs() << "\t\t Load instructions: \t# " << NumLoadInst << "\n";
	errs() << "\t\t Store instructions: \t# " << NumStoreInst << "\n";
	errs() << "\t Computation instructions: \t# " << NumComputation << "\n";
	errs() << "\t Branch instructions: \t# " << NumBranch << "\n";
	errs() << "\t Other instructions: \t# " << NumOthers << "\n";
	errs() << "\t -----------\n";
	errs() << "\t Total instructions: \t# " << NumTotalInst << "\n";
	errs() << "\t End\n";

	errs() << "\n\tFinish GetInstDistribution Pass!\n";
	errs() << "\n---------------------------------------\n";

	// Clear BBidFreq, BranchIDFreq
	BBidFreq.clear();
	BranchIDFreq.clear();
	
	return false;
}


char GetInstDistribution::ID = 0;

INITIALIZE_PASS(GetInstDistribution, "getInstDist",
				"This pass is used to get instruction distribution of kernels",
				false,
				false
				)


ModulePass *llvm::createGetInstDistributionPass() {
	return new 	GetInstDistribution();
}