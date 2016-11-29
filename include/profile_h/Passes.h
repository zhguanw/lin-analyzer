#ifndef PASSES_H
#define PASSES_H

namespace llvm{
	class ModulePass;
	class Pass;
	
	/// GetInstDistribution Pass
	ModulePass *createGetInstDistributionPass();
	void initializeGetInstDistributionPass(PassRegistry &Registry);

	/// AssignBasicBlockID Pass and QueryBasicBlockID Pass
	ModulePass *createAssignBasicBlockIDPass();
	void initializeAssignBasicBlockIDPass(PassRegistry &Registry);

	ModulePass *createQueryBasicBlockIDPass();
	void initializeQueryBasicBlockIDPass(PassRegistry &Registry);

	/// AssignLoadStoreID Pass and QueryLoadStoreID Pass
	ModulePass *createAssignLoadStoreIDPass();
	void initializeAssignLoadStoreIDPass(PassRegistry &Registry);

	ModulePass *createQueryLoadStoreIDPass();
	void initializeQueryLoadStoreIDPass(PassRegistry &Registry);

	/// ExtractLoopInfo Pass and LoopNumber Pass
	Pass* createExtractLoopInfoPass();
	void initializeExtractLoopInfoPass(PassRegistry &Registry);

	Pass* createLoopNumberPass();
	void initializeLoopNumberPass(PassRegistry &Registry);

	/// GetLoopBound pass
	ModulePass* createGetLoopBoundPass();
	void initializeGetLoopBoundPass(PassRegistry &Registry);

	/// CodeInstrument Pass
	ModulePass* createCodeInstrumentPass();
	void initializeCodeInstrumentPass(PassRegistry &Registry);

	/// AnalysisProfiling Pass
	ModulePass* createAnalysisProfilingPass();
	void initializeAnalysisProfilingPass(PassRegistry &Registry);

	/// InstrumentForDDDG Pass
	ModulePass* createInstrumentForDDDGPass();
	void initializeInstrumentForDDDGPass(PassRegistry &Registry);

	// Basic Block ID to frequency
	typedef std::map<uint64_t, uint64_t> BBidFreqMap;
	typedef std::map<BasicBlock*, uint64_t> BBFreqMap;
	// Basic Block ID to BBFreqMap
	typedef std::map<uint64_t, BBidFreqMap> BranchidFreqMap;
	typedef std::map<BasicBlock*, BBFreqMap> BranchFreqMap;

	extern BBidFreqMap BBidFreq;
	extern BranchidFreqMap BranchIDFreq;

	// Calculate number of arithmatic operations for loops inside functions except for the "main" function
	// The map used to trace binary operators: (functionName, loopName) -> std::vector<BB id> > fnlpNamePair2BinaryOpBBidMap.
	// After we get the basic block frequency, we can get the total number of binary operations executed by
	// multiplication.
	typedef std::map<std::pair<std::string, std::string>, std::vector<uint64_t> > fnlpNamePair2BinaryOpBBidMapTy;
	extern fnlpNamePair2BinaryOpBBidMapTy fnlpNamePair2BinaryOpBBidMap;

#ifdef _MSC_VER // If we use microsoft visual studio
	const uint32_t INITIAL_BIGVALUE = UINT_MAX;
#else 
	const uint32_t INITIAL_BIGVALUE = 0xffffffff;
#endif

}

#endif // End PASSES_H