#ifndef AUXILIARY_H
#define AUXILIARY_H

#include "llvm/Support/Timer.h"
#include "llvm/IR/BasicBlock.h"
#define ENABLE_TIMER

#include <list>
#include <fstream>

#ifdef ENABLE_BTREE_MAP
#include "profile_h/btree_map.h"
#endif // End of ENABLE_BTREE_MAP

#define READ_PORT_PER_PARTITION 2
#define WRITE_PORT_PER_PARTITION 1
#define MAX_WRITE_PORT_PER_PARTITION (WRITE_PORT_PER_PARTITION+1)

#define INFINITE_HARDWARE 999999999

extern std::string inputFileName;
extern std::string inputPath;
extern std::string outputPath;
extern std::string config_filename;
extern std::vector<std::string> kernel_names;
extern bool enable_profiling_time_only;
extern bool memory_trace_gen;
extern bool enable_no_trace;
extern bool show_cfg_detailed;
extern bool show_cfg_only;
extern bool show_dddg_bf_opt;
extern bool show_dddg_af_opt;
extern bool verbose_print;
extern bool enable_store_buffer;
extern bool enable_shared_load_removal;
extern bool disable_shared_load_removal;
extern bool enable_repeated_store_removal;
extern bool enable_memory_disambiguation;
extern bool enable_tree_height_reduction_float;
extern bool enable_tree_height_reduction_integer;
extern bool increase_load_latency;
extern bool disable_fp_unit_threshold;
extern bool enable_extra_scalar;
extern bool enable_rw_rw_memory;
extern bool target_vc707;
extern std::vector<std::string> target_loops;
extern std::ofstream summary;

typedef std::map<std::string, std::string> getElementPtrName2arrayNameMapTy;
extern getElementPtrName2arrayNameMapTy getElementPtrName2arrayNameMap;

#define VERBOSE_PRINT(X) \
	do { if (verbose_print) { X; } \
	} while (0)

namespace llvm {
	typedef std::pair<std::string, std::string> fn_lpNamePairTy;

	/// FunctionName --> number of loops inside it
	typedef std::map<std::string, uint64_t> funcName2loopNumMapTy;
	extern funcName2loopNumMapTy funcName2loopNumMap;

	/// Loop Basic Block ID --> (the corresponding Loop Header BB id, loop depth)
#ifdef ENABLE_BTREE_MAP
	typedef btree::btree_map<uint64_t, std::pair<uint64_t, unsigned> > loopBBid2HeaderAndloopDepthPairMapTy;
#else
	typedef std::map<uint64_t, std::pair<uint64_t, unsigned> > loopBBid2HeaderAndloopDepthPairMapTy;
#endif // End of ENABLE_BTREE_MAP
	extern loopBBid2HeaderAndloopDepthPairMapTy loopBBid2HeaderAndloopDepthPairMap;

	/// headerID_list_vec contains a vector of headerID_list
	/// headerID_list:  The header BB id of the innermost loop --> The header BB id of the outermost loop
	typedef std::vector<std::list<uint64_t> > headerID_list_vecTy;
	extern headerID_list_vecTy headerID_list_vec;

	typedef std::map<std::string, headerID_list_vecTy> funcName2headIDlistVecMapTy;
	extern funcName2headIDlistVecMapTy funcName2headIDlistVecMap;

	//typedef std::map<std::string, std::vector<std::pair<unsigned, unsigned> > > loop2subloopMapTy;
	/// loopName --> subloop <depth, bound> map
	typedef std::map<std::string, std::map<unsigned, unsigned> > loop2subloopMapTy;

	typedef std::map<std::string, loop2subloopMapTy> func2loopInfoMapTy;
	extern func2loopInfoMapTy func2loopInfoMap;

	typedef std::map<std::string, std::vector<std::string> > funcName2loopsMapTy;
	extern funcName2loopsMapTy funcName2loopsMap;

	/// loopName2BoundList: the whole loop name -> loop bounds list
	/// Loop bounds sequence: 
	/// The top-level loop -> lower-level loop -> ... -> the innermost-level loop
	typedef std::map<std::string, std::list<unsigned> > loopName2BoundListTy;
	typedef std::map<std::string, loopName2BoundListTy> funcName2loopBoundsTy;
	extern funcName2loopBoundsTy funcName2loopBounds;

	/// Used to recognize the relationship between basic block ids and loop name
	typedef std::map<uint64_t, std::string> BBids2loopNameMapTy;
	extern BBids2loopNameMapTy BBids2loopNameMap;

	typedef std::map<BasicBlock*, std::string> BB2loopNameMapTy;
	extern BB2loopNameMapTy BB2loopNameMap;

	typedef std::map<std::string, BBids2loopNameMapTy> funcName2loopInfoMapTy;
	extern funcName2loopInfoMapTy funcName2loopInfoMap;

#ifdef ENABLE_BTREE_MAP
	typedef btree::btree_map<uint64_t, std::list<unsigned> > lpBBid2EvalIterIndxMapTy;
#else
	typedef std::map<uint64_t, std::list<unsigned> > lpBBid2EvalIterIndxMapTy;
#endif // End of ENABLE_BTREE_MAP
	extern lpBBid2EvalIterIndxMapTy lpBBid2EvalIterIndxMap;

	/// Header Basic Block ID --> loop bound Map
	typedef std::map<uint64_t, uint64_t> BBheaderID2loopBoundMapTy;
	extern BBheaderID2loopBoundMapTy BBheaderID2loopBoundMap;

	typedef std::pair<std::string, std::string> bbFuncNamePairTy;
	typedef std::pair<std::string, unsigned> lpNameLevelPairTy;
	typedef std::map<std::string, unsigned > LpName2numLevelMapTy;
	typedef std::map<bbFuncNamePairTy, lpNameLevelPairTy> bbFuncNamePair2lpNameLevelPairMapTy;

	extern bbFuncNamePair2lpNameLevelPairMapTy bbFuncNamePair2lpNameLevelPairMap;
	extern bbFuncNamePair2lpNameLevelPairMapTy headerBBFuncnamePair2lpNameLevelPairMap;
	extern bbFuncNamePair2lpNameLevelPairMapTy exitBBFuncnamePair2lpNameLevelPairMap;
	extern LpName2numLevelMapTy LpName2numLevelMap;

}

#endif // End of AUXILIARY_H