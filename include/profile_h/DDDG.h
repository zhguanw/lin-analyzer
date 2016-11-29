#ifndef __DDDG_H__
#define __DDDG_H__

#include <stack>
#include <string>
#include <zlib.h>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <set>
#include <map>
#include <fstream>

#include "profile_h/file_func.h"
#include "profile_h/opcode_func.h"
#include "profile_h/auxiliary.h"
/*#define HANDLE_INST(num, opc, clas) case num: return opc;*/
using namespace std;

struct edge_node_info{
  unsigned sink_node;
  int par_id;
};

//data structure used to tract dependency
typedef unordered_map<std::string, unsigned int> string_to_uint;
typedef unordered_map<long long int, unsigned int> uint_to_uint;
typedef unordered_multimap<unsigned int, edge_node_info> multi_uint_to_node_info;
typedef std::unordered_map<std::string, std::string> instName2bbNameMapTy;
typedef std::map<std::pair<std::string, std::string>, std::string> headerBBFuncNamePair2lastInstMapTy;
typedef std::pair<std::string, std::string> lpNameLevelStrPairTy;
typedef std::map<lpNameLevelStrPairTy, std::string> lpNameLevelPair2headBBnameMapTy;
/// (Function Name, Basic Block Name) --> number of instructions in a BB Map
typedef std::map<std::pair<std::string, std::string>, unsigned> funcBBNmPair2numInstInBBMapTy;
//typedef std::map<std::string, std::vector<uint64_t> > loopName2levelLpBoundVecMapTy;
typedef std::map<std::string, std::vector<unsigned> > loopName2levelUnrollVecMapTy;
typedef std::map<std::string, std::pair<std::string, unsigned> > headerBBlastInst2loopNameLevelPairMapTy;
typedef std::map<std::string, unsigned> wholeloopName2loopBoundMapTy;
typedef std::map<std::string, bool> wholeloopName2perfectOrNotMapTy;

extern instName2bbNameMapTy instName2bbNameMap;
extern headerBBFuncNamePair2lastInstMapTy headerBBFuncNamePair2lastInstMap;
extern headerBBFuncNamePair2lastInstMapTy exitingBBFuncNamePair2lastInstMap;
extern lpNameLevelPair2headBBnameMapTy lpNameLevelPair2headBBnameMap;
extern lpNameLevelPair2headBBnameMapTy lpNameLevelPair2exitingBBnameMap;
extern funcBBNmPair2numInstInBBMapTy funcBBNmPair2numInstInBBMap;
extern loopName2levelUnrollVecMapTy loopName2levelUnrollVecMap;
extern wholeloopName2loopBoundMapTy wholeloopName2loopBoundMap;
extern wholeloopName2perfectOrNotMapTy wholeloopName2perfectOrNotMap;
//extern loopName2levelLpBoundVecMapTy loopName2levelLpBoundVecMap;

namespace llvm {
	//funcBBNmPair2numInstInBBMapTy funcBBNmPair2numInstInBBMap;
}

typedef std::pair<uint64_t, uint64_t> line_from_to_Ty;

class BaseDatapath;

class DDDG
{

private:
  BaseDatapath *datapath;
	std::string inputPath;

public:
  DDDG(BaseDatapath *_datapath, std::string _trace_name, std::string input_path);
  int num_edges();
  int num_nodes();
  int num_of_register_dependency();
  int num_of_memory_dependency();
  void output_method_call_graph(std::string bench);
  void output_dddg();
  bool build_initial_dddg();

private:
	line_from_to_Ty getTraceLineFromTo(std::string loopName, unsigned loopLevel, unsigned unroll_factor);
	void extract_trace_file(gzFile& trace_file);
	void extract_trace_file(gzFile& trace_file, uint64_t from, uint64_t to);
  void parse_instruction_line(std::string line);
  void parse_parameter(std::string line, int param_tag);
  void parse_result(std::string line);
  void parse_forward(std::string line);
  void parse_call_parameter(std::string line, int param_tag);

  std::string trace_name;
  std::string curr_dynamic_function;

  uint8_t curr_microop;
  uint8_t prev_microop;
  std::string prev_bblock;
  std::string curr_bblock;

  std::string callee_function;
  std::string callee_dynamic_function;

  bool last_parameter;
  int num_of_parameters;
  //Used to track the instruction that initialize call function parameters
  int last_call_source;

  std::string curr_instid;
  std::vector<long long int> parameter_value_per_inst;
  std::vector<unsigned> parameter_size_per_inst;
  std::vector<std::string> parameter_label_per_inst;
  std::vector<std::string> method_call_graph;
  /*unordered_map<unsigned, bool> to_ignore_methodid;*/
  int num_of_instructions;
  int num_of_reg_dep;
  int num_of_mem_dep;

	//register dependency tracking table using hash_map(hash_map)
	//memory dependency tracking table
	//edge multimap
	multi_uint_to_node_info register_edge_table;
	multi_uint_to_node_info memory_edge_table;
	//keep track of currently executed methods
	stack<std::string> active_method;
	//manage methods
	string_to_uint function_counter;
  string_to_uint register_last_written;
	uint_to_uint address_last_written;
	//loopName2levelLpBoundVecMapTy loopName2levelLpBoundVecMap;
};

#endif
