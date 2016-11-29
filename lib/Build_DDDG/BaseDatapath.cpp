#include <sstream>
#include <fstream>
#include <boost/tokenizer.hpp>

#include "profile_h/opcode_func.h"
#include "profile_h/BaseDatapath.h"
#include "llvm/Support/GraphWriter.h"
//#include "profile_h/global_variables.h"

using namespace llvm;

const std::string output_title = "type,function name,loop name,loop level,loop latency,IL,II,Res_mem,Res_op,Rec_II,DSP,BRAM18K,FF,LUT,Fadd,Fsub,Fmul,Fdiv,En_sharedLoad,shared_loads,repeated_stores,enable_TreeHeightReduction";

bool compare_ready_queue_basedon_alapTime(const std::pair<unsigned, unsigned>& first, const std::pair<unsigned, unsigned>& second) {
	return (first.second < second.second);
}

bool compare_arrayName2resII(const std::pair<std::string, double>& first, const std::pair<std::string, double>& second) {
	return (first.second < second.second);
}

BaseDatapath::BaseDatapath(std::string bench, string trace_file, string config_file, string input_path, std::string target_loop, unsigned lp_level, unsigned target_unroll_factor, unsigned IL_asap)
{
  benchName = bench;
	inputPath = input_path;
	//configure_fileName = inputPath + "config_example";
	target_loop_name = target_loop;
	target_loop_level = lp_level;
	target_lp_level_unroll_factor = target_unroll_factor;
	IL_asap_ii = IL_asap;
  DDDG *dddg;
  dddg = new DDDG(this, trace_file, inputPath);
  /*Build Initial DDDG*/
  if (dddg->build_initial_dddg())
  {
		std::cout << "DEBUG-INFO: [trace-analysis_trace-generation] Cannot build Dynamic Data Dependence Graph\n";
    exit(0);
  }
  delete dddg;
	std::cout << "DEBUG-INFO: [DDDG-analysis] Analyzing DDDG\n";
  numTotalNodes = microop.size();

  BGL_FORALL_VERTICES(v, graph_, Graph)
    nameToVertex[get(boost::vertex_index, graph_, v)] = v;
  vertexToName = get(boost::vertex_index, graph_);

  std::vector<std::string> dynamic_methodid(numTotalNodes, "");
  initDynamicMethodID(dynamic_methodid);

  for (auto dynamic_func_it = dynamic_methodid.begin(), E = dynamic_methodid.end();
       dynamic_func_it != E; dynamic_func_it++)
  {
    char func_id[256];
    int count;
    sscanf((*dynamic_func_it).c_str(), "%[^-]-%d\n", func_id, &count);
		if (functionNames.find(func_id) == functionNames.end()) {
			functionNames.insert(func_id);
		}
  }
	
  //parse_config(bench, config_file);

	num_cycles = 0;

	///FIXME: We set numOfPortsPerPartition to 1000, so that we do not have memory port limitations. 
	/// 1000 ports are sufficient. 
	/// Later, we need to add memory port limitation below (struct will be better) to take read/write
	/// ports into consideration.
	numOfPortsPerPartition = 1000;
}

BaseDatapath::~BaseDatapath() {}

void BaseDatapath::addDddgEdge(unsigned int from, unsigned int to, uint8_t parid)
{
  if (from != to)
    add_edge(from, to, EdgeProperty(parid), graph_);
}

void BaseDatapath::insertMicroop(int node_microop) {
	microop.push_back(node_microop);
}

//optimizationFunctions
void BaseDatapath::setGlobalGraph()
{

  std::cerr << "=============================================" << std::endl;
  std::cerr << "      Optimizing...            " << benchName << std::endl;
  std::cerr << "=============================================" << std::endl;
  finalIsolated.assign(numTotalNodes, 1);
}

void BaseDatapath::initialize_graph() {
	// Perform initialization before analysis
	setGlobalGraph();
}

void BaseDatapath::memoryAmbiguation()
{
	std::cout << "DEBUG-INFO: [trace-analysis_memory-disambiguation] Perform memory disambiguation on DDDG\n";

  std::unordered_multimap<std::string, std::string> pair_per_load;
  std::unordered_set<std::string> paired_store;
  std::unordered_map<std::string, bool> store_load_pair;

  std::vector<std::string> instid(numTotalNodes, "");
  std::vector<std::string> dynamic_methodid(numTotalNodes, "");
  std::vector<std::string> prev_basic_block(numTotalNodes, "");

  initInstID(instid);
  initDynamicMethodID(dynamic_methodid);
  initPrevBasicBlock(prev_basic_block);
  std::vector< Vertex > topo_nodes;
  boost::topological_sort(graph_, std::back_inserter(topo_nodes));
  //nodes with no incoming edges to first
  for (auto vi = topo_nodes.rbegin(); vi != topo_nodes.rend(); ++vi)
  {
    unsigned node_id = vertexToName[*vi];

    int node_microop = microop.at(node_id);
    if (!is_store_op(node_microop))
      continue;
    //iterate its children to find a load op
    out_edge_iter out_edge_it, out_edge_end;
    for (tie(out_edge_it, out_edge_end) = out_edges(*vi, graph_); out_edge_it != out_edge_end; ++out_edge_it)
    {
      int child_id = vertexToName[target(*out_edge_it, graph_)];
      int child_microop = microop.at(child_id);
      if (!is_load_op(child_microop))
        continue;
      std::string node_dynamic_methodid = dynamic_methodid.at(node_id);
      std::string load_dynamic_methodid = dynamic_methodid.at(child_id);
      if (node_dynamic_methodid.compare(load_dynamic_methodid) != 0)
        continue;

      std::string store_unique_id (node_dynamic_methodid + "-" + instid.at(node_id) + "-" + prev_basic_block.at(node_id));
      std::string load_unique_id (load_dynamic_methodid+ "-" + instid.at(child_id) + "-" + prev_basic_block.at(child_id));

      if (store_load_pair.find(store_unique_id + "-" + load_unique_id ) != store_load_pair.end())
        continue;
      //add to the pair
      store_load_pair[store_unique_id + "-" + load_unique_id] = 1;
      paired_store.insert(store_unique_id);
      auto load_range = pair_per_load.equal_range(load_unique_id);
      bool found_store = 0;
      for (auto store_it = load_range.first; store_it != load_range.second; store_it++)
      {
        if (store_unique_id.compare(store_it->second) == 0)
        {
          found_store = 1;
          break;
        }
      }
      if (!found_store)
      {
        pair_per_load.insert(make_pair(load_unique_id,store_unique_id));
      }
    }
  }
  if (store_load_pair.size() == 0)
    return;

  std::vector<newEdge> to_add_edges;
  std::unordered_map<std::string, unsigned> last_store;

  for (unsigned node_id = 0; node_id < numTotalNodes; node_id++)
  {
    int node_microop = microop.at(node_id);
    if (!is_memory_op(node_microop))
      continue;
    std::string unique_id (dynamic_methodid.at(node_id) + "-" + instid.at(node_id) + "-" + prev_basic_block.at(node_id));
    if (is_store_op(node_microop))
    {
      auto store_it = paired_store.find(unique_id);
      if (store_it == paired_store.end())
        continue;
      last_store[unique_id] = node_id;
    }
    else
    {
      assert(is_load_op(node_microop));
      auto load_range = pair_per_load.equal_range(unique_id);
      if (std::distance(load_range.first, load_range.second) == 1)
        continue;
      for (auto load_store_it = load_range.first; load_store_it != load_range.second; ++load_store_it)
      {
        assert(paired_store.find(load_store_it->second) != paired_store.end());
        auto prev_store_it = last_store.find(load_store_it->second);
        if (prev_store_it == last_store.end())
          continue;
        unsigned prev_store_id = prev_store_it->second;
        if (!doesEdgeExist(prev_store_id, node_id))
        {
          to_add_edges.push_back({prev_store_id, node_id, -1});
          dynamicMemoryOps.insert(load_store_it->second + "-" + prev_basic_block.at(prev_store_id));
          dynamicMemoryOps.insert(load_store_it->first + "-" + prev_basic_block.at(node_id));
        }
      }
    }
  }
  updateGraphWithNewEdges(to_add_edges);
}

/*
 * Read: graph_, microop
 * Modify: graph_
 */
void BaseDatapath::removePhiNodes()
{
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_PHINodes-removal] Remove PHI nodes on DDDG\n");

  EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

  std::set<Edge> to_remove_edges;
  std::vector<newEdge> to_add_edges;

  vertex_iter vi, vi_end;
  int removed_phi = 0;
  for (tie(vi, vi_end) = vertices(graph_); vi != vi_end; ++vi)
  {
    unsigned node_id = vertexToName[*vi];
    int node_microop = microop.at(node_id);
    if (node_microop != LLVM_IR_PHI && node_microop != LLVM_IR_BitCast)
      continue;
    //find its children

    std::vector< pair<unsigned, int> > phi_child;

    out_edge_iter out_edge_it, out_edge_end;
    for (tie(out_edge_it, out_edge_end) = out_edges(*vi, graph_);
      out_edge_it != out_edge_end; ++out_edge_it)
    {
      to_remove_edges.insert(*out_edge_it);
      phi_child.push_back(make_pair(vertexToName[target(*out_edge_it, graph_)],
                                     edge_to_parid[*out_edge_it]));
    }
    if (phi_child.size() == 0)
      continue;
    //find its parents
    in_edge_iter in_edge_it, in_edge_end;
    for (tie(in_edge_it, in_edge_end) = in_edges(*vi, graph_);
      in_edge_it != in_edge_end; ++in_edge_it)
    {
      unsigned parent_id = vertexToName[source(*in_edge_it, graph_)];
      to_remove_edges.insert(*in_edge_it);

      for (auto child_it = phi_child.begin(), chil_E = phi_child.end();
        child_it != chil_E; ++child_it)
        to_add_edges.push_back({parent_id, child_it->first, child_it->second});
    }
		// To release memory allocated for phi_child vector, if using std::vector::clear(),
		// it will not free up the memory. Therefore, the following code uses an anonymous
		// temporary.
    std::vector<pair<unsigned, int> >().swap(phi_child);
    removed_phi++;
  }

  updateGraphWithIsolatedEdges(to_remove_edges);
  updateGraphWithNewEdges(to_add_edges);
  //cleanLeafNodes();
}

/*
 * Read: lineNum.gz, flattenConfig, microop
 * Modify: graph_
 */
void BaseDatapath::loopFlatten()
{
  std::unordered_set<int> flatten_config;
  if (!readFlattenConfig(flatten_config))
    return;
  std::cerr << "-------------------------------" << std::endl;
  std::cerr << "         Loop Flatten          " << std::endl;
  std::cerr << "-------------------------------" << std::endl;
  std::vector<int> lineNum(numTotalNodes, -1);
  initLineNum(lineNum);

  std::vector<unsigned> to_remove_nodes;

  for(unsigned node_id = 0; node_id < numTotalNodes; node_id++)
  {
    int node_linenum = lineNum.at(node_id);
    auto it = flatten_config.find(node_linenum);
    if (it == flatten_config.end())
      continue;
    if (is_compute_op(microop.at(node_id)))
      microop.at(node_id) = LLVM_IR_Move;
    else if (is_branch_op(microop.at(node_id)))
      to_remove_nodes.push_back(node_id);
  }
  updateGraphWithIsolatedNodes(to_remove_nodes);
  //cleanLeafNodes();
}

void BaseDatapath::cleanLeafNodes()
{
  EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

  /*track the number of children each node has*/
  std::vector<int> num_of_children(numTotalNodes, 0);
  std::vector<unsigned> to_remove_nodes;

  std::vector< Vertex > topo_nodes;
  boost::topological_sort(graph_, std::back_inserter(topo_nodes));
  //bottom nodes first
  for (auto vi = topo_nodes.begin(); vi != topo_nodes.end(); ++vi)
  {
    Vertex node_vertex = *vi;
    if (boost::degree(node_vertex, graph_) == 0)
      continue;
    unsigned  node_id = vertexToName[node_vertex];
    int node_microop = microop.at(node_id);
    if (num_of_children.at(node_id) == boost::out_degree(node_vertex, graph_)
      && node_microop != LLVM_IR_SilentStore
      && node_microop != LLVM_IR_Store
      && node_microop != LLVM_IR_Ret
      && !is_branch_op(node_microop))
    {
      to_remove_nodes.push_back(node_id);
      //iterate its parents
      in_edge_iter in_edge_it, in_edge_end;
      for (tie(in_edge_it, in_edge_end) = in_edges(node_vertex, graph_); in_edge_it != in_edge_end; ++in_edge_it)
      {
        int parent_id = vertexToName[source(*in_edge_it, graph_)];
        num_of_children.at(parent_id)++;
      }
    }
    else if (is_branch_op(node_microop))
    {
      //iterate its parents
      in_edge_iter in_edge_it, in_edge_end;
      for (tie(in_edge_it, in_edge_end) = in_edges(node_vertex, graph_); in_edge_it != in_edge_end; ++in_edge_it)
      {
        if (edge_to_parid[*in_edge_it] == CONTROL_EDGE)
        {
          int parent_id = vertexToName[source(*in_edge_it, graph_)];
          num_of_children.at(parent_id)++;
        }
      }
    }
  }
  updateGraphWithIsolatedNodes(to_remove_nodes);
}

/*
 * Read: graph_, instid, microop
 * Modify: microop
 */
void BaseDatapath::removeInductionDependence()
{
  //set graph
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_induction-variables-removal] Remove induction variables on DDDG\n");

  std::vector<std::string> instid(numTotalNodes, "");
  initInstID(instid);

  std::vector< Vertex > topo_nodes;
  boost::topological_sort(graph_, std::back_inserter(topo_nodes));
  //nodes with no incoming edges to first
  for (auto vi = topo_nodes.rbegin(); vi != topo_nodes.rend(); ++vi)
  {
    unsigned node_id = vertexToName[*vi];
    std::string node_instid = instid.at(node_id);

		if (node_instid.find("indvars") != std::string::npos) {
			if (microop.at(node_id) == LLVM_IR_Add) {
				microop.at(node_id) = LLVM_IR_IndexAdd;
			}
		}
		else {
			in_edge_iter in_edge_it, in_edge_ie;
			for (tie(in_edge_it, in_edge_ie) = in_edges(*vi, graph_); in_edge_it != in_edge_ie; ++in_edge_it) {
				Vertex parent_vertex = source(*in_edge_it, graph_);
				unsigned parent_id = vertexToName[parent_vertex];
				std::string node_name = instid.at(parent_id);
				if (node_name.find("indvars") == std::string::npos) {
					continue;
				}
				if (microop.at(node_id) == LLVM_IR_Add) {
					microop.at(node_id) = LLVM_IR_IndexAdd;
				}

			}
		}
    
  }
}

//called in the end of the whole flow
void BaseDatapath::dumpStats()
{
  clearGraph();
  dumpGraph();
  writeMicroop(microop);
  writeFinalLevel();
  writeGlobalIsolated();
}

#ifdef ALADDIN_H
void BaseDatapath::loopPipelining()
{
  if (!readPipeliningConfig())
  {
    std::cerr << "Loop Pipelining is not ON." << std::endl;
    return ;
  }

  std::unordered_map<int, int > unrolling_config;
  if (!readUnrollingConfig(unrolling_config))
  {
    std::cerr << "Loop Unrolling is not defined. " << std::endl;
    std::cerr << "Loop pipelining is only applied to unrolled loops." << std::endl;
    return ;
  }

  if (loopBound.size() <= 2)
    return;
  std::cerr << "-------------------------------" << std::endl;
  std::cerr << "         Loop Pipelining        " << std::endl;
  std::cerr << "-------------------------------" << std::endl;

  EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

  vertex_iter vi, vi_end;
  std::set<Edge> to_remove_edges;
  std::vector<newEdge> to_add_edges;

  //After loop unrolling, we define strict control dependences between basic block,
  //where all the instructions in the following basic block depend on the prev branch instruction
  //To support loop pipelining, which allows the next iteration
  //starting without waiting until the prev iteration finish, we move the control dependences
  //between last branch node in the prev basic block and instructions in the next basic block
  //to first non isolated instruction in the prev basic block and instructions in the next basic block...
  std::map<unsigned, unsigned> first_non_isolated_node;
  auto bound_it = loopBound.begin();
  unsigned node_id = *bound_it;
  //skip first region
  bound_it++;
  while ( (unsigned)node_id < numTotalNodes)
  {
    assert(is_branch_op(microop.at(*bound_it)));
    while (node_id < *bound_it &&  (unsigned) node_id < numTotalNodes)
    {
      if (nameToVertex.find(node_id) == nameToVertex.end()
          || boost::degree(nameToVertex[node_id], graph_) == 0
          || is_branch_op(microop.at(node_id)) ) {
        node_id++;
        continue;
      }
      else {
        first_non_isolated_node[*bound_it] = node_id;
        node_id = *bound_it;
        break;
      }
    }
    if (first_non_isolated_node.find(*bound_it) == first_non_isolated_node.end())
      first_non_isolated_node[*bound_it] = *bound_it;
    bound_it++;
    if (bound_it == loopBound.end() - 1 )
      break;
  }
  int prev_branch = -1;
  int prev_first = -1;
  for(auto first_it = first_non_isolated_node.begin(), E = first_non_isolated_node.end(); first_it != E; ++first_it)
  {
    unsigned br_node = first_it->first;
    unsigned first_id = first_it->second;
    if (is_call_op(microop.at(br_node))) {
      prev_branch = -1;
      continue;
    }
    if (prev_branch != -1)
    {
      //adding dependence between prev_first and first_id
      if (!doesEdgeExist(prev_first, first_id))
        to_add_edges.push_back({(unsigned)prev_first, first_id, CONTROL_EDGE});
      //adding dependence between first_id and prev_branch's children
      out_edge_iter out_edge_it, out_edge_end;
      for (tie(out_edge_it, out_edge_end) = out_edges(nameToVertex[prev_branch], graph_); out_edge_it != out_edge_end; ++out_edge_it)
      {
        Vertex child_vertex = target(*out_edge_it, graph_);
        unsigned child_id = vertexToName[child_vertex];
        if (child_id <= first_id
            || edge_to_parid[*out_edge_it] != CONTROL_EDGE)
          continue;
        if (!doesEdgeExist(first_id, child_id))
          to_add_edges.push_back({first_id, child_id, 1});
      }
    }
    //update first_id's parents, dependence become strict control dependence
    in_edge_iter in_edge_it, in_edge_end;
    for (tie(in_edge_it, in_edge_end) = in_edges(nameToVertex[first_id], graph_); in_edge_it != in_edge_end; ++in_edge_it)
    {
      Vertex parent_vertex = source(*in_edge_it, graph_);
      unsigned parent_id = vertexToName[parent_vertex];
      if (is_branch_op(microop.at(parent_id)))
        continue;
      to_remove_edges.insert(*in_edge_it);
      to_add_edges.push_back({parent_id, first_id, CONTROL_EDGE});
    }
    //remove control dependence between br node to its children
    out_edge_iter out_edge_it, out_edge_end;
    for (tie(out_edge_it, out_edge_end) = out_edges(nameToVertex[br_node], graph_); out_edge_it != out_edge_end; ++out_edge_it) {
      if (is_call_op(microop.at(vertexToName[target(*out_edge_it, graph_)])))
        continue;
      if (edge_to_parid[*out_edge_it] != CONTROL_EDGE)
        continue;
      to_remove_edges.insert(*out_edge_it);
    }
    prev_branch = br_node;
    prev_first = first_id;
  }

  updateGraphWithIsolatedEdges(to_remove_edges);
  updateGraphWithNewEdges(to_add_edges);
  cleanLeafNodes();
}

#endif // End of ALADDIN_H
/*
 * Read: graph_, lineNum.gz, unrollingConfig, microop
 * Modify: graph_
 * Write: loop_bound
 */
void BaseDatapath::loopUnrolling()
{
  std::unordered_map<int, int > unrolling_config;
  readUnrollingConfig(unrolling_config);

  std::cerr << "-------------------------------" << std::endl;
  std::cerr << "         Loop Unrolling        " << std::endl;
  std::cerr << "-------------------------------" << std::endl;

  std::vector<unsigned> to_remove_nodes;
  std::unordered_map<std::string, unsigned> inst_dynamic_counts;
  std::vector<unsigned> nodes_between;
  std::vector<newEdge> to_add_edges;
  std::vector<int> lineNum(numTotalNodes, -1);
  initLineNum(lineNum);

  bool first = false;
  int iter_counts = 0;
  int prev_branch = -1;

  for(unsigned node_id = 0; node_id < numTotalNodes; node_id++)
  {
    if (nameToVertex.find(node_id) == nameToVertex.end())
      continue;
    Vertex node_vertex = nameToVertex[node_id];
    if (boost::degree(node_vertex, graph_) == 0
       && !is_call_op(microop.at(node_id)))
      continue;
    if (!first)
    {
      first = true;
      loopBound.push_back(node_id);
      prev_branch = node_id;
    }
    assert(prev_branch != -1);
    if (prev_branch != node_id &&
      !(is_dma_op(microop.at(prev_branch)) && is_dma_op(microop.at(node_id))) ) {
      to_add_edges.push_back({(unsigned)prev_branch, node_id, CONTROL_EDGE});
    }
    if (!is_branch_op(microop.at(node_id)))
      nodes_between.push_back(node_id);
    else
    {
      //for the case that the first non-isolated node is also a call node;
      if (is_call_op(microop.at(node_id)) && *loopBound.rbegin() != node_id)
      {
        loopBound.push_back(node_id);
        prev_branch = node_id;
      }

      int node_linenum = lineNum.at(node_id);
      auto unroll_it = unrolling_config.find(node_linenum);
      //not unrolling branch
      if (unroll_it == unrolling_config.end())
      {
        if (!is_call_op(microop.at(node_id))) {
          nodes_between.push_back(node_id);
          continue;
        }
        // Enforce dependences between branch nodes, including call nodes
        // Except for the case that both two branches are DMA operations.
        // (Two DMA operations can go in parallel.)
        if (!doesEdgeExist(prev_branch, node_id) &&
              !( is_dma_op(microop.at(prev_branch)) &&
                  is_dma_op(microop.at(node_id)) ) )
          to_add_edges.push_back({(unsigned)prev_branch, node_id, CONTROL_EDGE});
        for (auto prev_node_it = nodes_between.begin(), E = nodes_between.end();
                   prev_node_it != E; prev_node_it++)
        {
          if (!doesEdgeExist(*prev_node_it, node_id) &&
              !( is_dma_op(microop.at(*prev_node_it)) &&
                   is_dma_op(microop.at(node_id)) )  ) {
            to_add_edges.push_back({*prev_node_it, node_id, CONTROL_EDGE});
          }
        }
        nodes_between.clear();
        nodes_between.push_back(node_id);
        prev_branch = node_id;
      }
      else
      {
        int factor = unroll_it->second;
        int node_microop = microop.at(node_id);
        char unique_inst_id[256];
        sprintf(unique_inst_id, "%d-%d", node_microop, node_linenum);
        auto it = inst_dynamic_counts.find(unique_inst_id);
        if (it == inst_dynamic_counts.end())
        {
          inst_dynamic_counts[unique_inst_id] = 1;
          it = inst_dynamic_counts.find(unique_inst_id);
        }
        else
          it->second++;
        if (it->second % factor == 0)
        {
          loopBound.push_back(node_id);
          iter_counts++;
          for (auto prev_node_it = nodes_between.begin(), E = nodes_between.end();
                     prev_node_it != E; prev_node_it++)
          {
            if (!doesEdgeExist(*prev_node_it, node_id)) {
              to_add_edges.push_back({*prev_node_it, node_id, CONTROL_EDGE});
            }
          }
          nodes_between.clear();
          nodes_between.push_back(node_id);
          prev_branch = node_id;
        }
        else
          to_remove_nodes.push_back(node_id);
      }
    }
  }
  loopBound.push_back(numTotalNodes);

  if (iter_counts == 0 && unrolling_config.size() != 0 )
  {
    std::cerr << "-------------------------------" << std::endl;
    std::cerr << "Loop Unrolling Factor is Larger than the Loop Trip Count."
              << std::endl;
    std::cerr << "Loop Unrolling is NOT applied. Please choose a smaller "
              << "unrolling factor." << std::endl;
    std::cerr << "-------------------------------" << std::endl;
  }
  updateGraphWithNewEdges(to_add_edges);
  updateGraphWithIsolatedNodes(to_remove_nodes);
  cleanLeafNodes();
}


void BaseDatapath::removeSharedLoads() {
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_shared-load-removal] Remove shared loads inside DDDG" << std::endl);

  EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

  std::unordered_map<unsigned, long long int> address;
  initAddress(address);

  vertex_iter vi, vi_end;

  std::set<Edge> to_remove_edges;
  std::vector<newEdge> to_add_edges;

  shared_loads = 0;

  unsigned node_id = 0;
	std::unordered_map<unsigned, unsigned> address_loaded;
  while (node_id < numTotalNodes)
  {
    
    if (nameToVertex.find(node_id) == nameToVertex.end()){
			// If some vertices are deleted from DDDG, we need to check whether they are
			// still resided inside the DDDG
      node_id++;
      continue;
    }

    if (boost::degree(nameToVertex[node_id], graph_) == 0){
      node_id++;
      continue;
    }

    int node_microop = microop.at(node_id);
    long long int node_address = address[node_id];
    auto addr_it = address_loaded.find(node_address);
		if (is_store_op(node_microop) && addr_it != address_loaded.end()){
			address_loaded.erase(addr_it);
		}
    else if (is_load_op(node_microop)) {
			if (addr_it == address_loaded.end()){
				address_loaded.insert(std::make_pair(node_address, node_id));
			}
      else {
        shared_loads++;
        microop.at(node_id) = LLVM_IR_Move;
        unsigned prev_load = addr_it->second;
        //iterate through its children
        Vertex load_node = nameToVertex[node_id];
        out_edge_iter out_edge_it, out_edge_end;
        for (tie(out_edge_it, out_edge_end) = out_edges(load_node, graph_); out_edge_it != out_edge_end; ++out_edge_it)
        {
          Edge curr_edge = *out_edge_it;
          Vertex child_vertex = target(curr_edge, graph_);
          unsigned child_id = vertexToName[child_vertex];
          Vertex prev_load_vertex = nameToVertex[prev_load];
          if (!doesEdgeExistVertex(prev_load_vertex, child_vertex))
            to_add_edges.push_back({prev_load, child_id, edge_to_parid[curr_edge]});
          to_remove_edges.insert(*out_edge_it);
        }
        in_edge_iter in_edge_it, in_edge_end;
        for (tie(in_edge_it, in_edge_end) = in_edges(load_node, graph_); in_edge_it != in_edge_end; ++in_edge_it)
          to_remove_edges.insert(*in_edge_it);
      }
		}
		else {
			// Do nothing here.
		}
    node_id++;
  }

  updateGraphWithIsolatedEdges(to_remove_edges);
  updateGraphWithNewEdges(to_add_edges);
  //cleanLeafNodes();

	//VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_shared-load-removal] Finished" << std::endl);
}

/*
 * Read: loopBound, flattenConfig, graph_, instid, dynamicMethodID,
 *       prevBasicBlock
 * Modify: graph_
 */
void BaseDatapath::storeBuffer() {
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_store-buffer] Analyze store buffers on DDDG" << std::endl);

  EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

  std::vector<std::string> instid(numTotalNodes, "");
  std::vector<std::string> dynamic_methodid(numTotalNodes, "");
  std::vector<std::string> prev_basic_block(numTotalNodes, "");

  initInstID(instid);
  initDynamicMethodID(dynamic_methodid);
  initPrevBasicBlock(prev_basic_block);

  std::vector<newEdge> to_add_edges;
  std::vector<unsigned> to_remove_nodes;

  unsigned node_id = 0;
	for (unsigned node_id = 0; node_id < numTotalNodes; node_id++) {

    if (nameToVertex.find(node_id) == nameToVertex.end()
        || boost::degree(nameToVertex[node_id], graph_) == 0) {
      ++node_id;
      continue;
    }

    if (is_store_op(microop.at(node_id))) {
      //remove this store
      std::string store_unique_id (dynamic_methodid.at(node_id) + "-" + instid.at(node_id) + "-" + prev_basic_block.at(node_id));
      //dynamic stores, cannot disambiguated in the static time, cannot remove
      if (dynamicMemoryOps.find(store_unique_id) != dynamicMemoryOps.end()) {
        ++node_id;
        continue;
      }
      Vertex node = nameToVertex[node_id];
      out_edge_iter out_edge_it, out_edge_end;

      std::vector<Vertex> store_child;
      for (tie(out_edge_it, out_edge_end) = out_edges(node, graph_);
            out_edge_it != out_edge_end; ++out_edge_it) {
        Vertex child_vertex = target(*out_edge_it, graph_);
        int child_id = vertexToName[child_vertex];
        if (is_load_op(microop.at(child_id))) {
          std::string load_unique_id (dynamic_methodid.at(child_id) + "-"
                + instid.at(child_id) + "-" + prev_basic_block.at(child_id));
          if (dynamicMemoryOps.find(load_unique_id) != dynamicMemoryOps.end())
            continue;
          else
            store_child.push_back(child_vertex);
        }
      }

      if (store_child.size() > 0) {
        bool parent_found = false;
        Vertex store_parent;
        in_edge_iter in_edge_it, in_edge_end;
        for (tie(in_edge_it, in_edge_end) = in_edges(node, graph_);
          in_edge_it != in_edge_end; ++in_edge_it) {
          //parent node that generates value
          if (edge_to_parid[*in_edge_it] == 1) {
            parent_found = true;
            store_parent = source(*in_edge_it, graph_);
            break;
          }
        }

        if (parent_found) {
          for (auto load_it = store_child.begin(), E = store_child.end();
            load_it != E; ++load_it) {
            Vertex load_node = *load_it;
            to_remove_nodes.push_back(vertexToName[load_node]);

            out_edge_iter out_edge_it, out_edge_end;
						for (tie(out_edge_it, out_edge_end) = out_edges(load_node, graph_);
							out_edge_it != out_edge_end; ++out_edge_it) {
							to_add_edges.push_back({ (unsigned)vertexToName[store_parent],
								(unsigned)vertexToName[target(*out_edge_it, graph_)],
								edge_to_parid[*out_edge_it] });
						}
          }
        }

      }

    }

  }
  updateGraphWithNewEdges(to_add_edges);
  updateGraphWithIsolatedNodes(to_remove_nodes);
  //cleanLeafNodes();
}

/*
 * Read: loopBound, flattenConfig, graph_, address, instid, dynamicMethodID,
 *       prevBasicBlock
 * Modify: graph_
 */
void BaseDatapath::removeRepeatedStores() {
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_repeated-store-removal] Remove repeated stores on DDDG" << std::endl);

  std::unordered_map<unsigned, long long int> address;
  initAddress(address);

  std::vector<std::string> instid(numTotalNodes, "");
  std::vector<std::string> dynamic_methodid(numTotalNodes, "");
  std::vector<std::string> prev_basic_block(numTotalNodes, "");

  initInstID(instid);
  initDynamicMethodID(dynamic_methodid);
  initPrevBasicBlock(prev_basic_block);

  repeated_stores = 0;
  int node_id = numTotalNodes - 1;
	unordered_map<unsigned, int> address_store_map;
  while (node_id >=0 ) {

    if (nameToVertex.find(node_id) == nameToVertex.end()
        || boost::degree(nameToVertex[node_id], graph_) == 0
        || !is_store_op(microop.at(node_id))) {
      --node_id;
      continue;
    }
    long long int node_address = address[node_id];
    auto addr_it = address_store_map.find(node_address);

    if (addr_it == address_store_map.end())
      address_store_map[node_address] = node_id;
    else {
      //remove this store
      std::string store_unique_id (dynamic_methodid.at(node_id) + "-" + instid.at(node_id) + "-" + prev_basic_block.at(node_id));
      //dynamic stores, cannot disambiguated in the run time, cannot remove
      if (dynamicMemoryOps.find(store_unique_id) == dynamicMemoryOps.end()
          && boost::out_degree(nameToVertex[node_id], graph_)== 0) {
          microop.at(node_id) = LLVM_IR_SilentStore;
          repeated_stores++;
      }
    }
    --node_id;

  }
  //cleanLeafNodes();
}

/*
 * Read: loopBound, flattenConfig, graph_, microop
 * Modify: graph_
 */
void BaseDatapath::treeHeightReduction_integer()
{
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_tree-height-reduction-integer] Perform tree height reduction for integer operations on DDDG" << std::endl);

  EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

  std::vector<bool> updated(numTotalNodes, 0);
  std::vector<int> bound_region(numTotalNodes, 0);

	/*
  int region_id = 0;
  unsigned node_id = 0;
  auto bound_it = loopBound.begin();
  while (node_id < *bound_it)
  {
    bound_region.at(node_id) = region_id;
    node_id++;
    if (node_id == *bound_it)
    {
      region_id++;
      bound_it++;
      if (bound_it == loopBound.end())
        break;
    }
  }
	*/

  std::set<Edge> to_remove_edges;
  std::vector<newEdge> to_add_edges;

  //nodes with no outgoing edges to first (bottom nodes first)
  for(int node_id = numTotalNodes -1; node_id >= 0; node_id--)
  {
    if (nameToVertex.find(node_id) == nameToVertex.end()
       || boost::degree(nameToVertex[node_id], graph_) == 0
       || updated.at(node_id)
       || !is_associative(microop.at(node_id)) )
      continue;
    updated.at(node_id) = 1;
    int node_region = bound_region.at(node_id);

    std::list<unsigned> nodes;
    std::vector<Edge> tmp_remove_edges;
    std::vector<pair<int, bool> > leaves;
    std::vector<int> associative_chain;
    associative_chain.push_back(node_id);
    int chain_id = 0;
    while (chain_id < associative_chain.size())
    {
      int chain_node_id = associative_chain.at(chain_id);
      int chain_node_microop = microop.at(chain_node_id);
      if (is_associative(chain_node_microop))
      {
        updated.at(chain_node_id) = 1;
        int num_of_chain_parents = 0;
        in_edge_iter in_edge_it, in_edge_end;
        for (tie(in_edge_it, in_edge_end) = in_edges(nameToVertex[chain_node_id] , graph_); in_edge_it != in_edge_end; ++in_edge_it)
        {
          int parent_id = vertexToName[source(*in_edge_it, graph_)];
          if (is_branch_op(microop.at(parent_id)))
            continue;
          num_of_chain_parents++;
        }
        if (num_of_chain_parents == 2)
        {
          nodes.push_front(chain_node_id);
          for (tie(in_edge_it, in_edge_end) = in_edges(nameToVertex[chain_node_id] , graph_); in_edge_it != in_edge_end; ++in_edge_it)
          {
            Vertex parent_node = source(*in_edge_it, graph_);
            int parent_id = vertexToName[parent_node];
            assert(parent_id < chain_node_id);
            int parent_region = bound_region.at(parent_id);
            int parent_microop = microop.at(parent_id);
            if (is_branch_op(parent_microop))
              continue;
            Edge curr_edge = *in_edge_it;
            tmp_remove_edges.push_back(curr_edge);

            if (parent_region == node_region)
            {
              updated.at(parent_id) = 1;
              if (!is_associative(parent_microop))
                leaves.push_back(make_pair(parent_id, 0));
              else
              {
                out_edge_iter out_edge_it, out_edge_end;
                int num_of_children = 0;
                for (tie(out_edge_it, out_edge_end) = out_edges(parent_node, graph_); out_edge_it != out_edge_end; ++out_edge_it) {
                  if (edge_to_parid[*out_edge_it] != CONTROL_EDGE)
                    num_of_children++;
                }

                if (num_of_children == 1)
                  associative_chain.push_back(parent_id);
                else
                  leaves.push_back(make_pair(parent_id, 0));
              }
            }
            else
              leaves.push_back(make_pair(parent_id, 1));
          }
        }
        else
          leaves.push_back(make_pair(chain_node_id, 0));
      }
      else
        leaves.push_back(make_pair(chain_node_id, 0));
      chain_id++;
    }
    //build the tree
    if (nodes.size() < 3)
      continue;

    for(auto it = tmp_remove_edges.begin(), E = tmp_remove_edges.end(); it != E; it++)
      to_remove_edges.insert(*it);

    std::map<unsigned, unsigned> rank_map;
    auto leaf_it = leaves.begin();

    while (leaf_it != leaves.end())
    {
      if (leaf_it->second == 0)
        rank_map[leaf_it->first] = 0;
      else
        rank_map[leaf_it->first] = numTotalNodes;
      ++leaf_it;
    }
    //reconstruct the rest of the balanced tree
    auto node_it = nodes.begin();

    while (node_it != nodes.end())
    {
      unsigned node1, node2;
      if (rank_map.size() == 2)
      {
        node1 = rank_map.begin()->first;
        node2 = (++rank_map.begin())->first;
      }
      else
        findMinRankNodes(node1, node2, rank_map);
      assert((node1 != numTotalNodes) && (node2 != numTotalNodes));
      to_add_edges.push_back({node1, *node_it, 1});
      to_add_edges.push_back({node2, *node_it, 1});

      //place the new node in the map, remove the two old nodes
      rank_map[*node_it] = max(rank_map[node1], rank_map[node2]) + 1;
      rank_map.erase(node1);
      rank_map.erase(node2);
      ++node_it;
    }
  }
  updateGraphWithIsolatedEdges(to_remove_edges);
  updateGraphWithNewEdges(to_add_edges);
  //cleanLeafNodes();
}

void BaseDatapath::treeHeightReduction_float()
{
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_tree-height-reduction-float] Perform tree height reduction for floating point operations on DDDG" << std::endl);

	EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

	std::vector<bool> updated(numTotalNodes, 0);
	std::vector<int> bound_region(numTotalNodes, 0);

	/*
	int region_id = 0;
	unsigned node_id = 0;
	auto bound_it = loopBound.begin();
	while (node_id < *bound_it)
	{
	bound_region.at(node_id) = region_id;
	node_id++;
	if (node_id == *bound_it)
	{
	region_id++;
	bound_it++;
	if (bound_it == loopBound.end())
	break;
	}
	}
	*/

	std::set<Edge> to_remove_edges;
	std::vector<newEdge> to_add_edges;

	//nodes with no outgoing edges to first (bottom nodes first)
	for (int node_id = numTotalNodes - 1; node_id >= 0; node_id--)
	{
		if (nameToVertex.find(node_id) == nameToVertex.end()
			|| boost::degree(nameToVertex[node_id], graph_) == 0
			|| updated.at(node_id)
			|| !is_fassociative(microop.at(node_id)))
			continue;
		updated.at(node_id) = 1;
		int node_region = bound_region.at(node_id);

		std::list<unsigned> nodes;
		std::vector<Edge> tmp_remove_edges;
		std::vector<pair<int, bool> > leaves;
		std::vector<int> associative_chain;
		associative_chain.push_back(node_id);
		int chain_id = 0;
		while (chain_id < associative_chain.size())
		{
			int chain_node_id = associative_chain.at(chain_id);
			int chain_node_microop = microop.at(chain_node_id);
			if (is_fassociative(chain_node_microop))
			{
				updated.at(chain_node_id) = 1;
				int num_of_chain_parents = 0;
				in_edge_iter in_edge_it, in_edge_end;
				for (tie(in_edge_it, in_edge_end) = in_edges(nameToVertex[chain_node_id], graph_); in_edge_it != in_edge_end; ++in_edge_it)
				{
					int parent_id = vertexToName[source(*in_edge_it, graph_)];
					if (is_branch_op(microop.at(parent_id)))
						continue;
					num_of_chain_parents++;
				}
				if (num_of_chain_parents == 2)
				{
					nodes.push_front(chain_node_id);
					for (tie(in_edge_it, in_edge_end) = in_edges(nameToVertex[chain_node_id], graph_); in_edge_it != in_edge_end; ++in_edge_it)
					{
						Vertex parent_node = source(*in_edge_it, graph_);
						int parent_id = vertexToName[parent_node];
						assert(parent_id < chain_node_id);
						int parent_region = bound_region.at(parent_id);
						int parent_microop = microop.at(parent_id);
						if (is_branch_op(parent_microop))
							continue;
						Edge curr_edge = *in_edge_it;
						tmp_remove_edges.push_back(curr_edge);

						if (parent_region == node_region)
						{
							updated.at(parent_id) = 1;
							if (!is_fassociative(parent_microop))
								leaves.push_back(make_pair(parent_id, 0));
							else
							{
								out_edge_iter out_edge_it, out_edge_end;
								int num_of_children = 0;
								for (tie(out_edge_it, out_edge_end) = out_edges(parent_node, graph_); out_edge_it != out_edge_end; ++out_edge_it) {
									if (edge_to_parid[*out_edge_it] != CONTROL_EDGE)
										num_of_children++;
								}

								if (num_of_children == 1)
									associative_chain.push_back(parent_id);
								else
									leaves.push_back(make_pair(parent_id, 0));
							}
						}
						else
							leaves.push_back(make_pair(parent_id, 1));
					}
				}
				else
					leaves.push_back(make_pair(chain_node_id, 0));
			}
			else
				leaves.push_back(make_pair(chain_node_id, 0));
			chain_id++;
		}
		//build the tree
		if (nodes.size() < 3)
			continue;

		for (auto it = tmp_remove_edges.begin(), E = tmp_remove_edges.end(); it != E; it++)
			to_remove_edges.insert(*it);

		std::map<unsigned, unsigned> rank_map;
		auto leaf_it = leaves.begin();

		while (leaf_it != leaves.end())
		{
			if (leaf_it->second == 0)
				rank_map[leaf_it->first] = 0;
			else
				rank_map[leaf_it->first] = numTotalNodes;
			++leaf_it;
		}
		//reconstruct the rest of the balanced tree
		auto node_it = nodes.begin();

		while (node_it != nodes.end())
		{
			unsigned node1, node2;
			if (rank_map.size() == 2)
			{
				node1 = rank_map.begin()->first;
				node2 = (++rank_map.begin())->first;
			}
			else
				findMinRankNodes(node1, node2, rank_map);
			assert((node1 != numTotalNodes) && (node2 != numTotalNodes));
			to_add_edges.push_back({ node1, *node_it, 1 });
			to_add_edges.push_back({ node2, *node_it, 1 });

			//place the new node in the map, remove the two old nodes
			rank_map[*node_it] = max(rank_map[node1], rank_map[node2]) + 1;
			rank_map.erase(node1);
			rank_map.erase(node2);
			++node_it;
		}
	}
	updateGraphWithIsolatedEdges(to_remove_edges);
	updateGraphWithNewEdges(to_add_edges);
	//cleanLeafNodes();
}

void BaseDatapath::findMinRankNodes(unsigned &node1, unsigned &node2, std::map<unsigned, unsigned> &rank_map)
{
  unsigned min_rank = numTotalNodes;
  for (auto it = rank_map.begin(); it != rank_map.end(); ++it)
  {
    unsigned node_rank = it->second;
    if (node_rank < min_rank)
    {
      node1 = it->first;
      min_rank = node_rank;
    }
  }
  min_rank = numTotalNodes;
  for (auto it = rank_map.begin(); it != rank_map.end(); ++it)
  {
    unsigned node_rank = it->second;
    if ((it->first != node1) && (node_rank < min_rank))
    {
      node2 = it->first;
      min_rank = node_rank;
    }
  }
}

void BaseDatapath::updateGraphWithNewEdges(std::vector<newEdge> &to_add_edges)
{
  for(auto it = to_add_edges.begin(); it != to_add_edges.end(); ++it)
  {
    if (it->from != it->to && !doesEdgeExist(it->from, it->to))
      get(boost::edge_weight, graph_)[add_edge(it->from, it->to, graph_).first] = it->parid;
  }
}

void BaseDatapath::updateGraphWithIsolatedNodes(std::vector<unsigned> &to_remove_nodes)
{
  for(auto it = to_remove_nodes.begin(); it != to_remove_nodes.end(); ++it)
    clear_vertex(nameToVertex[*it], graph_);
}

void BaseDatapath::updateGraphWithIsolatedEdges(std::set<Edge> &to_remove_edges)
{
  for (auto it = to_remove_edges.begin(), E = to_remove_edges.end(); it!=E; ++it)
    remove_edge(*it, graph_);
}

#ifdef ALADDIN_H
/*
 * Write per cycle activity to bench_stats. The format is:
 * cycle_num,num-of-muls,num-of-adds,num-of-bitwise-ops,num-of-reg-reads,num-of-reg-writes
 * If it is called from ScratchpadDatapath, it also outputs per cycle memory
 * activity for each partitioned array.
 */
void BaseDatapath::writePerCycleActivity()
{
  std::string bn(benchName);

  activity_map mul_activity, add_activity, bit_activity;
  activity_map ld_activity, st_activity;

  max_activity_map max_mul_per_function;
  max_activity_map max_add_per_function;
  max_activity_map max_bit_per_function;

  std::vector<std::string> comp_partition_names;
  std::vector<std::string> mem_partition_names;
  //registers.getRegisterNames(comp_partition_names);
  getMemoryBlocks(mem_partition_names);


  initPerCycleActivity(comp_partition_names, mem_partition_names,
                       ld_activity, st_activity,
                       mul_activity, add_activity, bit_activity,
                       max_mul_per_function, max_add_per_function,
                       max_bit_per_function,
                       num_cycles);

  updatePerCycleActivity(ld_activity, st_activity,
                         mul_activity, add_activity, bit_activity,
                         max_mul_per_function, max_add_per_function,
                         max_bit_per_function);

  outputPerCycleActivity(comp_partition_names, mem_partition_names,
                         ld_activity, st_activity,
                         mul_activity, add_activity, bit_activity,
                         max_mul_per_function, max_add_per_function,
                         max_bit_per_function);

}

void BaseDatapath::initPerCycleActivity(
     std::vector<std::string> &comp_partition_names,
     std::vector<std::string> &mem_partition_names,
     activity_map &ld_activity, activity_map &st_activity,
     activity_map &mul_activity, activity_map &add_activity,
     activity_map &bit_activity,
     max_activity_map &max_mul_per_function,
     max_activity_map &max_add_per_function,
     max_activity_map &max_bit_per_function,
     int num_cycles)
{
  for (auto it = comp_partition_names.begin(); it != comp_partition_names.end() ; ++it)
  {
    ld_activity.insert({*it, make_vector(num_cycles)});
    st_activity.insert({*it, make_vector(num_cycles)});
  }
  for (auto it = mem_partition_names.begin(); it != mem_partition_names.end() ; ++it)
  {
    ld_activity.insert({*it, make_vector(num_cycles)});
    st_activity.insert({*it, make_vector(num_cycles)});
  }
  for (auto it = functionNames.begin(); it != functionNames.end() ; ++it)
  {
    mul_activity.insert({*it, make_vector(num_cycles)});
    add_activity.insert({*it, make_vector(num_cycles)});
    bit_activity.insert({*it, make_vector(num_cycles)});
    max_mul_per_function.insert({*it, 0});
    max_add_per_function.insert({*it, 0});
    max_bit_per_function.insert({*it, 0});
  }
}

void BaseDatapath::updatePerCycleActivity(
     activity_map &ld_activity, activity_map &st_activity,
     activity_map &mul_activity, activity_map &add_activity,
     activity_map &bit_activity,
     max_activity_map &max_mul_per_function,
     max_activity_map &max_add_per_function,
     max_activity_map &max_bit_per_function)
{
  /*We use two ways to count the number of functional units in accelerators: one
   * assumes that functional units can be reused in the same region; the other
   * assumes no reuse of functional units. The advantage of reusing is that it
   * elimates the cost of duplicating functional units which can lead to high
   * leakage power and area. However, additional wires and muxes may need to be
   * added for reusing.
   * In the current model, we assume that multipliers can be reused, since the
   * leakage power and area of multipliers are relatively significant, and no
   * reuse for adders. This way of modeling is consistent with our observation
   * of accelerators generated with Vivado.*/
  std::vector<std::string> dynamic_methodid(numTotalNodes, "");
  initDynamicMethodID(dynamic_methodid);

  int num_adds_so_far = 0, num_bits_so_far = 0;
  auto bound_it = loopBound.begin();
  for(unsigned node_id = 0; node_id < numTotalNodes; ++node_id)
  {
    char func_id[256];
    int count;

    sscanf(dynamic_methodid.at(node_id).c_str(), "%[^-]-%d\n", func_id, &count);
    if (node_id == *bound_it) {
      if (max_add_per_function[func_id] < num_adds_so_far)
        max_add_per_function[func_id] = num_adds_so_far;
      if (max_bit_per_function[func_id] < num_bits_so_far)
        max_bit_per_function[func_id] = num_bits_so_far;
      num_adds_so_far = 0;
      num_bits_so_far = 0;
      bound_it++;
    }
    if (finalIsolated.at(node_id))
      continue;
    int node_level = newLevel.at(node_id);
    int node_microop = microop.at(node_id);

    if (is_mul_op(node_microop))
      mul_activity[func_id].at(node_level) +=1;
    else if (is_add_op(node_microop)) {
      add_activity[func_id].at(node_level) +=1;
      num_adds_so_far +=1;
    }
    else if (is_bit_op(node_microop)) {
      bit_activity[func_id].at(node_level) +=1;
      num_bits_so_far +=1;
    }
    else if (is_load_op(node_microop)) {
      std::string base_addr = baseAddress[node_id].first;
      if (ld_activity.find(base_addr) != ld_activity.end())
        ld_activity[base_addr].at(node_level) += 1;
    }
    else if (is_store_op(node_microop)) {
      std::string base_addr = baseAddress[node_id].first;
      if (st_activity.find(base_addr) != st_activity.end())
        st_activity[base_addr].at(node_level) += 1;
    }
  }
  for (auto it = functionNames.begin(); it != functionNames.end() ; ++it)
    max_mul_per_function[*it] = *(std::max_element(mul_activity[*it].begin(),
                                                 mul_activity[*it].end()));
}


void BaseDatapath::outputPerCycleActivity(
     std::vector<std::string> &comp_partition_names,
     std::vector<std::string> &mem_partition_names,
     activity_map &ld_activity, activity_map &st_activity,
     activity_map &mul_activity, activity_map &add_activity,
     activity_map &bit_activity,
     max_activity_map &max_mul_per_function,
     max_activity_map &max_add_per_function,
     max_activity_map &max_bit_per_function)
{
  ofstream stats, power_stats;
  std::string bn(benchName);
  std::string file_name = bn + "_stats";
  stats.open(file_name.c_str());
  file_name += "_power";
  power_stats.open(file_name.c_str());

  stats << "cycles," << num_cycles << "," << numTotalNodes << std::endl;
  power_stats << "cycles," << num_cycles << "," << numTotalNodes << std::endl;
  stats << num_cycles << "," ;
  power_stats << num_cycles << "," ;

  /*Start writing the second line*/
  for (auto it = functionNames.begin(); it != functionNames.end() ; ++it)
  {
    stats << *it << "-mul," << *it << "-add," << *it << "-bit,";
    power_stats << *it << "-mul," << *it << "-add," << *it << "-bit,";
  }
  stats << "reg,";
  power_stats << "reg,";
  for (auto it = mem_partition_names.begin();
       it != mem_partition_names.end() ; ++it) {
    stats << *it << "-read," << *it << "-write,";
  }
  stats << std::endl;
  power_stats << std::endl;
  /*Finish writing the second line*/

  /*Caculating the number of FUs and leakage power*/
  int max_reg_read =  0, max_reg_write = 0;
  for (unsigned level_id = 0; ((int) level_id) < num_cycles; ++level_id)
  {
    if (max_reg_read < regStats.at(level_id).reads )
      max_reg_read = regStats.at(level_id).reads ;
    if (max_reg_write < regStats.at(level_id).writes )
      max_reg_write = regStats.at(level_id).writes ;
  }
  int max_reg = max_reg_read + max_reg_write;
  int max_add = 0, max_bit = 0, max_mul = 0;
  for (auto it = functionNames.begin(); it != functionNames.end() ; ++it)
  {
    max_bit += max_bit_per_function[*it];
    max_add += max_add_per_function[*it];
    max_mul += max_mul_per_function[*it];
  }

  float add_leakage_power = ADD_leak_power * max_add;
  float mul_leakage_power = MUL_leak_power * max_mul;
  float reg_leakage_power = registers.getTotalLeakagePower()
                            + REG_leak_power * 32 * max_reg;
  float fu_leakage_power = mul_leakage_power
                           + add_leakage_power
                           + reg_leakage_power;
  /*Finish caculating the number of FUs and leakage power*/

  float fu_dynamic_energy = 0;

  /*Start writing per cycle activity */
  for (unsigned curr_level = 0; ((int)curr_level) < num_cycles ; ++curr_level)
  {
    stats << curr_level << "," ;
    power_stats << curr_level << ",";
    //For FUs
    for (auto it = functionNames.begin(); it != functionNames.end() ; ++it)
    {
      float curr_mul_dynamic_power = (MUL_switch_power + MUL_int_power)
                                       * mul_activity[*it].at(curr_level);
      float curr_add_dynamic_power = (ADD_switch_power + ADD_int_power)
                                       * add_activity[*it].at(curr_level);
      fu_dynamic_energy += ( curr_mul_dynamic_power + curr_add_dynamic_power )
                            * cycleTime;

      stats       << mul_activity[*it].at(curr_level) << ","
                  << add_activity[*it].at(curr_level) << ","
                  << bit_activity[*it].at(curr_level) << ",";
      power_stats << curr_mul_dynamic_power + mul_leakage_power << ","
                  << curr_add_dynamic_power + add_leakage_power << ","
                  << "0," ;
    }
    //For regs
    int curr_reg_reads = regStats.at(curr_level).reads;
    int curr_reg_writes = regStats.at(curr_level).writes;
    float curr_reg_dynamic_energy = (REG_int_power + REG_sw_power)
                                  *(curr_reg_reads + curr_reg_writes)
                                  * 32 * cycleTime;
    for (auto it = comp_partition_names.begin();
         it != comp_partition_names.end() ; ++it)
    {
      curr_reg_reads += ld_activity.at(*it).at(curr_level);
      curr_reg_writes += st_activity.at(*it).at(curr_level);
      curr_reg_dynamic_energy += registers.getReadEnergy(*it)
                                   * ld_activity.at(*it).at(curr_level)
                                 + registers.getWriteEnergy(*it)
                                   * st_activity.at(*it).at(curr_level);
    }
    fu_dynamic_energy += curr_reg_dynamic_energy;

    stats << curr_reg_reads << "," << curr_reg_writes << "," ;
    power_stats << curr_reg_dynamic_energy / cycleTime + reg_leakage_power;

    for (auto it = mem_partition_names.begin();
         it != mem_partition_names.end() ; ++it)
      stats << ld_activity.at(*it).at(curr_level) << ","
            << st_activity.at(*it).at(curr_level) << ",";
    stats << std::endl;
    power_stats << std::endl;
  }
  stats.close();
  power_stats.close();

  float avg_mem_power =0, avg_mem_dynamic_power = 0, mem_leakage_power = 0;

  getAverageMemPower(num_cycles, &avg_mem_power,
                     &avg_mem_dynamic_power, &mem_leakage_power);

  float avg_fu_dynamic_power = fu_dynamic_energy / (cycleTime * num_cycles);
  float avg_fu_power = avg_fu_dynamic_power + fu_leakage_power;
  float avg_power = avg_fu_power + avg_mem_power;

  float mem_area = getTotalMemArea();
  unsigned mem_size = getTotalMemSize();
  float fu_area = registers.getTotalArea()
                  + ADD_area * max_add
                  + MUL_area * max_mul
                  + REG_area * 32 * max_reg;
  float total_area = mem_area + fu_area;
  //Summary output:
  //Cycle, Avg Power, Avg FU Power, Avg MEM Power, Total Area, FU Area, MEM Area
  std::cerr << "===============================" << std::endl;
  std::cerr << "        Aladdin Results        " << std::endl;
  std::cerr << "===============================" << std::endl;
  std::cerr << "Running : " << benchName << std::endl;
  std::cerr << "Cycle : " << num_cycles << " cycles" << std::endl;
  std::cerr << "Avg Power: " << avg_power << " mW" << std::endl;
  std::cerr << "Avg FU Power: " << avg_fu_power << " mW" << std::endl;
  std::cerr << "Avg FU Dynamic Power: " << avg_fu_dynamic_power << " mW" << std::endl;
  std::cerr << "Avg FU leakage Power: " << fu_leakage_power << " mW" << std::endl;
  std::cerr << "Avg SRAM Power: " << avg_mem_power << " mW" << std::endl;
  std::cerr << "Avg SRAM Dynamic Power: " << avg_mem_dynamic_power << " mW" << std::endl;
  std::cerr << "Avg SRAM Leakage Power: " << mem_leakage_power << " mW" << std::endl;
  std::cerr << "Total Area: " << total_area << " uM^2" << std::endl;
  std::cerr << "FU Area: " << fu_area << " uM^2" << std::endl;
  std::cerr << "SRAM Area: " << mem_area << " uM^2" << std::endl;
  std::cerr << "SRAM size: " << mem_size / 1024 << " KB" << std::endl;
  std::cerr << "Num of Multipliers (32-bit): " << max_mul  << std::endl;
  std::cerr << "Num of Adders (32-bit): " << max_add << std::endl;
  std::cerr << "===============================" << std::endl;
  std::cerr << "        Aladdin Results        " << std::endl;
  std::cerr << "===============================" << std::endl;

  ofstream summary;
  file_name = bn + "_summary";
  summary.open(file_name.c_str());
  summary << "===============================" << std::endl;
  summary << "        Aladdin Results        " << std::endl;
  summary << "===============================" << std::endl;
  summary << "Running : " << benchName << std::endl;
  summary << "Cycle : " << num_cycles << " cycles" << std::endl;
  summary << "Avg Power: " << avg_power << " mW" << std::endl;
  summary << "Avg FU Power: " << avg_fu_power << " mW" << std::endl;
  summary << "Avg FU Dynamic Power: " << avg_fu_dynamic_power << " mW" << std::endl;
  summary << "Avg FU leakage Power: " << fu_leakage_power << " mW" << std::endl;
  summary << "Avg SRAM Power: " << avg_mem_power << " mW" << std::endl;
  summary << "Avg SRAM Dynamic Power: " << avg_mem_dynamic_power << " mW" << std::endl;
  summary << "Avg SRAM Leakage Power: " << mem_leakage_power << " mW" << std::endl;
  summary << "Total Area: " << total_area << " uM^2" << std::endl;
  summary << "FU Area: " << fu_area << " uM^2" << std::endl;
  summary << "SRAM Area: " << mem_area << " uM^2" << std::endl;
  summary << "SRAM size: " << mem_size << " B" << std::endl;
  summary << "Num of Multipliers (32-bit): " << max_mul  << std::endl;
  summary << "Num of Adders (32-bit): " << max_add << std::endl;
  summary << "===============================" << std::endl;
  summary << "        Aladdin Results        " << std::endl;
  summary << "===============================" << std::endl;
  summary.close();
}
#endif // End of ALADDIN_H

void BaseDatapath::writeGlobalIsolated()
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_isolated.gz";
  write_gzip_bool_file(file_name, finalIsolated.size(), finalIsolated);
}

/*
* Read: graph, getElementPtr.gz, completePartitionConfig, PartitionConfig
* Modify: baseAddress
*/
void BaseDatapath::initBaseAddress()
{
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [trace-analysis_base-address-initialization] Initializing Base Address\n");

	std::unordered_map<std::string, unsigned> comp_part_config;
	readCompletePartitionConfig(comp_part_config);
	std::unordered_map<std::string, partitionEntry> part_config;
	readPartitionConfig(part_config);

	std::unordered_map<unsigned, pair<std::string, long long int> > getElementPtr;
	initGetElementPtr(getElementPtr);

	edgeToParid = get(boost::edge_weight, graph_);

	vertex_iter vi, vi_end;
	for (tie(vi, vi_end) = vertices(graph_); vi != vi_end; ++vi)
	{
		if (boost::degree(*vi, graph_) == 0)
			continue;
		Vertex curr_node = *vi;
		unsigned node_id = vertexToName[curr_node];
		int node_microop = microop.at(node_id);
		if (!is_memory_op(node_microop))
			continue;
		bool modified = 0;
		//iterate its parents, until it finds the root parent
		while (true) {
			bool found_parent = 0;
			in_edge_iter in_edge_it, in_edge_end;

			for (tie(in_edge_it, in_edge_end) = in_edges(curr_node, graph_);
				in_edge_it != in_edge_end; ++in_edge_it) {
				int edge_parid = edgeToParid[*in_edge_it];
				if ((node_microop == LLVM_IR_Load && edge_parid != 1)
					|| (node_microop == LLVM_IR_GetElementPtr && edge_parid != 1)
					|| (node_microop == LLVM_IR_Store && edge_parid != 2))
					continue;

				unsigned parent_id = vertexToName[source(*in_edge_it, graph_)];
				int parent_microop = microop.at(parent_id);
				if (parent_microop == LLVM_IR_GetElementPtr
					|| parent_microop == LLVM_IR_Load) {
					//remove address calculation directly
					baseAddress[node_id] = getElementPtr[parent_id];
					curr_node = source(*in_edge_it, graph_);
					node_microop = parent_microop;
					found_parent = 1;
					modified = 1;
					break;
				}
				else if (parent_microop == LLVM_IR_Alloca) {
					std::string part_name = getElementPtr[parent_id].first;
					baseAddress[node_id] = getElementPtr[parent_id];
					modified = 1;
					break;
				}
			}
			if (!found_parent)
				break;
		}
		if (!modified)
			baseAddress[node_id] = getElementPtr[node_id];

		std::string part_name = baseAddress[node_id].first;
		if (partitionArrayName.size() !=0 ) {
			///FIXME: Bug here, need to check users want to partition which arrays, if they only partition part of input/out arrays.
			///Status: Fixed.
			std::set<std::string>::iterator it_ptName;
			it_ptName = partitionArrayName.find(part_name);
			if (it_ptName != partitionArrayName.end()) {
				if ((part_config.size() != 0 || comp_part_config.size() != 0) && part_config.find(part_name) == part_config.end() &&
					comp_part_config.find(part_name) == comp_part_config.end()) {
					std::cerr << "Unknown partition : " << part_name << "@inst: "
						<< node_id << std::endl;
					exit(-1);
				}
			}
			else {
				noPartitionArrayName.insert(part_name);
			}
		}
		else {
			noPartitionArrayName.insert(part_name);
		}
	}
	writeBaseAddress();
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [trace-analysis_base-address-initialization] Finished initializing Base Address\n");
}

void BaseDatapath::writeBaseAddress()
{
  //ostringstream file_name;
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  //file_name << benchName << "_baseAddr.gz";
	file_name += "_baseAddr.gz";
  gzFile gzip_file;
  gzip_file = gzopen(file_name.c_str(), "w");
  for (auto it = baseAddress.begin(), E = baseAddress.end(); it != E; ++it)
    gzprintf(gzip_file, "node:%u,part:%s,base:%lld\n", it->first, it->second.first.c_str(), it->second.second);
  gzclose(gzip_file);
}

bool BaseDatapath::completePartition() {
	std::unordered_map<std::string, unsigned> comp_part_config;
	if (!readCompletePartitionConfig(comp_part_config))
		return 0;

	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_convert-memory2register] Convert memory into registers based on partitioning configuration file\n" << std::endl);

	for (auto it = comp_part_config.begin(); it != comp_part_config.end(); ++it)
	{
		std::string base_addr = it->first;
		unsigned size = it->second;

		registers.createRegister(base_addr, size);
	}

	return 1;
}

void BaseDatapath::scratchpadPartition() {
	//read the partition config file to get the address range
	// <base addr, <type, part_factor> >
	
	/// Create scratchpad for arrays without partition configuration,
	/// but size and wordsize of each array are set to "0".
	std::set<std::string>::iterator it_noParArray, ie_noParArray;
	it_noParArray = noPartitionArrayName.begin();
	ie_noParArray = noPartitionArrayName.end();
	for (; it_noParArray != ie_noParArray; it_noParArray++) {
		std::string noParArray_name = *it_noParArray;
		setScratchpad(noParArray_name, 0, 0);
	}

	std::unordered_map<std::string, partitionEntry> part_config;
	if (!readPartitionConfig(part_config))
		return;

	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [optimization_array-partition] Partition arrays based on partitioning configuration file\n" << std::endl);
	std::string bn(benchName);

	std::unordered_map<unsigned, pair<long long int, unsigned> > address;
	initAddressAndSize(address);
	//set scratchpad
	for (auto it = part_config.begin(); it != part_config.end(); ++it)
	{
		std::string base_addr = it->first;
		unsigned size = it->second.array_size; //num of bytes
		unsigned p_factor = it->second.part_factor;
		if (p_factor == 1) {
			continue;
		}
		unsigned wordsize = it->second.wordsize; //in bytes
		unsigned per_size = ceil(((float)size) / p_factor);

		for (unsigned i = 0; i < p_factor; i++)
		{
			ostringstream oss;
			oss << base_addr << "-" << i;
			setScratchpad(oss.str(), per_size, wordsize);
		}
	}
	///FIXME: It seems that base address of a partition array is not updated
	for (unsigned node_id = 0; node_id < numTotalNodes; node_id++)
	{
		if (!is_memory_op(microop.at(node_id)))
			continue;

		if (baseAddress.find(node_id) == baseAddress.end())
			continue;
		std::string base_label = baseAddress[node_id].first;
		long long int base_addr = baseAddress[node_id].second;

		auto part_it = part_config.find(base_label);
		if (part_it != part_config.end())
		{
			std::string p_type = part_it->second.type;
			assert((!p_type.compare("block")) || (!p_type.compare("cyclic")));

			unsigned num_of_elements = part_it->second.array_size;
			unsigned p_factor = part_it->second.part_factor;
			if (p_factor == 1) {
				continue;
			}
			long long int abs_addr = address[node_id].first;
			unsigned data_size = address[node_id].second / 8; //in bytes
			unsigned rel_addr = (abs_addr - base_addr) / data_size;
			if (!p_type.compare("block"))  //block partition
			{
				ostringstream oss;
				unsigned num_of_elements_in_2 = next_power_of_two(num_of_elements);
				oss << base_label << "-"
					<< (int)(rel_addr / ceil(num_of_elements_in_2 / p_factor));
				baseAddress[node_id].first = oss.str();
			}
			else // cyclic partition
			{
				ostringstream oss;
				oss << base_label << "-" << (rel_addr) % p_factor;
				baseAddress[node_id].first = oss.str();
			}
		}
	}
	writeBaseAddress();
}

void BaseDatapath::setScratchpad(std::string baseName, unsigned num_of_bytes, unsigned word_size){
	assert(!partitionExist(baseName));
	//size: number of words
	unsigned new_id = baseToPartitionID.size();
	baseToPartitionID[baseName] = new_id;

	OccupiedMemBWPerPartition bwPerPartition = { 0, 0 };
	occupiedMemPerPartition.push_back(bwPerPartition);
	//occupiedBWPerPartition.push_back(0);
	sizePerPartition.push_back(num_of_bytes);
}

bool BaseDatapath::partitionExist(std::string baseName) {
	std::unordered_map<std::string, unsigned>::iterator part_it;
	part_it = baseToPartitionID.find(baseName);
	return (part_it != baseToPartitionID.end()) ? 1 : 0;
}

void BaseDatapath::writeFinalLevel()
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_level.gz";
  write_gzip_file(file_name, newLevel.size(), newLevel);
}

void BaseDatapath::writeMicroop(std::vector<int> &microop)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_microop.gz";
  write_gzip_file(file_name, microop.size(), microop);
}

void BaseDatapath::initPrevBasicBlock(std::vector<std::string> &prevBasicBlock)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_prevBasicBlock.gz";
  read_gzip_string_file(file_name, prevBasicBlock.size(), prevBasicBlock);
}

void BaseDatapath::initCurBasicBlock(std::vector<std::string> &curBasicBlock)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
	file_name += "_curBasicBlock.gz";
	read_gzip_string_file(file_name, curBasicBlock.size(), curBasicBlock);
}

void BaseDatapath::initDynamicMethodID(std::vector<std::string> &methodid)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_dynamic_funcid.gz";
  read_gzip_string_file(file_name, methodid.size(), methodid);
}

void BaseDatapath::initMethodID(std::vector<int> &methodid)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_methodid.gz";
  read_gzip_file(file_name, methodid.size(), methodid);
}

void BaseDatapath::initInstID(std::vector<std::string> &instid)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_instid.gz";
  read_gzip_string_file(file_name, instid.size(), instid);
}

void BaseDatapath::initAddress(std::unordered_map<unsigned, long long int> &address)
{
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_memaddr.gz";
  gzFile gzip_file;
  gzip_file = gzopen(file_name.c_str(), "r");
  while (!gzeof(gzip_file))
  {
    char buffer[256];
    if (gzgets(gzip_file, buffer, 256) == NULL)
      break;
    unsigned node_id, size;
    long long int addr;
    sscanf(buffer, "%d,%lld,%d\n", &node_id, &addr, &size);
    address[node_id] = addr;
  }
  gzclose(gzip_file);
}

void BaseDatapath::initAddressAndSize(std::unordered_map<unsigned, pair<long long int, unsigned> > &address)
{
  //std::string file_name(inputPath+benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_memaddr.gz";
  gzFile gzip_file;
  gzip_file = gzopen(file_name.c_str(), "r");
  while (!gzeof(gzip_file))
  {
    char buffer[256];
    if (gzgets(gzip_file, buffer, 256) == NULL)
      break;
    unsigned node_id, size;
    long long int addr;
    sscanf(buffer, "%d,%lld,%d\n", &node_id, &addr, &size);
    address[node_id] = make_pair(addr, size);
  }
  gzclose(gzip_file);
}

void BaseDatapath::initLineNum(std::vector<int> &line_num)
{
  //ostringstream file_name;
	//std::string file_name(inputPath+benchName);
	std::string file_name(outputPath + benchName);
	file_name += "_linenum.gz";
  //file_name << benchName << "_linenum.gz";
  read_gzip_file(file_name, line_num.size(), line_num);
}

void BaseDatapath::initGetElementPtr(std::unordered_map<unsigned, pair<std::string, long long int> > &get_element_ptr)
{
  //ostringstream file_name;
  //file_name << benchName << "_getElementPtr.gz";
	//std::string file_name(inputPath+benchName);
	std::string file_name(outputPath + benchName);
	file_name += "_getElementPtr.gz";
  gzFile gzip_file;
  gzip_file = gzopen(file_name.c_str(), "r");
  while (!gzeof(gzip_file))
  {
    char buffer[256];
    if (gzgets(gzip_file, buffer, 256) == NULL)
      break;
    unsigned node_id;
    long long int address;
    char label[256];
    sscanf(buffer, "%d,%[^,],%lld\n", &node_id, label, &address);
    get_element_ptr[node_id] = make_pair(label, address);
  }
  gzclose(gzip_file);

	/// change intermedia array name to original array name
	std::unordered_map<unsigned, pair<std::string, long long int> >::iterator it = get_element_ptr.begin();
	std::unordered_map<unsigned, pair<std::string, long long int> >::iterator ie = get_element_ptr.end();
	for (; it != ie; ++it) {
		std::string fake_arrayName = it->second.first;
		getElementPtrName2arrayNameMapTy::iterator it_found = getElementPtrName2arrayNameMap.find(fake_arrayName);
		if (it_found != getElementPtrName2arrayNameMap.end()) {
			it->second.first = it_found->second;
		}
	}
}


// For parallel simulation
void BaseDatapath::setGraphForStepping()
{
	std::cerr << "=============================================" << std::endl;
	std::cerr << "      Scheduling...            " << benchName << std::endl;
	std::cerr << "=============================================" << std::endl;

	newLevel.assign(numTotalNodes, 0);
	maxParentOpLatency.assign(numTotalNodes,0);

	edgeToParid = get(boost::edge_weight, graph_);

	numTotalEdges = boost::num_edges(graph_);
	numParents.assign(numTotalNodes, 0);
	latestParents.assign(numTotalNodes, 0);
	executedNodes = 0;
	totalConnectedNodes = 0;
	for (unsigned node_id = 0; node_id < numTotalNodes; node_id++)
	{
		Vertex node = nameToVertex[node_id];
		if (boost::degree(node, graph_) != 0 || is_dma_op(microop.at(node_id)))
		{
			finalIsolated.at(node_id) = 0;
			numParents.at(node_id) = boost::in_degree(node, graph_);
			totalConnectedNodes++;
		}
	}
	executingQueue.clear();
	readyToExecuteQueue.clear();
	initExecutingQueue();
}

void BaseDatapath::dumpGraph()
{
	std::string bn(benchName + "_graph.dot");
  std::ofstream out(bn);

	EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);

	NameVecTy cur_basic_block(numTotalNodes, "");
	initCurBasicBlock(cur_basic_block);

	NameVecTy funcNamesVec(numTotalNodes, "");
	initDynamicMethodID(funcNamesVec);
	NameVecTy::iterator it, ie;
	for (it = funcNamesVec.begin(), ie = funcNamesVec.end(); it != ie; ++it) {
		std::size_t funName_len = it->find("-");
		std::string funName = it->substr(0, funName_len);
		it->swap(funName);
	}

  //write_graphviz(out, graph_);
	write_graphviz(out, graph_, color_writer(graph_, vertexToName, cur_basic_block, funcNamesVec, microop, bbFuncNamePair2lpNameLevelPairMap), edge_color_writer(graph_, edge_to_parid));

	// Display the graph
	GraphProgram::Name Program = GraphProgram::DOT;
	DisplayGraph(bn, true, Program);
	out.close();
}

int BaseDatapath::clearGraph()
{
  std::vector< Vertex > topo_nodes;
  boost::topological_sort(graph_, std::back_inserter(topo_nodes));
  //bottom nodes first
  std::vector<int> earliest_child(numTotalNodes, num_cycles);
  for (auto vi = topo_nodes.begin(); vi != topo_nodes.end(); ++vi)
  {
    unsigned node_id = vertexToName[*vi];
    if (finalIsolated.at(node_id))
      continue;
    unsigned node_microop = microop.at(node_id);
    if (!is_memory_op(node_microop) && ! is_branch_op(node_microop))
      if ((earliest_child.at(node_id) - 1 ) > newLevel.at(node_id))
        newLevel.at(node_id) = earliest_child.at(node_id) - 1;

    in_edge_iter in_i, in_end;
    for (tie(in_i, in_end) = in_edges(*vi , graph_); in_i != in_end; ++in_i)
    {
      int parent_id = vertexToName[source(*in_i, graph_)];
      if (earliest_child.at(parent_id) > newLevel.at(node_id))
        earliest_child.at(parent_id) = newLevel.at(node_id);
    }
  }
  updateRegStats();
  return num_cycles;
}

void BaseDatapath::updateRegStats()
{
  regStats.assign(num_cycles, {0, 0, 0});
  for(unsigned node_id = 0; node_id < numTotalNodes; node_id++)
  {
    if (finalIsolated.at(node_id))
      continue;
    if (is_control_op(microop.at(node_id)) ||
        is_index_op(microop.at(node_id)))
      continue;
    int node_level = newLevel.at(node_id);
    int max_children_level 		= node_level;

    Vertex node = nameToVertex[node_id];
    out_edge_iter out_edge_it, out_edge_end;
    std::set<int> children_levels;
    for (tie(out_edge_it, out_edge_end) = out_edges(node, graph_); out_edge_it != out_edge_end; ++out_edge_it)
    {
      int child_id = vertexToName[target(*out_edge_it, graph_)];
      int child_microop = microop.at(child_id);
      if (is_control_op(child_microop))
        continue;

      if (is_load_op(child_microop))
        continue;

      int child_level = newLevel.at(child_id);
      if (child_level > max_children_level)
        max_children_level = child_level;
      if (child_level > node_level  && child_level != num_cycles - 1)
        children_levels.insert(child_level);

    }
    for (auto it = children_levels.begin(); it != children_levels.end(); it++)
        regStats.at(*it).reads++;

    if (max_children_level > node_level && node_level != 0 )
      regStats.at(node_level).writes++;
  }
}

void BaseDatapath::calculateInstructionDistribution() {
	// Statistical counters
	uint64_t num_of_storeInst = 0;
	uint64_t num_of_loadInst = 0;
	uint64_t num_of_faddInst = 0;
	uint64_t num_of_fsubInst = 0;
	uint64_t num_of_fmulInst = 0;
	uint64_t num_of_fdivInst = 0;
	uint64_t num_of_fcmpInst = 0;

	uint64_t num_of_bitInst = 0;
	uint64_t num_of_memInst = 0;
	uint64_t num_of_computeInst = 0;
	uint64_t num_of_branchInst = 0;
	uint64_t num_of_otherInst = 0;
	uint64_t num_of_indexOpInst = 0;
	uint64_t num_of_controlInst = 0; // Phi instruction
	uint64_t num_of_callInst = 0;

	subTraceInst = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	for (unsigned i = 0; i < numTotalNodes; i++) {
		unsigned opcode = microop.at(i);

		if (is_memory_op(opcode)) {
			if (is_store_op(opcode)) {
				num_of_storeInst++;
				subTraceInst.num_st_inst++;
			}
			else if (is_load_op(opcode)) {
				num_of_loadInst++;
				subTraceInst.num_ld_inst++;
			}
			else {
				assert(false && "DEBUG-INFO: [Instruction_Distribution] Neither store nor load instruction, ERROR!\n");
			}
			num_of_memInst++;
		}
		else if (is_compute_op(opcode)) {
			if (is_bit_op(opcode)) {
				num_of_bitInst++;
				subTraceInst.num_bitwise_inst++;
			}
			else if (is_fadd_op(opcode)) {
				num_of_faddInst++;
				subTraceInst.num_fadd_inst++;
			}
			else if (is_fsub_op(opcode)) {
				num_of_fsubInst++;
				subTraceInst.num_fsub_inst++;
			}
			else if (is_fmul_op(opcode)) {
				num_of_fmulInst++;
				subTraceInst.num_fmul_inst++;
			}
			else if (is_fdiv_op(opcode)) {
				num_of_fdivInst++;
				subTraceInst.num_fdiv_inst++;
			}
			else if (is_fcmp_op(opcode)) {
				num_of_fcmpInst++;
				subTraceInst.num_fcmp_inst++;
			}
			else {
				// Integer operations
				subTraceInst.num_integer_inst++;
			}
			num_of_computeInst++;
		}
		else if (is_branch_op(opcode)) {
			num_of_branchInst++;
			subTraceInst.num_br_inst++;
		}
		else {
			if (is_index_op(opcode)) {
				num_of_indexOpInst++;
				///FIXME: In the current implementation, we treat it as the control instruction.
				subTraceInst.num_control_inst++;
			}
			if (is_control_op(opcode)){
				// Phi instruction
				num_of_controlInst++;
				subTraceInst.num_control_inst++;
			}
			if (is_call_op(opcode)) {
				num_of_callInst++;
			}
			num_of_otherInst++;
		}
	}
	/*
	/// The following part of code is used to test whether PHI instructions
	/// are removed or not.
	std::list<Vertex> topo_nodes;
	/// 1. Topologically sort the DDDG graph
	//boost::topological_sort(graph_arg, std::back_inserter(topo_nodes));
	uint64_t num_of_phiInst = 0;
	uint64_t num_of_phiIndegree = 0;
	uint64_t num_of_phiOutdegree = 0;
	uint64_t nodes_in_graph = 0;
	uint64_t num_vertice = 0;
	nodes_in_graph = num_vertices(graph_);
	boost::topological_sort(graph_, std::back_inserter(topo_nodes));
	num_vertice = topo_nodes.size();
	if (nodes_in_graph != num_vertice) {
		cout << "nodes_in_graph is not equal to num_vertice" << endl;
	}
	std::list<Vertex>::iterator it, ie;
	for (it = topo_nodes.begin(), ie = topo_nodes.end(); it != ie; ++it) {
		Vertex tmp = *it;
		unsigned tmp_opcode = getvertexOpcode(tmp);
		if (tmp_opcode == LLVM_IR_PHI) {
			if (boost::degree(tmp, graph_) != 0) {
				num_of_phiInst++;
			}
		}
	}
	*/

	// Write instruction distribution into a file
	//std::string file_name = inputPath + benchName + "_instDistribution.csv";
	std::string file_name = outputPath + benchName + "_instDistribution.csv";

	ofstream instDistribution_file(file_name);
	if (instDistribution_file.is_open()) {
		//instDistribution_file << "kernel name,store,load,bit-wise,indexOp,control,call,memory,compute,branch,other" << std::endl;

		instDistribution_file << benchName << "," << num_of_storeInst << "," << num_of_loadInst << ","; 
		instDistribution_file << num_of_bitInst << "," << num_of_indexOpInst << "," << num_of_controlInst << ",";
		instDistribution_file << num_of_callInst << "," << num_of_memInst << "," << num_of_computeInst << ","; 
		instDistribution_file << num_of_branchInst << "," << num_of_otherInst << std::endl;

		instDistribution_file.close();
	}
	else {
		assert(false && "DEBUG-INFO: [Instruction_Distribution] Error: Could not open instDistribution!\n");
	}
}

/// calculateTimestampDDDG function is used to calculate timestamp for all nodes
/// in Graph graph_arg.
/// Algorithm used in this function is adopted from Justin Holewinski's PLDI2012
/// paper, "Dynamic Trace-based analysis of vectorization potential of applications".
//void BaseDatapath::calculateTimestampDDDG(Graph graph_arg) {
void BaseDatapath::calculateTimestampDDDG() {
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [analysis_timestamp-calculation] Calculate timestamp of nodes on DDDG\n" << std::endl);

	assert(staticInstID2OpcodeMap.empty() != true && "staticInstID2OpcodeMap is empty, ERROR, please check");

	/// Change operation opcodes of related to induction variables from add/sub to indexAdd/indexSub
	staticInstID2OpcodeMapTy::iterator it_opcode = staticInstID2OpcodeMap.begin();
	staticInstID2OpcodeMapTy::iterator ie_opcode = staticInstID2OpcodeMap.end();
	for (; it_opcode != ie_opcode; ++it_opcode) {
		std::string static_instID_str = it_opcode->first;
		unsigned opcode_inst = it_opcode->second;
		if (static_instID_str.find("indvars") == std::string::npos){
			continue;
		}

		switch (opcode_inst) {
		case LLVM_IR_Add:
			it_opcode->second = LLVM_IR_IndexAdd;
			break;
		case LLVM_IR_Sub:
			it_opcode->second = LLVM_IR_IndexSub;
			break;
		default:
			// Do nothing
			break;
		}
	}

	/// Find interesting instructions: Currently, we only focus on floating point arithmetic operations
	unsigned size_static_inst = staticInstID2OpcodeMap.size();
	staticInstID2OpcodeMapTy interestedStaticInstID2OpcodeMap;
	staticInstID2OpcodeMapTy::iterator it_map = staticInstID2OpcodeMap.begin();
	staticInstID2OpcodeMapTy::iterator ie_map = staticInstID2OpcodeMap.end();
	for (; it_map != ie_map; ++it_map) {
		std::string instID_str = it_map->first;
		unsigned opcode = it_map->second;
		/*
		if (opcode == LLVM_IR_FAdd || opcode == LLVM_IR_FSub || opcode == LLVM_IR_FMul || opcode == LLVM_IR_FDiv) {
			interestedStaticInstID2OpcodeMap.insert(std::make_pair(instID_str, opcode));
		}*/
		switch (opcode) {
		//case LLVM_IR_Add:
		//case LLVM_IR_Sub:
		//case LLVM_IR_Mul:
		//case LLVM_IR_UDiv:
		//case LLVM_IR_SDiv:
		case LLVM_IR_FAdd:
		case LLVM_IR_FSub:
		case LLVM_IR_FMul:
		case LLVM_IR_FDiv:
		case LLVM_IR_FCmp:
			interestedStaticInstID2OpcodeMap.insert(std::make_pair(instID_str, opcode));
			break;
		default:
			// Do nothing
			break;
		}
	}
	unsigned size_interest_inst = interestedStaticInstID2OpcodeMap.size();

	std::list<Vertex> topo_nodes;
	/// 1. Topologically sort the DDDG graph
	//boost::topological_sort(graph_arg, std::back_inserter(topo_nodes));
	boost::topological_sort(graph_, std::back_inserter(topo_nodes));
	
	/// 2. Calculate timestamp for the DDDG graph
	
	// Get dynamic instruction id vector
	std::vector<std::string> dynInstIDvec(numTotalNodes, "");
	initInstID(dynInstIDvec);

	// Get vertex_id to memory address mapping first
	vertexID2AddressAndSizeMapTy vertexID2AddressAndSizeMap;
	initAddressAndSize(vertexID2AddressAndSizeMap);
	//std::unordered_map<unsigned, long long int> vertexID2AddressMap;
	//initAddress(vertexID2AddressMap);

	// Initialize timestamp of all vertices in the graph
	/*
	for (unsigned i = 0; i < size_static_inst; i++) {
		std::string temp_inst_id = staticInstIDvec.at(i);
		staticInstID2TimestampMap.insert(std::make_pair(temp_inst_id, 0));
	}*/
	Vertex2TimestampMap.clear();
	std::list<Vertex>::reverse_iterator vi, ve;
	for (vi = topo_nodes.rbegin(), ve = topo_nodes.rend(); vi != ve; ++vi) {
		Vertex tmp_vertex = *vi;
		Vertex2TimestampMap.insert(std::make_pair(tmp_vertex, 0));
	}

	// Calculate and update timestamp for all vertices in the graph
	staticInstID2TimestampInfoMap.clear();
	it_map = interestedStaticInstID2OpcodeMap.begin();
	ie_map = interestedStaticInstID2OpcodeMap.end();
	for (; it_map != ie_map; ++it_map) {
		std::string inst_id_str = it_map->first;
		timestamp2VertexMMapTy timestamp2VertexMMap;
		for (vi = topo_nodes.rbegin(), ve = topo_nodes.rend(); vi != ve; vi++) {
			unsigned current_vertex_id = vertexToName[*vi];
			std::string current_nodeID_str = dynInstIDvec.at(current_vertex_id);

			in_edge_itr_ty in_edge_i, in_edge_e;
			unsigned max_timestamp = 0;
			//unsigned max_timestamp = staticInstID2TimestampMap[current_instID];
			//for (boost::tie(in_edge_i, in_edge_e) = boost::in_edges(*vi, graph_arg); in_edge_i != in_edge_e; ++in_edge_i) {
				//Vertex parent_vertex = boost::source(*in_edge_i, graph_arg);

			for (boost::tie(in_edge_i, in_edge_e) = boost::in_edges(*vi, graph_); in_edge_i != in_edge_e; ++in_edge_i) {
				Vertex parent_vertex = boost::source(*in_edge_i, graph_);
				// Get the inst_id for this parent_vertex
				//unsigned pre_node_id = vertexToName[parent_vertex];
				//std::string pre_node_instID = dynInstIDvec.at(pre_node_id);
				unsigned tmp_timestamp = getvertexTimestamp(parent_vertex);
				if (tmp_timestamp > max_timestamp) {
					max_timestamp = tmp_timestamp;
				}
			} // Done for all parent vertices

			if (current_nodeID_str == inst_id_str) {
				max_timestamp++;
				timestamp2VertexMMap.insert(std::make_pair(max_timestamp, *vi));
			}
			//staticInstID2TimestampMap[current_instID] = max_timestamp;
			Vertex2TimestampMap[*vi] = max_timestamp;
		} // Explored all vertices in topologically sorted graph

		// Store timestamp for instances of an interesting static instruction
		staticInstID2TimestampInfoMap.insert(std::make_pair(inst_id_str, timestamp2VertexMMap));
	} // Calculate timestamps for all static instructions

	/// 3. Analyze consecutive memory accesses or non-unit constant stride accesses and 
	///		 refine 'staticInstID2TimestampInfoMap'
	///		 Method: For each interesting static instruction in 'interestedStaticInstID2OpcodeMap',
	///						 (a). We find whether its two predecessors contain load instructions by analyzing
	///						 their opcodes. If yes, we assign a pair (loadVertex, its address) into a 
	///						 vector, 'addressVec'. Otherwise, just assign (loadVertex, 0) as the pair.
	///						 (b). Then we find the successor for this interesting static instruction and check
	///						 whether it is a store instruction by analyzing its opcode. If yes, we also assign 
	///						 a pair (storeVertex, its address) into 'addressVec'. Otherwise, just assign the
	///						 pair with (storeVertex, 0).
	///						 (c). After finishing collecting all information in 'addressVec', we use std::sort
	///						 to sort addressVec based on address.
	///						 (d). Detect whether these addresses are unit-stride memory accesses or constant
	///						 non-unit-stride memory accesses, then based on this information, re-partition
	///						 timestamp in 'staticInstID2TimestampInfoMap'

	num_dynInst = 0;
	num_unitVectorInst = 0;
	num_constNonUnitVectorInst = 0;
	total_sizeOf_vectorInst = 0;
	num_singletonInst = 0;
	num_byte_transferred = 0;
	memSizeTy memSize = vertexID2AddressAndSizeMap.begin()->second.second >> BYTE_SIZE_SHIFT;
	//instancesToPartitionIDMap.clear();
	EdgeWeightMap edge_to_parid = get(boost::edge_weight, graph_);
	staticInstID2TimestampInfoMapTy::iterator it_mem = staticInstID2TimestampInfoMap.begin();
	staticInstID2TimestampInfoMapTy::iterator ie_mem = staticInstID2TimestampInfoMap.end();
	for (; it_mem != ie_mem; ++it_mem) {
		// Independent instruction instances.
		//timestamp2VertexMMapTy::iterator it_mmap = it_mem->second.begin();
		//timestamp2VertexMMapTy::iterator ie_mmap = it_mem->second.end();

		// For each timestamp, we analyze its vertices associated to this timestamp.
		unsigned maxTimestamp = it_mem->second.size() == 0 ? 0 : it_mem->second.rbegin()->first;
		num_dynInst += it_mem->second.size();
		for (unsigned time = 1; time < maxTimestamp + 1; time++) {
			unsigned partID = 1;
			VertexAddressTuplePairVecTy VertexAddressTuplePairForATimestampVec;
			//unsigned element_size = it_mem->second.count(time);

			std::pair<timestamp2VertexMMapTy::iterator, timestamp2VertexMMapTy::iterator> interestVertices;
			interestVertices = it_mem->second.equal_range(time);
			timestamp2VertexMMapTy::iterator it_vertex = interestVertices.first;
			timestamp2VertexMMapTy::iterator ie_vertex = interestVertices.second;
			/// FIXME: Currently, we assume that data of all memory accesses contains the same size. Thus, 
			///				 we just simply set memSize = vertexID2AddressAndSizeMap.begin()->second >> BYTE_SIZE_SHIFT.
			///				 In the future, memSize should supports and is related to each load/store instruction. This 
			///				 need another map to track the information.
			for (; it_vertex != ie_vertex; ++it_vertex) {
				Vertex tmpVertex = it_vertex->second;
				
				// Analyze the two predecessors
				in_edge_itr_ty in_edge_i, in_edge_e;
				boost::tie(in_edge_i, in_edge_e) = boost::in_edges(tmpVertex, graph_);
				int edge_parid = 0;
				AddressTuple addr_tuple;
				for (; in_edge_i != in_edge_e; ++in_edge_i) {
					Vertex src_vertex = boost::source(*in_edge_i, graph_);
					unsigned src_vertex_id = vertexToName[src_vertex];
					unsigned vertex_opcode = getvertexOpcode(src_vertex);
					//std::pair<long long int, unsigned> load_AddressAndSize;
					long long int load_address;
					edge_parid = edge_to_parid[*in_edge_i];
					if (edge_parid == LEFT_NODE || edge_parid == RIGHT_NODE) {
						if (vertex_opcode == LLVM_IR_Load) {
							//load_AddressAndSize = vertexID2AddressAndSizeMap[src_vertex_id];
							//load_address = vertexID2AddressMap[src_vertex_id];
							load_address = vertexID2AddressAndSizeMap[src_vertex_id].first;

							if (edge_parid == LEFT_NODE) {
								addr_tuple.raddr1 = load_address;
							}
							else if (edge_parid == RIGHT_NODE) {
								addr_tuple.raddr2 = load_address;
							}
							else {
								assert(false && "Error: edge_parid is neither LEFT_NODE nor RIGHT_NODE!\n");
							}

							// Calculate memory bytes transferred
							num_byte_transferred += (vertexID2AddressAndSizeMap[src_vertex_id].second >> BYTE_SIZE_SHIFT);
						}
						else {
							if (edge_parid == LEFT_NODE) {
								addr_tuple.raddr1 = 0;
							}
							else if (edge_parid == RIGHT_NODE) {
								addr_tuple.raddr2 = 0;
							}
							else {
								assert(false && "Error: edge_parid is neither LEFT_NODE nor RIGHT_NODE!\n");
							}
						}
					}

				} // Collect 1st and 2rd operands for a binaryOp instruction instance

				// Analyze the only one successor
				out_edge_itr_ty out_edge_i, out_edge_e;
				boost::tie(out_edge_i, out_edge_e) = boost::out_edges(tmpVertex, graph_);
				for (; out_edge_i != out_edge_e; ++out_edge_i) {
					Vertex target_vertex = boost::target(*out_edge_i, graph_);
					unsigned target_opcode = getvertexOpcode(target_vertex);
					unsigned target_vertex_id = vertexToName[target_vertex];
					//std::pair<long long int, unsigned> store_AddressAndSize;
					long long int store_address;
					if (target_opcode == LLVM_IR_Store) {
						//store_AddressAndSize = vertexID2AddressAndSizeMap[target_vertex_id];
						//store_address = vertexID2AddressMap[target_vertex_id];
						store_address = vertexID2AddressAndSizeMap[target_vertex_id].first;
						// Calculate memory bytes transferred
						num_byte_transferred += (vertexID2AddressAndSizeMap[target_vertex_id].second >> BYTE_SIZE_SHIFT);
					}
					else {
						//store_AddressAndSize = std::make_pair(0, 0);
						store_address = 0;
					}
					addr_tuple.waddr = store_address;

				} // Collect write address of a binaryOp instruction instance
				//std::string tmp_static_instID_str = it_mem->first;
				//unsigned tmp_timestamp = it_mmap->first;
				//Key key = { tmp_static_instID_str, tmpVertex, tmp_timestamp };
				VertexAddressTuplePairForATimestampVec.push_back(std::make_pair(tmpVertex, addr_tuple));
			} // Collect address tuple for a interesting vertex (an interesting instruction instance)

			// Sort VertexAddressTuplePairForATimestampVec based on AddressTuplePair
			VertexAddressTuplePairVecTy::iterator it_vertexAddrTuple = VertexAddressTuplePairForATimestampVec.begin();
			VertexAddressTuplePairVecTy::iterator ie_vertexAddrTuple = VertexAddressTuplePairForATimestampVec.end();
			std::sort(it_vertexAddrTuple, ie_vertexAddrTuple, pairCompare<pairVertexAddressTupleTy>);

			std::string tmp_static_instID_str = it_mem->first;

			vectorizationAnalysis(VertexAddressTuplePairForATimestampVec, tmp_static_instID_str, time, memSize, partID);
		} // Collect address tuples with the same timestamp for a static interesting instruction

	} // Collect address tuples for all static interesting instructions

//#ifdef NEED_TO_MODEFY_THIS_PART
	/// Write staticInstID2TimestampMap into a file
	//std::string file_name = inputPath + benchName + "_timestamp.csv";
	std::string file_name = outputPath + benchName + "_timestamp.csv";

	ofstream timestamp_file(file_name);
	if (timestamp_file.is_open()) {
			
		double vectorPercent = (double) total_sizeOf_vectorInst / (double) num_dynInst * 100.0;
		double averageLengthPerVector;
		double unitVecPercent;
		double nonUnitVecPercent;
		// Arithmetic Intensity is number of FLOPs per memory access
		double ArithmeticIntensity;

		if( num_unitVectorInst == 0 && num_constNonUnitVectorInst == 0 ) {
			averageLengthPerVector = 0.0;
			unitVecPercent = 0.0;
			nonUnitVecPercent = 0.0;
		}
		else {
			averageLengthPerVector = (double)total_sizeOf_vectorInst / (double)(num_unitVectorInst + num_constNonUnitVectorInst);
			unitVecPercent = (double)num_unitVectorInst / (double)(num_unitVectorInst + num_constNonUnitVectorInst) * 100.0;
			nonUnitVecPercent = (double)num_constNonUnitVectorInst / (double)(num_unitVectorInst + num_constNonUnitVectorInst) * 100.0;
		}

		assert(num_byte_transferred!=0 && "No any memory byte transferred!");
		ArithmeticIntensity = (double)num_dynInst / (double)num_byte_transferred;

		timestamp_file << "Total number of instances of interesting instructions, " << num_dynInst << std::endl;
		//timestamp_file << num_dynInst << ",";

		timestamp_file << "Total number of singleton instructions, " << num_singletonInst << std::endl;
		//timestamp_file << num_singletonInst << ",";

		timestamp_file << "Total number of unit-stride vector instructions of interesting instructions, " << num_unitVectorInst << std::endl;
		//timestamp_file << num_unitVectorInst << ",";

		timestamp_file << "Total number of constant non-unit-stride vector instructions of interesting instructions, " << num_constNonUnitVectorInst << std::endl;
		//timestamp_file << num_constNonUnitVectorInst << ",";

		// Total size of vector instructions: how many instruction instances that construct all these vector instructions
		timestamp_file << "Total size of vector instructions, " << total_sizeOf_vectorInst << std::endl;
		//timestamp_file << total_sizeOf_vectorInst << ",";

		timestamp_file << "Average length per vector instruction, " << averageLengthPerVector << std::endl;
		//timestamp_file << averageLengthPerVector << ",";

		timestamp_file << "Vector Percentage (%), " << vectorPercent << std::endl;
		//timestamp_file << vectorPercent << ",";

		timestamp_file << "Unit-stride vector instruction percentage (%), " << unitVecPercent << std::endl;
		//timestamp_file << unitVecPercent << ",";

		timestamp_file << "Constant non-unit-stride vector instruction percentage (%), " << nonUnitVecPercent << std::endl;
		//timestamp_file << nonUnitVecPercent << ",";

		timestamp_file << "Arithmetic Intensity of " << benchName << " kernel, " << ArithmeticIntensity << std::endl;
		//timestamp_file << ArithmeticIntensity << std::endl;

		timestamp_file.close();
	}
	else {
		assert(false && "Error: Could not open timestamp_file!\n");
	}
//#endif // End of NEED_TO_MODEFY_THIS_PART

	std::cout << "Finish timestamp calculation" << std::endl;
}

unsigned BaseDatapath::getvertexTimestamp(Vertex vertex_arg) {
	//return staticInstID2TimestampMap[vertex_instID];
	return Vertex2TimestampMap[vertex_arg];
}

unsigned BaseDatapath::getvertexOpcode(Vertex vertex_arg) const {
	unsigned vertex_id = vertexToName[vertex_arg];
	return microop.at(vertex_id);
}

void BaseDatapath::vectorizationAnalysis(VertexAddressTuplePairVecTy &pairVec, std::string& staticID, unsigned& timestamp_arg, memSizeTy memSize, unsigned& partID) {
	Vertex pre_vertex;
	AddressTuple pre_addrTuple;
	AddressTuple pre_stride = {0, 0, 0};
	AddressTuple init_tuple = { 0, 0, 0 };

	unsigned numInstancesInaPart = 0;
	bool init_flag = true;

	bool unitStrideFlag;
	bool constNonUnitStrideFlag;
	// Statistics
	//uint64_t NumUnitStrideVecInst = 0;
	//uint64_t NumConstantNonUnitStrideVecInst = 0;
	//uint64_t NumConstantVecInst = 0;
	//uint64_t NumInstance = 0;

	for (unsigned i = 0; i < pairVec.size(); i++) {
		
		Vertex cur_vertex = pairVec[i].first;
		AddressTuple cur_addrTuple = pairVec[i].second;
		AddressTuple cur_stride;

		//memSizeTy memSize = ver2addrAndSizeMap[cur_vertex].second >> BYTE_SIZE_SHIFT;
		KeyTy key = { staticID, cur_vertex, timestamp_arg };

		if (init_flag == true) {
			//partID = 1;
			init_flag = false;
			numInstancesInaPart = 1;
			cur_stride = init_tuple;
			unitStrideFlag = false;
			constNonUnitStrideFlag = false;
		}
		else {
			cur_stride = cur_addrTuple - pre_addrTuple;
			//if (numInstancesInaPart == 1 || numInstancesInaPart == 0) {
			if (numInstancesInaPart == 1) {
				assert(memSize!=0 && "memSize = 0, error!");
				unitStrideFlag = isUnitStride(cur_stride, memSize);
				constNonUnitStrideFlag = unitStrideFlag ? false : isConstantNonUnitStride(cur_stride, memSize);
				numInstancesInaPart++;

				if (i == (pairVec.size() - 1) ) {
					total_sizeOf_vectorInst += numInstancesInaPart;
					if (unitStrideFlag == true) {
						num_unitVectorInst++;
					}
					else if (constNonUnitStrideFlag == true) {
						num_constNonUnitVectorInst++;
					}
					else {
						assert(false && "numInstancesInaPart exceeds GPU architecture limitation, but unitStrideFlag ^ constNonUnitStrideFlag is false");
					}
				}

			}
			else {
				if (cur_stride == pre_stride) {
					// Stride of (A[i] - A[i-1]) is equal to that of (A[i+1] - A[i])
					numInstancesInaPart++;

					if ((i == (pairVec.size() - 1)) && numInstancesInaPart != GPU_VEC_ADD_MUL_LIMITATION) {
						total_sizeOf_vectorInst += numInstancesInaPart;
						if (unitStrideFlag == true) {
							num_unitVectorInst++;
						}
						else if (constNonUnitStrideFlag == true) {
							num_constNonUnitVectorInst++;
						}
						else {
							assert(false && "numInstancesInaPart exceeds GPU architecture limitation, but unitStrideFlag ^ constNonUnitStrideFlag is false");
						}
					}

				}
				else {
					// Stride of (A[i] - A[i-1]) is not equal to that of (A[i+1] - A[i]), need to add a new partition

					// This part is reserved for collecting statistics information
					if (unitStrideFlag == true) {
						total_sizeOf_vectorInst += numInstancesInaPart;
						num_unitVectorInst++;
						numInstancesInaPart = 1;
						cur_stride = init_tuple;
					}

					if (constNonUnitStrideFlag == true) {
						// If a partition only contains two instances and has a constant-non-unit stride (A[i] - A[i-1]), 
						// which is different from the stride (A[i+1] - A[i]) of the consecutive vertice. Then just split 
						// A[i-1] and A[i] into two seperate partitions.

						// If a partition contains more than two instances, then we need to increment num_constNonUnitVectorInst.
						// FIXME: Should we also increment total_sizeOf_vectorInst? This depends on architectures we target. If
						//				we target on vector machine, there are constant non-unit stride vector instructions, therefore,
						//				we need to increase total_sizeOf_vectorInst by numInstancesInaPart. I am not sure whether GPU in
						//				embedded domain can support constant non-unit stride vector instructions, for example Mali-T628
						//				GPU.
						//				In this version, I assume embedded-GPUs support constant non-unit stride vector instructions.

						if (numInstancesInaPart == 2) {
							KeyTy tmp_key = { staticID, pre_vertex, timestamp_arg };
							instancesToPartitionIDMap[tmp_key]++;
							num_singletonInst++;
							numInstancesInaPart = 2;
							// cur_stride should be assigned to pre_stride and not be initialized
							// If cur_vertex is the last instance in pairVec, then we need to seperate it with pre_vertex
							if ( i == (pairVec.size() - 1) ) {
								num_singletonInst++;
								partID++;
							}
						}
						else {
							num_constNonUnitVectorInst++;
							total_sizeOf_vectorInst += numInstancesInaPart;
							numInstancesInaPart = 1;
							cur_stride = init_tuple;
						}
					}
					assert( (unitStrideFlag^constNonUnitStrideFlag) && "Error! Both unitStrideFlag and constNonUnitStrideFlag are false or true");

					partID++;
				}

				/*
				if (numInstancesInaPart == GPU_VEC_ADD_MUL_LIMITATION) {
					total_sizeOf_vectorInst += numInstancesInaPart;

					// This part is reserved for collecting statistics information
					if (unitStrideFlag == true) {
						num_unitVectorInst++;
					}
					else if (constNonUnitStrideFlag == true) {
						num_constNonUnitVectorInst++;
					}
					else {
						// The program should not enter this part
						assert(false && "Neither unitStrideFlag nor constNonUnitStrideFlag");
					}

					numInstancesInaPart = 1;
					partID++;
					cur_stride = init_tuple;
				}*/
			}
		}
		
		instancesToPartitionIDMap.insert(std::make_pair(key, partID));
		//std::cout << "instancesToPartitionIDMap size = " << instancesToPartitionIDMap.size() << std::endl;

		// If number of instruction instances we grouped is exceeding the limitation of GPU computational
		// power, set a new partition for them
		if (numInstancesInaPart == GPU_VEC_ADD_MUL_LIMITATION) {
			total_sizeOf_vectorInst += numInstancesInaPart;
			if (unitStrideFlag == true) {
				num_unitVectorInst++;
			}else if (constNonUnitStrideFlag == true) {
				num_constNonUnitVectorInst++;
			}
			else {
				assert(false && "numInstancesInaPart exceeds GPU architecture limitation, but unitStrideFlag ^ constNonUnitStrideFlag is false");
			}

			partID++;
			init_flag = true;
			//numInstancesInaPart = 0;
			//cur_stride = init_tuple;
		}

		pre_vertex = cur_vertex;
		pre_addrTuple = cur_addrTuple;
		pre_stride = cur_stride;
	}

}

bool BaseDatapath::isUnitStride(const AddressTuple& strideTuple, memSizeTy size) {
	unsigned r1;
	unsigned r2;
	unsigned waddr;

	unsigned raddr1 = (strideTuple.raddr1 < 0) ? (-strideTuple.raddr1) : strideTuple.raddr1;
	unsigned raddr2 = (strideTuple.raddr2 < 0) ? (-strideTuple.raddr2) : strideTuple.raddr2;
	unsigned uwaddr = (strideTuple.waddr  < 0) ? (-strideTuple.waddr) : strideTuple.waddr;

	if (size == 1) {
		r1 = raddr1;
		r2 = raddr2;
		waddr = uwaddr;
	}
	else if (size == 2) {
		r1 = raddr1 >> 1;
		r2 = raddr2 >> 1;
		waddr = uwaddr >> 1;
	}
	else if (size == 4) {
		r1 = raddr1 >> 2;
		r2 = raddr2 >> 2;
		waddr = uwaddr >> 2;
	}
	else if (size == 8) {
		r1 = raddr1 >> 3;
		r2 = raddr2 >> 3;
		waddr = uwaddr >> 3;
	}
	else {
		r1 = raddr1 / size;
		r2 = raddr2 / size;
		waddr = uwaddr / size;
	}
	
	return (r1 == 0 || r1 == 1) && (r2 == 0 || r2 == 1) && (waddr == 0 || waddr == 1);
}

bool BaseDatapath::isConstantNonUnitStride(const AddressTuple& strideTuple, memSizeTy size) {
	unsigned r1;
	unsigned r2;
	unsigned waddr;

	unsigned raddr1 = (strideTuple.raddr1 < 0) ? (-strideTuple.raddr1) : strideTuple.raddr1;
	unsigned raddr2 = (strideTuple.raddr2 < 0) ? (-strideTuple.raddr2) : strideTuple.raddr2;
	unsigned uwaddr = (strideTuple.waddr  < 0) ? (-strideTuple.waddr)  : strideTuple.waddr;

	if (size == 1) {
		r1 = raddr1;
		r2 = raddr2;
		waddr = uwaddr;
	}
	else if (size == 2) {
		r1 = raddr1 >> 1;
		r2 = raddr2 >> 1;
		waddr = uwaddr >> 1;
	}
	else if (size == 4) {
		r1 = raddr1 >> 2;
		r2 = raddr2 >> 2;
		waddr = uwaddr >> 2;
	}
	else if (size == 8) {
		r1 = raddr1 >> 3;
		r2 = raddr2 >> 3;
		waddr = uwaddr >> 3;
	}
	else {
		r1 = raddr1 / size;
		r2 = raddr2 / size;
		waddr = uwaddr / size;
	}

	return (r1 > 1) || (r2 > 1) || (waddr > 1);
}

void BaseDatapath::parallelismProfileDDDG() {
	VERBOSE_PRINT(std::cout << "DEBUG-INFO: [analysis_parallelism-profile] Analyze parallelism of nodes on DDDG\n" << std::endl);

	assert(staticInstID2OpcodeMap.empty() != true && "staticInstID2OpcodeMap is empty, ERROR, please check!\n");

	std::list<Vertex> topo_nodes;
	/// 1. Topologically sort the DDDG graph
	boost::topological_sort(graph_, std::back_inserter(topo_nodes));

	/// 2. Calculate level for the DDDG graph

	// Get dynamic instruction id vector
	std::vector<std::string> dynInstIDvec(numTotalNodes, "");
	initInstID(dynInstIDvec);

	// Initialize level of all vertices in the graph

	Vertex2TimestampMap.clear();
	std::list<Vertex>::reverse_iterator vi, ve;
	for (vi = topo_nodes.rbegin(), ve = topo_nodes.rend(); vi != ve; ++vi) {
		Vertex tmp_vertex = *vi;
		Vertex2LevelMap.insert(std::make_pair(tmp_vertex, 0));
	}

	// Calculate and update timestamp for all vertices in the graph
	timestamp2VertexMMapTy timestamp2VertexMMap;
	for (vi = topo_nodes.rbegin(), ve = topo_nodes.rend(); vi != ve; vi++) {
		unsigned current_vertex_id = vertexToName[*vi];
		std::string current_nodeID_str = dynInstIDvec.at(current_vertex_id);
		unsigned opcode = staticInstID2OpcodeMap[current_nodeID_str];

		unsigned vertex_level = 0;
		bool has_parent_vertex = false;

		in_edge_itr_ty in_edge_i, in_edge_e;
		for (boost::tie(in_edge_i, in_edge_e) = boost::in_edges(*vi, graph_); in_edge_i != in_edge_e; ++in_edge_i) {
			Vertex parent_vertex = boost::source(*in_edge_i, graph_);
			// Get the inst_id for this parent_vertex
			//unsigned pre_node_id = vertexToName[parent_vertex];
			//std::string pre_node_instID = dynInstIDvec.at(pre_node_id);

			unsigned tmp_level = getvertexLevel(parent_vertex);
			vertex_level = (tmp_level > vertex_level) ? tmp_level : vertex_level;
			has_parent_vertex = true;
		} // Done for all parent vertices

		if (opcode == LLVM_IR_PHI || opcode == LLVM_IR_GetElementPtr) {
			Vertex2LevelMap[*vi] = vertex_level;
		}
		else {
			Vertex2LevelMap[*vi] = has_parent_vertex == true ? (vertex_level + 1) : vertex_level;
		}

	} // Explored all vertices in topologically sorted graph

	/// 3. Calculate parallel vertices per level
	///		 We need to ignore instructions like LLVM_IR_Br, LLVM_IR_BitCast, LLVM_IR_Phi etc.
	Vertex2LevelMapTy::const_iterator it_level_max, template_it_level_max;
	//it_level_max = std::max_element(Vertex2LevelMap.begin(), Vertex2LevelMap.end(), Vertex2LevelMapCmpFunc);
	it_level_max = map_max_element(Vertex2LevelMap);
	unsigned max_level = it_level_max->second;
	parallelism_profile.clear();
	parallelism_profile.resize(max_level+1, 0);
#ifdef CHECK_INST_IN_EACH_LEVEL
	instInlevels.clear();
#endif // End of CHECK_INST_IN_EACH_LEVEL

	Vertex2LevelMapTy::iterator it_level = Vertex2LevelMap.begin();
	Vertex2LevelMapTy::iterator ie_level = Vertex2LevelMap.end();
	bool skip_flag = false;
	for (; it_level != ie_level; ++it_level) {
		Vertex vertex_tmp = it_level->first;
		unsigned current_vertex_id = vertexToName[vertex_tmp];
		std::string current_nodeID_str = dynInstIDvec.at(current_vertex_id);
		unsigned opcode_tmp = staticInstID2OpcodeMap[current_nodeID_str];
		skip_flag = ignoreDynamicInst(opcode_tmp);
		if (skip_flag != true) {
			parallelism_profile[it_level->second]++;
#ifdef CHECK_INST_IN_EACH_LEVEL
			instInlevels[it_level->second].push_back(current_nodeID_str);
#endif // End of CHECK_INST_IN_EACH_LEVEL
		}
	}

	/// Calculate average parallelism
	uint64_t total_ops = 0;
	ave_parallelism = 0.0;
	for (unsigned i = 0; i < parallelism_profile.size(); i++) {
		total_ops += parallelism_profile[i];
	}

	ave_parallelism = (float)total_ops / (float)max_level;

	//#ifdef NEED_TO_MODEFY_THIS_PART
	/// Write staticInstID2TimestampMap into a file
	//std::string file_name = inputPath + benchName + "_par_prof.csv";
	std::string file_name = outputPath + benchName + "_par_prof.csv";

	ofstream timestamp_file(file_name);
	if (timestamp_file.is_open()) {
		timestamp_file << "Level, Num of operations" << std::endl;
		for (unsigned i = 0; i < parallelism_profile.size(); i++) {
			timestamp_file << i << "," << parallelism_profile[i];
#ifdef CHECK_INST_IN_EACH_LEVEL
			for (unsigned j = 0; j < parallelism_profile[i]; j++) {
				timestamp_file << ", " << instInlevels[i][j];
			}
#endif
			timestamp_file << std::endl;
		}
		timestamp_file.close();
	}
	else {
		assert(false && "Error: Could not open parallelism_profile_file!\n");
	}
	//#endif // End of NEED_TO_MODEFY_THIS_PART

	std::cout << "Finish parallelism profile calculation" << std::endl;
}

unsigned BaseDatapath::getvertexLevel(Vertex vertex_arg) {
	return Vertex2LevelMap[vertex_arg];
}

bool BaseDatapath::ignoreDynamicInst(unsigned opcode_arg) {
	bool result;
	switch (opcode_arg) {
	case LLVM_IR_Move:
	//case LLVM_IR_Ret:
	//case LLVM_IR_Br:
	//case LLVM_IR_Switch:
	//case LLVM_IR_IndirectBr:
	//case LLVM_IR_Invoke:
	case LLVM_IR_Resume:
	case LLVM_IR_Unreachable:
	case LLVM_IR_Alloca:
	case LLVM_IR_GetElementPtr:
	case LLVM_IR_Trunc:
	case LLVM_IR_ZExt:
	case LLVM_IR_SExt:
	case LLVM_IR_FPToUI:
	case LLVM_IR_FPToSI:
	case LLVM_IR_UIToFP:
	case LLVM_IR_SIToFP:
	case LLVM_IR_FPTrunc:
	case LLVM_IR_FPExt:
	case LLVM_IR_PtrToInt:
	case LLVM_IR_IntToPtr:
	case LLVM_IR_BitCast:
	case LLVM_IR_AddrSpaceCast:
	//case LLVM_IR_PHI:
	//case LLVM_IR_Call:
	case LLVM_IR_Select:
	case LLVM_IR_VAArg:
	case LLVM_IR_ExtractElement:
	case LLVM_IR_InsertElement:
	case LLVM_IR_ShuffleVector:
	case LLVM_IR_ExtractValue:
	case LLVM_IR_InsertValue:
	case LLVM_IR_LandingPad:
	//case LLVM_IR_IndexAdd:
	//case LLVM_IR_IndexSub:
		result = true;
		break;
	default:
		result = false;
		break;
	}
	return result;
}

/*
bool Vertex2LevelMapCmpFunc(const pair<Vertex, unsigned> &p1, const pair<Vertex, unsigned> &p2) {
	return (p1.second < p2.second);
}*/

/// For parallel simulation
void BaseDatapath::copyToExecutingQueue()
{
  auto it = readyToExecuteQueue.begin();
  while (it != readyToExecuteQueue.end())
  {
    executingQueue.push_back(*it);
    it = readyToExecuteQueue.erase(it);
  }
}

/// For parallel simulation
bool BaseDatapath::step()
{

  stepExecutingQueue();
  copyToExecutingQueue();
	num_cycles++;
#ifdef WRITE_EXECUTION_RECORD
	execution_record << "\n";
	execution_record << num_cycles << ", ";
#endif // End of WRITE_EXECUTION_RECORD
	updateDelayForNodeID();

	if (executedNodes == totalConnectedNodes) {
		std::vector<unsigned>().swap(maxParentOpLatency);
		return 1;
	}
  return 0;
}

void BaseDatapath::updateDelayForNodeID() {
	bool constraint_flag = false;
	std::vector<std::map<node_idTy, unsigned>::iterator> removeNodeIDit;
	std::map<node_idTy, unsigned>::iterator it, ie;
	for (it = nodeID2delay.begin(), ie = nodeID2delay.end(); it != ie; ++it) {
		/// Constraints should be applied here. If a node_id fails to meet constraints,
		/// we do not decrease delay value of a node
		unsigned node_id = it->first;
		unsigned delay = it->second;
		constraint_flag = calculateConstraint(node_id, delay);
		if (constraint_flag == false) {
			//assert(it->second != 0 && "Delay of this node id is already 0, can not be deducted again! Error\n");
			if (it->second != 0) {
				it->second--;
#ifdef WRITE_EXECUTION_RECORD
				execution_record << it->first << "-rd, ";
#endif // End of WRITE_EXECUTION_RECORD
			}
		}
		else {
#ifdef WRITE_EXECUTION_RECORD
			execution_record << it->first << "-cd, ";
#endif // End of WRITE_EXECUTION_RECORD
		}
	}

	for (it = nodeID2delay.begin(), ie = nodeID2delay.end(); it != ie; ++it) {
		if (it->second == 0) {
			unsigned node_id = it->first;
			//readyToExecuteQueue.push_back(node_id);
			executingQueue.push_back(node_id);
			removeNodeIDit.push_back(it);
			removeConstraint(node_id);
		}
	}

	for (int i = 0; i < removeNodeIDit.size(); i++) {
		nodeID2delay.erase(removeNodeIDit.at(i));
	}
}

bool BaseDatapath::calculateConstraint(unsigned node_id, unsigned node_delay) {

	unsigned opcode = microop.at(node_id);
	unsigned op_latency = fpga_node_latency(opcode);
	if (op_latency != node_delay) {
		// This node has been counted before, no need to consider constraint again
		if (is_fadd_op(opcode)) {
			fadd_engine.erase(node_id);
		}
		else if (is_fsub_op(opcode)) {
			fsub_engine.erase(node_id);
		}
		else if (is_fmul_op(opcode)) {
			fmul_engine.erase(node_id);
		}
		else if (is_fdiv_op(opcode)) {
			fdiv_engine.erase(node_id);
		}
		else {
			// Do nothing here.
		}
		return false;
	}

	engineSlotTy::iterator it;
	unsigned size = 0;
	if (is_fadd_op(opcode)) {
		size = fadd_engine.size();
		if ( size < fpga_constraints.get_fadd_num() ) {
			// No instruction uses fadd_engine, thus no constraint
			fadd_engine.insert(node_id);
			return false;
		}
		else {
			// fadd_engine is occupied at this moment
			return true;
		}
	}
	else if (is_fsub_op(opcode)) {
		size = fsub_engine.size();
		if (size < fpga_constraints.get_fsub_num()) {
			// No instruction uses fsub_engine, thus no constraint
			fsub_engine.insert(node_id);
			return false;
		}
		else {
			// fsub_engine is occupied at this moment
			return true;
		}
	}
	else if (is_fmul_op(opcode)) {
		size = fmul_engine.size();
		if (size < fpga_constraints.get_fmul_num()) {
			// No instruction uses fmul_engine, thus no constraint
			fmul_engine.insert(node_id);
			return false;
		}
		else {
			// fmul_engine is occupied at this moment
			return true;
		}
	}
	else if (is_fdiv_op(opcode)) {
		size = fdiv_engine.size();
		if (size < fpga_constraints.get_fdiv_num()) {
			// No instruction uses fdiv_engine, thus no constraint
			fdiv_engine.insert(node_id);
			return false;
		}
		else {
			// fdiv_engine is occupied at this moment
			return true;
		}
	}
	else if (is_load_op(opcode)) {
		std::string node_part = baseAddress[node_id].first;
		return memoryBWconstraint(node_part, true);
	}
	else if (is_store_op(opcode)) {
		std::string node_part = baseAddress[node_id].first;
		return memoryBWconstraint(node_part, false);
	}
	else {
		// No constraint
		return false;
	}

}


void BaseDatapath::removeConstraint(unsigned node_id) {
	unsigned opcode = microop.at(node_id);
	engineSlotTy::iterator it;
	unsigned size = 0;
	/*
	// floating point unit constraints will be removed after executing 1 cycle.
  // This is because of its pipeline design.
	if (is_fadd_op(opcode)) {
		fadd_engine.erase(node_id);
	}
	else if (is_fsub_op(opcode)) {
		fsub_engine.erase(node_id);
	}
	else if (is_fmul_op(opcode)) {
		fmul_engine.erase(node_id);
	}
	else if (is_fdiv_op(opcode)) {
		fdiv_engine.erase(node_id);
	}*/

	if (is_load_op(opcode)) {
		std::string node_part = baseAddress[node_id].first;
		unsigned partition_id = findPartitionID(node_part);
		assert( (occupiedMemPerPartition.at(partition_id).readPort_num != 0) && "Error: occupied read port num is 0 already, but report read port constraint!\n");
		occupiedMemPerPartition.at(partition_id).readPort_num--;
	}
	else if (is_store_op(opcode)) {
		std::string node_part = baseAddress[node_id].first;
		unsigned partition_id = findPartitionID(node_part);
		assert((occupiedMemPerPartition.at(partition_id).writePort_num != 0) && "Error: occupied write port num is 0 already, but report write port constraint!\n");
		occupiedMemPerPartition.at(partition_id).writePort_num--;
	}
	else {
		// No constraint, do nothing here
	}
}

/// For parallel simulation
void BaseDatapath::stepExecutingQueue() {

	auto it = executingQueue.begin();
	int index = 0;
	while (it != executingQueue.end()) {
		unsigned node_id = *it;
		/*
		/// After loopFlattening, all the compute operations will be replaced with LLVM_IR_Move microop,
		/// This is a bug inside loopFlattening, we can not replace the opcodes.
		if ( is_compute_op(microop.at(node_id)) ) {
		cout << "This is a compute instruction: " << microop.at(node_id) << endl;
		}*/
		unsigned opcode = microop.at(node_id);
		if (is_memory_op(opcode))
		{
			std::string node_part = baseAddress[node_id].first;
			//if (registers.has(node_part) || canServicePartition(node_part))
			if (registers.has(node_part))
			{
				if (is_load_op(microop.at(node_id)))
					registers.getRegister(node_part)->increment_loads();
				else
					registers.getRegister(node_part)->increment_stores();

			}
			else {
				if (is_load_op(microop.at(node_id)))
					increment_loads(node_part);
				else
					increment_stores(node_part);
			}
		}

		executedNodes++;
#ifdef WRITE_EXECUTION_RECORD
		execution_record << node_id << "-f, ";
#endif // End of WRITE_EXECUTION_RECORD
		newLevel.at(node_id) = num_cycles;
		executingQueue.erase(it);
		updateChildren(node_id);
		it = executingQueue.begin();
		std::advance(it, index);
	}

}

/* BaseDatapath::memoryBWconstraint:
     Return: false --> no constraint
		         true  --> has constraints
		 Input arguments:
       load_or_store: true --> load;  false --> store */
bool BaseDatapath::memoryBWconstraint(std::string baseName, bool load_or_store) {
	unsigned partition_id = findPartitionID(baseName);
	if (load_or_store) {
		// Load operation
		if (occupiedMemPerPartition.at(partition_id).readPort_num < fpga_constraints.get_freadP_num()) {
			// No constraint
			occupiedMemPerPartition.at(partition_id).readPort_num++;
			return false;
		}
	}
	else {
		// Store operation
		if (occupiedMemPerPartition.at(partition_id).writePort_num < fpga_constraints.get_fwriteP_num()){
			// No constraint
			occupiedMemPerPartition.at(partition_id).writePort_num++;
			return false;
		}
	}

	return true;
}

/*
bool BaseDatapath::canServicePartition(std::string baseName) {
	unsigned partition_id = findPartitionID(baseName);
	return (occupiedBWPerPartition.at(partition_id) < numOfPortsPerPartition) ? 1 : 0;
}

bool BaseDatapath::addressRequest(std::string baseName){
	if (canServicePartition(baseName))
	{
		unsigned partition_id = findPartitionID(baseName);
		assert(occupiedBWPerPartition.at(partition_id) < numOfPortsPerPartition);
		occupiedBWPerPartition.at(partition_id)++;
		return true;
	}
	else
		return false;
}
*/

unsigned BaseDatapath::findPartitionID(std::string baseName) {
	auto partition_it = baseToPartitionID.find(baseName);
	if (partition_it != baseToPartitionID.end())
		return partition_it->second;
	else
	{
		std::cerr << "Unknown Partition Name:" << baseName << std::endl;
		std::cerr << "Need to Explicitly Declare How to Partition Each Array in the Config File" << std::endl;
		exit(0);
	}
}

void BaseDatapath::increment_loads(std::string partition) {
	partition_loads[partition]++;
}

void BaseDatapath::increment_stores(std::string partition) {
	partition_stores[partition]++;
}

/// For parallel simulation
void BaseDatapath::updateChildren(unsigned node_id)
{
  Vertex node = nameToVertex[node_id];
  out_edge_iter out_edge_it, out_edge_end;
  for (tie(out_edge_it, out_edge_end) = out_edges(node, graph_); out_edge_it != out_edge_end; ++out_edge_it)
  {
    unsigned child_id = vertexToName[target(*out_edge_it, graph_)];
    int edge_parid = edgeToParid[*out_edge_it];
    if (numParents[child_id] > 0)
    {
			unsigned curr_microop = microop.at(node_id);
			unsigned curr_op_latency = fpga_node_latency(curr_microop);
      numParents[child_id]--;
			unsigned maxParentLatency = maxParentOpLatency.at(child_id);
			if (curr_op_latency > maxParentLatency) {
				maxParentOpLatency.at(child_id) = curr_op_latency;
				maxParentLatency = curr_op_latency;
			}

			if (numParents[child_id] == 0)
      {
        unsigned child_microop = microop.at(child_id);
				unsigned op_latency = fpga_node_latency(child_microop);
				//if ((op_latency == 0 || fpga_node_latency(curr_microop) == 0)
				//	&& edge_parid != CONTROL_EDGE) {
				//if ( (op_latency == 0 && edge_parid != CONTROL_EDGE) || (child_microop == LLVM_IR_Br && edge_parid == CONTROL_EDGE) ) {
				if ((maxParentLatency == 0) || (child_microop == LLVM_IR_Br && edge_parid == CONTROL_EDGE)) {
					executingQueue.push_back(child_id);
				}
				else {
					std::map<node_idTy, unsigned>::iterator itDelay = nodeID2delay.find(child_id);
					/*
					if ( itDelay != nodeID2delay.end() ) {
						if (itDelay->second == 0) {
							readyToExecuteQueue.push_back(child_id);
							nodeID2delay.erase(itDelay);
						}
					}
					else
					{
						nodeID2delay.insert( std::make_pair(child_id, fpga_node_latency(child_microop)) );
					}*/
					if (itDelay == nodeID2delay.end()) {
						nodeID2delay.insert(std::make_pair(child_id, op_latency));
					}
				}
        numParents[child_id] = -1;
      }
    }
  }
}


unsigned BaseDatapath::getIthLoop(std::string bbName, std::string funcName) {

	bbFuncNamePairTy bb_func_pair = std::make_pair(bbName, funcName);
	bbFuncNamePair2lpNameLevelPairMapTy::iterator it = bbFuncNamePair2lpNameLevelPairMap.find(bb_func_pair);
	bbFuncNamePair2lpNameLevelPairMapTy::iterator ie = bbFuncNamePair2lpNameLevelPairMap.end();
	assert((it != ie) && "Error: There is no bb_func_pair key in bbFuncNamePair2lpNameLevelPairMap!\n");
	std::string loop_name = bbFuncNamePair2lpNameLevelPairMap[bb_func_pair].first;
	std::size_t pos = loop_name.find("-");
	std::string sub_str = loop_name.substr(pos + 1);
	unsigned ith_loop = (unsigned)std::stoi(sub_str);

	return ith_loop;
}

/// For parallel simulation
void BaseDatapath::initExecutingQueue() {
	for (unsigned i = 0; i < numTotalNodes; i++)
	{
		if (numParents[i] == 0 && finalIsolated[i] != 1) {
			executingQueue.push_back(i);
		}
	}
}

/*
void BaseDatapath::clearOccupiedBWPerPartition() {
	std::vector<unsigned>::iterator it, ie;
	for (it = occupiedBWPerPartition.begin(), ie = occupiedBWPerPartition.end(); it != ie; ++it) {
		*it = 0;
	}
}
*/
void BaseDatapath::clearOccupiedBWPerPartition() {
	std::vector<OccupiedMemBWPerPartition>::iterator it, ie;
	for (it = occupiedMemPerPartition.begin(), ie = occupiedMemPerPartition.end(); it != ie; ++it) {
		OccupiedMemBWPerPartition tmp = { 0, 0 };
		*it = tmp;
	}
}


uint64_t BaseDatapath::run_parallel_simulation() {
	num_cycles = 0;

#ifdef WRITE_EXECUTION_RECORD
	std::string func_name = kernel_names.at(0);
	//execution_record.open(inputPath + func_name + "_execution_record.txt");
	execution_record.open(outputPath + func_name + "_execution_record.txt");
	execution_record << "cycle, executed instructions\n";
#endif // End of WRITE_EXECUTION_RECORD

	/// Initialize executingQueue, readyToExecuteQueue
	setGraphForStepping();

#ifdef WRITE_EXECUTION_RECORD
	execution_record << num_cycles << ", ";
#endif // End of WRITE_EXECUTION_RECORD

	/// Run simulation
	while (!step()) {
		clearOccupiedBWPerPartition();
	}

#ifdef WRITE_EXECUTION_RECORD
	execution_record.close();
#endif // End of WRITE_EXECUTION_RECORD

	return num_cycles;
}

uint64_t BaseDatapath::fpga_estimation_one_more_subtrace_for_recII_calculation() {
	/// This function will be invoked if user applys loop pipelining. It is used to calculate recurrence II.
	/// In this version, Lin-analyzer will first sample one more sub-trace and calculate its iteration latency, 
	/// IL_ii, via ASAP scheduling without resource constraints. Then Lin-analyzer will sample the original 
	/// sub-trace again based on loop unrolling configurations and calculate its iteration latency, IL_asap, 
	/// via ASAP. Then the recurrence II is obtained by
	///						rec_II = IL_ii - IL_asap
	/// NOTE: Cannot invoke fpga_estimation_one_more_subtrace_for_recII_calculation() and fpga_estimation() 
	/// in the same BaseDatapath instance.

	removeInductionDependence();
	removePhiNodes();
	//dumpGraph();
	if (enable_store_buffer == true) {
		storeBuffer();
	}
	//dumpGraph();
	// Add fpga node latency as edge weights in the graph
	edge_iter edge_it, edge_ie;
	EdgeWeightMap edge2weightMap = boost::get(boost::edge_weight, graph_);
	for (boost::tie(edge_it, edge_ie) = edges(graph_); edge_it != edge_ie; ++edge_it) {
		uint8_t weight = edge2weightMap[*edge_it];
		if (weight == CONTROL_EDGE) {
			boost::put(boost::edge_weight, graph_, *edge_it, 0);
		}
		else {
			Vertex source_node = boost::source(*edge_it, graph_);
			unsigned source_node_id = vertexToName[source_node];
			unsigned source_opcode = microop.at(source_node_id);
			unsigned fpga_latency = fpga_node_latency(source_opcode);
			boost::put(boost::edge_weight, graph_, *edge_it, fpga_latency);
		}
	}

	/// Before performing scheduling, read array information of this application first
	arrayN2sizeWordsizeBytePr.clear();
	readArrayInfo(arrayN2sizeWordsizeBytePr);

	/// Perform ASAP scheduling on the graph
	max_asapTy max_asap = asap_scheduling(graph_, asapSchedTime);

	/// Perform ALAP scheduling on the graph to get the lower bound of memory bandwidth, which 
	/// releases Res_II_mem constraint.
	alap_scheduling(graph_, alapSchedTime, max_asap.max_cycles, max_asap.max_level);

	return max_asap.max_cycles;
}

uint64_t BaseDatapath::fpga_estimation() {
	std::cout << "DEBUG-INFO: [fpga_estimation] Estimating iteration latency (IL) and initiation interval (II)" << std::endl;

	num_cycles = 0;
	max_II = 1;
	removeInductionDependence();
	calculateInstructionDistribution();
	parallelismProfileDDDG();
	removePhiNodes();
	if (enable_store_buffer == true) {
		storeBuffer();
	}

	// Add fpga node latency as edge weights in the graph
	edge_iter edge_it, edge_ie;
	EdgeWeightMap edge2weightMap = boost::get(boost::edge_weight, graph_);
	for (boost::tie(edge_it, edge_ie) = edges(graph_); edge_it != edge_ie; ++edge_it) {
		uint8_t weight = edge2weightMap[*edge_it];
		if (weight == CONTROL_EDGE) {
			boost::put(boost::edge_weight, graph_, *edge_it, 0);
		}
		else {
			Vertex source_node = boost::source(*edge_it, graph_);
			unsigned source_node_id = vertexToName[source_node];
			unsigned source_opcode = microop.at(source_node_id);
			unsigned fpga_latency = fpga_node_latency(source_opcode);
			boost::put(boost::edge_weight, graph_, *edge_it, fpga_latency);
		}
	}

	/// Print the graph to test whether we successfully modify weight of edges in the graph.
	//dumpGraph();

	/// Before performing scheduling, read array information of this application first
	arrayN2sizeWordsizeBytePr.clear();
	readArrayInfo(arrayN2sizeWordsizeBytePr);

	/// Check status of loop pipelining
	enable_pipeline = readPipeliningConfig();

	/// Perform ASAP scheduling on the graph
	max_asapTy max_asap = asap_scheduling(graph_, asapSchedTime);
	uint64_t max_asap_il = max_asap.max_cycles;

	/// Perform ALAP scheduling on the graph
	alap_scheduling(graph_, alapSchedTime, max_asap_il, max_asap.max_level);

	// Claculate critical path
	critical_path_extraction(asapSchedTime, alapSchedTime, cPathNodes);

	/// Perform Resource Constrained List Scheduling
	uint64_t max_rc_list_il = rc_list_scheduling(graph_, asapSchedTime, alapSchedTime, rcListSchedTime, cPathNodes, enable_pipeline);

	LoopLatencyTy lp_latency = getLoopTotalLatency(max_rc_list_il, max_II, target_loop_name, target_loop_level, enable_pipeline);

	if (lp_latency.enable_pipeline == true) {
		num_cycles = lp_latency.with_pipeline;
	}
	else {
		num_cycles = lp_latency.no_pipeline;
	}

	loopInfoTy loop_info = { target_loop_name, target_loop_level, target_lp_level_unroll_factor, enable_pipeline, num_cycles, IL_asap, rc_list_il, max_II, ResII_mem, limited_mem_name, ResII_op, limited_op_name, RecII };
	resourceTy fpga_rs = { dsp_used, bram18k_used, ff_used, lut_used, fadd_used, fsub_used, fmul_used, fdiv_used };
	sharedMemTy sharedMem = {shared_loads, repeated_stores};
	dump_summary(summary, loop_info, fpga_rs, sharedMem, arrayName2memeff, limited_fop_unit_types, subTraceInst, ave_parallelism, arrayName2maxMemOpNum_subtrace, arrayName2aveLoadAccessPerBank_subtrace, arrayName2aveStoreAccessPerBank_subtrace);
	std::cout << "DEBUG-INFO: [fpga_estimation] Finished" << std::endl;
	return num_cycles;
}

max_asapTy BaseDatapath::asap_scheduling(Graph& graph_tmp, schedTimeTy& schedTime) {
	std::cout << "DEBUG-INFO: [fpga_estimation-ASAP] Without FPGA Resource Constraints" << std::endl;
	max_asapTy max_asap = {0, 0};

	// Initialize schedTime
	schedTime.assign(numTotalNodes, 0);

	// cStep2nodeIDMap
	cStep2nodeIDMapTy cStep2nodeIDMap_tmp;

	// Get edge weight map
	EdgeWeightMap edge2weightMap = boost::get(boost::edge_weight, graph_tmp);

	// Update schedTime without any resource constraints
	for (unsigned i = 0; i < numTotalNodes; i++) {
		/*
		if (i < __BEGIN__ || i > __END__) {
			///FIXME: This section is only used to analyze specific code region, need to remove it later.
			///       Numbers here only work for convolution2d
			continue;
		}*/

		unsigned curr_node_id = i;
		Vertex curr_vertex = nameToVertex[curr_node_id];
		in_edge_iter edge_it, edge_ie;
		// Iterate all parent nodes of this current node, and assign new weight with the largest value
		// of (parent's weight + its edge weight) among all parent nodes
		if (boost::in_degree(curr_vertex, graph_tmp) == 0) {
			schedTime[curr_node_id] = 0;
			continue;
		}

		unsigned maxCurrStartTime = 0;
		for (boost::tie(edge_it, edge_ie) = boost::in_edges(curr_vertex, graph_tmp); edge_it != edge_ie; ++edge_it) {
			Vertex parent_vertex = boost::source(*edge_it, graph_tmp);
			unsigned parent_node_id = vertexToName[parent_vertex];
			unsigned parent_start_time = schedTime[parent_node_id];
			unsigned weight = edge2weightMap[*edge_it];
			unsigned curr_vertex_start_time = parent_start_time + weight;
			maxCurrStartTime = (maxCurrStartTime < curr_vertex_start_time) ? curr_vertex_start_time : maxCurrStartTime;
		}
		schedTime[curr_node_id] = maxCurrStartTime;
		
		// Store node ids into a cStep2nodeIDMap_tmp map, so that we can calculate its FPGA resource requirement
		cStep2nodeIDMap_tmp[maxCurrStartTime].push_back(curr_node_id);
	}

	calculateFPGAResRequired(cStep2nodeIDMap_tmp, asap_FPGA_resources);

	schedTimeTy::iterator it_max = std::max_element(schedTime.begin(), schedTime.end());
	max_asap.max_level = *it_max;
	unsigned max_latency = 0;
	std::vector<unsigned> maxLevel_nodesVec = cStep2nodeIDMap_tmp[max_asap.max_level];
	for (unsigned i = 0; i < maxLevel_nodesVec.size(); i++) {
		unsigned node_id_at_maxLevel = maxLevel_nodesVec.at(i);
		unsigned opcode = microop.at(node_id_at_maxLevel);
		unsigned latency = fpga_node_latency(opcode);
		max_latency = (max_latency > latency) ? max_latency : latency;
	}

	if (enable_extra_scalar == true) {
		max_asap.max_cycles = max_asap.max_level + max_latency;
	}
	else {
		max_asap.max_cycles = max_asap.max_level + max_latency - 1;
	}
	
	IL_asap = max_asap.max_cycles;

	arrayName2NumPartTy arrayName2NumPart = asap_FPGA_resources.getArrayName2NumPartMap();
	arrayName2NumPartTy::iterator array2part_it = arrayName2NumPart.begin();
	arrayName2NumPartTy::iterator array2part_ie = arrayName2NumPart.end();

	std::cout << "DEBUG-INFO: [fpga_estimation-ASAP] IL of this loop level: " << max_asap.max_cycles << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ASAP] Resource utilization: DSP=" << asap_FPGA_resources.getDSP() << "; ";
	std::cout << "FF=" << asap_FPGA_resources.getFF() << "; " << "LUT=" << asap_FPGA_resources.getLUT() << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ASAP] Number of Fadd=" << asap_FPGA_resources.getNumFadd() << "; ";
	std::cout << "Fsub=" << asap_FPGA_resources.getNumFsub() << "; " << "Fmul=" << asap_FPGA_resources.getNumFmul() << "; ";
	std::cout << "Fdiv=" << asap_FPGA_resources.getNumFdiv() << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ASAP] Number of array partitions of ";
	for (; array2part_it != array2part_ie; ++array2part_it) {
		std::cout << array2part_it->first << "=" << array2part_it->second << "; ";
	}
	std::cout << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ASAP] Finished" << std::endl;
	std::cout << "##########" << std::endl;
	cStep2nodeIDMapTy().swap(cStep2nodeIDMap_tmp);
	return max_asap;
}

void BaseDatapath::alap_scheduling(Graph& graph_tmp, schedTimeTy& schedTime, unsigned maxCstep, unsigned max_level) {
	std::cout << "DEBUG-INFO: [fpga_estimation-ALAP] Without FPGA Resource Constraints" << std::endl;
	// Initialize schedTime
	schedTime.assign(numTotalNodes, 0);

	// Create cStep2nodeIDMap_tmp for FPGA resource calculation
	cStep2nodeIDMapTy cStep2nodeIDMap_tmp;

	// Get edge weight map
	EdgeWeightMap edge2weightMap = boost::get(boost::edge_weight, graph_tmp);

	// Update schedTime without any resource constraints
	for (int i = numTotalNodes - 1; i >= 0; i--) {
		/*
		if (i < __BEGIN__ || i > __END__) {
			///FIXME: This section is only used to analyze specific code region, need to remove it later.
			///       Numbers here only work for convolution2d
			continue;
		}*/

		unsigned curr_node_id = i;
		Vertex curr_vertex = nameToVertex[curr_node_id];
		out_edge_iter edge_it, edge_ie;
		// Iterate all parent nodes of this current node, and assign new weight with the largest value
		// of (parent's weight + its edge weight) among all parent nodes
		if (boost::out_degree(curr_vertex, graph_tmp) == 0) {
			schedTime[curr_node_id] = max_level;
			continue;
		}

		unsigned minCurrStartTime = max_level;
		for (boost::tie(edge_it, edge_ie) = boost::out_edges(curr_vertex, graph_tmp); edge_it != edge_ie; ++edge_it) {
			Vertex child_vertex = boost::target(*edge_it, graph_tmp);
			unsigned child_node_id = vertexToName[child_vertex];
			unsigned child_start_time = schedTime[child_node_id];
			unsigned weight = edge2weightMap[*edge_it];
			unsigned curr_vertex_start_time = child_start_time - weight;
			/*
			if (child_node_id < __BEGIN__ || child_node_id > __END__) {
				///FIXME: This section is only used to analyze specific code region, need to remove it later.
				///       Numbers here only work for convolution2d
				continue;
			}*/
			minCurrStartTime = (minCurrStartTime > curr_vertex_start_time) ? curr_vertex_start_time : minCurrStartTime;
		}
		schedTime[curr_node_id] = minCurrStartTime;

		// Store node ids into a cStep2nodeIDMap_tmp map, so that we can calculate its FPGA resource requirement
		cStep2nodeIDMap_tmp[minCurrStartTime].push_back(curr_node_id);
	}

	calculateFPGAResRequired(cStep2nodeIDMap_tmp, alap_FPGA_resources);

	//uint64_t max_level = *std::max_element(schedTime.begin(), schedTime.end());

	arrayName2NumPartTy arrayName2NumPart = alap_FPGA_resources.getArrayName2NumPartMap();
	arrayName2NumPartTy::iterator array2part_it = arrayName2NumPart.begin();
	arrayName2NumPartTy::iterator array2part_ie = arrayName2NumPart.end();

	alap_fadd_used = alap_FPGA_resources.getNumFadd();
	alap_fsub_used = alap_FPGA_resources.getNumFsub();
	alap_fmul_used = alap_FPGA_resources.getNumFmul();
	alap_fdiv_used = alap_FPGA_resources.getNumFdiv();

	std::cout << "DEBUG-INFO: [fpga_estimation-ALAP] IL of this loop level: " << maxCstep << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ALAP] Resource utilization: DSP=" << alap_FPGA_resources.getDSP() << "; ";
	std::cout << "FF=" << alap_FPGA_resources.getFF() << "; " << "LUT=" << alap_FPGA_resources.getLUT() << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ALAP] Number of Fadd=" << alap_fadd_used << "; ";
	std::cout << "Fsub=" << alap_fsub_used << "; " << "Fmul=" << alap_fmul_used << "; ";
	std::cout << "Fdiv=" << alap_fdiv_used << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ALAP] Number of array partitions of ";
	for (; array2part_it != array2part_ie; ++array2part_it) {
		std::cout << array2part_it->first << "=" << array2part_it->second << "; ";
	}
	std::cout << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-ALAP] Finished" << std::endl;
	std::cout << "##########" << std::endl;
	cStep2nodeIDMapTy().swap(cStep2nodeIDMap_tmp);
}

uint64_t BaseDatapath::rc_list_scheduling(Graph& graph_tmp, schedTimeTy& asap_time, schedTimeTy& alap_time, schedTimeTy& rsList_time, cPathNodeTy& cp_nodes, bool pipelining_enable) {
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] With FPGA Resource Constraints" << std::endl;

	completePartition();
	scratchpadPartition();
	pre_optimization();

	if (show_dddg_af_opt == true) {
		// After optimization, check the graph
		dumpGraph();
	}

	rsList_time.assign(numTotalNodes, 0);

	/// Read array partition configuration
	arrayName2arrayConfigTy arrayName2arrayConfig;
	initialize_array_partition_config(arrayName2arrayConfig);

	/// Initialize cFPGA_constraint class
	UsedFPGA_Res_With_Constraint cFPGA_constraint(arrayName2arrayConfig);
	/// Set upper bounds of different floating point units in cFPGA_constraint.
	if (disable_fp_unit_threshold == false) {
		cFPGA_constraint.set_floating_unit_threshold(alap_fadd_used, alap_fsub_used, alap_fmul_used, alap_fdiv_used);
	}
	/// Get limited floating point operation unit types
	limited_fop_unit_types.clear();
	limited_fop_unit_types = cFPGA_constraint.get_constrained_fop_unit_type();

	/// Get graph edge weight map
	EdgeWeightMap edge2weightMap = boost::get(boost::edge_weight, graph_tmp);

	unsigned scheduledNodeCount = 0;
	unsigned cycle_tick = 0;

	/*
	if (i < __BEGIN__ || i > __END__) {
	///FIXME: This section is only used to analyze specific code region, need to remove it later.
	///       Numbers here only work for convolution2d
	continue;
	}*/
	/// Initialize readyFadd_subList, readyFmulList, readyFdivList, readyMemList and readyOthersList;
	initialize_starting_queue(graph_tmp, alap_time, asap_time);

	while (scheduledNodeCount != totalConnectedNodes) {
		/// Determine readyOperationLists
		determineReadyAndSelectOperationList(graph_tmp, alap_time, rsList_time, cycle_tick, cFPGA_constraint, scheduledNodeCount);
		determineExecutingMap(graph_tmp, alap_time, rsList_time, cycle_tick, cFPGA_constraint, scheduledNodeCount);
		updateExecutingMapAndReadyQueue(graph_tmp, alap_time, rsList_time, cycle_tick, cFPGA_constraint, scheduledNodeCount);
		cycle_tick++;
	}

	uint64_t iteration_latency = 1;
	if (enable_extra_scalar == true) {
		iteration_latency = cycle_tick + 1;
	}
	else {
		iteration_latency = cycle_tick - 1;
	}

	rc_list_il = iteration_latency;

	ResII_mem = getMemResII_based_on_rcList_scheduling(rsList_time, iteration_latency, arrayName2arrayConfig, cFPGA_constraint);
	ResII_op = getOpResII_based_on_rcList_scheduling(cFPGA_constraint);
	unsigned ResII = (ResII_mem > ResII_op) ? ResII_mem : ResII_op;

	RecII = getRecII(pipelining_enable, asap_time, cp_nodes);
	
	max_II = (ResII > RecII) ? ResII : RecII;

	arrayName2NumPartTy arrayName2NumPart = cFPGA_constraint.getArrayName2NumPartMap();
	arrayName2NumPartTy::iterator array2part_it = arrayName2NumPart.begin();
	arrayName2NumPartTy::iterator array2part_ie = arrayName2NumPart.end();

	arrayName2NumPartTy arrayName2bram18k_used = cFPGA_constraint.getArrayName2Bram18k_used();
	arrayName2NumPartTy::iterator array2bramUsed_it = arrayName2bram18k_used.begin();
	arrayName2NumPartTy::iterator array2bramUsed_ie = arrayName2bram18k_used.end();

	arrayName2memeff.clear();
	arrayName2memeff = cFPGA_constraint.getMemoryEfficiency();
	arrayName2memEfficiencyTy::iterator array2memeff_it = arrayName2memeff.begin();
	arrayName2memEfficiencyTy::iterator array2memeff_ie = arrayName2memeff.end();

	/// Calculate shared-load numbers
	if (disable_shared_load_removal == true) {
		shared_loads = 0;
		removeSharedLoads();
	}
	else {
		if (enable_sharedLoadRemoval == false || enable_shared_load_removal == false) {
			shared_loads = 0;
			removeSharedLoads();
		}
	}

	dsp_used = cFPGA_constraint.getDSP();
	bram18k_used = cFPGA_constraint.getBRAM18K();
	ff_used = cFPGA_constraint.getFF();
	lut_used = cFPGA_constraint.getLUT();
	fadd_used = cFPGA_constraint.getNumFadd();
	fsub_used = cFPGA_constraint.getNumFsub();
	fmul_used = cFPGA_constraint.getNumFmul();
	fdiv_used = cFPGA_constraint.getNumFdiv();

	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] IL of this loop level: " << iteration_latency << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Resource utilization: DSP=" << dsp_used << "; ";
	std::cout << "BRAM18K=" << bram18k_used << "; FF=" << ff_used << "; LUT=" << lut_used << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Number of Fadd=" << fadd_used << "; ";
	std::cout << "Fsub=" << fsub_used << "; " << "Fmul=" << fmul_used << "; ";
	std::cout << "Fdiv=" << fdiv_used << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Number of array partitions of ";
	for (; array2part_it != array2part_ie; ++array2part_it) {
		std::cout << array2part_it->first << "=" << array2part_it->second << "; ";
	}
	std::cout << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Bram18k used of array ";
	for (; array2bramUsed_it != array2bramUsed_ie; ++array2bramUsed_it) {
		std::cout << array2bramUsed_it->first << "=" << array2bramUsed_it->second << "; ";
	}
	std::cout << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Memory efficiency of array ";
	for (; array2memeff_it != array2memeff_ie; ++array2memeff_it) {
		std::cout << array2memeff_it->first << "=" << array2memeff_it->second << "; ";
	}
	std::cout << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] II of this loop level: " << max_II << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] RecII of this loop level: " << RecII << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] ResII_mem of this loop level: " << ResII_mem << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] ResII_op of this loop level: " << ResII_op << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Potential number of shared loads in this loop level: " << shared_loads << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Potential number of repeated stores in this loop level: " << repeated_stores << std::endl;
	std::cout << "DEBUG-INFO: [fpga_estimation-rc_list_scheduling] Finished" << std::endl;
	std::cout << "##########" << std::endl;
	return iteration_latency;
}

unsigned BaseDatapath::getRecII(bool enable_lpPipelining, schedTimeTy& asap_time, cPathNodeTy& cp_nodes) {
	unsigned rec_ii;
	if (enable_lpPipelining == true) {
		int sub = (int)(IL_asap_ii - IL_asap);
		int sub_ori = sub;
		assert((sub >= 0) && "Error: IL_asap_ii is smaller than IL_asap!\n");
		
		if (sub > 1) {
			/// Registers between floating point units will be removed if enabling loop pipelining. To make prediction
			/// accurate, we need to substract rec_ii with number of floating point units consisting of the last "rec_ii"
			/// cycles in the critical path.
			cStep2nodeIDMapTy schedTime2nodesInCpathVecMap;
			for (unsigned i = 0; i < cp_nodes.size(); i++) {
				unsigned node_id = cp_nodes.at(i);
				unsigned sched_time = asap_time.at(node_id);
				schedTime2nodesInCpathVecMap[sched_time].push_back(node_id);
			}
			/*
			std::vector<unsigned> cp_last_nodes = schedTime2nodesInCpathVecMap.rbegin()->second;
			bool set_time_start_flag = false;
			for (int i = 0; i < cp_last_nodes.size(); i++) {
				unsigned node_id = cp_last_nodes.at(i);
				unsigned cpNode_opcode = microop.at(node_id);
				if (is_memory_op(cpNode_opcode)) {
					set_time_start_flag = true;
					break;
				}
			}

			unsigned time_start = 0;
			if (set_time_start_flag == true) {
				time_start = 1 + (IL_asap + 1) - sub;
			}
			else {
				time_start = 1 + IL_asap - sub;
			}*/

			//unsigned time_start = 1 + IL_asap - sub;
			unsigned time_start = IL_asap - sub;
			
			unsigned time_end = IL_asap;
			unsigned time_stamp = time_start;
			unsigned counter = 0;
			while (time_stamp <= time_end) {
				cStep2nodeIDMapTy::iterator found = schedTime2nodesInCpathVecMap.find(time_stamp);
				//assert(found != schedTime2nodesInCpathVecMap.end() && "Error: cannot find time stamp in schedTime2nodesInCpathVecMap!\n");
				while (found == schedTime2nodesInCpathVecMap.end()) {
					time_stamp++;
					assert(time_stamp < time_end && "Error: cannot find time stamp in schedTime2nodesInCpathVecMap and time_stamp > time_end!\n");
					found = schedTime2nodesInCpathVecMap.find(time_stamp);
				}
				std::vector<unsigned> nodesVec = found->second;
				unsigned max_latency = 1;
				bool count_flag = false;
				for (unsigned i = 0; i < nodesVec.size(); i++) {
					unsigned node_id = nodesVec.at(i);
					unsigned opcode = microop.at(node_id);
					if (!is_float_op(opcode)) {
						continue;
					}
					unsigned latency = fpga_node_latency(opcode);
					max_latency = (max_latency > latency) ? max_latency : latency;
					count_flag = true;
				}
				if (count_flag == true) {
					counter++;
				}
				time_stamp += max_latency;
			}

			sub -= counter;
			assert(sub>=0 && "Error: sub is less than 0!\n");
		}

		rec_ii = (sub_ori == 1) ? 1 : (unsigned) sub + 1;
	}
	else {
		rec_ii = 1;
	}
	return rec_ii;
}

void BaseDatapath::pre_optimization() {
	
	loopName2levelUnrollVecMapTy::iterator it_lpUnr = loopName2levelUnrollVecMap.find(target_loop_name);
	bool found_it_lpUnr = it_lpUnr != loopName2levelUnrollVecMap.end();
	assert(found_it_lpUnr && "Error: Can not find loop name inside loopName2levelUnrollVecMap!\n");
	std::vector<unsigned> unroll_factors = it_lpUnr->second;
	unsigned size = unroll_factors.size();
	unsigned unroll_factor_innermost = unroll_factors.back();
	std::string whole_lp_name = target_loop_name + "_" + std::to_string(size);
	wholeloopName2loopBoundMapTy::iterator it_lpBound = wholeloopName2loopBoundMap.find(whole_lp_name);
	bool found_it_lpBound = it_lpBound != wholeloopName2loopBoundMap.end();
	assert(found_it_lpBound && "Error: Can not find loop name inside wholeloopName2loopBoundMap!\n");
	unsigned innermost_bound = wholeloopName2loopBoundMap[whole_lp_name];
	enable_sharedLoadRemoval = false;
	if (unroll_factor_innermost == innermost_bound) {
		enable_sharedLoadRemoval = true;
	}

	if (enable_memory_disambiguation == true) {
		memoryAmbiguation();
	}

	if (disable_shared_load_removal == false) {
		if (enable_sharedLoadRemoval == true || enable_shared_load_removal == true) {
			removeSharedLoads();
		}
	}
	
	/*
	if (enable_store_buffer == true) {
		storeBuffer();
	}*/

	if (enable_repeated_store_removal == true) {
		removeRepeatedStores();
	}

	if (enable_tree_height_reduction_integer == true) {
		treeHeightReduction_integer();
	}
	
	if (enable_tree_height_reduction_float == true) {
		treeHeightReduction_float();
	}
	
}

void BaseDatapath::determineReadyAndSelectOperationList(Graph& graph_tmp, schedTimeTy& alapTime, schedTimeTy& rsListTime, unsigned& clock_tick, UsedFPGA_Res_With_Constraint& cFPGA_constraints, unsigned& scheduledCounters) {

	/// Make root nodes of the graph only starts their execution at the cycles of alap scheduling time.
	if (startingNodesList.size() != 0) {
		startingNodesList.sort(compare_ready_queue_basedon_alapTime);

		while (startingNodesList.size() != 0) {
			unsigned curr_node_id = startingNodesList.front().first;
			unsigned alap_time = startingNodesList.front().second;
			if (clock_tick == alap_time) {
				initialize_ready_queue(alap_time, curr_node_id);
				startingNodesList.pop_front();
			}
			else {
				// Exit the while loop
				break;
			}
		}

	}

	/// Schedule nodes with fpga latency 0 first until it is completely empty at this cycle_tick
	while (readyOthersList.size() != 0) {
		unsigned curr_node_id = readyOthersList.front().first;
		Vertex curr_vertex = nameToVertex[curr_node_id];
		readyOthersList.pop_front();
		rsListTime[curr_node_id] = clock_tick;
		scheduledCounters++;
		out_edge_iter outEg_it, outEg_ie;
		for (boost::tie(outEg_it, outEg_ie) = boost::out_edges(curr_vertex, graph_tmp); outEg_it != outEg_ie; ++outEg_it) {
			Vertex child_vertex = boost::target(*outEg_it, graph_tmp);
			unsigned child_node_id = vertexToName[child_vertex];
			numParents[child_node_id]--;
			if (numParents[child_node_id] == 0 && finalIsolated[child_node_id] != 1) {
				unsigned opcode = microop.at(child_node_id);
				updateReadyQueue(child_node_id, opcode, alapTime);
			}

		}
	}

	/// Select floating point operations from ready queues
	if (readyFaddList.size() != 0) {
		// Clear selectedQueue
		selectedFaddList.clear();

		// Sort the ready queue based on alap_schedTime
		readyFaddList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyFaddList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			bool success_or_not = cFPGA_constraints.try_to_occupy_one_fadd();
			// success_or_not: true --> it means that FPGA successfully allocates one fadd unit for it.
			if (success_or_not == true) {
				unsigned rd_node_id = readyFaddList.front().first;
				selectedFaddList.push_back(rd_node_id);
				readyFaddList.pop_front();
				rsListTime[rd_node_id] = clock_tick;
			}
			else {
				// We do not have enough FPGA resources already, nodes inside ready queue need to wait for resources to be released;
				// exit this loop
				break;
			}
		}
	} // End of determining readyFaddList and selectedFaddList

	if (readyFsubList.size() != 0) {
		// Clear selectedQueue
		selectedFsubList.clear();

		// Sort the ready queue based on alap_schedTime
		readyFsubList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyFsubList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			bool success_or_not = cFPGA_constraints.try_to_occupy_one_fsub();
			// success_or_not: true --> it means that FPGA successfully allocates one fsub unit for it.
			if (success_or_not == true) {
				unsigned rd_node_id = readyFsubList.front().first;
				selectedFsubList.push_back(rd_node_id);
				readyFsubList.pop_front();
				rsListTime[rd_node_id] = clock_tick;
			}
			else {
				// We do not have enough FPGA resources already, nodes inside ready queue need to wait for resources to be released;
				// exit this loop
				break;
			}
		}
	} // End of determining readyFsubList and selectedFsubList

	if (readyFmulList.size() != 0) {
		// Clear selectedQueue
		selectedFmulList.clear();

		// Sort the ready queue based on alap_schedTime
		readyFmulList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyFmulList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			bool success_or_not = cFPGA_constraints.try_to_occupy_one_fmul();
			// success_or_not: true --> it means that FPGA successfully allocates one fmul unit for it.
			if (success_or_not == true) {
				unsigned rd_node_id = readyFmulList.front().first;
				selectedFmulList.push_back(rd_node_id);
				readyFmulList.pop_front();
				rsListTime[rd_node_id] = clock_tick;
			}
			else {
				// We do not have enough FPGA resources already, nodes inside ready queue need to wait for resources to be released;
				// exit this loop
				break;
			}
		}
	} // End of determining readyFmulList and selectedFmulList

	if (readyFdivList.size() != 0) {
		// Clear selectedQueue
		selectedFdivList.clear();

		// Sort the ready queue based on alap_schedTime
		readyFdivList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyFdivList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			bool success_or_not = cFPGA_constraints.try_to_occupy_one_fdiv();
			// success_or_not: true --> it means that FPGA successfully allocates one fdiv unit for it.
			if (success_or_not == true) {
				unsigned rd_node_id = readyFdivList.front().first;
				selectedFdivList.push_back(rd_node_id);
				readyFdivList.pop_front();
				rsListTime[rd_node_id] = clock_tick;
			}
			else {
				// We do not have enough FPGA resources already, nodes inside ready queue need to wait for resources to be released;
				// exit this loop
				break;
			}
		}
	} // End of determining readyFdivList and selectedFdivList

	if (readyFcmpList.size() != 0) {
		/// In current implementation, we do not consider resources used inside Fcmp operations.
		// Clear selectedQueue
		selectedFcmpList.clear();

		// Sort the ready queue based on alap_schedTime
		readyFcmpList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyFcmpList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			unsigned rd_node_id = readyFcmpList.front().first;
			selectedFcmpList.push_back(rd_node_id);
			readyFcmpList.pop_front();
			rsListTime[rd_node_id] = clock_tick;
		}
	} // End of determining readyFcmpList and selectedFcmpList

	/// Select memory operation from ready queues
	if (readyLoadList.size() != 0) {
		// Clear selectedQueue
		selectedLoadList.clear();

		// Sort the ready queue based on alap_schedTime
		readyLoadList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyLoadList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			unsigned rd_node_id = readyLoadList.front().first;
			std::string array_name = baseAddress[rd_node_id].first;
			bool success_or_not = cFPGA_constraints.try_to_occupy_one_read_port(array_name);
			// success_or_not: true --> it means that FPGA successfully allocates one memory read port for it.
			if (success_or_not == true) {
				selectedLoadList.push_back(rd_node_id);
				readyLoadList.pop_front();
				rsListTime[rd_node_id] = clock_tick;
			}
			else {
				// We do not have enough FPGA resources already, nodes inside ready queue need to wait for resources to be released;
				// exit this loop
				break;
			}
		}
	} // End of determining readyLoadList and selectedLoadList

	if (readyStoreList.size() != 0) {
		// Clear selectedQueue
		selectedStoreList.clear();

		// Sort the ready queue based on alap_schedTime
		readyStoreList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyStoreList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			unsigned rd_node_id = readyStoreList.front().first;
			std::string array_name = baseAddress[rd_node_id].first;
			bool success_or_not = cFPGA_constraints.try_to_occupy_one_write_port(array_name);
			// success_or_not: true --> it means that FPGA successfully allocates one memory write port for it.
			if (success_or_not == true) {
				selectedStoreList.push_back(rd_node_id);
				readyStoreList.pop_front();
				rsListTime[rd_node_id] = clock_tick;
			}
			else {
				// We do not have enough FPGA resources already, nodes inside ready queue need to wait for resources to be released;
				// exit this loop
				break;
			}
		}
	} // End of determining readyStoreList and selectedStoreList

	/// Select integer operations from ready queues
	if (readyIntegerOpList.size() != 0) {
		/// In current implementation, we do not consider resources used inside integer operations.
		// Clear selectedQueue
		selectedIntegerOpList.clear();

		// Sort the ready queue based on alap_schedTime
		readyIntegerOpList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyIntegerOpList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			unsigned rd_node_id = readyIntegerOpList.front().first;
			selectedIntegerOpList.push_back(rd_node_id);
			readyIntegerOpList.pop_front();
			rsListTime[rd_node_id] = clock_tick;
		}
	} // End of determining readyIntegerOpList and selectedIntegerOpList

	/// Select other nodes with one latency cycle from ready queues
	if (readyCallOpList.size() != 0) {
		// Clear selectedQueue
		selectedCallOpList.clear();

		// Sort the ready queue based on alap_schedTime
		readyCallOpList.sort(compare_ready_queue_basedon_alapTime);
		unsigned size_readyQueue = readyCallOpList.size();
		for (unsigned i = 0; i < size_readyQueue; i++) {
			unsigned rd_node_id = readyCallOpList.front().first;
			selectedCallOpList.push_back(rd_node_id);
			readyCallOpList.pop_front();
			rsListTime[rd_node_id] = clock_tick;
		}
	} // End of determining readyCallOpList and selectedCallOpList
}

void BaseDatapath::determineExecutingMap(Graph& graph_tmp, schedTimeTy& alapTime, schedTimeTy& rsListTime, unsigned& clock_tick, UsedFPGA_Res_With_Constraint& cFPGA_constraints, unsigned& scheduledCounter) {
	/// Select floating point operations from selected queues and determine 
	/// executing maps

	if (selectedFaddList.size() != 0) {
		while (selectedFaddList.size() != 0) {
			unsigned st_node_id = selectedFaddList.front();
			unsigned opcode = microop.at(st_node_id);
			unsigned node_latency = fpga_node_latency(opcode);
			executingFaddMap.insert(std::make_pair(st_node_id, node_latency));
			selectedFaddList.pop_front();
			// Donot need to wait for this hardware unit to be released 
			// because of its pipeline design
			cFPGA_constraints.release_one_fadd();
		}
	} // End of determining executingFaddMap and update readyFaddList

	if (selectedFsubList.size() != 0) {
		while (selectedFsubList.size() != 0) {
			unsigned st_node_id = selectedFsubList.front();
			unsigned opcode = microop.at(st_node_id);
			unsigned node_latency = fpga_node_latency(opcode);
			executingFsubMap.insert(std::make_pair(st_node_id, node_latency));
			selectedFsubList.pop_front();
			// Donot need to wait for this hardware unit to be released 
			// because of its pipeline design
			cFPGA_constraints.release_one_fsub();
		}
	} // End of determining executingFsubMap and update readyFsubList

	if (selectedFmulList.size() != 0) {
		while (selectedFmulList.size() != 0) {
			unsigned st_node_id = selectedFmulList.front();
			unsigned opcode = microop.at(st_node_id);
			unsigned node_latency = fpga_node_latency(opcode);
			executingFmulMap.insert(std::make_pair(st_node_id, node_latency));
			selectedFmulList.pop_front();
			// Donot need to wait for this hardware unit to be released 
			// because of its pipeline design
			cFPGA_constraints.release_one_fmul();
		}
	} // End of determining executingFmulMap and update readyFmulList

	if (selectedFdivList.size() != 0) {
		while (selectedFdivList.size() != 0) {
			unsigned st_node_id = selectedFdivList.front();
			unsigned opcode = microop.at(st_node_id);
			unsigned node_latency = fpga_node_latency(opcode);
			executingFdivMap.insert(std::make_pair(st_node_id, node_latency));
			selectedFdivList.pop_front();
			// Donot need to wait for this hardware unit to be released 
			// because of its pipeline design
			cFPGA_constraints.release_one_fdiv();
		}
	} // End of determining executingFdivMap and update readyFdivList

	if (selectedFcmpList.size() != 0) {
		// Since FPGA node latency of Fcmp is only one cycle, we do not need
		// to put it inside executingFcmpMap, just delete all elements inside
		// selectedFcmpList
		while (selectedFcmpList.size() != 0) {
			unsigned st_node_id = selectedFcmpList.front();
			// rsListTime[st_node_id] has been update in determineReadyAndSelectOperationList(...)
			// function
			//rsListTime[st_node_id] = clock_tick;
			scheduledCounter++;
			selectedFcmpList.pop_front();

			// Update children nodes
			updateChildrenNodes(graph_tmp, st_node_id, alapTime);
		}
	} // End of determining executingFcmpMap and update readyFcmpList

	if (selectedLoadList.size() != 0) {
		// Since FPGA node latency of Load is only one cycle (Actually two cycles,
		// but the second cycle is just used to write into registers and its value
		// at second cycle is available already. Thus, we only consider one cycle 
		// FPGA load operation in this simulator), we do not need to put it inside 
		// executingLoadMap
		while (selectedLoadList.size() != 0) {
			unsigned st_node_id = selectedLoadList.front();
			std::string array_name = baseAddress[st_node_id].first;
			//rsListTime[st_node_id] = clock_tick;
			scheduledCounter++;
			// Donot need to wait for this hardware unit to be released 
			// because of its pipeline design
			cFPGA_constraints.release_one_read_port(array_name);

			selectedLoadList.pop_front();

			// Update children nodes
			updateChildrenNodes(graph_tmp, st_node_id, alapTime);
		}
	} // End of determining executingLoadMap and update readyLoadList

	if (selectedStoreList.size() != 0) {
		// Since FPGA node latency of write is only one cycle we do not need to put 
		// it inside executingStoreMap
		while (selectedStoreList.size() != 0) {
			unsigned st_node_id = selectedStoreList.front();
			std::string array_name = baseAddress[st_node_id].first;
			//rsListTime[st_node_id] = clock_tick;
			scheduledCounter++;
			// Donot need to wait for this hardware unit to be released 
			// because of its pipeline design
			cFPGA_constraints.release_one_write_port(array_name);

			selectedStoreList.pop_front();

			// Update children nodes
			updateChildrenNodes(graph_tmp, st_node_id, alapTime);
		}
	} // End of determining executingStoreMap and update readyStoreList

	if (selectedIntegerOpList.size() != 0) {
		// Since FPGA node latency of Call operation is only one cycle, we do not 
		// need to put it inside executingFcmpMap;
		while (selectedIntegerOpList.size() != 0) {
			unsigned st_node_id = selectedIntegerOpList.front();
			unsigned opcode = microop.at(st_node_id);
			unsigned node_latency = fpga_node_latency(opcode);
			if (node_latency == 1) {
				//rsListTime[st_node_id] = clock_tick;
				scheduledCounter++;
				// Update children nodes
				updateChildrenNodes(graph_tmp, st_node_id, alapTime);
			}
			else {
				executingIntegerOpMap.insert(std::make_pair(st_node_id, node_latency));
			}
			selectedIntegerOpList.pop_front();
		}

	} // End of determining executingIntegerOpMap and update readyIntegerOpList

	if (selectedCallOpList.size() != 0) {
		// Since FPGA node latency of Call operation is only one cycle, we do not 
		// need to put it inside executingFcmpMap, just delete all elements inside
		// selectedCallOpList
		while (selectedCallOpList.size() != 0) {
			unsigned st_node_id = selectedCallOpList.front();
			//rsListTime[st_node_id] = clock_tick;
			scheduledCounter++;
			selectedCallOpList.pop_front();

			// Update children nodes
			updateChildrenNodes(graph_tmp, st_node_id, alapTime);
		}
	} // End of determining executingCallOpMap and update readyCallOpList
}

void BaseDatapath::updateExecutingMapAndReadyQueue(Graph& graph_tmp, schedTimeTy& alapTime, schedTimeTy& rsListTime, unsigned& clock_tick, UsedFPGA_Res_With_Constraint& cFPGA_constraints, unsigned& scheduledCounter) {
	/// Deduct one cycle from all delay associated to node id inside executingOperationMap
	if (executingFaddMap.size() != 0) {
		std::vector<executingQueueMapTy::iterator> delete_vec;
		executingQueueMapTy::iterator it = executingFaddMap.begin();
		executingQueueMapTy::iterator ie = executingFaddMap.end();
		for (; it != ie; ++it) {
			unsigned em_node_id = it->first;
			// Decrease delay of node "em_node_id" by one cycle
			it->second--;
			unsigned delay = it->second;
			if (delay == 0) {
				// rsListTime[em_node_id] has been update in determineReadyAndSelectOperationList(...)
				// function
				//rsListTime[em_node_id] = clock_tick;
				scheduledCounter++;
				// The FPGA Fadd unit has been released already at determineExecutingMap(...) function,
				// no need to release again

				// Update children nodes
				updateChildrenNodes(graph_tmp, em_node_id, alapTime);

				// Need to delete this entry inside executingOperationMap
				delete_vec.push_back(it);
			}
		}

		for (int i = 0; i < delete_vec.size(); i++) {
			executingFaddMap.erase(delete_vec[i]);
		}
	} // End of executingFaddMap

	if (executingFsubMap.size() != 0) {
		std::vector<executingQueueMapTy::iterator> delete_vec;
		executingQueueMapTy::iterator it = executingFsubMap.begin();
		executingQueueMapTy::iterator ie = executingFsubMap.end();
		for (; it != ie; ++it) {
			unsigned em_node_id = it->first;
			// Decrease delay of node "em_node_id" by one cycle
			it->second--;
			unsigned delay = it->second;
			if (delay == 0) {
				//rsListTime[em_node_id] = clock_tick;
				scheduledCounter++;
				// The FPGA Fsub unit has been released already at determineExecutingMap(...) function,
				// no need to release again

				// Update children nodes
				updateChildrenNodes(graph_tmp, em_node_id, alapTime);

				// Need to delete this entry inside executingOperationMap
				delete_vec.push_back(it);
			}
		}

		for (int i = 0; i < delete_vec.size(); i++) {
			executingFsubMap.erase(delete_vec[i]);
		}
	} // End of executingFsubMap

	if (executingFmulMap.size() != 0) {
		std::vector<executingQueueMapTy::iterator> delete_vec;
		executingQueueMapTy::iterator it = executingFmulMap.begin();
		executingQueueMapTy::iterator ie = executingFmulMap.end();
		for (; it != ie; ++it) {
			unsigned em_node_id = it->first;
			// Decrease delay of node "em_node_id" by one cycle
			it->second--;
			unsigned delay = it->second;
			if (delay == 0) {
				//rsListTime[em_node_id] = clock_tick;
				scheduledCounter++;
				// The FPGA Fmul unit has been released already at determineExecutingMap(...) function,
				// no need to release again

				// Update children nodes
				updateChildrenNodes(graph_tmp, em_node_id, alapTime);

				// Need to delete this entry inside executingOperationMap
				delete_vec.push_back(it);
			}
		}

		for (int i = 0; i < delete_vec.size(); i++) {
			executingFmulMap.erase(delete_vec[i]);
		}
	} // End of executingFmulMap

	if (executingFdivMap.size() != 0) {
		std::vector<executingQueueMapTy::iterator> delete_vec;
		executingQueueMapTy::iterator it = executingFdivMap.begin();
		executingQueueMapTy::iterator ie = executingFdivMap.end();
		for (; it != ie; ++it) {
			unsigned em_node_id = it->first;
			// Decrease delay of node "em_node_id" by one cycle
			it->second--;
			unsigned delay = it->second;
			if (delay == 0) {
				//rsListTime[em_node_id] = clock_tick;
				scheduledCounter++;
				// The FPGA Fdiv unit has been released already at determineExecutingMap(...) function,
				// no need to release again

				// Update children nodes
				updateChildrenNodes(graph_tmp, em_node_id, alapTime);

				// Need to delete this entry inside executingOperationMap
				delete_vec.push_back(it);
			}
		}

		for (int i = 0; i < delete_vec.size(); i++) {
			executingFdivMap.erase(delete_vec[i]);
		}
	} // End of executingFdivMap

	if (executingIntegerOpMap.size() != 0) {
		std::vector<executingQueueMapTy::iterator> delete_vec;
		executingQueueMapTy::iterator it = executingIntegerOpMap.begin();
		executingQueueMapTy::iterator ie = executingIntegerOpMap.end();
		for (; it != ie; ++it) {
			unsigned em_node_id = it->first;
			// Decrease delay of node "em_node_id" by one cycle
			it->second--;
			unsigned delay = it->second;
			if (delay == 0) {
				//rsListTime[em_node_id] = clock_tick;
				scheduledCounter++;
				//FIXME: In the current implementation, we assume that FPGA has enough resources for 
				//			 Integer Operation units and no resource limitation. Later if we want to add
				//			 resource limitation for integer operations, we can add it inside 
				//			 determineExecutingMap(...) function and release FPGA resources below.

				// Update children nodes
				updateChildrenNodes(graph_tmp, em_node_id, alapTime);

				// Need to delete this entry inside executingOperationMap
				delete_vec.push_back(it);
			}
		}

		for (int i = 0; i < delete_vec.size(); i++) {
			executingIntegerOpMap.erase(delete_vec[i]);
		}
	} // End of executingIntegerOpMap
}

void BaseDatapath::calculateFPGAResRequired(cStep2nodeIDMapTy& cStepMap, UsedFPGA_ResClass& fpga_resources) {

	/// In the current implementation, we only consider floating point operations
	unsigned sizeFadd_subEngine = 0;
	unsigned sizeFmulEngine = 0;
	unsigned sizeFdivEngine = 0;

	/// Initialize memory bandwidth
	std::map<std::string, unsigned> arrayName2numReadPort;
	std::map<std::string, unsigned> arrayName2numWritePort;
	arrayInfoMapTy::iterator arrayName_it = arrayN2sizeWordsizeBytePr.begin();
	arrayInfoMapTy::iterator arrayName_ie = arrayN2sizeWordsizeBytePr.end();
	for (; arrayName_it != arrayName_ie; ++arrayName_it) {
		std::string arrayName = arrayName_it->first;
		arrayName2numReadPort.insert(std::make_pair(arrayName, 0));
		arrayName2numWritePort.insert(std::make_pair(arrayName, 0));
		fpga_resources.increaseArrayPart(arrayName);
	}

	cStep2nodeIDMapTy::iterator it, ie;
	for (it = cStepMap.begin(), ie = cStepMap.end(); it != ie; ++it) {
		unsigned fadd_subCount = 0;
		unsigned fmulCount = 0;
		unsigned fdivCount = 0;

		std::map<std::string, unsigned>::iterator rd_it = arrayName2numReadPort.begin();
		std::map<std::string, unsigned>::iterator rd_ie = arrayName2numReadPort.end();
		std::map<std::string, unsigned>::iterator wr_it = arrayName2numWritePort.begin();
		std::map<std::string, unsigned>::iterator wr_ie = arrayName2numWritePort.end();
		for (; rd_it != rd_ie; ++rd_it) {
			rd_it->second = 0;
		}

		for (; wr_it != wr_ie; ++wr_it) {
			wr_it->second = 0;
		}

		unsigned size = it->second.size();
		for (unsigned i = 0; i < size; i++) {
			unsigned node_id = it->second.at(i);
			unsigned opcode = microop.at(node_id);
			if (is_fadd_op(opcode) || is_fsub_op(opcode)) {
				fadd_subCount++;
				if (fadd_subCount > sizeFadd_subEngine) {
					sizeFadd_subEngine++;
					switch (opcode) {
					case LLVM_IR_FAdd:
						fpga_resources.increaseFadd();
						break;
					case LLVM_IR_FSub:
						fpga_resources.increaseFsub();
						break;
					default:
						// Do nothing here
						assert(false && "Shouldn't jump to this branch!\n");
						break;
					}
				}
			}

			if (is_fmul_op(opcode)) {
				fmulCount++;
				if (fmulCount > sizeFmulEngine) {
					sizeFmulEngine++;
					fpga_resources.increaseFmul();
				}
			}

			if (is_fdiv_op(opcode)) {
				fdivCount++;
				if (fdivCount > sizeFdivEngine) {
					sizeFdivEngine++;
					fpga_resources.increaseFdiv();
				}
			}

			// Calculate memory partition
			if (is_load_op(opcode)) {
				///FIXME: This part will be failed if we call scratchpadPartition() function before this function
				std::string array_name = baseAddress[node_id].first;
				arrayName2numReadPort[array_name]++;
				unsigned RdportNum = READ_PORT_PER_PARTITION * fpga_resources.getArrayPartitionNum(array_name);
				if (arrayName2numReadPort[array_name] > RdportNum) {
					fpga_resources.increaseArrayPart(array_name);
				}
			}
			
			if (is_store_op(opcode)) {
				///FIXME: This part will be failed if we call scratchpadPartition() function before this function
				std::string array_name = baseAddress[node_id].first;
				arrayName2numWritePort[array_name]++;
				unsigned WrportNum = WRITE_PORT_PER_PARTITION * fpga_resources.getArrayPartitionNum(array_name);
				if (arrayName2numWritePort[array_name] > WrportNum) {
					fpga_resources.increaseArrayPart(array_name);
				}
			}

		}
	}

}

void BaseDatapath::critical_path_extraction(schedTimeTy& asap_schedTime, schedTimeTy& alap_schedTime, cPathNodeTy& cpNodeVec) {
	cpNodeVec.clear();
	unsigned size_asap = asap_schedTime.size();
	unsigned size_alap = alap_schedTime.size();
	assert( (size_asap==size_alap) && "Error: Mismatched size - asap_schedTime & alap_schedTime\n" );

	for (unsigned i = 0; i < numTotalNodes; i++) {
		/*
		if (i < __BEGIN__ || i > __END__) {
			///FIXME: This section is only used to analyze specific code region, need to remove it later.
			///       Numbers here only work for convolution2d
			continue;
		}*/
		unsigned asapTime = asap_schedTime[i];
		unsigned alapTime = alap_schedTime[i];
		if (asapTime == alapTime) {
			cpNodeVec.push_back(i);
		}
	}
}

void BaseDatapath::initialize_array_partition_config(arrayName2arrayConfigTy& arrayName2Config) {
	std::unordered_map<std::string, unsigned> comp_part_config;
	readCompletePartitionConfig(comp_part_config);
	std::unordered_map<std::string, partitionEntry> part_config;
	readPartitionConfig(part_config);
	arrayInfoMapTy::iterator arrayN_it = arrayN2sizeWordsizeBytePr.begin();
	arrayInfoMapTy::iterator arrayN_ie = arrayN2sizeWordsizeBytePr.end();
	for (; arrayN_it != arrayN_ie; ++arrayN_it) {
		std::string arrayName = arrayN_it->first;
		uint64_t size_in_byte = arrayN_it->second.first;
		unsigned wordsize_in_byte = arrayN_it->second.second;
		std::unordered_map<std::string, unsigned>::iterator it_comp_part = comp_part_config.find(arrayName);
		std::unordered_map<std::string, partitionEntry>::iterator it_part = part_config.find(arrayName);
		if (it_comp_part != comp_part_config.end()) {
			arrayName2Config.insert(std::make_pair(arrayName, std::make_pair(0, std::make_pair(comp_part_config[arrayName], wordsize_in_byte))));
		}
		else if (it_part != part_config.end()) {
			arrayName2Config.insert(std::make_pair(arrayName, std::make_pair(part_config[arrayName].part_factor, std::make_pair(part_config[arrayName].array_size, part_config[arrayName].wordsize))));
		}
		else {
			// No partition configuration
			arrayName2Config.insert(std::make_pair(arrayName, std::make_pair(1, std::make_pair(size_in_byte, wordsize_in_byte))));
		}
	}
}

void BaseDatapath::initialize_starting_queue(Graph& graph_tmp, schedTimeTy& alapTime, schedTimeTy& asapTime) {
	numParents.assign(numTotalNodes, 0);
	finalIsolated.assign(numTotalNodes, 1);
	totalConnectedNodes = 0;

	startingNodesList.clear();

	readyFaddList.clear();
	readyFsubList.clear();
	readyFmulList.clear();
	readyFdivList.clear();
	readyFcmpList.clear();
	readyLoadList.clear();
	readyStoreList.clear();
	readyIntegerOpList.clear();
	readyCallOpList.clear();
	readyOthersList.clear();

	selectedFaddList.clear();
	selectedFsubList.clear();
	selectedFmulList.clear();
	selectedFdivList.clear();
	selectedFcmpList.clear();
	selectedLoadList.clear();
	selectedStoreList.clear();
	selectedIntegerOpList.clear();
	selectedCallOpList.clear();

	// Executing queues are only used when nodes with latency larger than 1
	executingFaddMap.clear();
	executingFsubMap.clear();
	executingFmulMap.clear();
	executingFdivMap.clear();
	executingIntegerOpMap.clear();

	vertex_iter it, ie;
	for (boost::tie(it, ie) = boost::vertices(graph_tmp); it != ie; ++it) {
		unsigned curr_node_id = vertexToName[*it];
		/*
		if (curr_node_id > __END__) {
			///FIXME: This section is only used to analyze specific code region, need to remove it later.
			///       Numbers here only work for convolution2d
			continue;
		}*/
		unsigned degree = boost::degree(*it, graph_tmp);
		if (degree > 0) {
			unsigned inDegree = boost::in_degree(*it, graph_tmp);
			numParents[curr_node_id] = inDegree;
			totalConnectedNodes++;
			finalIsolated[curr_node_id] = 0;
			unsigned asap_time = asapTime[curr_node_id];
			unsigned alap_time = alapTime[curr_node_id];
			if (inDegree == 0) {
				startingNodesList.push_back(std::make_pair(curr_node_id, alap_time));
				/*
				if (alap_time == asap_time) {
					startingNodesList.push_back(std::make_pair(curr_node_id, alap_time));
				}
				else {
					unsigned opcode = microop.at(curr_node_id);
					switch (opcode) {
					case LLVM_IR_FAdd:
						readyFaddList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_FSub:
						readyFsubList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_FMul:
						readyFmulList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_FDiv:
						readyFdivList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_FCmp:
						readyFcmpList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_Load:
						readyLoadList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_Store:
						readyStoreList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_Add:
					case LLVM_IR_Sub:
					case LLVM_IR_Mul:
					case LLVM_IR_UDiv:
					case LLVM_IR_SDiv:
						readyIntegerOpList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					case LLVM_IR_Call:
						readyCallOpList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					default:
						readyOthersList.push_back(std::make_pair(curr_node_id, alap_time));
						break;
					}
				}
				*/
			}
		}
	}
}

void BaseDatapath::initialize_ready_queue(unsigned alap_time_tick, unsigned cur_node_id) {
	unsigned opcode = microop.at(cur_node_id);
	switch (opcode) {
	case LLVM_IR_FAdd:
		readyFaddList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_FSub:
		readyFsubList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_FMul:
		readyFmulList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_FDiv:
		readyFdivList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_FCmp:
		readyFcmpList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_Load:
		readyLoadList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_Store:
		readyStoreList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_Add:
	case LLVM_IR_Sub:
	case LLVM_IR_Mul:
	case LLVM_IR_UDiv:
	case LLVM_IR_SDiv:
		readyIntegerOpList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	case LLVM_IR_Call:
		readyCallOpList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	default:
		readyOthersList.push_back(std::make_pair(cur_node_id, alap_time_tick));
		break;
	}

}

void BaseDatapath::updateReadyQueue(unsigned child_nodeID, unsigned opcode, schedTimeTy& alapTime) {
	switch (opcode) {
	case LLVM_IR_FAdd:
		readyFaddList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_FSub:
		readyFsubList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_FMul:
		readyFmulList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_FDiv:
		readyFdivList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_FCmp:
		readyFcmpList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_Load:
		readyLoadList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_Store:
		readyStoreList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_Add:
	case LLVM_IR_Sub:
	case LLVM_IR_Mul:
	case LLVM_IR_UDiv:
	case LLVM_IR_SDiv:
		readyIntegerOpList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	case LLVM_IR_Call:
		readyCallOpList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	default:
		readyOthersList.push_back(std::make_pair(child_nodeID, alapTime[child_nodeID]));
		break;
	}
}

void BaseDatapath::updateChildrenNodes(Graph& graph_tmp, unsigned curr_node_Id, schedTimeTy& alap_time) {
	Vertex st_vertex = nameToVertex[curr_node_Id];
	out_edge_iter outEg_it, outEg_ie;
	for (boost::tie(outEg_it, outEg_ie) = boost::out_edges(st_vertex, graph_tmp); outEg_it != outEg_ie; ++outEg_it) {
		Vertex child_vertex = boost::target(*outEg_it, graph_tmp);
		unsigned child_node_id = vertexToName[child_vertex];
		unsigned opcode = microop.at(child_node_id);
		numParents[child_node_id]--;
		if (numParents[child_node_id] == 0 && finalIsolated[child_node_id] != 1) {
			updateReadyQueue(child_node_id, opcode, alap_time);
		}
	}
}

unsigned BaseDatapath::getMemResII_based_on_rcList_scheduling(schedTimeTy& rsList_time, unsigned max_level, arrayName2arrayConfigTy& arrayName2arrayConfig, UsedFPGA_Res_With_Constraint& cFPGA_constraints) {
	unsigned max_memResII = 1;
	limited_mem_name.assign("NO_LIMITED_ARRAY");

	std::map<std::string, int64_t> arrayName2prevRStime_read;
	std::map<std::string, int64_t> arrayName2prevRStime_write;
	//std::map<std::string, uint64_t> arrayName2numRead;
	//std::map<std::string, uint64_t> arrayName2numWrite;
	std::map<std::string, std::vector<int64_t> > arrayName2Iss_read;
	std::map<std::string, std::vector<int64_t> > arrayName2Iss_write;
	std::map<std::string, double> arrayName2ResII;
	arrayName2arrayConfigTy::iterator it = arrayName2arrayConfig.begin();
	arrayName2arrayConfigTy::iterator ie = arrayName2arrayConfig.end();
	for (; it != ie; ++it) {
		std::string array_name = it->first;
		unsigned part_factor = it->second.first;
		if (part_factor > 1) {
			for (unsigned i = 0; i < part_factor; i++) {
				std::string new_array_name = array_name + "-" + std::to_string(i);
				arrayName2prevRStime_read.insert(std::make_pair(new_array_name, -1));
				arrayName2prevRStime_write.insert(std::make_pair(new_array_name, -1));
				arrayName2ResII.insert(std::make_pair(new_array_name, 0.0));
			}
		}
		else {
			arrayName2prevRStime_read.insert(std::make_pair(array_name, -1));
			arrayName2prevRStime_write.insert(std::make_pair(array_name, -1));
			arrayName2ResII.insert(std::make_pair(array_name, 0.0));
		}
	}

	std::map<uint64_t, std::vector<unsigned> > cStep2nodeIDs;
	for (unsigned i = 0; i < numTotalNodes; i++) {
		unsigned node_id = i;
		uint64_t sched_time = rsList_time[node_id];
		cStep2nodeIDs[sched_time].push_back(node_id);
	}

	std::map<uint64_t, std::vector<unsigned> >::iterator it_cstep = cStep2nodeIDs.begin();
	std::map<uint64_t, std::vector<unsigned> >::iterator ie_cstep = cStep2nodeIDs.end();
	for (; it_cstep != ie_cstep; ++it_cstep) {
		uint64_t cstep_val = it_cstep->first;
		unsigned size = it_cstep->second.size();
		for (unsigned i = 0; i < size; i++) {
			unsigned node_id = it_cstep->second.at(i);
			unsigned opcode = microop.at(node_id);
			if (opcode != LLVM_IR_Load && opcode != LLVM_IR_Store) {
				continue;
			}

			std::string arrayName = baseAddress[node_id].first;
			std::string old_array_name;

			std::size_t found = arrayName.find("-");
			if (found != string::npos) {
				old_array_name = arrayName.substr(0, found);
			}
			else {
				old_array_name = arrayName;
			}
			

			if (opcode == LLVM_IR_Load) {
				if (arrayName2arrayConfig[old_array_name].first == 0) {
					/// Complete partitioning, no need to analyze its RecII
					continue;
				}

				arrayName2numRead[arrayName]++;

				if (arrayName2prevRStime_read[arrayName] == -1) {
					arrayName2prevRStime_read[arrayName] = cstep_val;
				}

				if (arrayName2prevRStime_read[arrayName] != cstep_val) {
					arrayName2Iss_read[arrayName].push_back(cstep_val - arrayName2prevRStime_read[arrayName]);
					arrayName2prevRStime_read[arrayName] = cstep_val;
				}
			}

			if (opcode == LLVM_IR_Store) {
				if (arrayName2arrayConfig[old_array_name].first == 0) {
					/// Complete partitioning, no need to analyze its RecII
					continue;
				}

				arrayName2numWrite[arrayName]++;

				if (arrayName2prevRStime_write[arrayName] == -1) {
					arrayName2prevRStime_write[arrayName] = cstep_val;
				}

				if (arrayName2prevRStime_write[arrayName] != cstep_val) {
					arrayName2Iss_write[arrayName].push_back(cstep_val - arrayName2prevRStime_write[arrayName]);
					arrayName2prevRStime_write[arrayName] = cstep_val;
				}

			}

		}
	}

	/// Analyze ResII for memory read
	std::map<std::string, double>::iterator it_arrayN = arrayName2ResII.begin();
	std::map<std::string, double>::iterator ie_arrayN = arrayName2ResII.end();
	for (; it_arrayN != ie_arrayN; ++it_arrayN) {

		double tmp_read_II = 0.0;
		double tmp_write_II = 0.0;
		std::string array_name = it_arrayN->first;

		/// Analyze ResII for memory read
		std::map<std::string, uint64_t>::iterator it_numRead = arrayName2numRead.find(array_name);
		if (it_numRead != arrayName2numRead.end()) {
			uint64_t numRead = it_numRead->second;
			std::map<std::string, std::vector<int64_t> >::iterator it_issRead = arrayName2Iss_read.find(array_name);
			if (it_issRead != arrayName2Iss_read.end()) {
				uint64_t Iss_read = *std::min_element(arrayName2Iss_read[array_name].begin(), arrayName2Iss_read[array_name].end());
				//unsigned total_ports = cFPGA_constraints.getArrayPartitionNum(array_name) * READ_PORT_PER_PARTITION;
				/// In this version, we calculate tmp_read_II based on each partition array, no need to consider how many partitions for an array.
				unsigned total_ports = READ_PORT_PER_PARTITION;
				//double tmp_rd_II = std::ceil(double(numRead*Iss_read) / double(total_ports));
				//errs() << "tmp_rd_II" << tmp_rd_II << "\n";
				tmp_read_II = std::ceil(double(numRead) / double(total_ports));
			}
		}

		/// Analyze ResII for memory write
		std::map<std::string, uint64_t>::iterator it_numWrite = arrayName2numWrite.find(array_name);
		if (it_numWrite != arrayName2numWrite.end()) {
			uint64_t numWrite = it_numWrite->second;
			unsigned total_ports = cFPGA_constraints.getArrayPartNameWritePortNum(array_name);
			std::map<std::string, std::vector<int64_t> >::iterator it_issWrite = arrayName2Iss_write.find(array_name);
			if (it_issWrite != arrayName2Iss_write.end()) {
				uint64_t Iss_write = *std::min_element(arrayName2Iss_write[array_name].begin(), arrayName2Iss_write[array_name].end());
				//unsigned total_ports = cFPGA_constraints.getArrayPartitionNum(array_name) * cFPGA_constraints.getArrayWritePortNumPerPart(array_name);
				/// In this version, we calculate tmp_read_II based on each partition array, no need to consider how many partitions for an array.
				//unsigned total_ports = cFPGA_constraints.getArrayWritePortNumPerPart(array_name);
				tmp_write_II = std::ceil(double(numWrite*Iss_write) / double(total_ports));
			}
			else {
				tmp_write_II = std::ceil(double(numWrite) / double(total_ports));
			}
		}

		arrayName2ResII[array_name] = (tmp_read_II > tmp_write_II) ? tmp_read_II : tmp_write_II;

	}
	std::map<std::string, double>::iterator max_it = std::max_element(arrayName2ResII.begin(), arrayName2ResII.end(), compare_arrayName2resII);
	std::string array_name_max = max_it->first;
	double mem_ResII = max_it->second;

	if (max_memResII < (unsigned)mem_ResII) {
		max_memResII = mem_ResII;
		limited_mem_name.assign(array_name_max);
	}
	return max_memResII;
}

unsigned BaseDatapath::getOpResII_based_on_rcList_scheduling(UsedFPGA_Res_With_Constraint& cFPGA_constraints) {
	unsigned max_opResII = 1;
	limited_op_name.assign("NO_LIMITED_FOP_UNIT");

	unsigned alap_fadd_num = alap_FPGA_resources.getNumFadd();
	unsigned alap_fsub_num = alap_FPGA_resources.getNumFsub();
	unsigned alap_fmul_num = alap_FPGA_resources.getNumFmul();
	unsigned alap_fdiv_num = alap_FPGA_resources.getNumFdiv();

	unsigned used_fadd_num = cFPGA_constraints.getNumFadd();
	unsigned used_fsub_num = cFPGA_constraints.getNumFsub();
	unsigned used_fmul_num = cFPGA_constraints.getNumFmul();
	unsigned used_fdiv_num = cFPGA_constraints.getNumFdiv();
	
	float opResII = 0.0;
	float tmp_opResII = 0.0;
	if (used_fadd_num != 0) {
		opResII = std::ceil(float(alap_fadd_num) / float(used_fadd_num));
		limited_op_name.assign("Fadd");
	}

	if (used_fsub_num != 0) {
		tmp_opResII = std::ceil(float(alap_fsub_num) / float(used_fsub_num));
		if (tmp_opResII > opResII) {
			opResII = tmp_opResII;
			limited_op_name.assign("Fsub");
		}
		//opResII = (opResII > tmp_opResII) ? opResII : tmp_opResII;
	}

	if (used_fmul_num != 0) {
		tmp_opResII = std::ceil(float(alap_fmul_num) / float(used_fmul_num));
		if (tmp_opResII > opResII) {
			opResII = tmp_opResII;
			limited_op_name.assign("Fmul");
		}
		//opResII = (opResII > tmp_opResII) ? opResII : tmp_opResII;
	}

	if (used_fdiv_num != 0) {
		tmp_opResII = std::ceil(float(alap_fdiv_num) / float(used_fdiv_num));
		if (tmp_opResII > opResII) {
			opResII = tmp_opResII;
			limited_op_name.assign("Fdiv");
		}
		//opResII = (opResII > tmp_opResII) ? opResII : tmp_opResII;
	}

	return ((max_opResII > (unsigned)opResII) ? max_opResII : (unsigned)opResII);
}

LoopLatencyTy BaseDatapath::getLoopTotalLatency(uint64_t iterationLat, unsigned max_ii, std::string lp_name, unsigned lp_level, bool pipelining_or_not) {
	LoopLatencyTy lp_lat = {0, 0, false};
	unsigned extra_cost = 0;

	/// Analyze maximum number of read/write of memory banks of arrays
	getTargetArrayName2maxReadWrite();

	//wholeloopName2loopBoundMap;
	//loopName2levelUnrollVecMap;
	loopName2levelUnrollVecMapTy::iterator it_lpUnr = loopName2levelUnrollVecMap.find(lp_name);
	assert( (it_lpUnr!=loopName2levelUnrollVecMap.end()) && "Error: can not find loop name in loopName2levelUnrollVecMap!\n" );
	std::string whole_lpName = lp_name + "_" + std::to_string(lp_level);
	std::vector<unsigned> target_lp_unr = it_lpUnr->second;
	wholeloopName2loopBoundMapTy::iterator it_wlpBound = wholeloopName2loopBoundMap.find(whole_lpName);
	assert( (it_wlpBound != wholeloopName2loopBoundMap.end()) && "Error: can not find loop name in wholeloopName2loopBoundMap!\n");

	unsigned size_lpLevel = lp_level;
	/// FIXME: In current implementation, iterationLat is the iteration latency of the innermost level loop. We can extend it to
	/// more general case in the future.
	for (int i = (int) (size_lpLevel-1); i >= 0; i--) {
		unsigned unroll_factor = target_lp_unr.at(i);
		std::string cur_whole_lpName = lp_name + "_" + std::to_string(i+1);
		unsigned cur_bound = wholeloopName2loopBoundMap[cur_whole_lpName];
		assert( (cur_bound != 0) && "Error: cur_bound is equal to 0!\n" );
		if (i == size_lpLevel-1) {
			// Extra cost for entering and exiting a loop
			unsigned extra_cost = 3;
			lp_lat.no_pipeline = iterationLat * (cur_bound / unroll_factor) + extra_cost;
			if (!pipelining_or_not) {
				calculateArrayName2maxReadWrite(cur_bound / unroll_factor);
			}
		}
		else {
			lp_lat.no_pipeline *= (cur_bound / unroll_factor);
			if (!pipelining_or_not) {
				calculateArrayName2maxReadWrite(cur_bound / unroll_factor);
			}
		}

		if (i != 0) {
			//lp_lat.no_pipeline += extra_cost;
			unsigned unroll_factor_upperLoop = target_lp_unr.at(i - 1);
			lp_lat.no_pipeline *= unroll_factor_upperLoop;
			if (!pipelining_or_not) {
				calculateArrayName2maxReadWrite(unroll_factor_upperLoop);
			}
			if (extra_cost != 0) {
				lp_lat.no_pipeline += unroll_factor_upperLoop + 1;
			}
		}
	}

	if (pipelining_or_not) {
		lp_lat.enable_pipeline = true;

		std::string cur_whole_lpName = lp_name + "_" + std::to_string(lp_level);
		unsigned cur_bound = wholeloopName2loopBoundMap[cur_whole_lpName];
		unsigned unroll_factor = target_lp_unr.at(lp_level - 1);
		unsigned cur_iterations = cur_bound / unroll_factor;
		unsigned total_iteration = 1;
		
		int start_level = 0;
		int i = 0;
		for (i = (lp_level - 2); i >= 0; i--) {
			std::string tmp_whole_lpName = lp_name + "_" + std::to_string(i + 1);
			bool perfectOrNot = wholeloopName2perfectOrNotMap[tmp_whole_lpName];
			unsigned tmp_bound = wholeloopName2loopBoundMap[tmp_whole_lpName];
			if (perfectOrNot == true) {
				cur_iterations *= tmp_bound;
			}
			else {
				start_level = i;
				break;
			}
		}
		start_level = i;

		for (i = start_level; i >= 0; i--) {
			std::string tmp_whole_lpName = lp_name + "_" + std::to_string(i + 1);
			unsigned tmp_bound = wholeloopName2loopBoundMap[tmp_whole_lpName];
			total_iteration *= tmp_bound;
		}

		/*
		for (i = (lp_level - 2); i >= 0; i--) {
			std::string tmp_whole_lpName = lp_name + "_" + std::to_string(i+1);
			unsigned tmp_bound = wholeloopName2loopBoundMap[tmp_whole_lpName];
			total_iteration *= tmp_bound;
		}
		total_iteration -= 1;
		*/

		//lp_lat.with_pipeline = max_ii * total_iteration + iterationLat;
		lp_lat.with_pipeline = (max_ii * (cur_iterations - 1) + iterationLat + 2) * total_iteration;
		calculateArrayName2maxReadWrite(cur_iterations * total_iteration);
	}

	writeLogOfArrayName2maxReadWrite();
	return lp_lat;
}

void BaseDatapath::getTargetArrayName2maxReadWrite() {

	arrayName2maxMemOpNumTy arrayName2numMemOp;
	arrayName2maxMemOpNum_subtrace.clear();
	arrayName2aveLoadAccessPerBank_subtrace.clear();
	arrayName2aveStoreAccessPerBank_subtrace.clear();

	/// Read array partition configuration
	arrayName2arrayConfigTy arrayName2arrayConfig;
	initialize_array_partition_config(arrayName2arrayConfig);

	arrayName2arrayConfigTy::iterator it_conf, ie_conf;
	it_conf = arrayName2arrayConfig.begin();
	ie_conf = arrayName2arrayConfig.end();
	for (; it_conf != ie_conf; ++it_conf) {
		std::string array_name = it_conf->first;
		arrayName2aveLoadAccessPerBank_subtrace.insert(std::make_pair(array_name, 0.0));
		arrayName2aveStoreAccessPerBank_subtrace.insert(std::make_pair(array_name, 0.0));
	}

	/// Calculate number of memory operations of each memory bank of array A
	arrayName2maxMemOpNumTy::iterator r_it, r_ie;
	r_it = arrayName2numRead.begin();
	r_ie = arrayName2numRead.end();
	for (; r_it != r_ie; ++r_it) {
		arrayName2numMemOp.insert(std::make_pair(r_it->first, r_it->second));
		
		std::string array_name = r_it->first;
		std::string old_array_name;
		std::size_t found = array_name.find("-");
		if (found != string::npos) {
			old_array_name = array_name.substr(0, found);
		}
		else {
			old_array_name = array_name;
		}
		arrayName2MemBankStatisTy::iterator found_ave = arrayName2aveLoadAccessPerBank_subtrace.find(old_array_name);
		if (found_ave != arrayName2aveLoadAccessPerBank_subtrace.end()) {
			uint64_t num_read_access = r_it->second;
			found_ave->second += (float)num_read_access;
		}
		else {
			assert(false && "DEBUG-INFO: [getTargetArrayName2maxReadWrite] Can not find array name in arrayName2aveLoadAccessPerBank_subtrace!\n");
		}

	}

	arrayName2maxMemOpNumTy::iterator w_it, w_ie;
	w_it = arrayName2numWrite.begin();
	w_ie = arrayName2numWrite.end();
	for (; w_it != w_ie; ++w_it) {
		arrayName2maxMemOpNumTy::iterator found_it = arrayName2numMemOp.find(w_it->first);
		if (found_it != arrayName2numMemOp.end()) {
			uint64_t numMemOp = found_it->second;
			found_it->second = numMemOp + w_it->second;
		}
		else {
			arrayName2numMemOp.insert(std::make_pair(w_it->first, w_it->second));
		}

		std::string array_name = w_it->first;
		std::string old_array_name;
		std::size_t found = array_name.find("-");
		if (found != string::npos) {
			old_array_name = array_name.substr(0, found);
		}
		else {
			old_array_name = array_name;
		}
		arrayName2MemBankStatisTy::iterator found_ave = arrayName2aveStoreAccessPerBank_subtrace.find(old_array_name);
		if (found_ave != arrayName2aveStoreAccessPerBank_subtrace.end()) {
			uint64_t num_write_access = w_it->second;
			found_ave->second += (float)num_write_access;
		}
		else {
			assert(false && "DEBUG-INFO: [getTargetArrayName2maxReadWrite] Can not find array name in arrayName2aveLoadAccessPerBank_subtrace!\n");
		}

	}

	/// Get average memory load access per memory bank
	arrayName2MemBankStatisTy::iterator it_ld_ave, ie_ld_ave;
	it_ld_ave = arrayName2aveLoadAccessPerBank_subtrace.begin();
	ie_ld_ave = arrayName2aveLoadAccessPerBank_subtrace.end();
	for (; it_ld_ave != ie_ld_ave; ++it_ld_ave) {
		std::string array_name = it_ld_ave->first;
		unsigned part_num = arrayName2arrayConfig[array_name].first;
		if (part_num == 0) {
			/// Complete partitioning, no need to worry memory bank conflict problems
			it_ld_ave->second = 0.0;
		}
		else {
			float ave_ld_acc = it_ld_ave->second;
			it_ld_ave->second = ave_ld_acc / (float)part_num;
		}
	}

	/// Get average memory write access per memory bank
	arrayName2MemBankStatisTy::iterator it_st_ave, ie_st_ave;
	it_st_ave = arrayName2aveStoreAccessPerBank_subtrace.begin();
	ie_st_ave = arrayName2aveStoreAccessPerBank_subtrace.end();
	for (; it_st_ave != ie_st_ave; ++it_st_ave) {
		std::string array_name = it_st_ave->first;
		unsigned part_num = arrayName2arrayConfig[array_name].first;
		if (part_num == 0) {
			/// Complete partitioning, no need to worry memory bank conflict problems
			it_st_ave->second = 0.0;
		}
		else {
			float ave_st_acc = it_st_ave->second;
			it_st_ave->second = ave_st_acc / (float)part_num;
		}
	}

	std::map<std::string, uint64_t>::iterator it, ie;
	it = arrayName2numMemOp.begin();
	ie = arrayName2numMemOp.end();
	for (; it != ie; ++it) {
		std::string array_name = it->first;
		std::size_t found = array_name.find("-");
		uint64_t numMemOp = it->second;
		if (found != std::string::npos) {
			std::string original_name = array_name.substr(0, found);
			arrayName2maxMemOpNumTy::iterator found_it = arrayName2maxMemOpNum.find(original_name);
			if (found_it == arrayName2maxMemOpNum.end()) {
				arrayName2maxMemOpNum.insert(std::make_pair(original_name, numMemOp));
				arrayName2maxMemOpNum_subtrace.insert(std::make_pair(original_name, numMemOp));
			}
			else {
				uint64_t max_memOp = found_it->second;
				found_it->second = numMemOp > max_memOp ? numMemOp : max_memOp;
				///FIXME: If we want to record which memory bank, then we need to add codes below. Currently,
				///       we do not consider which memory bank has the maximum number of read operation.
			}
		}
		else {
			arrayName2maxMemOpNum.insert(std::make_pair(array_name, numMemOp));
			arrayName2maxMemOpNum_subtrace.insert(std::make_pair(array_name, numMemOp));
		}
	}



}

void BaseDatapath::calculateArrayName2maxReadWrite(unsigned int iterations) {
	
	arrayName2maxMemOpNumTy::iterator it, ie;
	it = arrayName2maxMemOpNum.begin();
	ie = arrayName2maxMemOpNum.end();
	for (; it != ie; ++it) {
		std::string array_name = it->first;
		uint64_t maxMemOp = it->second;
		it->second = maxMemOp * iterations;
	}

}

void BaseDatapath::writeLogOfArrayName2maxReadWrite() {
	std::string file_output = outputPath + benchName + "_array.log";
	std::ofstream array_file;
	array_file.open(file_output);
	if (array_file.is_open()) {
		arrayName2maxMemOpNumTy::iterator it, ie;
		it = arrayName2maxMemOpNum.begin();
		ie = arrayName2maxMemOpNum.end();
		for (; it != ie; ++it) {
			array_file << it->first << ": " << it->second << std::endl;
		}
	}
	else {
		assert(false && "Error: Cannot open array file!\n");
	}

	array_file.close();
}

std::string BaseDatapath::getTargetLoopName() const {
	return target_loop_name;
}

unsigned BaseDatapath::getTargetLoopLevel() const {
	return target_loop_level;
}

unsigned BaseDatapath::getTargetLoopLevelUnrollFactor() const {
	return target_lp_level_unroll_factor;
}

uint64_t BaseDatapath::run_fpga_simulation() {

	num_cycles = 0;

	unsigned faddNum = 4;
	unsigned fsubNum = 1;
	unsigned fmulNum = 6;
	unsigned fdivNum = 1;
	unsigned bramNum = INFINITE_HARDWARE;
	unsigned readPNum = 2;
	unsigned writePNum = 1;
	set_fpga_constraints(faddNum, fsubNum, fmulNum, fdivNum, bramNum, readPNum, writePNum);

	std::string func_name = kernel_names.at(0);

	assert((funcName2loopNumMap.find(func_name) != funcName2loopNumMap.end()) && "Error! Function name is inside funcName2loopNumMap");
	unsigned numLoopInaFunc = funcName2loopNumMap.find(func_name)->second;

#ifdef WRITE_EXECUTION_RECORD
	//execution_record.open(inputPath + func_name + "_execution_record.txt");
	execution_record.open(outputPath + func_name + "_execution_record.txt");
	execution_record << "cycle, executed instructions\n";
#endif // End of WRITE_EXECUTION_RECORD

	if (numLoopInaFunc < 2) {
		/// Initialize executingQueue, readyToExecuteQueue
		setGraphForStepping();

#ifdef WRITE_EXECUTION_RECORD
		execution_record << num_cycles << ", ";
#endif // End of WRITE_EXECUTION_RECORD

		/// Run simulation
		while (!step()) {
			//clearOccupiedBWPerPartition();
		}
	}
	else {
		NameVecTy cur_basic_block(numTotalNodes, "");
		initCurBasicBlock(cur_basic_block);

		/// Initialize executingQueue, readyToExecuteQueue
		setGraphForStepping(cur_basic_block);

#ifdef WRITE_EXECUTION_RECORD
		execution_record << num_cycles << ", ";
#endif // End of WRITE_EXECUTION_RECORD

		/// Run simulation
		while (!step(cur_basic_block)) {
			//clearOccupiedBWPerPartition();
		}
	}
	
#ifdef WRITE_EXECUTION_RECORD
	execution_record.close();
#endif // End of WRITE_EXECUTION_RECORD

	/// Return FPGA
	return num_cycles;
}

//stepFunctions
//multiple function, each function is a separate graph
// For FPGA simulation
void BaseDatapath::setGraphForStepping(NameVecTy& curBBname)
{
	std::cerr << "=============================================" << std::endl;
	std::cerr << "      Scheduling...            " << benchName << std::endl;
	std::cerr << "=============================================" << std::endl;

	newLevel.assign(numTotalNodes, 0);

	edgeToParid = get(boost::edge_weight, graph_);

	numTotalEdges = boost::num_edges(graph_);
	numParents.assign(numTotalNodes, 0);
	latestParents.assign(numTotalNodes, 0);
	executedNodes = 0;
	totalConnectedNodes = 0;
	for (unsigned node_id = 0; node_id < numTotalNodes; node_id++)
	{
		Vertex node = nameToVertex[node_id];
		if (boost::degree(node, graph_) != 0 || is_dma_op(microop.at(node_id)))
		{
			finalIsolated.at(node_id) = 0;
			numParents.at(node_id) = boost::in_degree(node, graph_);
			totalConnectedNodes++;
		}
	}

	loopExecutionTracer.clear();

	executingQueue.clear();
	readyToExecuteQueue.clear();
	waitForLoopQueue.clear();
	initExecutingQueue(curBBname);
}

///FIXME: need to transfer node_id2funcNameMap, node_id2bbNameMap as arguments to trace execution
/// For FPGA simulation
void BaseDatapath::initExecutingQueue(NameVecTy& curBBname) {
	/// For multiple function, each function is a separate graph. Thus, we do not need funcName2loopExecutionTracerMap.
	/// Instead, we only need local loopExecutionTracerMap.

	///FIXME: func_name should be related to the graph, but in the first step, since all our testbenches do not
	///       have multiple functions, thus we just use the first elements inside kernel_names. For multiple functions
	///       we need to update the following statement.
	std::string func_name = kernel_names.at(0);

	assert((funcName2loopNumMap.find(func_name) != funcName2loopNumMap.end()) && "Error! Function name is inside funcName2loopNumMap");
	unsigned numLoopInaFunc = funcName2loopNumMap.find(func_name)->second;

	for (unsigned i = 0; i < numLoopInaFunc; i++) {
		if (i == 0) {
			loopExecutionTracer.push_back(Loop_Status::EXECUTING);
		}
		else {
			loopExecutionTracer.push_back(Loop_Status::INIT);
		}
	}

	/// Initialize executingQueue.
	for (unsigned i = 0; i < numTotalNodes; i++)
	{
		if (numParents[i] == 0 && finalIsolated[i] != 1) {
			assert((curBBname.size() > i) && "Error: size of curBBname is less than i!\n");
			std::string bb_name = curBBname.at(i);
			bbFuncNamePairTy bb_func_pair = std::make_pair(bb_name, func_name);
			bbFuncNamePair2lpNameLevelPairMapTy::iterator it = bbFuncNamePair2lpNameLevelPairMap.find(bb_func_pair);
			bbFuncNamePair2lpNameLevelPairMapTy::iterator ie = bbFuncNamePair2lpNameLevelPairMap.end();
			if (it != ie) {
				std::string loop_name = bbFuncNamePair2lpNameLevelPairMap[bb_func_pair].first;
				std::size_t pos = loop_name.find("-");
				std::string sub_str = loop_name.substr(pos + 1);
				unsigned ith_loop = (unsigned)std::stoi(sub_str);

				assert((loopExecutionTracer.size() > ith_loop) && "Error: ith loop is not inside loopExecutionTracer!\n");
				if (loopExecutionTracer.at(ith_loop) == Loop_Status::EXECUTING) {
					executingQueue.push_back(i);
				}
				else {
					waitForLoopQueue.push_back(i);
					/*
					unsigned opcode = microop.at(i);
					unsigned op_latency = fpga_node_latency(opcode);
					if ( op_latency == 0) {
						readyToExecuteQueue.push_back(i);
					}
					else {
						nodeID2delay.insert( std::make_pair(i, op_latency) );
					}
					*/
				}
			}
			else {
				/// Statements here are outside loops, just put into executingQueue
				executingQueue.push_back(i);
			}
		}
	}

}

/// For FPGA simulation
bool BaseDatapath::step(NameVecTy& bbnames)
{
	stepExecutingQueue(bbnames);
	copyToExecutingQueue();
	num_cycles++;
#ifdef WRITE_EXECUTION_RECORD
	execution_record << "\n";
	execution_record << num_cycles << ", ";
#endif // End of WRITE_EXECUTION_RECORD
	updateDelayForNodeID();
	updateWaitForLoopQueue(bbnames);

	if (executedNodes == totalConnectedNodes)
		return 1;
	return 0;
}

/// For FPGA simulation
void BaseDatapath::stepExecutingQueue(NameVecTy& bbnames) {

	auto it = executingQueue.begin();
	int index = 0;
	while (it != executingQueue.end()) {
		unsigned node_id = *it;
		/*
		/// After loopFlattening, all the compute operations will be replaced with LLVM_IR_Move microop,
		/// This is a bug inside loopFlattening, we can not replace the opcodes.
		if ( is_compute_op(microop.at(node_id)) ) {
		cout << "This is a compute instruction: " << microop.at(node_id) << endl;
		}*/
		unsigned opcode = microop.at(node_id);
		if (is_memory_op(opcode))
		{
			std::string node_part = baseAddress[node_id].first;
			//if (registers.has(node_part) || canServicePartition(node_part))
			if ( registers.has(node_part) )
			{
				if (is_load_op(microop.at(node_id)))
					registers.getRegister(node_part)->increment_loads();
				else
					registers.getRegister(node_part)->increment_stores();

			}
			else {
				if (is_load_op(microop.at(node_id)))
					increment_loads(node_part);
				else
					increment_stores(node_part);
			}
		}

		executedNodes++;
#ifdef WRITE_EXECUTION_RECORD
		execution_record << node_id << "-f, ";
#endif // End of WRITE_EXECUTION_RECORD
		newLevel.at(node_id) = num_cycles;
		executingQueue.erase(it);
		updateChildren(node_id, bbnames);
		it = executingQueue.begin();
		std::advance(it, index);
	}

}

/// For FPGA simulation
void BaseDatapath::updateChildren(unsigned node_id, NameVecTy& bb_names)
{
	std::string func_name = kernel_names.at(0);

	Vertex node = nameToVertex[node_id];
	out_edge_iter out_edge_it, out_edge_end;
	for (tie(out_edge_it, out_edge_end) = out_edges(node, graph_); out_edge_it != out_edge_end; ++out_edge_it)
	{
		unsigned child_id = vertexToName[target(*out_edge_it, graph_)];
		int edge_parid = edgeToParid[*out_edge_it];
		if (numParents[child_id] > 0)
		{
			numParents[child_id]--;
			if (numParents[child_id] == 0)
			{
				std::string bb_name = bb_names.at(child_id);
				bbFuncNamePairTy bb_func_pair = std::make_pair(bb_name, func_name);
				bbFuncNamePair2lpNameLevelPairMapTy::iterator it = bbFuncNamePair2lpNameLevelPairMap.find(bb_func_pair);
				bbFuncNamePair2lpNameLevelPairMapTy::iterator ie = bbFuncNamePair2lpNameLevelPairMap.end();
				if (it != ie) {
					unsigned ith_loop = getIthLoop(bb_name, func_name);
					if (loopExecutionTracer.at(ith_loop) == Loop_Status::EXECUTING) {
						unsigned child_microop = microop.at(child_id);
						unsigned curr_microop = microop.at(node_id);
						unsigned op_latency = fpga_node_latency(child_microop);
						/*
						if ( op_latency == 0 && edge_parid != CONTROL_EDGE) {
							executingQueue.push_back(child_id);
						}*/
						if ((op_latency == 0 && edge_parid != CONTROL_EDGE) || (child_microop == LLVM_IR_Br && edge_parid == CONTROL_EDGE)) {
							executingQueue.push_back(child_id);
						}
						else {
							nodeID2delay.insert(std::make_pair(child_id, op_latency));
						}
					}
					else {
						assert( (loopExecutionTracer.at(ith_loop) != Loop_Status::FINISHED) && "Error: There exists a node that belongs to a loop with FINISHED status!\n");
						waitForLoopQueue.push_back(child_id);
					}
				}
				else {
					// this child node is not inside any loops
					executingQueue.push_back(child_id);
				}
				numParents[child_id] = -1;
			}
		}
	}
}

/// For FPGA simulation
///FIXME: potential bug: if there exists a node id belongs to a loop with "FINISH" status, it will make the whole
///       program to run in while loop forever. Maybe it is the graph's problem. Be careful.
/*
void BaseDatapath::copyToExecutingQueue(NameVecTy& bb_names)
{
	std::string func_name = kernel_names.at(0);
	auto it = readyToExecuteQueue.begin();
	while (it != readyToExecuteQueue.end())
	{
		std::string bb_name = bb_names.at(*it);
		bbFuncNamePairTy bb_func_pair = std::make_pair(bb_name, func_name);
		bbFuncNamePair2lpNameLevelPairMapTy::iterator it_map = bbFuncNamePair2lpNameLevelPairMap.find(bb_func_pair);
		bbFuncNamePair2lpNameLevelPairMapTy::iterator ie_map = bbFuncNamePair2lpNameLevelPairMap.end();
		// This node id belongs to which loop
		unsigned ith_loop = getIthLoop(bb_name, func_name);

		if (it_map != ie_map ) {
			if (loopExecutionTracer.at(ith_loop) == Loop_Status::EXECUTING) {
				// This node id belongs to a loop with "EXECUTING" status
				executingQueue.push_back(*it);
				it = readyToExecuteQueue.erase(it);
			}
			else {
				// This node id belongs to a loop with "INIT" status. We need to check whether there exists an element
				// that belongs to a loop with "EXECUTING" status in executingQueue, readyToExecuteQueue and 
				// nodeID2delay. If no node belongs to a loop with "EXECUTING" status, then
				// 1. Set the loop with "EXECUTING" status to "FINISH" status;
				// 2. Set the next loop of the "FINISH" loop to status "EXECUTING" (from "INIT");
				// 3. If this node id belongs to this loop, then push back to executingQueue, otherwise keep it inside
				//    readyToExecuteQueue and increases iterator it.
				bool executingLoop_found = false;
				unsigned size_loopTracer = loopExecutionTracer.size();
				unsigned executing_loop = 0;
				for (unsigned i = 0; i < size_loopTracer; i++) {
					if (loopExecutionTracer.at(i) == Loop_Status::EXECUTING) {
						executing_loop = i;
						break;
					}
				}

				unsigned size_exeQueue = executingQueue.size();
				for (unsigned i = 0; i < size_exeQueue; i++) {
					unsigned node_id = executingQueue.at(i);
					std::string bbName = bb_names.at(node_id);
					unsigned ithLoop = getIthLoop(bbName, func_name);
					if (ithLoop == executing_loop) {
						executingLoop_found = true;
						break;
					}
				}

				if (executingLoop_found == false) {
					unsigned size_readyQueue = readyToExecuteQueue.size();
					for (unsigned i = 0; i < size_readyQueue; i++) {
						unsigned node_id = readyToExecuteQueue.at(i);
						std::string bbName = bb_names.at(node_id);
						unsigned ithLoop = getIthLoop(bbName, func_name);
						if (ithLoop == executing_loop) {
							executingLoop_found = true;
							break;
						}
					}

					if (executingLoop_found == false) {
						std::map<node_idTy, unsigned>::iterator it_delay = nodeID2delay.begin();
						std::map<node_idTy, unsigned>::iterator ie_delay = nodeID2delay.end();
						for (; it_delay != ie_delay; it_delay++) {
							unsigned node_id = it_delay->first;
							std::string bbName = bb_names.at(node_id);
							unsigned ithLoop = getIthLoop(bbName, func_name);
							if (ithLoop == executing_loop) {
								executingLoop_found = true;
								break;
							}
						}

						if (executingLoop_found == false) {
							if ( executing_loop < size_loopTracer-1 ) {
								loopExecutionTracer.at(executing_loop+1) = Loop_Status::EXECUTING;
							}
							loopExecutionTracer.at(executing_loop) = Loop_Status::FINISHED;

							if (ith_loop == executing_loop) {
								// This node id belongs to a new loop with "EXECUTING" status
								executingQueue.push_back(*it);
								it = readyToExecuteQueue.erase(it);
							}
							else {
								// This node id doesn't belong to a new loop with "EXECUTING" status
								it++;
							}
						}
						else {
							it++;
							executingLoop_found = false;
						}

					}
					else {
						it++;
						executingLoop_found = false;
					}
				}
				else {
					it++;
					executingLoop_found = false;
				}

			}
		}
		else {
			// This node id is not inside any loops
			executingQueue.push_back(*it);
			it = readyToExecuteQueue.erase(it);
		}
	}
}
*/

/// BaseDatapath::updateWaitForLoopQueue(NameVecTy& bb_names):
///		Check sizes of executingQueue, readyToExecuteQueue and nodeID2delay; 
///		If all 0, then need to update loops' status and copy data to executingQueue or nodeID2delay;
///		Otherwise, just ignore and do nothing.
void BaseDatapath::updateWaitForLoopQueue(NameVecTy& bb_names) {
	std::string func_name = kernel_names.at(0);

	/// Get ith loop with "EXECUTING" status
	int executingLoop = -1;
	unsigned size_loopTracer = loopExecutionTracer.size();
	for (unsigned i = 0; i<size_loopTracer; i++) {
		if (loopExecutionTracer.at(i) == Loop_Status::EXECUTING) {
			executingLoop = i;
			break;
		}
	}

	unsigned size_Executing = executingQueue.size();
	unsigned size_Ready = readyToExecuteQueue.size();
	unsigned size_Delay = nodeID2delay.size();
	if ( size_Executing == 0 && size_Ready == 0 && size_Delay == 0 ) {
		// Change loop status
		loopExecutionTracer.at(executingLoop) = Loop_Status::FINISHED;
		if (executingLoop < (int) (size_loopTracer - 1) ) {
			executingLoop += 1;
			loopExecutionTracer.at(executingLoop) = Loop_Status::EXECUTING;
			unsigned size_wait = waitForLoopQueue.size();
			for (unsigned i = 0; i < size_wait; i++) {
				unsigned node_id = waitForLoopQueue.at(i);
				std::string bb_name = bb_names.at(node_id);
				unsigned ith_loop = getIthLoop(bb_name, func_name);
				assert( (loopExecutionTracer.at(ith_loop) != Loop_Status::FINISHED) && "Error: There exists a node that belongs to a loop with FINISHED status!\n" );
				if ( (int) ith_loop == executingLoop){
					// Only push back node ids to executingQueue or nodeID2delay when ith_loop has status "EXECUTING"
					unsigned op_latency = fpga_node_latency(microop.at(node_id));
					if (op_latency == 0) {
						executingQueue.push_back(node_id);
					}
					else {
						nodeID2delay.insert(std::make_pair(node_id, op_latency));
					}
				}
			}
		}
	}

}

int BaseDatapath::shortestDistanceBetweenNodes(unsigned int from, unsigned int to)
{
  std::list<pair<unsigned int, unsigned int> > queue;
  queue.push_back({from, 0});
  while(queue.size() != 0)
  {
    unsigned int curr_node = queue.front().first;
    unsigned int curr_dist = queue.front().second;
    out_edge_iter out_edge_it, out_edge_end;
    for (tie(out_edge_it, out_edge_end) = out_edges(nameToVertex[curr_node], graph_); out_edge_it != out_edge_end; ++out_edge_it)
    {
      if (get(boost::edge_weight, graph_, *out_edge_it) != CONTROL_EDGE)
      {
        int child_id = vertexToName[target(*out_edge_it, graph_)];
        if (child_id == to)
          return curr_dist + 1;
        queue.push_back({child_id, curr_dist + 1});
      }
    }
    queue.pop_front();
  }
  return -1;
}

//readConfigs
void BaseDatapath::readArrayInfo(arrayInfoMapTy& arrayInfoMap) {
	ifstream config_file;
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
	file_name += "_array_info";
	config_file.open(file_name.c_str());
	if (!config_file.is_open()){
		assert(false && "DEBUG-INFO: [configuration-array_info] Please provide array information first!\n");
	}
		
	while (!config_file.eof())
	{
		std::string wholeline;
		getline(config_file, wholeline);
		if (wholeline.size() == 0)
			break;
		char type[256];
		char array_name[256];
		uint64_t size_in_byte;
		unsigned wordsize_in_byte;
		sscanf(wholeline.c_str(), "%[^,],%[^,],%lld,%d\n", type, array_name, &size_in_byte, &wordsize_in_byte);
		arrayInfoMap.insert(std::make_pair(array_name, std::make_pair(size_in_byte, wordsize_in_byte)));
	}
	config_file.close();
}

bool BaseDatapath::readPipeliningConfig()
{
	bool enable_pipelining = false;
  ifstream config_file;
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_pipelining_config";
  config_file.open(file_name.c_str());
  if (!config_file.is_open())
    return 0;

	while (!config_file.eof())
	{
		std::string wholeline;
		getline(config_file, wholeline);
		if (wholeline.size() == 0)
			break;

		char config_type[256];
		char func[256];
		unsigned ith_loop;
		unsigned loop_level;
		sscanf(wholeline.c_str(), "%[^,],%[^,],%d,%d,\n", config_type, func, &ith_loop, &loop_level);
		std::string conf_type(config_type);
		std::string func_name(func);
		if (conf_type.compare("pipeline")) {
			continue;
		}
		std::string lp_name = func_name + "_loop-" + std::to_string(ith_loop);
		if (lp_name == target_loop_name && loop_level == target_loop_level) {
			enable_pipelining = true;
		}
	}

	config_file.close();
	return enable_pipelining;
}

bool BaseDatapath::readUnrollingConfig(std::unordered_map<int, int > &unrolling_config)
{
  ifstream config_file;
  //std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_unrolling_config";
  config_file.open(file_name.c_str());
  if (!config_file.is_open())
    return 0;
  while(!config_file.eof())
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
    unrolling_config[line_num] =factor;
  }
  config_file.close();
  return 1;
}

bool BaseDatapath::readFlattenConfig(std::unordered_set<int> &flatten_config)
{
  ifstream config_file;
  std::string file_name(benchName);
  file_name += "_flatten_config";
  config_file.open(file_name.c_str());
  if (!config_file.is_open())
    return 0;
  while(!config_file.eof())
  {
    std::string wholeline;
    getline(config_file, wholeline);
    if (wholeline.size() == 0)
      break;
    char func[256];
    int line_num;
    sscanf(wholeline.c_str(), "%[^,],%d\n", func, &line_num);
    flatten_config.insert(line_num);
  }
  config_file.close();
  return 1;
}

bool BaseDatapath::readCompletePartitionConfig(std::unordered_map<std::string, unsigned> &config)
{
	//std::string comp_partition_file(inputPath + benchName);
	std::string comp_partition_file(outputPath + benchName);
	comp_partition_file += "_complete_partition_config";

	if (!fileExists(comp_partition_file))
		return 0;

	ifstream config_file;
	config_file.open(comp_partition_file);
	if (!config_file.is_open()) {
		return 0;
	}

	std::string wholeline;
	while (!config_file.eof())
	{
		getline(config_file, wholeline);
		if (wholeline.size() == 0)
			break;

		unsigned size;
		char config_type[256];
		char type[256];
		char base_addr[256];
		sscanf(wholeline.c_str(), "%[^,],%[^,],%[^,],%d\n", config_type, type, base_addr, &size);
		config[base_addr] = size;
		partitionArrayName.insert(base_addr);
	}
	config_file.close();
	
  return 1;
}

bool BaseDatapath::readPartitionConfig(std::unordered_map<std::string, partitionEntry> & partition_config)
{
  ifstream config_file;
	//std::string file_name(inputPath + benchName);
	std::string file_name(outputPath + benchName);
  file_name += "_partition_config";
  if (!fileExists(file_name))
    return 0;

  config_file.open(file_name.c_str());
	if (!config_file.is_open()) {
		return 0;
	}

  std::string wholeline;
  while (!config_file.eof())
	{
    getline(config_file, wholeline);
    if (wholeline.size() == 0) 
			break;

    unsigned size, p_factor, wordsize;
		char config_type[256];
    char type[256];
    char base_addr[256];
    sscanf(wholeline.c_str(), "%[^,],%[^,],%[^,],%d,%d,%d\n", config_type, type, base_addr, &size, &wordsize, &p_factor);
    std::string p_type(type);
		partitionEntry partEntry(p_type, size, wordsize, p_factor);
    partition_config[base_addr] = partEntry;
		partitionArrayName.insert(base_addr);
  }
  config_file.close();
  return 1;
}

void dump_summary(ofstream& summary_file, loopInfoTy loop_info, resourceTy fpga_rs, sharedMemTy sharedMem, arrayName2memEfficiencyTy& Memefficiency, std::vector<std::string>& limited_op_types, subTraceInstTy& sub_traceInst, float aveParallelism, arrayName2maxMemOpNumTy& arrayN2maxMemOpNum, arrayName2MemBankStatisTy arrayN2aveLdPerBank, arrayName2MemBankStatisTy arrayN2aveStPerBank) {
	unsigned resII_mem = loop_info.ResII_mem;
	unsigned resII_op = loop_info.ResII_op;
	unsigned recII = loop_info.RecII;
	summary_file << "==========================" << std::endl;
	summary_file << "loop name: " << loop_info.loop_name << std::endl;
	summary_file << "loop level: " << loop_info.lp_level << std::endl;
	summary_file << "loop unrolling factor: " << loop_info.lp_level_unroll_factor << std::endl;
	summary_file << "loop pipelining enabled: " << loop_info.enable_pipelining << std::endl;
	summary_file << "total FPGA cycles: " << loop_info.num_cycles << std::endl;
	summary_file << "--------------------------" << std::endl;
	summary_file << "number of shared loads: " << sharedMem.shared_loads << std::endl;
	summary_file << "number of repeated stores: " << sharedMem.repeated_stores << std::endl;
	summary_file << "--------------------------" << std::endl;
	summary_file << "ideal iteration latency: " << loop_info.IL_asap << std::endl;
	summary_file << "Iteration Latency: " << loop_info.rc_list_il << std::endl;
	summary_file << "Initiation Interval (if applicable): " << loop_info.max_II << std::endl;
	summary_file << "resII_mem: " << resII_mem << std::endl;
	summary_file << "resII_op: " << resII_op << std::endl;
	summary_file << "recII: " << recII << std::endl;
	summary_file << "limited factor: ";
	if (resII_mem > resII_op && resII_mem > recII && resII_mem > 1) {
		summary_file << "Memory, array name = " << loop_info.limited_mem_name << std::endl;
	}
	else if (resII_op > resII_mem && resII_op > recII && resII_op > 1) {
		summary_file << "Floating point operation, " << loop_info.limited_op_name << std::endl;
	}
	else if (recII > resII_mem && recII > resII_op && recII > 1) {
		summary_file << "Loop carried dependency" << std::endl;
	}
	else {
		summary_file << "No" << std::endl;
	}
	summary_file << "--------------------------" << std::endl;
	summary_file << "limited by DSP resource: ";
	unsigned size_vec = limited_op_types.size();
	if (size_vec > 1) {
		summary_file << "Yes" << std::endl;
		for (unsigned i = 0; i < size_vec; i++) {
			std::string fop_type = limited_op_types.at(i);
			if (fop_type == "FADD_UNIT_CONSTRAINT") {
				summary_file << "\t\tLimited floating point operation type: Fadd" << std::endl;
			}
			else if (fop_type == "FSUB_UNIT_CONSTRAINT") {
				summary_file << "\t\tLimited floating point operation type: Fsub" << std::endl;
			}
			else if (fop_type == "FMUL_UNIT_CONSTRAINT") {
				summary_file << "\t\tLimited floating point operation type: Fmul" << std::endl;
			}
			else if (fop_type == "FDIV_UNIT_CONSTRAINT") {
				summary_file << "\t\tLimited floating point operation type: Fdiv" << std::endl;
			}
			else {
				assert(false && "Error: Do not have this floating point operation type!\n");
			}
		}
	}
	else {
		summary_file << "No" << std::endl;
	}

	summary_file << "--------------------------" << std::endl;
	summary_file << "DSP used: " << fpga_rs.dsp_used << std::endl;
	summary_file << "BRAM18K used: " << fpga_rs.bram18k_used << std::endl;
	summary_file << "FF used: " << fpga_rs.ff_used << std::endl;
	summary_file << "LUT used: " << fpga_rs.lut_used << std::endl;
	summary_file << "Fadd used: " << fpga_rs.fadd_used << std::endl;
	summary_file << "Fsub used: " << fpga_rs.fsub_used << std::endl;
	summary_file << "Fmul used: " << fpga_rs.fmul_used << std::endl;
	summary_file << "Fdiv used: " << fpga_rs.fdiv_used << std::endl;
	summary_file << "--------------------------" << std::endl;
	summary_file << "memory efficiency: " << std::endl;
	arrayName2memEfficiencyTy::iterator it = Memefficiency.begin();
	arrayName2memEfficiencyTy::iterator ie = Memefficiency.end();
	for (; it != ie; ++it) {
		summary_file << it->first << ": " << it->second << std::endl;
	}
	summary_file << "--------------------------" << std::endl;
	summary_file << "  Sub-trace Information   " << std::endl;
	/// Instructions
	summary_file << "load instruction num: " << sub_traceInst.num_ld_inst << std::endl;
	summary_file << "store instruction num: " << sub_traceInst.num_st_inst << std::endl;
	summary_file << "fadd instruction num: " << sub_traceInst.num_fadd_inst << std::endl;
	summary_file << "fsub instruction num: " << sub_traceInst.num_fsub_inst << std::endl;
	summary_file << "fmul instruction num: " << sub_traceInst.num_fmul_inst << std::endl;
	summary_file << "fdiv instruction num: " << sub_traceInst.num_fdiv_inst << std::endl;
	summary_file << "fcmp instruction num: " << sub_traceInst.num_fcmp_inst << std::endl;
	summary_file << "integer instruction num: " << sub_traceInst.num_integer_inst << std::endl;
	summary_file << "bitwise instruction num: " << sub_traceInst.num_bitwise_inst << std::endl;
	// Control instruction includes PHI instructions and index addition/substraction
	summary_file << "control instruction num: " << sub_traceInst.num_control_inst << std::endl;
	summary_file << "branch instruction num: " << sub_traceInst.num_br_inst << std::endl;
	/// average parallelism
	summary_file << "parallelism profile -- average parallelism: " << aveParallelism << std::endl;
	/// Memory bank information
	arrayName2maxMemOpNumTy::iterator it_maxMem, ie_maxMem;
	it_maxMem = arrayN2maxMemOpNum.begin();
	ie_maxMem = arrayN2maxMemOpNum.end();
	uint64_t maxMemOpNum = 0;
	for (; it_maxMem != ie_maxMem; ++it_maxMem) {
		uint64_t max_memop_num = it_maxMem->second;
		maxMemOpNum = maxMemOpNum > max_memop_num ? max_memop_num : max_memop_num;
	}
	summary_file << "maximum number of memory access to a memory bank: " << maxMemOpNum << std::endl;
	arrayName2MemBankStatisTy::iterator it_aveLdPBk, ie_aveLdPBk;
	it_aveLdPBk = arrayN2aveLdPerBank.begin();
	ie_aveLdPBk = arrayN2aveLdPerBank.end();
	float max_aveLdperbank = 0.0;
	for (; it_aveLdPBk != ie_aveLdPBk; ++it_aveLdPBk) {
		float aveLd = it_aveLdPBk->second;
		max_aveLdperbank = max_aveLdperbank > aveLd ? max_aveLdperbank : aveLd;
	}
	summary_file << "maximum average number of memory bank load access: " << max_aveLdperbank << std::endl;
	arrayName2MemBankStatisTy::iterator it_aveStPBk, ie_aveStPBk;
	it_aveStPBk = arrayN2aveStPerBank.begin();
	ie_aveStPBk = arrayN2aveStPerBank.end();
	float max_aveStperbank = 0.0;
	for (; it_aveStPBk != ie_aveStPBk; ++it_aveStPBk) {
		float aveSt = it_aveStPBk->second;
		max_aveStperbank = max_aveStperbank > aveSt ? max_aveStperbank : aveSt;
	}
	summary_file << "maximum average number of memory bank store access: " << max_aveStperbank << std::endl;
	summary_file << "--------------------------" << std::endl;
	summary_file << "Finished " << loop_info.loop_name + "_" + std::to_string(loop_info.lp_level) << std::endl;
	summary_file << "--------------------------" << std::endl;
	
}

void BaseDatapath::set_fpga_constraints(unsigned faddNum, unsigned fsubNum, unsigned fmulNum, unsigned fdivNum, unsigned bramNum, unsigned readPNum, unsigned writePNum) {
	fpga_constraints.setConstraints(faddNum, fsubNum, fmulNum, fdivNum, bramNum, readPNum, writePNum);
}

/* Tokenizes an input string and returns a vector. */
void BaseDatapath::tokenizeString(std::string input,
                                  std::vector<int>& tokenized_list)
{
  using namespace boost;
  tokenizer<> tok(input);
  for(tokenizer<>::iterator beg = tok.begin(); beg != tok.end(); ++beg)
  {
    int value;
    istringstream(*beg) >> value;
    tokenized_list.push_back(value);
  }
}