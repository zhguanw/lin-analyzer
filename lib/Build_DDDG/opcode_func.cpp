#include "profile_h/opcode_func.h"
#include "profile_h/InstrumentForDDDGPass.h"
#include "profile_h/auxiliary.h"

extern bool increase_load_latency;

bool is_associative(unsigned microop)
{
  if (microop == LLVM_IR_Add )
    return true;
  return false;
}

bool is_fassociative(unsigned microop)
{
	if (microop == LLVM_IR_FAdd)
		return true;
	return false;
}

bool is_memory_op(unsigned microop)
{
  if (microop == LLVM_IR_Load || microop == LLVM_IR_Store)
    return true;
  return false;
}

bool is_compute_op(unsigned microop)
{
  switch(microop)
  {
    case LLVM_IR_Add : case LLVM_IR_FAdd : case LLVM_IR_Sub : case LLVM_IR_FSub :
    case LLVM_IR_Mul : case LLVM_IR_FMul : case LLVM_IR_UDiv : case LLVM_IR_SDiv :
    case LLVM_IR_FDiv : case LLVM_IR_URem : case LLVM_IR_SRem : case LLVM_IR_FRem :
    case LLVM_IR_Shl : case LLVM_IR_LShr : case LLVM_IR_AShr : case LLVM_IR_And :
		case LLVM_IR_Or: case LLVM_IR_Xor: case LLVM_IR_ICmp : case LLVM_IR_FCmp :
      return true;
    default:
      return false;
  }
}

bool is_store_op(unsigned microop)
{
  if (microop == LLVM_IR_Store)
    return true;
  return false;
}

bool is_load_op(unsigned microop)
{
  if (microop == LLVM_IR_Load)
    return true;
  return false;
}

bool is_bit_op(unsigned microop)
{
  switch (microop)
  {
    case LLVM_IR_Shl: case LLVM_IR_LShr: case LLVM_IR_AShr :
      case LLVM_IR_And : case LLVM_IR_Or: case LLVM_IR_Xor :
      return true;
    default:
      return false;
  }
}

bool is_control_op (unsigned microop)
{
  if (microop == LLVM_IR_PHI)
    return true;
  return is_branch_op(microop);
}

bool is_phi_op (unsigned microop) {
	return (microop == LLVM_IR_PHI) ? true : false;
}

bool is_branch_op (unsigned microop)
{
  if (microop == LLVM_IR_Br || microop == LLVM_IR_Switch)
    return true;
  return is_call_op(microop);
}

bool is_call_op(unsigned microop)
{
  if (microop == LLVM_IR_Call)
    return true;
  return is_dma_op(microop);
}

bool is_index_op (unsigned microop)
{
	if (microop == LLVM_IR_IndexAdd || microop == LLVM_IR_IndexSub)
    return true;
  return false;
}

bool is_dma_load(unsigned microop)
{
  return microop == LLVM_IR_DMALoad ;
}

bool is_dma_store(unsigned microop)
{
  return microop == LLVM_IR_DMAStore ;
}

bool is_dma_op(unsigned microop)
{
  return is_dma_load(microop) || is_dma_store(microop);
}

bool is_mul_op(unsigned microop)
{
  switch(microop)
  {
    case LLVM_IR_Mul: case LLVM_IR_UDiv: case LLVM_IR_FMul : case LLVM_IR_SDiv :
    case LLVM_IR_FDiv : case LLVM_IR_URem : case LLVM_IR_SRem : case LLVM_IR_FRem:
      return true;
    default:
      return false;
  }
}

bool is_add_op(unsigned microop)
{
  if (microop == LLVM_IR_Add || microop == LLVM_IR_Sub)
    return true;
  return false;
}

bool is_fadd_op(unsigned microop){
	return ( (microop == LLVM_IR_FAdd) ? 1 : 0 );
}

bool is_fsub_op(unsigned microop){
	return ((microop == LLVM_IR_FSub) ? 1 : 0);
}

bool is_fmul_op(unsigned microop){
	return ((microop == LLVM_IR_FMul) ? 1 : 0);
}

bool is_fdiv_op(unsigned microop){
	return ((microop == LLVM_IR_FDiv) ? 1 : 0);
}

bool is_fcmp_op(unsigned microop){
	return ((microop == LLVM_IR_FCmp) ? 1 : 0);
}

bool is_float_op(unsigned microop) {
	bool flag = (microop == LLVM_IR_FAdd || microop == LLVM_IR_FSub || microop == LLVM_IR_FMul || microop == LLVM_IR_FDiv);
	return (flag ? 1 : 0);
}

unsigned fpga_node_latency (unsigned  microop)
{
	switch (microop) {
		case LLVM_IR_Shl:
		case LLVM_IR_LShr:
		case LLVM_IR_AShr:
		case LLVM_IR_And:
		case LLVM_IR_Or:
		case LLVM_IR_Xor:
		case LLVM_IR_ICmp:
		case LLVM_IR_Br:
		case LLVM_IR_IndexAdd:
		case LLVM_IR_IndexSub:
			return 0;
		case LLVM_IR_Add:
		case LLVM_IR_Sub: 
			return ADDSUB_LATENCY;
		case LLVM_IR_Call:
			return 0;
		case LLVM_IR_Store:
			return STORE_LATENCY;
		case LLVM_IR_SilentStore:
			return 0;
		case LLVM_IR_Load:
			// Observation from Vivado HLS: 
			//   1. When enabling partitioning without pipelining, a load operation 
			//      from BRAM needs 2 access latency
			//   2. When no partitioning, a load operation from  BRAM needs 1 access 
			//      latency
			// FIXME: In current implementation, we apply 2 load latency to all arrays 
			//        for simplicity. But to further improving accuracy, it is better
			//        to ONLY associate 2 load latency with partitioned arrays and use
			//        1 load latency for normal arrays.
			if (increase_load_latency == true) {
				return LOAD_LATENCY;
			}
			else {
				return LOAD_LATENCY - 1;
			}
		case LLVM_IR_Mul:
			// 64 bits -- 18 cycles
			// 50 bits -- 11 cycles
			// 32 bits -- 6 cycles
			// 24 bits -- 3 cycles
			// 20 bits -- 3 cycles
			// 18 bits -- 1 cycles
			// 16 bits -- 1 cycles
			// 10 bits -- 1 cycles
			// 8  bits -- 1 cycles
			// Currently, we only consider 32-bit applications
			return MUL32BIT_LATENCY;
		case LLVM_IR_UDiv:
		case LLVM_IR_SDiv:
			// 64 bits -- 68 cycles
			// 50 bits -- 54 cycles
			// 32 bits -- 36 cycles
			// 24 bits -- 28 cycles
			// 16 bits -- 20 cycles
			// 10 bits -- 14 cycles
			// 8  bits -- 12 cycles
			// Currently, we only consider 32-bit applications
			return DIV32BIT_LATENCY;
		case LLVM_IR_FAdd:
		case LLVM_IR_FSub:
			// 32/64 bits -- 5 cycles
			return FADDSUB32BIT_LATENCY;
		case LLVM_IR_FMul:
			// 64 bits -- 6 cycles
			// 32 bits -- 4 cycles
			return FMUL32BIT_LATENCY;
		case LLVM_IR_FDiv:
			// 64 bits -- 31 cycles
			// 32 bits -- 16 cycles
			return FDIV32BIT_LATENCY;
		case LLVM_IR_FCmp:
			return FCMP_LATENCY;
		default: 
			return 0;
	}
}
