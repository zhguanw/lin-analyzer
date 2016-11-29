#ifndef __XILINX_NODE_LATENCY_H__
#define __XILINX_NODE_LATENCY_H__

/*
#define	ADD_LATENCY       1
#define	SUB_LATENCY				1
#define MUL_LATENCY				2
#define DIV_LATENCY				2
#define AND_LATENCY				1
#define OR_LATENCY				1
#define XOR_LATENCY				1
#define SHL_LATENCY				1
#define BR_LATENCY				1
#define ICMP_LATENCY			1
#define IndexAdd_LATENCY	1
#define IndexSub_LATENCY	1

#define FADD_LATENCY			5
#define FSUB_LATENCY			5
#define FMUL_LATENCY			4
#define FDIV_LATENCY      6
#define FCMP_LATENCY      1
*/

#define LOAD_LATENCY	    2
#define STORE_LATENCY			1

#define ADDSUB_LATENCY		1
#define MUL32BIT_LATENCY	6
#define DIV32BIT_LATENCY  36

#define FADDSUB32BIT_LATENCY 5
#define FMUL32BIT_LATENCY 4
#define FDIV32BIT_LATENCY 16
#define FCMP_LATENCY			1

#define BRAM_R_PORTS      2
#define	BRAM_W_PORTS			1

#endif
