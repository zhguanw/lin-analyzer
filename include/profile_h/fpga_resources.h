#ifndef FPGA_RESOURCES_H
#define FPGA_RESOURCES_H

#include <map>
#include <cmath>
#include <string>
#include <cstdint>
#include <iostream>
//#include "profile_h/lin-profile.h"
#include "profile_h/generic_func.h"
#include "profile_h/auxiliary.h"

/// Define FPGA chip we use
//#define __XILINX_ZYNQ7000_ZC702__
//#define __XILINX_VIRTEX7_VC707__

/// Array name --> Array partition number
typedef std::map<std::string, unsigned> arrayName2NumPartTy;
typedef std::map<std::string, unsigned> arrayName2NumPortPerPartTy;
/// Array name --> (Array partition number, Total Array Size in Byte, Word size in Byte) Map
/// If an array is completely partitioned, then its Array partition number is 0.
typedef std::map<std::string, std::pair<unsigned, std::pair<uint64_t, unsigned> > > arrayName2arrayConfigTy;
/// Array name --> memory inefficiency per partition
typedef std::map<std::string, float> arrayName2memEfficiencyTy;

class FPGA_ResClass {
public:
	FPGA_ResClass(): dsp(0), bram18k(0), ff(0), lut(0) {}
	FPGA_ResClass(unsigned n_dsp, unsigned n_bram, unsigned n_ff, unsigned n_lut) : dsp(n_dsp), bram18k(n_bram), ff(n_ff), lut(n_lut) {}
	virtual ~FPGA_ResClass() {}

	unsigned getDSP() const { return dsp; }
	unsigned getBRAM18K() const { return bram18k; }
	unsigned getFF() const { return ff; }
	unsigned getLUT() const { return lut; }

	void setDSP(unsigned num) { dsp = num; }
	void setBRAM18K(unsigned num) { bram18k = num; }
	void setFF(unsigned num) { ff = num; }
	void setLUT(unsigned num) { lut = num; }

private:
	unsigned dsp;
	unsigned bram18k;
	unsigned ff;
	unsigned lut;
};

class UsedFPGA_ResClass : public FPGA_ResClass {
public:
	///FIXME: The resources used by fadd/sub, fmul, fdiv here are also related to its bitwidth,
	///				But we assume all are 32-bit floating point operation, later we need to change here.
	UsedFPGA_ResClass() : num_fadd(0), num_fsub(0), num_fmul(0), num_fdiv(0) {}

	void increaseFadd() { updateUsedFadd_sub(); num_fadd++; }
	void increaseFsub() { updateUsedFadd_sub(); num_fsub++; }
	void increaseFmul() { updateUsedFmul(); num_fmul++; }
	void increaseFdiv() { updateUsedFdiv(); num_fdiv++; }
	void increaseArrayPart(std::string arrayName);

	void setArrayPartNum(std::string arrayName, unsigned partNum);

	unsigned getNumFadd() const { return num_fadd; }
	unsigned getNumFsub() const { return num_fsub; }
	unsigned getNumFmul() const { return num_fmul; }
	unsigned getNumFdiv() const { return num_fdiv; }
	unsigned getArrayPartitionNum(std::string arrayName);
	arrayName2NumPartTy getArrayName2NumPartMap() const;

	const unsigned faddsub_DSP = 2;
	const unsigned faddsub_FF = 205;
	const unsigned faddsub_LUT = 390;
	const unsigned fmul_DSP = 3;
	const unsigned fmul_FF = 143;
	const unsigned fmul_LUT = 321; 
	const unsigned fdiv_DSP = 0;
	const unsigned fdiv_FF = 761;
	const unsigned fdiv_LUT = 994;

private:
	void updateUsedFadd_sub();
	void updateUsedFmul();
	void updateUsedFdiv();

	unsigned num_fadd;
	unsigned num_fsub;
	unsigned num_fmul;
	unsigned num_fdiv;
	arrayName2NumPartTy array2numPart;
};

/// Currently, we focus on a FPGA chip on Xilinx ZC702 development board. The resource limitations are
/// defined inside the following class.
class UsedFPGA_Res_With_Constraint : public UsedFPGA_ResClass {
public:
	explicit UsedFPGA_Res_With_Constraint(arrayName2arrayConfigTy& arrayName2arrayCfg);
	// If an array partition only contains bits less than or equal to bram_sizeInbit_threshold,
	// then it will be implemented with distributed RAM instead of BRAM.
	const float bram_sizeInbit_threshold = 512;

	// Default is the available resource of Xilinx Zedboard or ZC702
	unsigned max_dsp = 220;
	unsigned max_bram18k = 280;
	unsigned max_ff = 106400;
	unsigned max_lut = 53200;

/*
#if defined(__XILINX_ZYNQ7000_ZC702__)
	const unsigned max_dsp = 220;
	const unsigned max_bram18k = 280;
	const unsigned max_ff = 106400;
	const unsigned max_lut = 53200;
#elif defined(__XILINX_VIRTEX7_VC707__)
	const unsigned max_dsp = 2800;
	const unsigned max_bram18k = 2060;
	const unsigned max_ff = 607200;
	const unsigned max_lut = 303600;
#else
	const unsigned max_dsp = 220;
	const unsigned max_bram18k = 280;
	const unsigned max_ff = 106400;
	const unsigned max_lut = 53200;
#endif
*/

	const unsigned size_BRAM18K_in_bit = 18432;
	bool check_resources_limitation(unsigned bram18k_usage, unsigned dsp_usage, unsigned ff_usage, unsigned lut_usage);
	arrayName2memEfficiencyTy getMemoryEfficiency() const;
	arrayName2NumPartTy getArrayName2Bram18k_used() const;
	unsigned getTotalBram18k_usage() const;
	unsigned getArrayWritePortNumPerPart(std::string original_array_name);
	unsigned getArrayPartNameWritePortNum(std::string array_name);

	/// try_to_occupy_one_*():
	///		return: false	-- resource constraint, need to wait;
	///						true	-- user successfully occupies one compute unit or memory read/write port.
	bool try_to_occupy_one_fadd();
	bool try_to_occupy_one_fsub();
	bool try_to_occupy_one_fmul();
	bool try_to_occupy_one_fdiv();
	bool try_to_occupy_one_read_port(std::string arrayName);
	bool try_to_occupy_one_write_port(std::string arrayName);

	void release_one_fadd();
	void release_one_fsub();
	void release_one_fmul();
	void release_one_fdiv();
	void release_one_read_port(std::string arrayName);
	void release_one_write_port(std::string arrayName);

	void set_floating_unit_threshold(unsigned fadd_limited, unsigned fsub_limited, unsigned fmul_limited, unsigned fdiv_limited);
	std::vector<std::string> get_constrained_fop_unit_type();

private:
	void setBRAM18K_usage();
	bool askForMore_Fadd();
	bool askForMore_Fsub();
	bool askForMore_Fmul();
	bool askForMore_Fdiv();

	arrayName2arrayConfigTy arrayName2arrayConfig;
	arrayName2NumPartTy arrayName2bram18k_used;
	arrayName2memEfficiencyTy arrayName2memEfficiency;
	arrayName2NumPartTy arrayPartName2readPort;
	arrayName2NumPartTy arrayPartName2writePort;
	unsigned bram18k_total_used;
	unsigned fadd_used;
	unsigned fsub_used;
	unsigned fmul_used;
	unsigned fdiv_used;
	arrayName2NumPartTy arrayPartName2readPort_used;
	arrayName2NumPartTy arrayPartName2writePort_used;
	arrayName2NumPortPerPartTy arrayName2NumWritePortPerPart;

	bool set_unit_threshold_flag;
	unsigned threshold_fadd_num;
	unsigned threshold_fsub_num;
	unsigned threshold_fmul_num;
	unsigned threshold_fdiv_num;

	bool dsp_resource_limited;
	bool fadd_limited_flag;
	bool fsub_limited_flag;
	bool fmul_limited_flag;
	bool fdiv_limited_flag;
};
#endif // End of FPGA_RESOURCES_H