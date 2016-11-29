#include "profile_h/fpga_resources.h"

void UsedFPGA_ResClass::updateUsedFadd_sub() {
	unsigned old_dsp = getDSP(); 
	unsigned old_ff = getFF(); 
	unsigned old_lut = getLUT(); 
	setDSP(old_dsp + faddsub_DSP);
	setFF(old_ff + faddsub_FF);
	setLUT(old_lut + faddsub_LUT);
}

void UsedFPGA_ResClass::updateUsedFmul() {
	unsigned old_dsp = getDSP();
	unsigned old_ff = getFF();
	unsigned old_lut = getLUT();
	setDSP(old_dsp + fmul_DSP);
	setFF(old_ff + fmul_FF);
	setLUT(old_lut + fmul_LUT);
}

void UsedFPGA_ResClass::updateUsedFdiv() {
	unsigned old_dsp = getDSP();
	unsigned old_ff = getFF();
	unsigned old_lut = getLUT();
	setDSP(old_dsp + fdiv_DSP);
	setFF(old_ff + fdiv_FF);
	setLUT(old_lut + fdiv_LUT);
}

void UsedFPGA_ResClass::increaseArrayPart(std::string arrayName) {
	arrayName2NumPartTy::iterator it;
	it = array2numPart.find(arrayName);
	if ( it == array2numPart.end() ) {
		array2numPart.insert(std::make_pair(arrayName, 1));
	}
	else {
		array2numPart[arrayName]++;
	}
}

unsigned UsedFPGA_ResClass::getArrayPartitionNum(std::string arrayName) { 
	return array2numPart[arrayName];
}

void UsedFPGA_ResClass::setArrayPartNum(std::string arrayName, unsigned partNum) {
	arrayName2NumPartTy::iterator it;
	it = array2numPart.find(arrayName);
	if (it == array2numPart.end()) {
		array2numPart.insert(std::make_pair(arrayName, partNum));
	}
	else {
		array2numPart[arrayName] = partNum;
	}
}

arrayName2NumPartTy UsedFPGA_ResClass::getArrayName2NumPartMap() const {
	return array2numPart;
}

UsedFPGA_Res_With_Constraint::UsedFPGA_Res_With_Constraint(arrayName2arrayConfigTy& arrayName2arrayCfg) {
	arrayName2arrayConfig = arrayName2arrayCfg;
	bram18k_total_used = 0;
	fadd_used = 0;
	fsub_used = 0;
	fmul_used = 0;
	fdiv_used = 0;
	set_unit_threshold_flag = false;
	threshold_fadd_num = INFINITE_HARDWARE;
	threshold_fsub_num = INFINITE_HARDWARE;
	threshold_fmul_num = INFINITE_HARDWARE;
	threshold_fdiv_num = INFINITE_HARDWARE;

	dsp_resource_limited = false;
	fadd_limited_flag = false;
	fsub_limited_flag = false;
	fmul_limited_flag = false;
	fdiv_limited_flag = false;

	if (disable_fp_unit_threshold == true) {
		max_dsp = INFINITE_HARDWARE;
		max_bram18k = INFINITE_HARDWARE;
		max_ff = INFINITE_HARDWARE;
		max_lut = INFINITE_HARDWARE;
	}

	if (target_vc707 == true && disable_fp_unit_threshold == false) {
		max_dsp = 2800;
		max_bram18k = 2060;
		max_ff = 607200;
		max_lut = 303600;
	}

	setBRAM18K_usage();
}

/// check_resources_limitation():
///		return: false -- has no constraint
///						true  -- has constraints
bool UsedFPGA_Res_With_Constraint::check_resources_limitation(unsigned bram18k_usage, unsigned dsp_usage, unsigned ff_usage, unsigned lut_usage) {
	bool dsp_flag = false;
	bool bram18k_flag = false;
	bool ff_flag = false;
	bool lut_flag = false;

	if (dsp_usage > max_dsp) {
#ifdef DEBUG_INFORMATION_OUTPUT
		std::cout << "DEBUG-INFO: [fpga_resources-utilization] DSP usage exceeds!" << std::endl;
#endif
		dsp_flag = true;
	}

	if (bram18k_usage > max_bram18k) {
#ifdef DEBUG_INFORMATION_OUTPUT
		std::cout << "DEBUG-INFO: [fpga_resources-utilization] BRAM18K usage exceeds!" << std::endl;
#endif
		bram18k_flag = true;
	}
	
	if (ff_usage > max_ff) {
#ifdef DEBUG_INFORMATION_OUTPUT
		std::cout << "DEBUG-INFO: [fpga_resources-utilization] FF usage exceeds!" << std::endl;
#endif
		ff_flag = true;
	}
	
	if (lut_usage > max_lut) {
#ifdef DEBUG_INFORMATION_OUTPUT
		std::cout << "DEBUG-INFO: [fpga_resources-utilization] LUT usage exceeds!" << std::endl;
#endif
		lut_flag = true;
	}

	return (dsp_flag | bram18k_flag | ff_flag | lut_flag);
}

/// askForMore_Fadd():
///		return: false -- it refuses to add more fadd resources; 
///						true	-- it succeeds in adding 1 more fadd resources.
bool UsedFPGA_Res_With_Constraint::askForMore_Fadd() {
	unsigned new_bram18k = getBRAM18K();
	unsigned new_dsp = getDSP() + faddsub_DSP;
	unsigned new_ff = getFF() + faddsub_FF;
	unsigned new_lut = getLUT() + faddsub_LUT;
	bool constraint_flag = check_resources_limitation(new_bram18k, new_dsp, new_ff, new_lut);
	/// If constraint_flag = false, it means no constraint.
	if (constraint_flag == false) {
		increaseFadd();
	}

	/// If it exceeds fpga resources, we do nothing.

	return !constraint_flag;
}

/// askForMore_Fsub():
///		return: false -- it refuses to add more fsub resources; 
///						true	-- it succeeds in adding 1 more fsub resources.
bool UsedFPGA_Res_With_Constraint::askForMore_Fsub() {
	unsigned new_bram18k = getBRAM18K();
	unsigned new_dsp = getDSP() + faddsub_DSP;
	unsigned new_ff = getFF() + faddsub_FF;
	unsigned new_lut = getLUT() + faddsub_LUT;
	bool constraint_flag = check_resources_limitation(new_bram18k, new_dsp, new_ff, new_lut);
	/// If constraint_flag = false, it means no constraint.
	if (constraint_flag == false) {
		increaseFsub();
	}

	/// If it exceeds fpga resources, we do nothing.

	return !constraint_flag;
}

/// askForMore_Fmul():
///		return: false -- it refuses to add more fmul resources; 
///						true	-- it succeeds in adding 1 more fmul resources.
bool UsedFPGA_Res_With_Constraint::askForMore_Fmul() {
	unsigned new_bram18k = getBRAM18K();
	unsigned new_dsp = getDSP() + fmul_DSP;
	unsigned new_ff = getFF() + fmul_FF;
	unsigned new_lut = getLUT() + fmul_LUT;
	bool constraint_flag = check_resources_limitation(new_bram18k, new_dsp, new_ff, new_lut);
	/// If constraint_flag = false, it means no constraint.
	if (constraint_flag == false) {
		increaseFmul();
	}

	/// If it exceeds fpga resources, we do nothing.

	return !constraint_flag;
}

/// askForMore_Fdiv():
///		return: false -- it refuses to add more fdiv resources; 
///						true	-- it succeeds in adding 1 more fdiv resources.
bool UsedFPGA_Res_With_Constraint::askForMore_Fdiv() {
	unsigned new_bram18k = getBRAM18K();
	unsigned new_dsp = getDSP() + fdiv_DSP;
	unsigned new_ff = getFF() + fdiv_FF;
	unsigned new_lut = getLUT() + fdiv_LUT;
	bool constraint_flag = check_resources_limitation(new_bram18k, new_dsp, new_ff, new_lut);
	/// If constraint_flag = false, it means no constraint.
	if (constraint_flag == false) {
		increaseFdiv();
	}

	/// If it exceeds fpga resources, we do nothing.

	return !constraint_flag;
}

void UsedFPGA_Res_With_Constraint::setBRAM18K_usage() {
	arrayPartName2readPort.clear();
	arrayPartName2readPort_used.clear();
	arrayPartName2writePort.clear();
	arrayPartName2writePort_used.clear();
	arrayName2arrayConfigTy::iterator it = arrayName2arrayConfig.begin();
	arrayName2arrayConfigTy::iterator ie = arrayName2arrayConfig.end();

	for (; it != ie; ++it) {
		std::string arrayName = it->first;
		unsigned arrayPartitionNum = it->second.first;
		uint64_t arrayTotalSizeInByte = it->second.second.first;
		unsigned arrayWordSizeInByte = it->second.second.second;
		if (arrayPartitionNum !=0 ) {
			float sizeInBit_perPartition = (float)arrayTotalSizeInByte * 8.0 / (float)arrayPartitionNum;
			float numBram18k_perPartition_float = 0.0f;
			unsigned numBram18k_perPartition = 0;
			float memEfficiencyPerPartition = 0.0f;
			if (sizeInBit_perPartition <= bram_sizeInbit_threshold) {
				numBram18k_perPartition_float = 0.0f;
				numBram18k_perPartition = 0;
				memEfficiencyPerPartition = 100.0f;
			}
			else {
				numBram18k_perPartition_float = ceil(sizeInBit_perPartition / (float)size_BRAM18K_in_bit);
				numBram18k_perPartition = next_power_of_two((unsigned)numBram18k_perPartition_float);
				memEfficiencyPerPartition = sizeInBit_perPartition / ((float)numBram18k_perPartition * (float)size_BRAM18K_in_bit);
			}
			unsigned bram18k_used_per_array = arrayPartitionNum * numBram18k_perPartition;
			bram18k_total_used += bram18k_used_per_array;
			arrayName2bram18k_used.insert(std::make_pair(arrayName, bram18k_used_per_array));
			arrayName2memEfficiency.insert(std::make_pair(arrayName, memEfficiencyPerPartition));
		}
		else {
			/// Completely partition this array into all registers. No bram used
			/// FIXME: When completely partition an array, it will use more FF and LUT resources, later if we want to 
			///				 support other resources usage, we need to add it below.
			arrayName2bram18k_used.insert(std::make_pair(arrayName, 0));
			arrayName2memEfficiency.insert(std::make_pair(arrayName, 0));
		}

		setArrayPartNum(arrayName, arrayPartitionNum);
		arrayName2NumWritePortPerPart[arrayName] = WRITE_PORT_PER_PARTITION;
	}

	if (bram18k_total_used > max_bram18k) {
		std::cout << "DEBUG-INFO: [fpga_resources-utilization] BRAM18K resource is not enough, please use smaller partition factor or use other larger FPGA chips!" << std::endl;
		std::cout << "DEBUG-INFO: [fpga_resources-partition] Cannot partition arrays with their partition factor, use factor=1 as default" << std::endl;
		unsigned bram18k_total_without_part = 0;
		arrayName2bram18k_used.clear();
		arrayName2memEfficiency.clear();

		arrayName2arrayConfigTy::iterator it = arrayName2arrayConfig.begin();
		arrayName2arrayConfigTy::iterator ie = arrayName2arrayConfig.end();
		for (; it != ie; ++it) {
			std::string arrayName = it->first;
			unsigned arrayPartitionNum = it->second.first;
			uint64_t arrayTotalSizeInByte = it->second.second.second;
			if (arrayPartitionNum != 0) {
				float sizeInBit_perPartition = (float)arrayTotalSizeInByte * 8.0;
				float numBram18k_perPartition_float = ceil(sizeInBit_perPartition / (float)size_BRAM18K_in_bit);
				unsigned numBram18k_perPartition = next_power_of_two((unsigned)numBram18k_perPartition_float);
				float memEfficiencyPerPartition = sizeInBit_perPartition / ((float)numBram18k_perPartition * (float)size_BRAM18K_in_bit);
				unsigned bram18k_used_per_array = numBram18k_perPartition;
				bram18k_total_without_part += bram18k_used_per_array;
				arrayName2bram18k_used.insert(std::make_pair(arrayName, bram18k_used_per_array));
				arrayName2memEfficiency.insert(std::make_pair(arrayName, memEfficiencyPerPartition));
				arrayPartName2readPort.insert(std::make_pair(arrayName, READ_PORT_PER_PARTITION));
				arrayPartName2writePort.insert(std::make_pair(arrayName, WRITE_PORT_PER_PARTITION));
			}
			else {
				/// Completely partition this array into all registers. No bram used
				arrayName2bram18k_used.insert(std::make_pair(arrayName+"-register", 0));
				arrayName2memEfficiency.insert(std::make_pair(arrayName + "-register", 0));
				arrayPartName2readPort.insert(std::make_pair(arrayName, INFINITE_HARDWARE));
				arrayPartName2writePort.insert(std::make_pair(arrayName, INFINITE_HARDWARE));
			}

			arrayPartName2readPort_used.insert(std::make_pair(arrayName, 0));
			arrayPartName2writePort_used.insert(std::make_pair(arrayName, 0));
		}

		if (bram18k_total_without_part > max_bram18k) {
			std::cout << "DEBUG-INFO: [fpga_resources-utilization] Size of Arrays in this application are too large for this FPGA chip, please select larger FPGA chips!" << std::endl;
		}
		setBRAM18K(bram18k_total_without_part);
	}
	else {
		setBRAM18K(bram18k_total_used);
		arrayPartName2readPort.clear();
		arrayPartName2writePort.clear();
		arrayName2arrayConfigTy::iterator it = arrayName2arrayConfig.begin();
		arrayName2arrayConfigTy::iterator ie = arrayName2arrayConfig.end();
		for (; it != ie; ++it) {
			std::string arrayName = it->first;
			unsigned arrayPartitionNum = it->second.first;
			uint64_t arrayTotalSizeInByte = it->second.second.second;
			if (arrayPartitionNum > 1) {
				for (unsigned i = 0; i < arrayPartitionNum; i++) {
					std::string array_part_name = arrayName + "-" +std::to_string(i);
					arrayPartName2readPort.insert(std::make_pair(array_part_name, READ_PORT_PER_PARTITION));
					arrayPartName2writePort.insert(std::make_pair(array_part_name, WRITE_PORT_PER_PARTITION));
					arrayPartName2readPort_used.insert(std::make_pair(array_part_name, 0));
					arrayPartName2writePort_used.insert(std::make_pair(array_part_name, 0));
				}
			}
			else if (arrayPartitionNum == 0) {
				arrayPartName2readPort.insert(std::make_pair(arrayName, INFINITE_HARDWARE));
				arrayPartName2writePort.insert(std::make_pair(arrayName, INFINITE_HARDWARE));
				arrayPartName2readPort_used.insert(std::make_pair(arrayName, 0));
				arrayPartName2writePort_used.insert(std::make_pair(arrayName, 0));
			}
			else {
				arrayPartName2readPort.insert(std::make_pair(arrayName, READ_PORT_PER_PARTITION));
				arrayPartName2writePort.insert(std::make_pair(arrayName, WRITE_PORT_PER_PARTITION));
				arrayPartName2readPort_used.insert(std::make_pair(arrayName, 0));
				arrayPartName2writePort_used.insert(std::make_pair(arrayName, 0));
			}

		}
	}
}

arrayName2memEfficiencyTy UsedFPGA_Res_With_Constraint::getMemoryEfficiency() const {
	return arrayName2memEfficiency;
}

arrayName2NumPartTy UsedFPGA_Res_With_Constraint::getArrayName2Bram18k_used() const {
	return arrayName2bram18k_used;
}

unsigned UsedFPGA_Res_With_Constraint::getTotalBram18k_usage() const {
	return getBRAM18K();
}

unsigned UsedFPGA_Res_With_Constraint::getArrayWritePortNumPerPart(std::string original_array_name) {
	arrayName2NumPortPerPartTy::iterator it = arrayName2NumWritePortPerPart.find(original_array_name);
	assert(it != arrayName2NumWritePortPerPart.end() && "Error: Input of getArrayWritePortNumPerPart() is original array name!\n");
	return arrayName2NumWritePortPerPart.at(original_array_name);
}

unsigned UsedFPGA_Res_With_Constraint::getArrayPartNameWritePortNum(std::string array_name) {
	arrayName2NumPartTy::iterator it = arrayPartName2writePort.find(array_name);
	assert(it != arrayPartName2writePort.end() && "Error: Input of getArrayPartNameWritePortNum() cannot be found in arrayPartName2writePort!\n");
	return arrayPartName2writePort.at(array_name);
}

/// try_to_occupy_one_fadd():
///		return: false	-- resource constraint, need to wait;
///						true	-- user successfully occupies one compute unit.
bool UsedFPGA_Res_With_Constraint::try_to_occupy_one_fadd() {
	unsigned num_fadd_available = getNumFadd();
	if (fadd_used < num_fadd_available) {
		fadd_used++;
		return true;
	}
	else {
		if ((set_unit_threshold_flag == true) && (num_fadd_available != 0)) {
			if (num_fadd_available >= threshold_fadd_num) {
				return false;
			}
		}

		// can_add_more_compute_unit:
		//		false -- cannot add more compute unit, user has to wait for compute units;
		//		true  -- it succeeds in adding 1 more compute unit.
		bool can_add_more_compute_unit = askForMore_Fadd();
		if (can_add_more_compute_unit == true) {
			fadd_used++;
		}
		return can_add_more_compute_unit;
	}
}

/// try_to_occupy_one_fsub():
///		return: false	-- resource constraint, need to wait;
///						true	-- user successfully occupies one compute unit.
bool UsedFPGA_Res_With_Constraint::try_to_occupy_one_fsub() {
	unsigned num_fsub_available = getNumFsub();
	if (fsub_used < num_fsub_available) {
		fsub_used++;
		return true;
	}
	else {
		if ((set_unit_threshold_flag == true) && (num_fsub_available != 0)) {
			if (num_fsub_available >= threshold_fsub_num) {
				return false;
			}
		}

		// can_add_more_compute_unit:
		//		false -- cannot add more compute unit, user has to wait for compute units;
		//		true  -- it succeeds in adding 1 more compute unit.
		bool can_add_more_compute_unit = askForMore_Fsub();
		if (can_add_more_compute_unit == true) {
			fsub_used++;
		}
		return can_add_more_compute_unit;
	}
}

/// try_to_occupy_one_fmul():
///		return: false	-- resource constraint, need to wait;
///						true	-- user successfully occupies one compute unit.
bool UsedFPGA_Res_With_Constraint::try_to_occupy_one_fmul() {
	unsigned num_fmul_available = getNumFmul();
	if (fmul_used < num_fmul_available) {
		fmul_used++;
		return true;
	}
	else {
		if ((set_unit_threshold_flag == true) && (num_fmul_available != 0)) {
			if (num_fmul_available >= threshold_fmul_num) {
				return false;
			}
		}

		// can_add_more_compute_unit:
		//		false -- cannot add more compute unit, user has to wait for compute units;
		//		true  -- it succeeds in adding 1 more compute unit.
		bool can_add_more_compute_unit = askForMore_Fmul();
		if (can_add_more_compute_unit == true) {
			fmul_used++;
		}
		return can_add_more_compute_unit;
	}
}

/// try_to_occupy_one_fdiv():
///		return: false	-- resource constraint, need to wait;
///						true	-- user successfully occupies one compute unit.
bool UsedFPGA_Res_With_Constraint::try_to_occupy_one_fdiv() {
	unsigned num_fdiv_available = getNumFdiv();
	if (fdiv_used < num_fdiv_available) {
		fdiv_used++;
		return true;
	}
	else {
		if ((set_unit_threshold_flag == true) && (num_fdiv_available != 0)) {
			if (num_fdiv_available >= threshold_fdiv_num) {
				return false;
			}
		}

		// can_add_more_compute_unit:
		//		false -- cannot add more compute unit, user has to wait for compute units;
		//		true  -- it succeeds in adding 1 more compute unit.
		bool can_add_more_compute_unit = askForMore_Fdiv();
		if (can_add_more_compute_unit == true) {
			fdiv_used++;
		}
		return can_add_more_compute_unit;
	}
}

/// try_to_occupy_one_read_port():
///		return: false	-- memory read port constraint, need to wait;
///						true	-- user successfully occupies one memory read port.
bool UsedFPGA_Res_With_Constraint::try_to_occupy_one_read_port(std::string arrayName) {
	arrayName2NumPartTy::iterator it = arrayPartName2readPort.find(arrayName);
	arrayName2NumPartTy::iterator it_used = arrayPartName2readPort_used.find(arrayName);
	assert((it != arrayPartName2readPort.end()) && "Error: Array name is not inside arrayPartName2readPort!\n");
	assert((it_used != arrayPartName2readPort_used.end()) && "Error: Array name is not inside arrayPartName2readPort_used!\n");
	unsigned read_ports_available = it->second;
	unsigned read_ports_used = it_used->second;
	if (read_ports_used >= read_ports_available) {
		return false;
	}
	else {
		it_used->second++;
		return true;
	}
}

/// try_to_occupy_one_write_port():
///		return: false	-- memory write port constraint, need to wait;
///						true	-- user successfully occupies one memory write port.
bool UsedFPGA_Res_With_Constraint::try_to_occupy_one_write_port(std::string arrayName) {
	arrayName2NumPartTy::iterator it = arrayPartName2writePort.find(arrayName);
	arrayName2NumPartTy::iterator it_used = arrayPartName2writePort_used.find(arrayName);
	assert((it != arrayPartName2writePort.end()) && "Error: Array name is not inside arrayPartName2writePort!\n");
	assert((it_used != arrayPartName2writePort_used.end()) && "Error: Array name is not inside arrayPartName2writePort_used!\n");
	unsigned write_ports_available = it->second;
	unsigned write_ports_used = it_used->second;
	if (write_ports_used >= write_ports_available) {
		if ((write_ports_available<MAX_WRITE_PORT_PER_PARTITION) && (enable_rw_rw_memory == true)) {
			it->second++;
			it_used->second++;
			std::size_t found = arrayName.find("-");
			if (found == string::npos) {
				arrayName2NumWritePortPerPart[arrayName]++;
			}
			else {
				std::string original_arrayName = arrayName.substr(0, found);
				arrayName2NumWritePortPerPart[original_arrayName]++;
			}
			
			return true;
		}
		else {
			return false;
		}
	}
	else {
		it_used->second++;
		return true;
	}
}

void UsedFPGA_Res_With_Constraint::release_one_fadd() {
	assert( (fadd_used!=0) && "Error: fadd_used is already 0, cannot be deducted again!\n" );
	fadd_used--;
}

void UsedFPGA_Res_With_Constraint::release_one_fsub() {
	assert((fsub_used != 0) && "Error: fsub_used is already 0, cannot be deducted again!\n");
	fsub_used--;
}

void UsedFPGA_Res_With_Constraint::release_one_fmul() {
	assert((fmul_used != 0) && "Error: fmul_used is already 0, cannot be deducted again!\n");
	fmul_used--;
}

void UsedFPGA_Res_With_Constraint::release_one_fdiv() {
	assert((fdiv_used != 0) && "Error: fdiv_used is already 0, cannot be deducted again!\n");
	fdiv_used--;
}

void UsedFPGA_Res_With_Constraint::release_one_read_port(std::string arrayName) {
	arrayName2NumPartTy::iterator it = arrayPartName2readPort_used.find(arrayName);
	assert((it != arrayPartName2readPort_used.end()) && "Error: Array name is not inside arrayPartName2readPort_used!\n");
	unsigned read_ports_used = it->second;
	assert((read_ports_used != 0) && "Errors: Read ports of this array are not used at all, no need to deduct!\n");
	it->second--;
}

void UsedFPGA_Res_With_Constraint::release_one_write_port(std::string arrayName) {
	arrayName2NumPartTy::iterator it = arrayPartName2writePort_used.find(arrayName);
	assert((it != arrayPartName2writePort_used.end()) && "Error: Array name is not inside arrayPartName2writePort_used!\n");
	unsigned write_ports_used = it->second;
	assert((write_ports_used != 0) && "Errors: Write ports of this array are not used at all, no need to deduct!\n");
	it->second--;
}

void UsedFPGA_Res_With_Constraint::set_floating_unit_threshold(unsigned fadd_limited, unsigned fsub_limited, unsigned fmul_limited, unsigned fdiv_limited) {
	set_unit_threshold_flag = true;
	unsigned total_dsp_may = (fadd_limited + fsub_limited) * faddsub_DSP + fmul_limited * fmul_DSP + fdiv_limited * fdiv_DSP;
	if (total_dsp_may > max_dsp) {
		float scale = (float)total_dsp_may / (float)max_dsp;
		fadd_limited = (unsigned) std::ceil((float)fadd_limited / scale);
		fsub_limited = (unsigned)std::ceil((float)fsub_limited / scale);
		fmul_limited = (unsigned)std::ceil((float)fmul_limited / scale);
		fdiv_limited = (unsigned)std::ceil((float)fdiv_limited / scale);

		std::vector<unsigned> find_max;
		dsp_resource_limited = true;
		find_max.push_back(fadd_limited*faddsub_DSP);
		find_max.push_back(fsub_limited*faddsub_DSP);
		find_max.push_back(fmul_limited*fmul_DSP);
		find_max.push_back(fdiv_limited*fdiv_DSP);
		unsigned max_val = *std::max_element(find_max.begin(), find_max.end());
		if (max_val == (fadd_limited*faddsub_DSP)) {
			fadd_limited_flag = true;
		}

		if (max_val == (fsub_limited*faddsub_DSP)) {
			fsub_limited_flag = true;
		}

		if (max_val == (fmul_limited*fmul_DSP)) {
			fmul_limited_flag = true;
		}

		if (max_val == (fdiv_limited*fdiv_DSP)) {
			fdiv_limited_flag = true;
		}

	}
	else {
		fadd_limited = (fadd_limited == 0) ? INFINITE_HARDWARE : fadd_limited;
		fsub_limited = (fsub_limited == 0) ? INFINITE_HARDWARE : fsub_limited;
		fmul_limited = (fmul_limited == 0) ? INFINITE_HARDWARE : fmul_limited;
		fdiv_limited = (fdiv_limited == 0) ? INFINITE_HARDWARE : fdiv_limited;
	}

	threshold_fadd_num = fadd_limited;
	threshold_fsub_num = fsub_limited;
	threshold_fmul_num = fmul_limited;
	threshold_fdiv_num = fdiv_limited;
}

std::vector<std::string> UsedFPGA_Res_With_Constraint::get_constrained_fop_unit_type() {
	std::vector<std::string> constrained_types;
	if (dsp_resource_limited == true) {
		if (fadd_limited_flag == true) {
			constrained_types.push_back("FADD_UNIT_CONSTRAINT");
		}
		if (fsub_limited_flag == true) {
			constrained_types.push_back("FSUB_UNIT_CONSTRAINT");
		}
		if (fmul_limited_flag == true) {
			constrained_types.push_back("FMUL_UNIT_CONSTRAINT");
		}
		if (fdiv_limited_flag == true) {
			constrained_types.push_back("FDIV_UNIT_CONSTRAINT");
		}
	}
	else {
		constrained_types.push_back("NO_FOP_UNIT_CONSTRAINT");
	}

	return constrained_types;
}