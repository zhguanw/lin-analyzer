#ifndef DYNAMIC_DATAPATH_H
#define DYNAMIC_DATAPATH_H

#include <string>
#include "profile_h/BaseDatapath.h"
#include "profile_h/lin-profile.h"

class DynamicDatapath : public BaseDatapath {
public:
	DynamicDatapath(std::string kernel_name, std::string trace_file_name, std::string input_path, std::string lp_name, unsigned lp_level, unsigned lp_unroll_factor, bool enable_pipelining, unsigned IL_asap);
	~DynamicDatapath();

	unsigned getIL_asap_ii() const;
private:
	std::string target_loop_name;
	unsigned target_loop_level;
	unsigned target_lp_level_unroll_factor;
	unsigned IL_asap_ii;
};

#endif // End of DYNAMIC_DATAPATH_H