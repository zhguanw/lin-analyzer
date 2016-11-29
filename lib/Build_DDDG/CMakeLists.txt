include_directories( ${ZLIB_INCLUDE_DIRS} )

add_llvm_library(BuildDDDGlib
	DynamicDatapath.cpp
	BaseDatapath.cpp
	DDDG.cpp
	SlotTracker.cpp
	TraceFunctions.cpp
	file_func.cpp
	Registers.cpp
	opcode_func.cpp
	generic_func.cpp
	fpga_resources.cpp
	
	LINK_LIBS
	${ZLIB_LIBRARY}
	)
	
#target_link_libraries(BuildDDDGlib ${ZLIB_LIBRARY})