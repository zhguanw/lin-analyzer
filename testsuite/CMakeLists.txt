cmake_minimum_required(VERSION 2.8)

set(TEST_POLYBENCH_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TEST_POLYBENCH_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR})
MESSAGE(STATUS "TEST_POLYBENCH_SRC_DIR: ${TEST_POLYBENCH_SRC_DIR}")

SET(LLVM_BUILD_BIN_DIR ${LLVM_BINARY_DIR}/bin)
MESSAGE(STATUS "LLVM_BUILD_BIN_DIR: ${LLVM_BUILD_BIN_DIR}")

#set up LLVM binary
set(CLANG ${LLVM_BUILD_BIN_DIR}/clang)
#set(CLANG ${LLVM_BINARY_DIR}/bin/clang)
set(OPT ${LLVM_BUILD_BIN_DIR}/opt)
#set(OPT ${LLVM_BINARY_DIR}/bin/opt)
set(LLC ${LLVM_BUILD_BIN_DIR}/llc)
#set(LLVM ${LLVM_BINARY_DIR}/bin/llc)
set(LLI ${LLVM_BUILD_BIN_DIR}/lli)
#set(LLI ${LLVM_BINARY_DIR}/bin/lli)
set(LLVM_LINK ${LLVM_BUILD_BIN_DIR}/llvm-link)
#set(LLVM_LINK ${LLVM_BINARY_DIR}/bin/llvm-link)

#set up lin-analyzer binary in cmake/ folder
set(LIN_PROFILER ${LLVM_BUILD_BIN_DIR}/lin-analyzer)

#Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${TEST_POLYBENCH_SRC_DIR}/cmake")
MESSAGE(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

#set(BENCHSUITE_NAME "polybench")
add_custom_target(generating_llvm_ir 
	COMMENT "Compiling all polybench benchmarks into llvm ir"
	)

add_custom_target(generating_linked_llvm_ir
	COMMENT "Linking all llvm IRs into a single llvm IR"
	)

add_custom_target(clean_temp_bin_bc_files
	COMMENT "Clean all llvm temp IRs"
	)

add_custom_target(run_polybench
	COMMENT "Running polybench with lin-analyzer"
	)

add_custom_target(polybench
	COMMENT "Finish Polybench testing"
	)

add_dependencies(generating_llvm_ir clean_temp_bin_bc_files)
add_dependencies(generating_linked_llvm_ir generating_llvm_ir)
add_dependencies(run_polybench generating_linked_llvm_ir lin-analyzer)
add_dependencies(polybench run_polybench)

macro(compile_polybench test_file kernel_func postfix)
	set(TEST_NAME ${test_file})
	set(KERNEL_FUNC_NAME ${kernel_func})

	set(TESTSUITE_SRC_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
	set(TESTSUITE_BIN_ROOT ${CMAKE_CURRENT_BINARY_DIR})

	set(TESTSUITE_SRC_BC_DIR ${TESTSUITE_BIN_ROOT})

	set(INPUT_SRC "${TESTSUITE_SRC_ROOT}/${test_file}.${postfix}")
	set(OUTPUT_BC "${TESTSUITE_SRC_BC_DIR}/${test_file}.bc")
	
	#if (TESTSUITE_SRC_ROOT AND TESTSUITE_BIN_ROOT)
		#MESSAGE(STATUS "TESTSUITE_SRC_ROOT: " ${TESTSUITE_SRC_ROOT})
		#MESSAGE(STATUS "TESTSUITE_BIN_ROOT: " ${TESTSUITE_BIN_ROOT})
	#endif (TESTSUITE_SRC_ROOT AND TESTSUITE_BIN_ROOT)
	#MESSAGE(STATUS "INPUT_SRC: ${INPUT_SRC}")
	#MESSAGE(STATUS "OUTPUT_BC: ${OUTPUT_BC}")	
	add_custom_command(OUTPUT ${OUTPUT_BC}
		COMMAND ${CLANG} -g -O1 -emit-llvm
			-c ${INPUT_SRC} -o ${OUTPUT_BC}
		DEPENDS ${INPUT_SRC} ${CLANG}
		WORKING_DIRECTORY ${TESTSUITE_SRC_BC_DIR}
		COMMENT "Compiling ${test_file}.${postfix} to ${OUTPUT_BC}" 
	)
	add_custom_target(${test_file}_${KERNEL_FUNC_NAME}_ir DEPENDS ${OUTPUT_BC})
	add_dependencies(generating_llvm_ir ${test_file}_${KERNEL_FUNC_NAME}_ir)
	
endmacro(compile_polybench)

macro(link_polybench kernelName bc_files)
	get_filename_component(kernelNameWOpostfix ${kernelName} NAME_WE)
	set(KERNEL_NAME ${kernelNameWOpostfix})
	set(TESTSUITE_LINK_SRC_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
	set(TESTSUITE_LINK_BIN_ROOT ${CMAKE_CURRENT_BINARY_DIR})
	#MESSAGE(STATUS "kernel_func name: ${KERNEL_NAME}")
	set(KERNEL_BC "${KERNEL_NAME}.bc")
	set(ABSTRACT_KERNEL_BC "abstract_kernel.bc")
	set(MAIN_BC "main.bc")
	set(LINKED_BC "${TESTSUITE_LINK_BIN_ROOT}/${KERNEL_NAME}_linked.bc")
	set(OPT_BC "${TESTSUITE_LINK_BIN_ROOT}/${KERNEL_NAME}_linked_opt.bc")
	set(OPT_FLAGS -mem2reg -instnamer -lcssa -indvars)
	
	MESSAGE(STATUS "LINKED_BC: ${LINKED_BC}")
	add_custom_command(OUTPUT ${LINKED_BC}
		#COMMAND ${LLVM_LINK} ${bc_files} -o ${LINKED_BC}
		COMMAND ${LLVM_LINK} ${MAIN_BC} ${ABSTRACT_KERNEL_BC} ${KERNEL_BC} -o ${LINKED_BC}
		DEPENDS ${bc_files} ${LLVM_LINK}
		WORKING_DIRECTORY ${TESTSUITE_LINK_BIN_ROOT}
		COMMENT "linking all *.bc into ${KERNEL_NAME}_linked.bc"
	)
	add_custom_command(OUTPUT ${OPT_BC}
		COMMAND ${OPT} ${OPT_FLAGS} ${LINKED_BC} -o ${OPT_BC}
		DEPENDS ${LINKED_BC} ${OPT}
		WORKING_DIRECTORY ${TESTSUITE_LINK_BIN_ROOT}
		COMMENT "opting the ${KERNEL_NAME}_linked.bc into ${KERNEL_NAME}_linked_opt.bc"
	)
	add_custom_target(${KERNEL_NAME}_opt_ir DEPENDS ${OPT_BC})
	add_dependencies(generating_linked_llvm_ir ${KERNEL_NAME}_opt_ir generating_llvm_ir)
endmacro(link_polybench)

#macro(remove_bin_bc_files kernelName)
#	file(GLOB temp_bin_bc_files RELATIVE "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/*.bc")
#	foreach(temp_file ${temp_bin_bc_files})
#		file(REMOVE ${temp_file})
#	endforeach(temp_file)
#	add_custom_target(${kernelName}_clean ${temp_bin_bc_files})
#	add_dependencies(clean_temp_bin_bc_files ${kernelName}_clean)
#endmacro(remove_bin_bc_files)

macro(test_polybench kernelName)
        set(BC_PATH "${CMAKE_CURRENT_BINARY_DIR}/")
        get_filename_component(NameWOpostfix ${kernelName} NAME_WE)
        set(Kernel_Name "${NameWOpostfix}")
	set(INPUT_BC "${CMAKE_CURRENT_BINARY_DIR}/${Kernel_Name}_linked_opt.bc")

	set(DYNAMIC_TRACE "${CMAKE_CURRENT_BINARY_DIR}/dynamic_trace.gz")
	set(INSTRUMENTED_BC "${Kernel_Name}_linked_opt_trace.bc")
	
	add_custom_command(OUTPUT ${DYNAMIC_TRACE} ${INSTRUMENTED_BC}
		COMMAND ${LIN_PROFILER} ${INPUT_BC} ${BC_PATH} ${Kernel_Name}
		DEPENDS ${INPUT_BC} ${LIN_PROFILER}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Testing ${Kernel_Name} with lin-analyzer"
	)
	add_custom_target(${Kernel_Name}_lin-analyzer Depends ${DYNAMIC_TRACE} ${INSTRUMENTED_BC})
	add_dependencies(run_polybench ${Kernel_Name}_lin-analyzer)
endmacro(test_polybench)

add_subdirectory(Ecobench)

