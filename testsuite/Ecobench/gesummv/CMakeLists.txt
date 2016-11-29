#compile_polybench("gesummv" "gesummv" "c")
set(kernel_name "gesummv")
#remove_bin_bc_files(${kernel_name})
#file(GLOB test_files *.c)
#MESSAGE(STATUS "CMAKE_CURRENT_SOURCE_DIR in gesummv: ${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB test_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
MESSAGE(STATUS "test_files in compile_polybench: ${test_files}")
foreach(test_file ${test_files})
	get_filename_component(file_name ${test_file} NAME_WE)
	MESSAGE(STATUS "file_name in compile_polybench: ${file_name}")
	compile_polybench(${file_name} ${kernel_name} "c")
	list(APPEND bc_files ${file_name}.bc)
endforeach(test_file)

MESSAGE(STATUS "bc_files: ${bc_files}")
link_polybench(${kernel_name} ${bc_files})

#Test polybench using Lin-profiler
#test_polybench(${kernel_name})
