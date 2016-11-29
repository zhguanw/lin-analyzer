#add_subdirectory(Giri)
if (NOT TARGET LLVMLinProfiler)
	add_subdirectory(Profile)
endif (NOT TARGET LLVMLinProfiler)

#add_subdirectory(Profile)
add_subdirectory(Build_DDDG)
