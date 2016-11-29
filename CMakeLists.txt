#Modified from Shang cmakelist.txt
cmake_minimum_required(VERSION 2.8)
project(lin-profile)

IF (MSVC)
   # Activate C++ exception handling
   IF (NOT CMAKE_CXX_FLAGS MATCHES "/EHsc")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc") 
   ENDIF ()
ENDIF (MSVC)

set(LIN_PROF_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LIN_PROF_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

#if(MSVC) # Work around for the fucking msvc bug on rtti.
  set(LLVM_REQUIRES_RTTI 1)
  set(LLVM_ENABLE_EH 1)
#endif (MSVC)

# Set your project compile flags.
# E.g. if using the C++ header files
# you will need to enable C++11 support
# for your compiler.

#Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${LIN_PROF_SOURCE_DIR}/cmake")
#MESSAGE(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

FIND_PACKAGE(BOOST REQUIRED)
FIND_PACKAGE(ZLIB REQUIRED)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_BINARY_DIR}/include
)

IF (BOOST_FOUND)
	include_directories( ${BOOST_INCLUDE_DIR} )
	MESSAGE(STATUS "ADD BOOST directory: ${BOOST_INCLUDE_DIR}")
	ADD_DEFINITIONS(-DBOOST_EXCEPTION_DISABLE)
ENDIF (BOOST_FOUND) 

IF (ZLIB_FOUND)
	include_directories( ${ZLIB_INCLUDE_DIRS} )
	MESSAGE(STATUS "ADD ZLIB directory: ${ZLIB_INCLUDE_DIRS}")
ENDIF (ZLIB_FOUND) 
#MESSAGE(STATUS "LLVM_TOOLS_BINARY_DIR: ${LLVM_TOOLS_BINARY_DIR}")
#MESSAGE(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
#MESSAGE(STATUS "LLVM_BINARY_DIR: ${LLVM_BINARY_DIR}")
add_definitions( -D_GNU_SOURCE )

set(LLVM_LINK_COMPONENTS
	#$(LLVM_TARGETS_TO_BUILD)
	Support
	)

install(DIRECTORY include
	DESTINATION .
	PATTERN ".svn" EXCLUDE
	)

set(LLVM_LINK_COMPONENTS
  ${LLVM_TARGETS_TO_BUILD}
  Analysis
  BitWriter
  CodeGen
  Core
  IPA
  IPO
  IRReader
  InstCombine
  Instrumentation
  MC
  ObjCARCOpts
  ScalarOpts
  Support
  Target
  TransformUtils
  Vectorize
  )
	
#add a dummy file in lin-profile lib
add_llvm_tool(lin-analyzer
	lin-profile.cpp
	)
	
target_link_libraries(lin-analyzer
	LLVMLinProfiler
	BuildDDDGlib
	#${ZLIB_LIBRARY}
	)
	
add_subdirectory(include)
add_subdirectory(lib)
# add_subdirectory(tools)

OPTION(ENABLE_TESTSUITE "setup the testsuite environment for lin-analyzer" OFF)
IF (ENABLE_TESTSUITE)
	add_subdirectory(testsuite)
ENDIF (ENABLE_TESTSUITE) 
