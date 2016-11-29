#ifndef LIN_PROFILE_H
#define LIN_PROFILE_H

#define BUILD_DDDG_H

#include "profile_h/AssignBasicBlockIDPass.h"
#include "profile_h/AssignLoadStoreIDPass.h"
#include "profile_h/QueryBasicBlockIDPass.h"
#include "profile_h/QueryLoadStoreIDPass.h"
#include "profile_h/LoopNumberPass.h"
#include "profile_h/ExtractLoopInfoPass.h"

#ifndef BUILD_DDDG_H
//#include "profile_h/GetLoopBoundPass.h"
#include "profile_h/CodeInstrumentPass.h"
#include "profile_h/AnalysisProfilingPass.h"

#else // define BUILD_DDDG_H
#include "profile_h/InstrumentForDDDGPass.h"
#endif // End of ifndef BUILD_DDDG_H

//#include "profile_h/Passes.h"
#include "profile_h/auxiliary.h"

//#define ENABLE_INSTDISTRIBUTION 
#ifdef ENABLE_INSTDISTRIBUTION
#include "profile_h/GetInstDistributionPass.h"
#endif // End of ENABLE_INSTDISTRIBUTION

#define TWO_TIME_PROFILING
#if defined(TWO_TIME_PROFILING) || defined(ENABLE_INSTDISTRIBUTION)
#include "profile_h/GetLoopBoundPass.h"
#endif // End of TWO_TIME_PROFILING or ENABLE_INSTDISTRIBUTION

static llvm::cl::opt<std::string>
InputFilename(llvm::cl::Positional, llvm::cl::desc("<input bitcode file>"), llvm::cl::init("-"), llvm::cl::value_desc("filename"));

static llvm::cl::opt<std::string>
OutputFilename("o", llvm::cl::desc("<output bitcode file>"), llvm::cl::value_desc("filename"));

void parse_input_arguments(int argc, char **argv);
#endif // End LIN_PROFILE_H
