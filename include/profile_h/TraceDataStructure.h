#ifndef TRACEDATASTRUCTURE_H
#define TRACEDATASTRUCTURE_H

#include "llvm/Pass.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"

#include <list>
#include <forward_list>

//#include "profile_h/ExtractLoopInfoPass.h"
#include "profile_h/auxiliary.h"
#include "profile_h/Passes.h"

#define USE_LIST_FOR_TRACE_F
#define USE_LIST_FOR_RECORDLS
#define USE_LIST_FOR_DATADEPTYPE

namespace llvm {

	enum class RecordType : unsigned {
		NOType = 'NA',   // Not a type
		LOADType = 'L',	// Load record
		STOREType = 'S',	// Store record
	};

	enum class DepType : unsigned {
		INITType = 'IT',  // Initial Type
		RAWType  = 'RAW', // RAW Dependence; flow (aka true) dependences, leading from a store to a load
		WARType  = 'WAR', // WAR Dependence; anti dependences, leading from a load to a store
		WAWType  = 'WAW'  // WAW Dependence; output dependences, leading from stores to stores
	};

	/// The class below is used to record information of 
	/// load or store instructions at runtime
	class RecordLoadStore
	{
	public:
		//RecordLoadStore() { errs() << "ERROR! Designers should send arguments to RecordLoadStore class! \n"; }
		RecordLoadStore();
		RecordLoadStore(RecordType ls_type, uint64_t lsid, unsigned char* ptr, uint64_t size);

		RecordLoadStore& operator=(RecordLoadStore const &rhs);

		// Access function members for querying
		RecordType getlsRecordType() const;
		uint64_t getlsID() const;
		uintptr_t getls_address() const;
		uint64_t getls_size() const;
		unsigned char* getls_ptr() const;

		void setEval_Iteration_Index(std::list<unsigned> evalIterList);
		std::list<unsigned> getEval_Iteration_Index() const;
		//void setLSinst(Instruction* inst);
		//Instruction* getLSinst() const;

		~RecordLoadStore() {}

	private:
		//Instruction* Inst;
		RecordType record_type;
		uint64_t lsID;
		unsigned char* ls_ptr;
		uint64_t ls_size;
		uintptr_t ls_address;  // Load/store address 
		/// Evaluation Iteration Index sequence : The top-level -> lower-level -> the innermost-level loop
		std::list<unsigned> eval_Iteration_Index;
	};
	/// End of the record class declaration

	/// Record the trace:
	class RecordTrace
	{
	public:
		RecordTrace();
		RecordTrace(RecordLoadStore lsRecord, unsigned lineid, BasicBlock* BB, uint64_t traceID);

		// Set recordTrace
		void setRecordTrace(RecordLoadStore lsRecord, unsigned lineid, BasicBlock* BB, uint64_t BBid, uint64_t traceID);

		// Access function members for querying
		uint64_t getlsID() const;
		unsigned getlineID() const;
		uintptr_t getls_address() const;
		uint64_t getls_size() const;
		std::string getRecordType() const;
		RecordLoadStore getls_record() const;

		//Instruction* getRecordInst() const;
		BasicBlock* getBelongedBB() const;
		uint64_t getBelongedBBid() const;
		Function* getBelongedFunc() const;
		uint64_t getTraceID() const;
		fn_lpNamePairTy getfn_lpNamePair() const;

		void setEval_Iteration_Index(std::list<unsigned> eval_list);

		~RecordTrace() {}

	private:
		RecordLoadStore ls_record;
		uint64_t lineID;
		unsigned iteration;
		BasicBlock* belongedBB;
		uint64_t belongedBBid;
		uint64_t trace_id;
		std::string funcName;
		std::string loopName;
	}; // End of RecordTrace class

	/// Record dynamic data dependence
	class DynDataDep 
	{
	public:
		DynDataDep();
		DynDataDep(DepType depTy, uint64_t srclsid, uint64_t dstlsid, std::vector<int64_t> depDis, \
							 fn_lpNamePairTy src_fn_lpNamePair, fn_lpNamePairTy dst_fn_lpNamePair, bool inter_loop, bool inter_func);
		std::pair<uint64_t, uint64_t> getDepPair() const;
		std::vector<int64_t> getDepDistance() const;
		std::string getsrcFuncName() const;
		std::string getsrcLoopName() const;
		std::string getdstFuncName() const;
		std::string getdstLoopName() const;
		std::string getDepType() const;
		bool isLoopIndependent() const;
		bool isCrossLevel() const;
		bool isloopcarrieddependent() const;
		bool isinterloopdependent() const;
		bool isinterfuncdependent() const;
		void setSrcLSid(uint64_t src_lsid);
		void setDstLSid(uint64_t dst_lsid);
		void setDepType(DepType depTy);

	private:
		DepType depType;
		std::vector<int64_t> depDistance;
		uint64_t srcLSID;
		uint64_t dstLSID;
		unsigned loop_level;
		bool loop_independent;
		bool cross_level;
		bool loop_carried_dependency;
		bool inter_loop_dependency;
		bool inter_function_dependency;
		/*
		std::string funcName;
		std::string loopName;
		*/

		fn_lpNamePairTy srcfn_lpNamePair;
		fn_lpNamePairTy dstfn_lpNamePair;
	};

#ifndef USE_LIST_FOR_TRACE_F
	typedef std::vector<RecordTrace> TraceType;
#else
	typedef std::list<RecordTrace> TraceType;
#endif // End of USE_LIST_FOR_TRACE_F
	extern TraceType Trace_f;

	// Vector used to record load/store information at runtime
#ifndef USE_LIST_FOR_RECORDLS
	extern std::vector<RecordLoadStore> recordls;
#else
	extern std::list<RecordLoadStore> recordls;
#endif // End of USE_LIST_FOR_RECORDLS

	//typedef std::vector<uintptr_t> lsAddressVecTy;
	//extern lsAddressVecTy lsAddressVec;

#ifndef USE_LIST_FOR_DATADEPTYPE
	typedef std::vector<DynDataDep> DataDepType;
#else
	typedef std::list<DynDataDep> DataDepType;
#endif // End of USE_LIST_FOR_DATADEPTYPE
	extern DataDepType DynDataDepVec;
	extern DataDepType RefinedDynDataDepVec;
	extern DataDepType RefinedInterDynDataDepVec;
}
#endif // End of TRACEDATASTRUCTURE_H