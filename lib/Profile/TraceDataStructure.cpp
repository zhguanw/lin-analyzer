#include "profile_h/TraceDataStructure.h"

using namespace llvm;

RecordLoadStore::RecordLoadStore() : record_type(RecordType::NOType), lsID(0), ls_ptr(nullptr), ls_size(0), ls_address(0){}

RecordLoadStore::RecordLoadStore(RecordType ls_type, uint64_t lsid, unsigned char* ptr, uint64_t size) :
	 							 record_type(ls_type), lsID(lsid), ls_ptr(ptr), ls_size(size)
{
	ls_address = reinterpret_cast<uintptr_t>(ptr);
	
	// Initialize eval_Iteration_Index list
	eval_Iteration_Index.clear();
}

RecordLoadStore& RecordLoadStore::operator=(RecordLoadStore const &rhs){
	if (this != &rhs) {
		record_type = rhs.getlsRecordType();
		lsID = rhs.getlsID();
		ls_ptr = rhs.getls_ptr();
		ls_size = rhs.getls_size();
		ls_address = reinterpret_cast<uintptr_t>(ls_ptr);
		eval_Iteration_Index.assign(rhs.eval_Iteration_Index.begin(), rhs.eval_Iteration_Index.end());
	}
	return *this;
}

RecordType RecordLoadStore::getlsRecordType() const {
	return record_type;
}

uint64_t RecordLoadStore::getlsID() const {
	return lsID;
}

uintptr_t RecordLoadStore::getls_address() const {
	return ls_address;
}

uint64_t RecordLoadStore::getls_size() const {
	return ls_size;
}

unsigned char* RecordLoadStore::getls_ptr() const{
	return ls_ptr;
}

void RecordLoadStore::setEval_Iteration_Index(std::list<unsigned> evalIterList){
	eval_Iteration_Index.clear();
	eval_Iteration_Index.assign(evalIterList.begin(), evalIterList.end());
}

std::list<unsigned> RecordLoadStore::getEval_Iteration_Index() const {
	return eval_Iteration_Index;
}

/// RecordTrace member functions declaration
RecordTrace::RecordTrace() {
	RecordLoadStore init_lsRecord(RecordType::NOType, 0, NULL, 0);
	ls_record = init_lsRecord;
	lineID = 0;
	belongedBB = NULL;
	belongedBBid = 0;
	trace_id = 0;
	funcName = "";
	loopName = "";
}

RecordTrace::RecordTrace(RecordLoadStore lsRecord, unsigned lineid, BasicBlock* BB, uint64_t traceID) :
						 ls_record(lsRecord), lineID(lineid), belongedBB(BB), trace_id(traceID) {
	funcName = BB->getParent()->getName();
	loopName = BB2loopNameMap.at(BB);
}

void RecordTrace::setRecordTrace(RecordLoadStore lsRecord, unsigned lineid, BasicBlock* BB, uint64_t BBid, uint64_t traceID){
	ls_record = lsRecord;
	lineID = lineid;
	belongedBB = BB;
	belongedBBid = BBid;
	trace_id = traceID;
	funcName = BB->getParent()->getName();
	loopName = BB2loopNameMap.at(BB);
}

uint64_t RecordTrace::getlsID() const {
	return ls_record.getlsID();
}

unsigned RecordTrace::getlineID() const {
	return lineID;
}

uintptr_t RecordTrace::getls_address() const{
	return ls_record.getls_address();
}

uint64_t RecordTrace::getls_size() const{
	return ls_record.getls_size();
}

std::string RecordTrace::getRecordType() const {
	std::string str;
	RecordType record_type = ls_record.getlsRecordType();

	switch (record_type) {
	case RecordType::LOADType:
		str = "LOADType";
		break;
	case RecordType::STOREType:
		str = "STOREType";
		break;
	default:
		str = "Error, Not a type";
	}

	return str;
}

RecordLoadStore RecordTrace::getls_record() const {
	return ls_record;
}

BasicBlock* RecordTrace::getBelongedBB() const {
	return belongedBB;
}

uint64_t RecordTrace::getBelongedBBid() const {
	return belongedBBid;
}

Function* RecordTrace::getBelongedFunc() const {
	BasicBlock* bb = getBelongedBB();
	return bb->getParent();
}

uint64_t RecordTrace::getTraceID() const {
	return trace_id;
}

void RecordTrace::setEval_Iteration_Index(std::list<unsigned> eval_list) {
	ls_record.setEval_Iteration_Index(eval_list);
}

fn_lpNamePairTy RecordTrace::getfn_lpNamePair() const {
	return std::make_pair(funcName, loopName);
}

DynDataDep::DynDataDep() {
	srcLSID = 0;
	dstLSID = 0;
	depType = DepType::INITType;
	depDistance.clear();
	loop_level = 0;
	loop_independent = true;
	cross_level = false;
	srcfn_lpNamePair = std::make_pair("", "");
	dstfn_lpNamePair = std::make_pair("", "");
}

DynDataDep::DynDataDep(DepType depTy, uint64_t srclsid, uint64_t dstlsid, std::vector<int64_t> depDis, \
											 fn_lpNamePairTy src_fn_lpNamePair, fn_lpNamePairTy dst_fn_lpNamePair, bool inter_loop, bool inter_func) :
											 depType(depTy), srcLSID(srclsid), dstLSID(dstlsid), srcfn_lpNamePair(src_fn_lpNamePair), dstfn_lpNamePair(dst_fn_lpNamePair),
											 inter_loop_dependency(inter_loop), inter_function_dependency(inter_func)
{
	if (inter_loop || inter_func) {
		loop_carried_dependency = false;
		loop_independent = true;
		cross_level = false;
	}
	else {
		depDistance.assign(depDis.begin(), depDis.end());
		loop_independent = true;
		loop_carried_dependency = false;

		unsigned size = depDistance.size();
		loop_level = size;
		for (unsigned i = 0; i < depDistance.size(); i++) {
			if (depDistance.at(i) == INITIAL_BIGVALUE) {
				loop_level--;
			}

			if (depDistance.at(i) != 0) {
				loop_independent = false;
				loop_carried_dependency = true;
				break;
			}
		}

		if (loop_level != size) {
			cross_level = true;
		}
		else {
			cross_level = false;
		}
	}
}

std::pair<uint64_t, uint64_t> DynDataDep::getDepPair() const {
	return std::make_pair(srcLSID, dstLSID);
}

std::vector<int64_t> DynDataDep::getDepDistance() const {
	return depDistance;
}

std::string DynDataDep::getsrcFuncName() const {
	return srcfn_lpNamePair.first;
}

std::string DynDataDep::getsrcLoopName() const {
	return srcfn_lpNamePair.second;
}

std::string DynDataDep::getdstFuncName() const {
	return dstfn_lpNamePair.first;
}

std::string DynDataDep::getdstLoopName() const {
	return dstfn_lpNamePair.second;
}

std::string DynDataDep::getDepType() const {
	std::string str;
	DepType dependence_type = depType;

	switch (dependence_type) {
	case DepType::INITType:
		str = "INITType";
		break;
	case DepType::RAWType:
		str = "RAWType";
		break;
	case DepType::WARType:
		str = "WARType";
		break;
	case DepType::WAWType:
		str = "WAWType";
		break;
	default:
		str = "Error, Not a type";
	}

	return str;
}

bool DynDataDep::isLoopIndependent() const {
	return loop_independent;
}

bool DynDataDep::isCrossLevel() const {
	return cross_level;
}

bool DynDataDep::isloopcarrieddependent() const {
	return loop_carried_dependency;
}

bool DynDataDep::isinterloopdependent() const {
	return inter_loop_dependency;
}

bool DynDataDep::isinterfuncdependent() const {
	return inter_function_dependency;
}

void DynDataDep::setDepType(DepType depTy) {
	depType = DepType::INITType;
}

void DynDataDep::setSrcLSid(uint64_t src_lsid) {
	srcLSID = src_lsid;
}

void DynDataDep::setDstLSid(uint64_t dst_lsid) {
	dstLSID = dst_lsid;
}