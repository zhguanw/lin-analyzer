#ifndef SLOTTRACKER_H
#define SLOTTRACKER_H

#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallString.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/DebugInfo.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/ValueSymbolTable.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/CallingConv.h"
#include "llvm/Support/Debug.h"
#include <algorithm>
#include <cctype>

namespace llvm {

	/// MDNode map iterators.
	typedef DenseMap<const MDNode*, unsigned>::iterator mdn_iterator;
	/// AttributeSet map iterators.
	typedef DenseMap<AttributeSet, unsigned>::iterator as_iterator;

	class SlotTracker {
	public:

		/// ValueMap - A mapping of Values to slot numbers.
		typedef DenseMap<const Value*, unsigned> ValueMap;

		/// Construct from a module
		explicit SlotTracker(const Module *M);
		/// Construct from a function, starting out in incorp state.
		explicit SlotTracker(const Function *F);

		/// Return the slot number of the specified value in it's type
		/// plane.  If something is not in the SlotTracker, return -1.
		int getLocalSlot(const Value *V);
		int getGlobalSlot(const GlobalValue *V);
		int getMetadataSlot(const MDNode *N);
		int getAttributeGroupSlot(AttributeSet AS);

		/// If you'd like to deal with a function instead of just a module, use
		/// this method to get its data into the SlotTracker.
		void incorporateFunction(const Function *F);

		/// After calling incorporateFunction, use this method to remove the
		/// most recently incorporated function from the SlotTracker. This
		/// will reset the state of the machine back to just the module contents.
		void purgeFunction();

		/// MDNode map iterators.
		mdn_iterator mdn_begin();
		mdn_iterator mdn_end();
		unsigned mdn_size() const;
		bool mdn_empty() const;

		/// AttributeSet map iterators.
		as_iterator as_begin();
		as_iterator as_end();
		unsigned as_size() const;
		bool as_empty() const;

		/// This function does the actual initialization.
		inline void initialize();

		// Implementation Details
	private:
		/// TheModule - The module for which we are holding slot numbers.
		const Module* TheModule;

		/// TheFunction - The function for which we are holding slot numbers.
		const Function* TheFunction;
		bool FunctionProcessed;

		/// mMap - The slot map for the module level data.
		ValueMap mMap;
		unsigned mNext;

		/// fMap - The slot map for the function level data.
		ValueMap fMap;
		unsigned fNext;

		/// mdnMap - Map for MDNodes.
		DenseMap<const MDNode*, unsigned> mdnMap;
		unsigned mdnNext;

		/// asMap - The slot map for attribute sets.
		DenseMap<AttributeSet, unsigned> asMap;
		unsigned asNext;

		/// CreateModuleSlot - Insert the specified GlobalValue* into the slot table.
		void CreateModuleSlot(const GlobalValue *V);

		/// CreateMetadataSlot - Insert the specified MDNode* into the slot table.
		void CreateMetadataSlot(const MDNode *N);

		/// CreateFunctionSlot - Insert the specified Value* into the slot table.
		void CreateFunctionSlot(const Value *V);

		/// \brief Insert the specified AttributeSet into the slot table.
		void CreateAttributeSetSlot(AttributeSet AS);

		/// Add all of the module level global variables (and their initializers)
		/// and function declarations, but not the contents of those functions.
		void processModule();

		/// Add all of the functions arguments, basic blocks, and instructions.
		void processFunction();

		SlotTracker(const SlotTracker &) LLVM_DELETED_FUNCTION;
		void operator=(const SlotTracker &)LLVM_DELETED_FUNCTION;
	};

	SlotTracker *createSlotTracker(const Module *M);
	static SlotTracker *createSlotTracker(const Value *V);

} // End of llvm namespace

#endif // End of SLOTTRACKER_H