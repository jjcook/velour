//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// atomics.h
//
// wrappers for atomic intrinsics
//

#ifndef _ATOMICS_H_
#define _ATOMICS_H_

#ifdef __GNUC__
#  define GCC_VERSION ( __GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

// **************************************
// *** non-atomic compiler intrinsics ***
// **************************************

#ifndef _mm_prefetch
#define _MM_HINT_NTA 0
#define _MM_HINT_T0 0
// #define _mm_prefetch(addr, hint)		asm volatile ("prefetcht0 %0" :: "m" (*(unsigned long *)addr));
#define _mm_prefetch(addr, hint)		// do nothing
#endif



// **************************************
// *** atomic compiler intrinsics *******
// **************************************
/*
namespace atomics {

template<typename T>
static inline void
fetch_and_add(T *target, T adjustment) {
#ifdef VELOUR_TBB
#  ifdef __GNUC__
#    if GCC_VERSION >= 40102
       __sync_fetch_and_add(target, adjustment);
#    else
#      error "g++ 4.1.2 or later required."
#    endif
#  elif defined __INTEL_COMPILER
#    if __INTEL_COMPILER >= 1010
#      error "TODO: icc compiler intrinsic!"
#    else
#      error "icc 10.1 or later required."
#    endif
#  else
#    error "Unrecognized compiler -- please use gcc or icc."
#  endif
#else // standard, non-atomic
  *target += adjustment;
#endif
}

template<typename T>
static inline void
fetch_and_or(T *target, T adjustment) {
#ifdef VELOUR_TBB
#  ifdef __GNUC__
#    if GCC_VERSION >= 40102
       __sync_fetch_and_or(target, adjustment);
#    else
#      error "g++ 4.1.2 or later required."
#    endif
#  elif defined __INTEL_COMPILER
#    if __INTEL_COMPILER >= 1010
#      error "TODO: icc compiler intrinsic!"
#    else
#      error "icc 10.1 or later required."
#    endif
#  else
#    error "Unrecognized compiler -- please use gcc or icc."
#  endif
#else // standard, non-atomic
  *target |= adjustment;
#endif
}

// FIXME: Ensures atomic load, should be better 
template<typename T>
static inline T
fetch_and_nop(T *target) {
#ifdef VELOUR_TBB
#  ifdef __GNUC__
#    if GCC_VERSION >= 40102
       T retVal = __sync_fetch_and_or(target, 0);
#    else
#      error "g++ 4.1.2 or later required."
#    endif
#  elif defined __INTEL_COMPILER
#    if __INTEL_COMPILER >= 1010
#      error "TODO: icc compiler intrinsic!"
#    else
#      error "icc 10.1 or later required."
#    endif
#  else
#    error "Unrecognized compiler -- please use gcc or icc."
#  endif
#else // standard, non-atomic
  T retVal = *target;
#endif
  return retVal;
}

template<typename T>
static inline T 
compare_and_swap(T* target, T original, T new_val) {
  T old_val = *target;
#ifdef VELOUR_TBB
#  ifdef __GNUC__
#    if GCC_VERSION >= 40102
       old_val = __sync_val_compare_and_swap(target, original, new_val);
       return old_val;
#    else
#      error "g++ 4.1.2 or later required."
#    endif
#  elif defined __INTEL_COMPILER
#    if __INTEL_COMPILER >= 1010
#      error "TODO: icc compiler intrinsic!"
#    else
#      error "icc 10.1 or later required."
#    endif
#  else
#    error "Unrecognized compiler -- please use gcc or icc."
#  endif
#else // standard, non-atomic
  *target = new_val;
  return old_val;
#endif
}

} // namespace atomics
*/
#endif /* _ATOMICS_H_ */

