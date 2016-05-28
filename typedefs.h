#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <stdint.h>
#include <inttypes.h>

typedef uint_fast64_t num_type;
#define NUM_TYPE_MAX UINT_FAST64_MAX
#define PRI_NUM_TYPE "%" PRIuFAST64
#define SCN_NUM_TYPE "%" SCNuFAST64

#define MAX_EXP 63
typedef uint_fast8_t exp_type;
#define EXP_TYPE_MAX UINT_FAST8_MAX
#define PRI_EXP_TYPE "%" PRIuFAST8
#define EXP_TYPE_MAX_MASK (((exp_type)1)<<7)

// f(t) = (2+1)*(3+1)*...*(pt+1) - 2*3*...*pt, t >= 1
// log2(f(3))  <  8 < log2(f(4))
// log2(f(5))  < 16 < log2(f(6))
// log2(f(9))  < 32 < log2(f(10))
// log2(f(15)) < 64 < log2(f(16))
#define MAX_POW_COUNT 15
typedef uint_fast8_t pow_idx_type;
#define PRI_POW_IDX_TYPE "%" PRIuFAST8

// num_type x
#define round_sqrt(x) llrint(sqrt((double)(x)))

#endif/*TYPEDEFS_H*/

