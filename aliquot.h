#ifndef ALIQUOT_H
#define ALIQUOT_H

#include <stdint.h>
#include <inttypes.h>

typedef uint_fast64_t num_type;
#define NUM_TYPE_MAX UINT64_MAX
#define PRI_NUM_TYPE "%" PRIuFAST64
#define SCN_NUM_TYPE "%" SCNuFAST64

void aliquot_init(const char *primes_filename);
void aliquot_terminate();
num_type calc_sigma(num_type num);
num_type calc_aliquot_sum(num_type num);

#endif/*ALIQUOT_H*/

