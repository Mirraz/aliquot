#ifndef INVERSE_H
#define INVERSE_H

#include <stdbool.h>
#include "typedefs.h"

typedef struct {
	num_type prime;
	exp_type exp;
} prime_pow_struct;

// returns: whether to stop the search algorithm?
typedef bool (*aliquot_inverse_cb_type)(prime_pow_struct prime_pow_list[], pow_idx_type pow_count);

void aliquot_inverse_init(const char *primes_filename);
void aliquot_inverse_terminate();
void aliquot_inverse_set_callback(aliquot_inverse_cb_type aliquot_inverse_cb);
void aliquot_inverse_run(num_type req_aliquot_sum);

num_type calc_pow(num_type base, exp_type exp);

#endif/*INVERSE_H*/

