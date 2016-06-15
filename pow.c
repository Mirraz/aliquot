#include <assert.h>
#include "typedefs.h"

// returns: base^exp
static num_type calc_pow_impl(num_type base, exp_type exp) __attribute__ ((unused));
static num_type calc_pow_impl(num_type base, exp_type exp) {
	assert(exp > 0);
	// here simple sequential multiplication faster than fast exponentiation by squaring
	num_type result = base;
	for (exp_type i=1; i<exp; ++i) {
		assert(result <= NUM_TYPE_MAX / base);
		result *= base;
	}
	return result;
}

// returns: base^exp + base^(exp-1) + ... + base + 1
static num_type calc_pow_sigma_impl(num_type base, exp_type exp) __attribute__ ((unused));
static num_type calc_pow_sigma_impl(num_type base, exp_type exp) {
	assert(exp > 0);
	num_type result = base + 1;
	for (exp_type i=1; i<exp; ++i) {
		assert(result <= NUM_TYPE_MAX / base);
		result *= base;
		assert(result < NUM_TYPE_MAX);
		++result;
	}
	return result;
}
