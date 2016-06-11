#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// exit, EXIT_FAILURE
#include <stddef.h>		// size_t
#include <stdbool.h>
#include <memory.h>
#include <math.h>
#include <fenv.h>

#include "typedefs.h"
#include "primes.h"
#include "inverse.h"

#define dprintf(...) fprintf(stderr, __VA_ARGS__)
#define assert_may_overflow assert
//#define assert_may_overflow(cond) ;

typedef enum {
	PRIME_STATUS_ERROR,
	PRIME_STATUS_UNKNOWN,
	PRIME_STATUS_NOT_PRIME,
	PRIME_STATUS_PRIME
} prime_status_enum;

typedef struct {
	num_type prime;        // === base
	num_type pow;          // = prime^exp
	num_type pow_sigma;    // = sigma(pow)
	num_type prefix_mul;
	num_type prefix_sigma; // = sigma(prefix_mul)
	size_t prime_idx;
	exp_type exp;
	prime_status_enum prime_status;
	bool is_prime_in_list;
} calc_struct;

const num_type first_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};

prime_type *primes;
size_t primes_count;
num_type primes_max;

aliquot_inverse_cb_type aliquot_inverse_external_cb;

#ifndef NDEBUG
void calc_list_dprint(calc_struct calc_list[], pow_idx_type pow_count) {
	pow_idx_type i;
	dprintf("exp:              "); for (i=0; i<pow_count; ++i) dprintf(PRI_EXP_TYPE "\t", calc_list[i].exp         ); dprintf("\n");
	dprintf("prime:            "); for (i=0; i<pow_count; ++i) dprintf(PRI_NUM_TYPE "\t", calc_list[i].prime       ); dprintf("\n");
	dprintf("pow:              "); for (i=0; i<pow_count; ++i) dprintf(PRI_NUM_TYPE "\t", calc_list[i].pow         ); dprintf("\n");
	dprintf("pow_sigma:        "); for (i=0; i<pow_count; ++i) dprintf(PRI_NUM_TYPE "\t", calc_list[i].pow_sigma   ); dprintf("\n");
	dprintf("prefix_mul:       "); for (i=0; i<pow_count; ++i) dprintf(PRI_NUM_TYPE "\t", calc_list[i].prefix_mul  ); dprintf("\n");
	dprintf("prefix_sigma:     "); for (i=0; i<pow_count; ++i) dprintf(PRI_NUM_TYPE "\t", calc_list[i].prefix_sigma); dprintf("\n");
	
	dprintf("is_prime_in_list: "); for (i=0; i<pow_count; ++i) dprintf("%s" "\t", calc_list[i].is_prime_in_list ? "true" : "false"); dprintf("\n");
	dprintf("prime_idx:        "); for (i=0; i<pow_count; ++i) dprintf("%zu" "\t", calc_list[i].prime_idx   ); dprintf("\n");
	dprintf("prime_status:     "); for (i=0; i<pow_count; ++i) dprintf( "%u" "\t", calc_list[i].prime_status); dprintf("\n");
}
#endif

static bool aliquot_inverse_cb(calc_struct prime_calc_list[], pow_idx_type pow_count) {
	static prime_pow_struct prime_pow_list[MAX_POW_COUNT];
	for (pow_idx_type i=0; i<pow_count; ++i) {
		prime_pow_list[i].prime = prime_calc_list[i].prime;
		prime_pow_list[i].exp   = prime_calc_list[i].exp;
	}
	return aliquot_inverse_external_cb(prime_pow_list, pow_count);
}

// returns: base^exp
num_type calc_pow(num_type base, exp_type exp) {
	assert(exp > 0);
	assert(base > 1);
	if (exp == 1) return base;
	// here simple sequential multiplication faster than fast exponentiation by squaring
	num_type result = 1;
	for (exp_type i=0; i<exp; ++i) {
		assert(result <= NUM_TYPE_MAX / base);
		result *= base;
	}
	return result;
}

// returns: base^exp + base^(exp-1) + ... + base + 1
num_type calc_pow_sigma(num_type base, exp_type exp) {
	assert(exp > 0);
	assert(base > 1);
	if (exp == 1) {
		assert(base < NUM_TYPE_MAX);
		return base + 1;
	}
	num_type result = 1;
	for (exp_type i=0; i<exp; ++i) {
		assert(result <= NUM_TYPE_MAX / base);
		result *= base;
		assert(result < NUM_TYPE_MAX);
		++result;
	}
	return result;
}

bool is_big_odd_number_prime(num_type n) {
	assert(n > primes_max);
	assert(n & 1);
	num_type n_sqrt = round_sqrt(n);
	num_type p = 3;
	size_t idx = 2;
	while (p <= n_sqrt) {
		if (!(n % p)) return false;
		if (idx < primes_count) {
			p = primes[idx++];
		} else {
			p += 2;
		}
	}
	return true;
}

bool is_prime(num_type n) {
	assert(n > 1);
	assert(primes_count > 0);
	if (!(n & 1)) return n == 2;
	if (n <= primes_max) {
		size_t left_idx = 0, right_idx = primes_count-1, mid_idx;
		do {
			mid_idx = left_idx + (right_idx - left_idx) / 2;
			if      (primes[mid_idx] < n) left_idx  = mid_idx + 1;
			else if (primes[mid_idx] > n) right_idx = mid_idx - 1;
			else return true;
		} while (left_idx <= right_idx);
		return false;
	} else {
		return is_big_odd_number_prime(n);
	}
}

// --- find last prime --- {

num_type calc_aliquot_for_pow(exp_type exp, num_type prime) {
	assert(exp > 1);
	return calc_pow_sigma(prime, exp-1);
}

// XXX [may be outdated] max req_aliquot_sum < 2^54.19906
// (without assert_may_overflow): n = 2^55, 3^35, 5^24, 7^21: find_base_for_pow(exp = 42): p_mid   = 3: 3^41 overflows
// (within  assert_may_overflow):                  s >= 2^42: find_base_for_pow(exp = 42): p_right = 3: 3^41 overflows

// returns: 0 if not found, value>1 if found
num_type find_base_for_pow(exp_type exp, num_type req_aliquot_sum) {
	assert(exp > 1);
	
	exp_type sum_log2 = floor_log2(req_aliquot_sum);
	num_type p_left;
	if (sum_log2 > (exp << 1)) {
		exp_type p_left_log2  = (sum_log2 - 2) / (exp - 1) - 1;
		assert(p_left_log2 >= 1);
		p_left  = ((num_type)1) << p_left_log2;
	} else {
		p_left = 2;
	}
	num_type p_right;
	assert(exp - 2 <= EXP_TYPE_MAX / 3);
	assert(sum_log2 <= EXP_TYPE_MAX / 2);
	// preventing premature overflow
	if (3 * (exp - 2) < 2 * sum_log2) {
		// p_right_log2 = (sum_log2 + 3) / (exp - 1)
		assert(sum_log2 < EXP_TYPE_MAX);
		assert(exp <= EXP_TYPE_MAX - (sum_log2 + 1));
		exp_type p_right_log2 = (sum_log2 + exp + 1) / (exp - 1);
		assert(p_right_log2 <= MAX_EXP);
		p_right = ((num_type)1) << p_right_log2;
	} else {
		p_right = 3;
	}
	
	assert(p_left >= 2);
	assert(p_left < p_right);
	assert             (calc_aliquot_for_pow(exp, p_left ) < req_aliquot_sum);
	assert_may_overflow(calc_aliquot_for_pow(exp, p_right) > req_aliquot_sum);
	
	++p_left; --p_right;
	
	while (p_left <= p_right) {
		num_type p_mid = p_left + (p_right - p_left) / 2;
		num_type aliquot_sum = calc_aliquot_for_pow(exp, p_mid);
		if      (aliquot_sum < req_aliquot_sum) p_left  = p_mid + 1;
		else if (aliquot_sum > req_aliquot_sum) p_right = p_mid - 1;
		else return p_mid;
	}
	return 0;
}

// returns: 0 if not found, value>1 if found
num_type calc_last_base_for_exp_1(num_type prefix_sigma, num_type prefix_aliquot, num_type req_aliquot_sum) {
	assert(req_aliquot_sum > prefix_sigma);
	num_type p_numerator = req_aliquot_sum - prefix_sigma;
	num_type p_denominator = prefix_aliquot;
	assert(p_denominator <= NUM_TYPE_MAX / 2);
	assert(p_numerator >= 2*p_denominator);
	if (p_numerator % p_denominator) return 0;
	return p_numerator / p_denominator;
}

num_type calc_aliquot(num_type prefix_sigma, num_type prefix_aliquot, exp_type exp, num_type prime) {
	assert(prefix_sigma > 1); assert(prefix_aliquot >= 1); assert(prefix_sigma > prefix_aliquot);
	assert(exp >= 2); assert(prime >= 2);
	num_type pow, pows_sum;
	//pow = calc_pow(prime, exp);
	//pows_sum = calc_pow_sigma(prime, exp-1);
	{
		assert(prime <= NUM_TYPE_MAX / prime);
		pow = prime * prime;
		pows_sum = prime + 1;
		for (exp_type i=0; i<exp-2; ++i) {
			assert(pow <= NUM_TYPE_MAX / prime);
			pow *= prime;
			assert(pows_sum <= NUM_TYPE_MAX / prime);
			pows_sum *= prime;
			assert(pows_sum < NUM_TYPE_MAX);
			++pows_sum;
		}
	}
	assert(pow <= NUM_TYPE_MAX / prefix_aliquot);
	num_type summand1 = prefix_aliquot * pow;
	assert(pows_sum <= NUM_TYPE_MAX / prefix_sigma);
	num_type summand2 = prefix_sigma * pows_sum;
	assert(summand1 <= NUM_TYPE_MAX - summand2);
	return summand1 + summand2;
}

#define min(x,y) ((x) < (y) ? (x) : (y))

// returns: 0 if not found, value>1 if found
num_type find_last_base(
		num_type prefix_sigma, num_type prefix_aliquot,
		exp_type exp,
		num_type min_base,
		num_type req_aliquot_sum
) {
	assert(exp >= 2);
	assert(min_base >= 5);
	assert(min_base % 2 == 1);
	
	exp_type sum_log2 = floor_log2(req_aliquot_sum);
	exp_type prefix_sigma_log2 = floor_log2(prefix_sigma);
	exp_type prefix_aliq_log2  = floor_log2(prefix_aliquot);
	num_type p_left;
	assert(exp <= EXP_TYPE_MAX/2 - 1);
	assert(prefix_aliq_log2  <= EXP_TYPE_MAX - 2*(exp - 1));
	assert(prefix_sigma_log2 <= EXP_TYPE_MAX - 2*(exp + 1));
	if (sum_log2 >= prefix_aliq_log2 + 2*(exp - 1) && sum_log2 >= prefix_sigma_log2 + 2*(exp + 1)) {
		exp_type p_left_log2_01 = (sum_log2 - prefix_aliq_log2  - 2) / exp;
		exp_type p_left_log2_02 = (sum_log2 - prefix_sigma_log2 - 4) / (exp - 1);
		exp_type p_left_log2 = min(p_left_log2_02, p_left_log2_01);
		assert(p_left_log2 > 0);
		--p_left_log2;
		assert(p_left_log2 <= MAX_EXP);
		p_left = ((num_type)1) << p_left_log2;
	} else {
		p_left = 2;
	}
	num_type p_right;
	{
		assert(exp <= EXP_TYPE_MAX - sum_log2);
		assert(sum_log2 + exp >= prefix_aliq_log2);
		// p_right_log2_01 = (sum_log2 - prefix_aliq_log2 + 1) / exp
		exp_type p_right_log2_01 = (sum_log2 + exp - prefix_aliq_log2) / exp;
		// p_right_log2_02 = (sum_log2 - prefix_sigma_log2 + 3) / (exp - 1)
		exp_type p_right_log2_02 = (sum_log2 + exp - prefix_sigma_log2 + 1) / (exp - 1);
		exp_type p_right_log2 = min(p_right_log2_02, p_right_log2_01);
		assert(p_right_log2 <= MAX_EXP);
		p_right = ((num_type)1) << p_right_log2;
	}
	
	assert(p_left >= 2);
	assert(p_left < p_right);
	assert             (calc_aliquot(prefix_sigma, prefix_aliquot, exp, p_left ) < req_aliquot_sum);
	assert_may_overflow(calc_aliquot(prefix_sigma, prefix_aliquot, exp, p_right) > req_aliquot_sum);
	
	++p_left; --p_right;
	if (p_left < min_base) p_left = min_base;
	
	while (p_left <= p_right) {
		num_type p_mid = p_left + (p_right - p_left) / 2;
		num_type aliquot_sum = calc_aliquot(prefix_sigma, prefix_aliquot, exp, p_mid);
		if      (aliquot_sum < req_aliquot_sum) p_left  = p_mid + 1;
		else if (aliquot_sum > req_aliquot_sum) p_right = p_mid - 1;
		else return p_mid;
	}
	return 0;
}

void find_last_prime_calc(calc_struct prime_calc_list[], pow_idx_type pow_count, num_type req_aliquot_sum) {
	assert(pow_count > 1);
	exp_type exp = prime_calc_list[pow_count-1].exp;
	
	num_type prime;
	num_type prefix_sigma = prime_calc_list[pow_count-2].prefix_sigma;
	num_type prefix_aliquot = prefix_sigma - prime_calc_list[pow_count-2].prefix_mul;
	if (exp == 1) {
		prime = calc_last_base_for_exp_1(prefix_sigma, prefix_aliquot, req_aliquot_sum);
		assert(prime == 0 || prime > prime_calc_list[pow_count-1].prime);
	} else {
		num_type min_base = prime_calc_list[pow_count-1].prime;
		assert(min_base > 2);
		assert(min_base <= NUM_TYPE_MAX - 2);
		assert(calc_aliquot(prefix_sigma, prefix_aliquot, exp, min_base) < req_aliquot_sum);
		prime = find_last_base(prefix_sigma, prefix_aliquot, exp, min_base+2, req_aliquot_sum);
	}
	
	if (prime != 0 && is_prime(prime)) {
		prime_calc_list[pow_count-1].prime = prime;
		aliquot_inverse_cb(prime_calc_list, pow_count);
	}
}

void find_last_prime_calc_for_pow(calc_struct prime_calc_list[], pow_idx_type pow_count, num_type req_aliquot_sum) {
	num_type prime = find_base_for_pow(prime_calc_list[pow_count-1].exp, req_aliquot_sum);
	if (prime != 0 && is_prime(prime)) {
		prime_calc_list[pow_count-1].prime = prime;
		aliquot_inverse_cb(prime_calc_list, pow_count);
	}
}

// --- find last prime --- }
// --- iterate primes --- {

// maybe out_prime_calc == in_prime_calc
void next_maybe_prime_calc(const calc_struct *in_prime_calc, calc_struct *out_prime_calc) {
	if (in_prime_calc->is_prime_in_list && in_prime_calc->prime_idx < primes_count-1) {
		size_t out_prime_idx = in_prime_calc->prime_idx + 1;
		out_prime_calc->prime = primes[out_prime_idx];
		out_prime_calc->prime_idx = out_prime_idx;
		out_prime_calc->is_prime_in_list = true;
		out_prime_calc->prime_status = PRIME_STATUS_PRIME;
	} else {
		assert(in_prime_calc->prime <= NUM_TYPE_MAX - 2);
		out_prime_calc->prime = in_prime_calc->prime + 2;
		out_prime_calc->is_prime_in_list = false;
		out_prime_calc->prime_status = PRIME_STATUS_UNKNOWN;
	}
}

bool is_prime_calc(calc_struct *cur_prime_calc) {
	assert(cur_prime_calc->prime_status != PRIME_STATUS_ERROR);
	if (cur_prime_calc->prime_status == PRIME_STATUS_PRIME) return true;
	assert(cur_prime_calc->prime_status == PRIME_STATUS_UNKNOWN);
	assert(!cur_prime_calc->is_prime_in_list);
	assert(cur_prime_calc->prime > primes_max);
	bool is_prime_res = is_big_odd_number_prime(cur_prime_calc->prime);
	cur_prime_calc->prime_status = is_prime_res ? PRIME_STATUS_PRIME : PRIME_STATUS_NOT_PRIME;
	return is_prime_res;
}

num_type prime_calc_recalc(calc_struct prime_calc_list[], pow_idx_type idx) {
	calc_struct *cur_prime_calc = &(prime_calc_list[idx]);
	cur_prime_calc->pow       = calc_pow      (cur_prime_calc->prime, cur_prime_calc->exp);
	cur_prime_calc->pow_sigma = calc_pow_sigma(cur_prime_calc->prime, cur_prime_calc->exp);
	if (idx != 0) {
		assert(cur_prime_calc->pow <= NUM_TYPE_MAX / prime_calc_list[idx-1].prefix_mul);
		cur_prime_calc->prefix_mul = prime_calc_list[idx-1].prefix_mul * cur_prime_calc->pow;
		// XXX if req_aliquot_sum is odd
		//     req_aliquot_sum < 8589934590 = 2^33-2
		//     for sum = 8589934590:
		//         prime_calc_list = [.0 = {.exp=1, .prime=4294967295}, .1 = {.exp=1, prime=4294967297}]
		//         idx = 1, cur_prime_calc->pow_sigma = 4294967298, prime_calc_list[idx-1].prefix_sigma = 4294967296
		assert(cur_prime_calc->pow_sigma <= NUM_TYPE_MAX / prime_calc_list[idx-1].prefix_sigma);
		cur_prime_calc->prefix_sigma = prime_calc_list[idx-1].prefix_sigma * cur_prime_calc->pow_sigma;
	} else {
		cur_prime_calc->prefix_mul = cur_prime_calc->pow;
		cur_prime_calc->prefix_sigma = cur_prime_calc->pow_sigma;
	}
	assert(cur_prime_calc->prefix_sigma > cur_prime_calc->prefix_mul);
	return cur_prime_calc->prefix_sigma - cur_prime_calc->prefix_mul;
}

// is_prime_in_list values: [true, ..., true, false, ..., false]
// begin_id <= idx < end_idx
bool is_prime_calc_list_slice(calc_struct prime_calc_list[], pow_idx_type begin_idx, pow_idx_type end_idx) {
	assert(begin_idx < end_idx);
	pow_idx_type idx = end_idx - 1;
	while (true) {
		if (prime_calc_list[idx].is_prime_in_list) {
			++idx;
			break;
		}
		if (idx == begin_idx) break;
		--idx;
	}
	while (idx < end_idx && is_prime_calc(&(prime_calc_list[idx]))) ++idx;
	return idx == end_idx;
}

// inc and fill: both use maybe-primes
// returns: true - success, false - failed
bool inc_and_fill_maybe_primes(
	calc_struct prime_calc_list[], pow_idx_type pow_count,
	pow_idx_type idx, num_type req_aliquot_sum
) {
	pow_idx_type initial_idx = idx;
	const calc_struct *prev_prime_calc = &(prime_calc_list[idx]);
	for (; idx<pow_count; ++idx) {
		calc_struct *cur_prime_calc = &(prime_calc_list[idx]);
		next_maybe_prime_calc(prev_prime_calc, cur_prime_calc);
		num_type prefix_aliquot_sum = prime_calc_recalc(prime_calc_list, idx);
		if (prefix_aliquot_sum > req_aliquot_sum) {
			return false;
		} else if (prefix_aliquot_sum == req_aliquot_sum) {
			if (idx == pow_count-1 && is_prime_calc_list_slice(prime_calc_list, initial_idx, pow_count)) {
				aliquot_inverse_cb(prime_calc_list, pow_count);
			}
			return false;
		}
		prev_prime_calc = cur_prime_calc;
	}
	return true;
}

// returns: false - next is calculated, true - end is reached
bool prime_calc_next(
	calc_struct prime_calc_list[], pow_idx_type pow_count,
	pow_idx_type min_inc_prime_idx, num_type req_aliquot_sum
) {
	assert(min_inc_prime_idx <= 1);
	assert(pow_count >= 2);
	assert(min_inc_prime_idx <= pow_count - 2);
	pow_idx_type idx = pow_count - 2;
	
	while (true) {
		bool fill_res = inc_and_fill_maybe_primes(prime_calc_list, pow_count, idx, req_aliquot_sum);
		if (fill_res) {
			if (is_prime_calc_list_slice(prime_calc_list, idx, pow_count-1)) return false;
		} else {
			if (idx == min_inc_prime_idx) return true;
			--idx;
		}
	}
}

void prime_calc(calc_struct exp_calc_list[], pow_idx_type pow_count, num_type req_aliquot_sum) {
	static calc_struct prime_calc_list[MAX_POW_COUNT];
	memcpy(prime_calc_list, exp_calc_list, sizeof(exp_calc_list[0])*pow_count);
	(void)exp_calc_list;
	
	if (pow_count == 1) {
		// consider parity
		if ((req_aliquot_sum & 1) != (prime_calc_list[0].exp & 1)) return;
		find_last_prime_calc_for_pow(prime_calc_list, pow_count, req_aliquot_sum);
		return;
	}
	
	for (pow_idx_type i=0; i<pow_count; ++i) {
		prime_calc_list[i].prime_idx = i;
		prime_calc_list[i].is_prime_in_list = true;
		prime_calc_list[i].prime_status = PRIME_STATUS_PRIME;
	}
	
	// consider parity
	pow_idx_type min_inc_prime_idx = 0;
	bool inc_first_prime = false;
	if (!(req_aliquot_sum & 1)) {
		pow_idx_type i = 1;
		while (i<pow_count && !(prime_calc_list[i].exp & 1)) ++i;
		if (i < pow_count) {
			min_inc_prime_idx = 1;
		} else {
			inc_first_prime = true;
		}
	} else {
		pow_idx_type i = 1;
		while (i<pow_count && !(prime_calc_list[i].exp & 1)) ++i;
		if (i < pow_count) {
			inc_first_prime = true;
		} else {
			if (!(prime_calc_list[0].exp & 1)) min_inc_prime_idx = 1;
		}
	}
	
	if (inc_first_prime) {
		// XXX if req_aliquot_sum is even: req_aliquot_sum < 7925355423235 = aliquot_sum(2^26*3^10) < 2^42.85
		bool fill_res = inc_and_fill_maybe_primes(prime_calc_list, pow_count, 0, req_aliquot_sum);
		if (!fill_res) return;
	}
	
	if (pow_count <= min_inc_prime_idx + 1) {
		find_last_prime_calc(prime_calc_list, pow_count, req_aliquot_sum);
		return;
	}
	
	bool is_end_reached;
	do {
#ifndef NDEBUG
		for (pow_idx_type i=0; i<pow_count-1; ++i) assert(prime_calc_list[i].prime_status == PRIME_STATUS_PRIME);
#endif
		find_last_prime_calc(prime_calc_list, pow_count, req_aliquot_sum);
		is_end_reached = prime_calc_next(prime_calc_list, pow_count, min_inc_prime_idx, req_aliquot_sum);
	} while (!is_end_reached);
}

// --- iterate primes --- }
// --- iterate exps --- {

num_type exp_calc_inc_exp(calc_struct exp_calc_list[], pow_idx_type idx) {
	calc_struct *cur_exp_calc = &(exp_calc_list[idx]);
	assert(cur_exp_calc->exp < MAX_EXP);
	++cur_exp_calc->exp;
	assert(cur_exp_calc->prefix_mul <= NUM_TYPE_MAX / cur_exp_calc->prime);
	cur_exp_calc->prefix_mul *= cur_exp_calc->prime;
	cur_exp_calc->pow *= cur_exp_calc->prime;
	assert(cur_exp_calc->pow <= NUM_TYPE_MAX - cur_exp_calc->pow_sigma);
	cur_exp_calc->pow_sigma += cur_exp_calc->pow;
	assert(idx == 0 || exp_calc_list[idx-1].prefix_sigma <= NUM_TYPE_MAX / cur_exp_calc->pow_sigma);
	cur_exp_calc->prefix_sigma = (
		idx != 0 ?
		exp_calc_list[idx-1].prefix_sigma * cur_exp_calc->pow_sigma :
		cur_exp_calc->pow_sigma
	);
	assert(cur_exp_calc->prefix_sigma > cur_exp_calc->prefix_mul);
	return cur_exp_calc->prefix_sigma - cur_exp_calc->prefix_mul;
}

num_type exp_calc_fill_exp(calc_struct exp_calc_list[], pow_idx_type idx) {
	calc_struct *cur_exp_calc = &(exp_calc_list[idx]);
	cur_exp_calc->exp = 1;
	cur_exp_calc->pow = cur_exp_calc->prime;
	assert(idx > 0);
	assert(exp_calc_list[idx-1].prefix_mul <= NUM_TYPE_MAX / cur_exp_calc->prime);
	cur_exp_calc->prefix_mul = exp_calc_list[idx-1].prefix_mul * cur_exp_calc->prime;
	assert(cur_exp_calc->prime < NUM_TYPE_MAX);
	cur_exp_calc->pow_sigma = cur_exp_calc->prime + 1;
	assert(exp_calc_list[idx-1].prefix_sigma <= NUM_TYPE_MAX / cur_exp_calc->pow_sigma);
	cur_exp_calc->prefix_sigma = exp_calc_list[idx-1].prefix_sigma * cur_exp_calc->pow_sigma;
	assert(cur_exp_calc->prefix_sigma > cur_exp_calc->prefix_mul);
	return cur_exp_calc->prefix_sigma - cur_exp_calc->prefix_mul;
}

num_type exp_calc_fill_any_exp(calc_struct exp_calc_list[], pow_idx_type idx) {
	if (idx != 0) {
		return exp_calc_fill_exp(exp_calc_list, idx);
	} else {
		calc_struct *cur_exp_calc = &(exp_calc_list[0]);
		cur_exp_calc->exp = 1;
		cur_exp_calc->pow = cur_exp_calc->prime;
		cur_exp_calc->prefix_mul = cur_exp_calc->prime;
		assert(cur_exp_calc->prime < NUM_TYPE_MAX);
		cur_exp_calc->pow_sigma = cur_exp_calc->prime + 1;
		cur_exp_calc->prefix_sigma = cur_exp_calc->pow_sigma;
		assert(cur_exp_calc->prefix_sigma > cur_exp_calc->prefix_mul);
		return cur_exp_calc->prefix_sigma - cur_exp_calc->prefix_mul;
	}
}

// XXX req_aliquot_sum < 43281303292150770 < 2^56

// returns: false - next is calculated, true - end is reached
bool exp_calc_next(calc_struct exp_calc_list[], pow_idx_type pow_count, num_type req_aliquot_sum) {
	assert(pow_count > 0);
	for (pow_idx_type idx = pow_count-1;; --idx) {
		num_type inc_prefix_aliquot_sum = exp_calc_inc_exp(exp_calc_list, idx);
		if (inc_prefix_aliquot_sum < req_aliquot_sum) {
			num_type fill_prefix_aliquot_sum;
			do {
				if (idx == pow_count-1) return false;
				++idx;
				fill_prefix_aliquot_sum = exp_calc_fill_exp(exp_calc_list, idx);
			} while (fill_prefix_aliquot_sum < req_aliquot_sum);
			if (fill_prefix_aliquot_sum == req_aliquot_sum && idx == pow_count-1) {
				aliquot_inverse_cb(exp_calc_list, pow_count);
			}
		} else if (inc_prefix_aliquot_sum == req_aliquot_sum && idx == pow_count-1) {
			aliquot_inverse_cb(exp_calc_list, pow_count);
		}
		if (idx == 0) return true;
	}
}

void exp_calc(num_type req_aliquot_sum) {
	assert(req_aliquot_sum > 1);
	
	static calc_struct exp_calc_list[MAX_POW_COUNT];
	for (pow_idx_type i=0; i<MAX_POW_COUNT; ++i) exp_calc_list[i].prime = first_primes[i];
	
	for (pow_idx_type pow_count=1;; ++pow_count) {
		num_type prefix_aliquot_sum;
		for (pow_idx_type fill_idx=0; fill_idx<pow_count; ++fill_idx) {
			prefix_aliquot_sum = exp_calc_fill_any_exp(exp_calc_list, fill_idx);
			assert(fill_idx == pow_count-1 || prefix_aliquot_sum < req_aliquot_sum);
		}
		if (prefix_aliquot_sum >= req_aliquot_sum) {
			if (prefix_aliquot_sum == req_aliquot_sum) aliquot_inverse_cb(exp_calc_list, pow_count);
			break;
		}
		
		if (pow_count == 1) {
			bool is_end_reached = exp_calc_next(exp_calc_list, pow_count, req_aliquot_sum);
			if (is_end_reached) continue;
		}
		
		bool is_end_reached;
		do {
			prime_calc(exp_calc_list, pow_count, req_aliquot_sum);
			is_end_reached = exp_calc_next(exp_calc_list, pow_count, req_aliquot_sum);
		} while (!is_end_reached);
	}
}

// --- iterate exps --- }
// --- main {

int default_rounding_direction;
primes_array_struct *primes_array;

void aliquot_inverse_init(const char *primes_filename) {
	assert(sizeof(first_primes)/sizeof(first_primes[0]) >= MAX_POW_COUNT);

	primes_array = primes_construct(primes_filename);
	primes = primes_get_array(primes_array);
	primes_count = primes_get_count(primes_array);
	assert(primes_count > 0);
	primes_max = primes[primes_count-1];
	
	default_rounding_direction = fegetround();
	fesetround(0);
}

void aliquot_inverse_terminate() {
	fesetround(default_rounding_direction);
	primes_destruct(primes_array);
}

void aliquot_inverse_set_callback(aliquot_inverse_cb_type cb) {
	aliquot_inverse_external_cb = cb;
}

void aliquot_inverse_run(num_type req_aliquot_sum) {
	if (req_aliquot_sum <= 1) {
		fprintf(stderr, "Wrong input: " PRI_NUM_TYPE "\n", req_aliquot_sum);
		exit(EXIT_FAILURE);
	}
	// for num_type === uint64_t
	// TODO: estimate max allowed value for even req_aliquot_sum
	if ((req_aliquot_sum & 1) && (req_aliquot_sum > 8589934589ull)) {
		fprintf(stderr, "Input number is too large\n");
		exit(EXIT_FAILURE);
	}
	exp_calc(req_aliquot_sum);
}

// --- main }

