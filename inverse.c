#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <memory.h>
#include <math.h>
#include <fenv.h>

#include "typedefs.h"
#include "primes.h"

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

void aliquot_inverse_cb(calc_struct calc_list[], pow_idx_type pow_count);

// returns: base^exp
num_type calc_pow(num_type base, exp_type exp) {
	assert(exp > 0);
	assert(base > 1);
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
	num_type result = 1;
	for (exp_type i=0; i<exp; ++i) {
		assert(result <= NUM_TYPE_MAX / base);
		result *= base;
		assert(result < NUM_TYPE_MAX);
		++result;
	}
	return result;
}

// TODO: test, which is faster
num_type calc_pow__old(num_type base, exp_type exp) {
	num_type result = 1;
	exp_type mask = EXP_TYPE_MAX_MASK;
	while (mask != 0 && !(exp & mask)) mask >>= 1;
	while (mask > 0) {
		assert(result <= NUM_TYPE_MAX / result);
		result *= result;
		if (exp & mask) {
			assert(result <= NUM_TYPE_MAX / base);
			result *= base;
		}
		mask >>= 1;
	}
	return result;
}

// TODO: test, which is faster
num_type calc_pow_sigma__old(num_type base, exp_type exp) {
	assert(exp > 0);
	assert(exp < EXP_TYPE_MAX);
	assert(base > 1);
	if      (exp == 1) return base + 1;
	else if (exp == 2) return (base + 1) * base + 1;
	else if (exp == 3) return ((base + 1) * base + 1) * base + 1;
	else {
		// preveinting premature overflow
		num_type pow = calc_pow(base, exp);
		num_type sum_suffix = (pow - 1) / (base - 1);
		assert(pow <= NUM_TYPE_MAX - sum_suffix);
		return pow + sum_suffix;
	}
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

// XXX max req_aliquot_sum < 2^54.19906
// (without assert_may_overflow): n = 2^55, 3^35, 5^24, 7^21: find_base_for_pow(exp = 42): p_mid   = 3: 3^41 overflows
// (within  assert_may_overflow):                  s >= 2^42: find_base_for_pow(exp = 42): p_right = 3: 3^41 overflows

// returns: 0 if not found, value>1 if found
num_type find_base_for_pow(exp_type exp, num_type req_aliquot_sum) {
	assert(exp > 1);
	// TODO: estimate rounding errors and detect overflow
	double dexp = 1.0 / ((double)(exp - 1));
	num_type p_left = llround(pow((double)req_aliquot_sum/2.0, dexp));
	num_type p_right = llrint(pow((double)req_aliquot_sum, dexp));
	assert(p_left >= 2);
	assert(p_right < NUM_TYPE_MAX);
	assert(p_left <= p_right);
	assert             (calc_aliquot_for_pow(exp, (p_left > 2 ? p_left-1 : 2)) < req_aliquot_sum);
	assert_may_overflow(calc_aliquot_for_pow(exp, p_right+1                  ) > req_aliquot_sum);
	if (p_left < 3) p_left = 3;
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
num_type calc_last_base(num_type prefix_sigma, num_type prefix_aliquot, exp_type exp, num_type req_aliquot_sum) {
	assert(exp <= 2);
	assert(req_aliquot_sum > prefix_sigma);
	num_type sub_req_sigma = req_aliquot_sum - prefix_sigma;
	num_type p_numerator;
	if (exp == 1) {
		p_numerator = sub_req_sigma;
	} else {
		assert(prefix_sigma <= NUM_TYPE_MAX / prefix_sigma);
		num_type summand1 = prefix_sigma*prefix_sigma;
		assert(prefix_aliquot <= NUM_TYPE_MAX / 4);
		assert(sub_req_sigma <= (NUM_TYPE_MAX / 4) / prefix_aliquot);
		num_type summand2 = 4*prefix_aliquot*sub_req_sigma;
		// XXX req_aliquot_sum <= 9126805519 < 2^33
		// examples:
		//     num = 2*3*38993^2, sum = 9123192222 => not overflows
		//     num = 2*3*39019^2, sum = 9135362406 => overflows
		assert(summand1 <= NUM_TYPE_MAX - summand2);
		num_type discr = summand1 + summand2;
		num_type discr_sqrt = round_sqrt(discr);
		if (discr_sqrt*discr_sqrt != discr) return 0;
		assert(discr_sqrt > prefix_sigma);
		p_numerator = discr_sqrt - prefix_sigma;
		if (p_numerator & 1) return 0;
		p_numerator >>= 1;
	}
	num_type p_denominator = prefix_aliquot;
	assert(p_numerator >= 2*p_denominator);
	if (p_numerator % p_denominator) return 0;
	return p_numerator / p_denominator;
}

num_type calc_aliquot(num_type prefix_sigma, num_type prefix_aliquot, exp_type exp, num_type prime) {
	assert(prefix_sigma > 1); assert(prefix_aliquot >= 1); assert(prefix_sigma > prefix_aliquot);
	assert(exp > 2); assert(prime >= 2);
	num_type pow = calc_pow(prime, exp);
	num_type pows_sum = (pow - 1) / (prime - 1);
	assert(pow <= NUM_TYPE_MAX / prefix_aliquot);
	num_type summand1 = prefix_aliquot * pow;
	assert(pows_sum <= NUM_TYPE_MAX / prefix_sigma);
	num_type summand2 = prefix_sigma * pows_sum;
	assert(summand1 <= NUM_TYPE_MAX - summand2);
	return summand1 + summand2;
}

// returns: 0 if not found, value>1 if found
num_type find_last_base(
		num_type prefix_sigma, num_type prefix_aliquot,
		exp_type exp,
		num_type min_base,
		num_type req_aliquot_sum
) {
	assert(exp > 2);
	assert(min_base >= 5);
	// TODO: estimate rounding errors and detect overflow
	double dexp = 1.0 / (double)exp;
	assert(prefix_sigma <= NUM_TYPE_MAX - prefix_aliquot);
	num_type sum_sigma_aliquot = prefix_sigma + prefix_aliquot;
	num_type p_left = llround(pow((double)req_aliquot_sum/(double)sum_sigma_aliquot, dexp));
	num_type p_right = llrint(pow((double)req_aliquot_sum/(double)prefix_aliquot   , dexp));
	assert(p_left >= 1);
	assert(p_right < NUM_TYPE_MAX);
	assert(p_left <= p_right);
	assert(p_left <= 2 || calc_aliquot(prefix_sigma, prefix_aliquot, exp, p_left-1 ) < req_aliquot_sum);
	assert_may_overflow(  calc_aliquot(prefix_sigma, prefix_aliquot, exp, p_right+1) > req_aliquot_sum);
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
	assert(pow_count > 0);
	exp_type exp = prime_calc_list[pow_count-1].exp;
	
	num_type prime;
	if (pow_count > 1) {
		num_type prefix_sigma = prime_calc_list[pow_count-2].prefix_sigma;
		num_type prefix_aliquot = prefix_sigma - prime_calc_list[pow_count-2].prefix_mul;
		if (exp <= 2) {
			prime = calc_last_base(prefix_sigma, prefix_aliquot, exp, req_aliquot_sum);
			assert(prime == 0 || prime > prime_calc_list[pow_count-1].prime);
		} else {
			num_type min_base = prime_calc_list[pow_count-1].prime;
			assert(min_base > 2);
			assert(min_base <= NUM_TYPE_MAX - 2);
			assert(calc_aliquot(prefix_sigma, prefix_aliquot, exp, min_base) < req_aliquot_sum);
			prime = find_last_base(prefix_sigma, prefix_aliquot, exp, min_base+2, req_aliquot_sum);
		}
	} else {
		assert(prime_calc_list[0].prime == 2);
		prime = find_base_for_pow(exp, req_aliquot_sum);
	}
	
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

// TODO
void calc_pow_and_pow_sigma(num_type prime, exp_type exp, num_type *p_pow, num_type *p_pow_sigma) {
	*p_pow       = calc_pow(prime, exp);
	*p_pow_sigma = calc_pow_sigma(prime, exp);
}

num_type prime_calc_recalc(calc_struct prime_calc_list[], pow_idx_type idx) {
	calc_struct *cur_prime_calc = &(prime_calc_list[idx]);
	calc_pow_and_pow_sigma(cur_prime_calc->prime, cur_prime_calc->exp,
		&(cur_prime_calc->pow), &(cur_prime_calc->pow_sigma));
	assert(idx == 0 || cur_prime_calc->pow <= NUM_TYPE_MAX / prime_calc_list[idx-1].prefix_mul);
	cur_prime_calc->prefix_mul =
		(idx > 0 ? prime_calc_list[idx-1].prefix_mul*cur_prime_calc->pow : cur_prime_calc->pow);
	assert(idx == 0 || cur_prime_calc->pow_sigma <= NUM_TYPE_MAX / prime_calc_list[idx-1].prefix_sigma);
	cur_prime_calc->prefix_sigma =
		(idx > 0 ? prime_calc_list[idx-1].prefix_sigma*cur_prime_calc->pow_sigma : cur_prime_calc->pow_sigma);
	assert(cur_prime_calc->prefix_sigma > cur_prime_calc->prefix_mul);
	return cur_prime_calc->prefix_sigma - cur_prime_calc->prefix_mul;
}

// inc and fill: both use maybe-primes
// returns: true - success, false - failed
bool inc_and_fill_maybe_primes(
	calc_struct prime_calc_list[], pow_idx_type pow_count,
	pow_idx_type idx, num_type req_aliquot_sum
) {
	const calc_struct *prev_prime_calc = &(prime_calc_list[idx]);
	for (; idx<pow_count; ++idx) {
		calc_struct *cur_prime_calc = &(prime_calc_list[idx]);
		next_maybe_prime_calc(prev_prime_calc, cur_prime_calc);
		num_type prefix_aliquot_sum = prime_calc_recalc(prime_calc_list, idx);
		if (prefix_aliquot_sum > req_aliquot_sum) {
			return false;
		} else if (prefix_aliquot_sum == req_aliquot_sum) {
			if (idx == pow_count-1 && is_prime_calc(cur_prime_calc)) aliquot_inverse_cb(prime_calc_list, pow_count);
			return false;
		}
		prev_prime_calc = cur_prime_calc;
	}
	return true;
}

// returns: false - next is calculated, true - end is reached
bool prime_calc_next(calc_struct prime_calc_list[], pow_idx_type pow_count, num_type req_aliquot_sum) {
	assert(pow_count >= 2);
	pow_idx_type idx = pow_count - 2;
	
	while (true) {
		bool fill_res = inc_and_fill_maybe_primes(prime_calc_list, pow_count, idx, req_aliquot_sum);
		if (fill_res) {
			while (idx < pow_count-1 && is_prime_calc( &(prime_calc_list[idx]) )) ++idx;
			if (idx == pow_count-1) return false;
		} else {
			if (idx == 0) return true;
			--idx;
		}
	}
}

void prime_calc(calc_struct exp_calc_list[], pow_idx_type pow_count, num_type req_aliquot_sum) {
//dprintf("[prime_calc][1]");
//for (pow_idx_type i=0; i<pow_count; ++i) dprintf(PRI_EXP_TYPE " ", exp_calc_list[i].exp);
//dprintf("\n");
	static calc_struct prime_calc_list[MAX_POW_COUNT];
	memcpy(prime_calc_list, exp_calc_list, sizeof(exp_calc_list[0])*pow_count);
	
	if (pow_count == 1) {
		find_last_prime_calc(prime_calc_list, pow_count, req_aliquot_sum);
		return;
	}
	
	for (pow_idx_type i=0; i<pow_count; ++i) {
		prime_calc_list[i].prime_idx = i;
		prime_calc_list[i].is_prime_in_list = true;
		prime_calc_list[i].prime_status = PRIME_STATUS_PRIME;
	}
	
	bool is_end_reached;
	do {
#ifndef NDEBUG
		for (pow_idx_type i=0; i<pow_count-1; ++i) assert(prime_calc_list[i].prime_status == PRIME_STATUS_PRIME);
#endif
//dprintf("[prime_calc][1] {\n");
//calc_list_dprint(prime_calc_list, pow_count);
//dprintf("prime_calc][1] }\n");
		find_last_prime_calc(prime_calc_list, pow_count, req_aliquot_sum);
		is_end_reached = prime_calc_next(prime_calc_list, pow_count, req_aliquot_sum);
	} while (!is_end_reached);
}

// --- iterate primes --- }
// --- iterate exps --- {

num_type exp_calc_inc_exp(calc_struct exp_calc_list[], pow_idx_type idx) {
	calc_struct *cur_exp_calc = &(exp_calc_list[idx]);
	++cur_exp_calc->exp;
	assert(cur_exp_calc->exp <= MAX_EXP);
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
//dprintf("[exp_calc_next][1] req_aliquot_sum = " PRI_NUM_TYPE "\n", req_aliquot_sum);
//dprintf("[exp_calc_next][2] pow_count = " PRI_POW_IDX_TYPE "\n", pow_count);
	pow_idx_type inc_idx = pow_count-1;
	while (true) {
		while (true) {
//dprintf("[exp_calc_next][3] inc_idx=" PRI_POW_IDX_TYPE "\n", inc_idx);
			num_type prefix_aliquot_sum = exp_calc_inc_exp(exp_calc_list, inc_idx);
//dprintf("[exp_calc_next][4] {\n");
//calc_list_dprint(exp_calc_list, inc_idx+1);
//dprintf("[exp_calc_next][4] }\n");
//dprintf("[exp_calc_next][5] prefix_aliquot_sum=" PRI_NUM_TYPE "\n", prefix_aliquot_sum);
			if (prefix_aliquot_sum < req_aliquot_sum) {
				break;
			} else if (prefix_aliquot_sum == req_aliquot_sum && inc_idx == pow_count-1) {
				aliquot_inverse_cb(exp_calc_list, pow_count);
			}
			if (inc_idx == 0) {
//dprintf("[exp_calc_next]return true[1]\n");
				return true;
			}
			--inc_idx;
		}
		pow_idx_type fill_idx;
		for (fill_idx=inc_idx+1; fill_idx<pow_count; ++fill_idx) {
//dprintf("[exp_calc_next][6] fill_idx=" PRI_POW_IDX_TYPE "\n", fill_idx);
			num_type prefix_aliquot_sum = exp_calc_fill_exp(exp_calc_list, fill_idx);
//dprintf("[exp_calc_next][7] {\n");
//calc_list_dprint(exp_calc_list, fill_idx+1);
//dprintf("[exp_calc_next][7] }\n");
//dprintf("[exp_calc_next][8] prefix_aliquot_sum=" PRI_NUM_TYPE "\n", prefix_aliquot_sum);
			if (prefix_aliquot_sum == req_aliquot_sum) {
				if (fill_idx == pow_count-1) aliquot_inverse_cb(exp_calc_list, pow_count);
				break;
			} else if (prefix_aliquot_sum > req_aliquot_sum) {
				break;
			}
		}
//dprintf("[exp_calc_next][9] fill_idx=" PRI_POW_IDX_TYPE "\n", fill_idx);
		if (fill_idx == pow_count) break;
		if (inc_idx == 0) {
//dprintf("[exp_calc_next]return true[2]\n");
			return true;
		}
//dprintf("[exp_calc_next][10]\n");
		--inc_idx;
	}
//dprintf("[exp_calc_next]return false\n");
	return false;
}

void exp_calc(num_type req_aliquot_sum) {
	assert(req_aliquot_sum > 1);
	
	static calc_struct exp_calc_list[MAX_POW_COUNT];
	for (pow_idx_type i=0; i<MAX_POW_COUNT; ++i) exp_calc_list[i].prime = first_primes[i];
	
	for (pow_idx_type pow_count=1;; ++pow_count) {
//dprintf("[exp_calc][1] pow_count=" PRI_POW_IDX_TYPE "\n", pow_count);
		num_type prefix_aliquot_sum;
		for (pow_idx_type fill_idx=0; fill_idx<pow_count; ++fill_idx) {
			prefix_aliquot_sum = exp_calc_fill_any_exp(exp_calc_list, fill_idx);
//dprintf("[exp_calc][2] {\n");
//calc_list_dprint(exp_calc_list, fill_idx+1);
//dprintf("[exp_calc][2] }\n");
//dprintf("[exp_calc][3] prefix_aliquot_sum=" PRI_NUM_TYPE "\n", prefix_aliquot_sum);
			assert(fill_idx == pow_count-1 || prefix_aliquot_sum < req_aliquot_sum);
		}
		if (prefix_aliquot_sum == req_aliquot_sum) aliquot_inverse_cb(exp_calc_list, pow_count);
		if (prefix_aliquot_sum >= req_aliquot_sum) break;
		
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

void aliquot_inverse_init() {
	assert(sizeof(first_primes)/sizeof(first_primes[0]) >= MAX_POW_COUNT);

	primes_array = primes_construct();
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

void aliquot_inverse(num_type req_aliquot_sum) {
	if (req_aliquot_sum <= 1) {
		fprintf(stderr, "Wrong input: " PRI_NUM_TYPE "\n", req_aliquot_sum);
		exit(EXIT_FAILURE);
	}
	exp_calc(req_aliquot_sum);
}

// --- main }

//==============================

num_type req_aliquot_sum_saved;

void aliquot_inverse_cb(calc_struct calc_list[], pow_idx_type pow_count) {
//dprintf("[aliquot_inverse_cb][1]");
//for (pow_idx_type i=0; i<pow_count; ++i) dprintf(PRI_NUM_TYPE "^" PRI_EXP_TYPE " ", calc_list[i].prime, calc_list[i].exp);
//dprintf("\n");
	
	num_type value = 1;
	for (pow_idx_type i=0; i<pow_count; ++i) {
		num_type pow = calc_pow(calc_list[i].prime, calc_list[i].exp);
		assert(value <= NUM_TYPE_MAX / pow);
		value *= pow;
	}
	printf(PRI_NUM_TYPE "\t" PRI_NUM_TYPE "\n", req_aliquot_sum_saved, value);
}

void run() {
	aliquot_inverse_init();
	
	num_type req_aliquot_sum;
	int scanf_res = scanf(SCN_NUM_TYPE, &req_aliquot_sum);
	assert(scanf_res == 1);
	req_aliquot_sum_saved = req_aliquot_sum;
	aliquot_inverse(req_aliquot_sum);
	
	aliquot_inverse_terminate();
}

int main() {
	run();
	return 0;
}

