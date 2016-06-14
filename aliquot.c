#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// exit, EXIT_FAILURE
#include <stddef.h>		// size_t
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#include "primes.h"

typedef uint_fast64_t num_type;
#define NUM_TYPE_MAX UINT64_MAX

typedef uint_fast8_t exp_type;

// returns: base^exp + base^(exp-1) + ... + base + 1
static num_type calc_pow_sigma(num_type base, exp_type exp) {
	num_type result = 1;
	for (exp_type i=0; i<exp; ++i) {
		assert(result <= NUM_TYPE_MAX / base);
		result *= base;
		assert(result < NUM_TYPE_MAX);
		++result;
	}
	return result;
}

#define round_sqrt(n) round(sqrt((double)n))

void factorize_cb(num_type prime, exp_type exp);

prime_type *factorize_primes;
size_t factorize_primes_count;

void factorize(num_type n) {
	assert(n > 0);
	if (n == 1) return;
	
	if (!(n & 1)) {
		exp_type exp = 0;
		do {
			n >>= 1;
			++exp;
		} while (!(n & 1));
		factorize_cb(2, exp);
	}
	
	num_type n_sqrt = round_sqrt(n);
	num_type p = 3;
	size_t idx = 2;
	
	while (p <= n_sqrt) {
		if (!(n % p)) {
			exp_type exp = 0;
			do {
				n /= p;
				++exp;
			} while (!(n % p));
			factorize_cb(p, exp);
			n_sqrt = round_sqrt(n);
		}
		if (idx < factorize_primes_count) {
			p = factorize_primes[idx++];
		} else {
			p += 2;
		}
	}
	if (n != 1) factorize_cb(n, 1);
}

num_type sigma;

void factorize_cb(num_type prime, exp_type exp) {
	assert(exp > 0);
	num_type mul = calc_pow_sigma(prime, exp);
	assert(sigma <= NUM_TYPE_MAX / mul);
	sigma *= mul;
}

num_type calc_sigma(num_type num) {
	sigma = 1;
	factorize(num);
	return sigma;
}

num_type calc_aliquot_sum(num_type num) {
	num_type sigma = calc_sigma(num);
	assert(sigma > num);
	return sigma - num;
}

primes_array_struct *primes_array;

void aliquot_init(const char *primes_filename) {
	primes_array = primes_construct(primes_filename);
	factorize_primes = primes_get_array(primes_array);
	factorize_primes_count = primes_get_count(primes_array);
}

void aliquot_terminate() {
	primes_destruct(primes_array);
}

// ----------------

void run_aliquot_sequence() {
	num_type n;
	int scanf_res = scanf("%" SCNuFAST64, &n);
	assert(scanf_res == 1);
	printf("%" PRIuFAST64 "\n", n);
	
	while (n != 1) {
		n = calc_aliquot_sum(n);
		printf("%" PRIuFAST64 "\n", n);
	}
}

void run_aliquot_sum_table() {
	num_type n;
	for (n=2; n<=16777216; ++n) {
		//printf("%" PRIuFAST64 "\t", n);
		printf("%" PRIuFAST64 "\n", calc_aliquot_sum(n));
	}
}

void run_aliquot_sum() {
	num_type n;
	int scanf_res = scanf("%" SCNuFAST64, &n);
	assert(scanf_res == 1);
	
	printf("%" PRIuFAST64 "\n", calc_aliquot_sum(n));
}

void run(const char *primes_filename) {
	aliquot_init(primes_filename);
	
	run_aliquot_sum();
	
	aliquot_terminate();
}

int main(int argc, char *argv[]) {
	if (argc != 2) {
		fprintf(stderr, "Too few arguments\n");
		exit(EXIT_FAILURE);
	}
	run(argv[1]);
	return 0;
}

