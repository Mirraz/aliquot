#include <assert.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

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

// ===================================

#define round_sqrt(n) round(sqrt((double)n))

size_t fill_primes(num_type primes[], size_t primes_size, num_type max_num) {
	if (primes_size == 0 || max_num < 2) return 0;
	primes[0] = 2;
	size_t primes_count = 1;
	num_type n = 3;
	num_type n_sqrt = 1;
	while (primes_count < primes_size && n <= max_num) {
		bool is_prime = true;
		for (size_t i=0; primes[i]<=n_sqrt; ++i) {
			if (!(n % primes[i])) {
				is_prime = false;
				break;
			}
		}
		assert(primes_count < primes_size);
		if (is_prime) primes[primes_count++] = n;
		n += 2;
		n_sqrt = round_sqrt(n);
	}
	return primes_count;
}

void factorize_cb(num_type prime, exp_type exp);

num_type *factorize_primes;
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

// ===================================

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

// ----------------

#define primes_size (1024*1024)
num_type primes[primes_size];

void primes_init(size_t primes_count) {
	assert(primes_count <= primes_size);
	factorize_primes_count = fill_primes(primes, primes_count, NUM_TYPE_MAX);
	assert(factorize_primes_count == primes_count);
	factorize_primes = primes;
}

// ----------------

void run_aliquot_sequence() {
	primes_init(1024*32);
	
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
	primes_init(1024*32);
	
	num_type n;
	for (n=2; n<=16777216; ++n) {
		printf("%" PRIuFAST64 "\t%"  PRIuFAST64 "\n", n, calc_aliquot_sum(n));
	}
}

void run_aliquot_sum() {
	primes_init(1024*32);
	
	num_type n;
	int scanf_res = scanf("%" SCNuFAST64, &n);
	assert(scanf_res == 1);
	
	printf("%" PRIuFAST64 "\n", calc_aliquot_sum(n));
}

int main() {
	run_aliquot_sum();
	return 0;
}

